"""2D PES for an upright H2 above Si atom 13 of Si(111): V(R, h).

R = H-H distance, h = height of the lower H above Si 13 (both in A); the upper H
is at h + R. Energies in eV relative to E(h = 8, R = 8), i.e. two H atoms far from
each other and from the surface.

    V = V0 + I_gp            (PBE)
      + D_gp                 (D3 dispersion, optional)
    V0 = V_HH(R) + V_H(h) + V_H(h + R)

V_HH is the h = 8 row (gas-phase H2), V_H the R = 8 column (one H against the
surface); I_gp is a GPR on the remaining 3-body part, which is zero on both cuts
and decays to the prior mean 0 off the grid. I_gp is switched off smoothly between
TAPER[0] and TAPER[1] in both R and h (the data show < 1 meV there), so the two
asymptotic channels are exactly the 1D cuts. Gradients are analytic.
"""
import numpy as np
from scipy.interpolate import CubicSpline
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

R_INF = 8.0
H_INF = 8.0
R_GAS = 0.75          # grid point used as the gas-phase H2 reference
SQ5 = np.sqrt(5.0)
TAPER = (5.5, 8.0)    # A; 3-body GP switched off over this range in R and h


class Curve1D:
    """Cubic spline through (x, y), constant above x[-1]. Below x[0] an exponential
    wall C + A exp(-k (x - x0)) matching the spline's value, slope and curvature."""

    def __init__(self, x, y):
        o = np.argsort(x)
        x, y = np.asarray(x, float)[o], np.asarray(y, float)[o]
        self.s = CubicSpline(x, y, bc_type=("not-a-knot", (1, 0.0)))
        self.x0, self.x1, self.y0, self.y1 = x[0], x[-1], y[0], y[-1]
        self.d0, c0 = float(self.s(x[0], 1)), float(self.s(x[0], 2))
        if self.d0 < 0 and c0 > 0:
            self.k = -c0 / self.d0
            self.A = -self.d0 / self.k
            self.C = self.y0 - self.A
        else:
            self.k = None

    def __call__(self, x, nu=0):
        x = np.asarray(x, float)
        out = np.asarray(self.s(np.clip(x, self.x0, self.x1), nu), float)
        lo, hi = x < self.x0, x > self.x1
        if self.k is not None:
            ex = self.A * np.exp(-self.k * (x[lo] - self.x0))
            out[lo] = self.C + ex if nu == 0 else -self.k * ex
        else:
            out[lo] = self.y0 + self.d0 * (x[lo] - self.x0) if nu == 0 else self.d0
        out[hi] = self.y1 if nu == 0 else 0.0
        return out


def _taper(x, nu=0):
    """1 below TAPER[0], 0 above TAPER[1], cosine in between (C1)."""
    a, b = TAPER
    t = np.clip((np.asarray(x, float) - a) / (b - a), 0.0, 1.0)
    if nu == 0:
        return 0.5 * (1 + np.cos(np.pi * t))
    return np.where((t > 0) & (t < 1), -0.5 * np.pi / (b - a) * np.sin(np.pi * t), 0.0)


def _kernel():
    # targets are not normalized, so the GP prior mean (and the far-field limit) is 0 eV
    return (ConstantKernel(1.0, (1e-6, 1e3)) * Matern([1.0, 1.0], (1e-2, 1e2), nu=2.5)
            + WhiteKernel(1e-4, (1e-10, 1e-1)))


def _gp_mean_grad(gp, x):
    """d(mean)/dx for a fitted GPR with kernel C * Matern(nu=2.5) + White."""
    c = gp.kernel_.k1.k1.constant_value
    ls = np.broadcast_to(gp.kernel_.k1.k2.length_scale, (x.shape[1],))
    a = gp.alpha_.ravel()
    out = np.empty_like(x)
    for i in range(0, len(x), 2000):
        diff = x[i:i + 2000, None, :] - gp.X_train_[None, :, :]
        d = np.sqrt(np.sum((diff / ls) ** 2, axis=-1))
        f = -(5.0 / 3.0) * c * (1 + SQ5 * d) * np.exp(-SQ5 * d)
        out[i:i + 2000] = np.einsum("nm,nmk,m->nk", f, diff / ls ** 2, a)
    return out * np.atleast_1d(gp._y_train_std)


class PES:
    """transform: 'none' -> GP inputs (R, h); 'exp' -> (exp(-R/lam), exp(-h/lam))."""

    def __init__(self, transform="none", lam=1.0, wall_ev=8.0, wall_sigma=0.05,
                 restarts=10, seed=0):
        self.transform, self.lam = transform, lam
        self.wall_ev, self.wall_sigma = wall_ev, wall_sigma
        self.restarts, self.seed = restarts, seed

    # ---------- coordinates
    def _x(self, R, h):
        X = np.column_stack([np.ravel(R), np.ravel(h)]).astype(float)
        return np.exp(-X / self.lam) if self.transform == "exp" else X

    def _dxdq(self, R, h):
        """d(feature)/d(R, h), per column."""
        X = np.column_stack([np.ravel(R), np.ravel(h)]).astype(float)
        return -np.exp(-X / self.lam) / self.lam if self.transform == "exp" else np.ones_like(X)

    # ---------- fitting
    def fit(self, df):
        """df: columns h, R, E_pbe, Edisp (absolute eV). Must contain the h = 8 row
        and the R = 8 column, which define the baseline."""
        ref = df[(np.isclose(df.h, H_INF)) & (np.isclose(df.R, R_INF))]
        if len(ref) != 1:
            raise ValueError("need the (h = 8, R = 8) point")
        self.e_ref, self.d_ref = float(ref.E_pbe.iloc[0]), float(ref.Edisp.iloc[0])
        e = df.E_pbe.values - self.e_ref
        dd = df.Edisp.values - self.d_ref

        row = df[np.isclose(df.h, H_INF)]
        col = df[np.isclose(df.R, R_INF)]
        self.v_hh = Curve1D(row.R.values, row.E_pbe.values - self.e_ref)
        self.v_h = Curve1D(col.h.values, col.E_pbe.values - self.e_ref)

        R, h = df.R.values, df.h.values
        y = e - self.v0(R, h)
        X = self._x(R, h)

        # larger noise on the high-energy wall so it does not set the hyperparameters
        etot_gas = (e + dd) - self.gas_ref(True, _df=df)
        sig = np.where(etot_gas > self.wall_ev, self.wall_sigma, 0.0)
        self.gp = GaussianProcessRegressor(_kernel(), alpha=sig ** 2 + 1e-10,
                                           normalize_y=False, random_state=self.seed,
                                           n_restarts_optimizer=self.restarts).fit(X, y)
        self.gp_d3 = GaussianProcessRegressor(_kernel(), alpha=1e-10, normalize_y=False,
                                              random_state=self.seed,
                                              n_restarts_optimizer=self.restarts).fit(X, dd)
        self._gas = {True: self.gas_ref(True, _df=df), False: self.gas_ref(False, _df=df)}
        return self

    def gas_ref(self, d3=True, _df=None):
        """Energy of gas-phase H2 at R = 0.75 A, h = 8 (relative to the atoms)."""
        if _df is None:
            return self._gas[d3]
        p = _df[np.isclose(_df.h, H_INF) & np.isclose(_df.R, R_GAS)]
        return float((p.E_pbe - self.e_ref + (p.Edisp - self.d_ref if d3 else 0)).iloc[0])

    # ---------- evaluation
    def v0(self, R, h, grad=False):
        R, h = np.ravel(R).astype(float), np.ravel(h).astype(float)
        if not grad:
            return self.v_hh(R) + self.v_h(h) + self.v_h(h + R)
        dR = self.v_hh(R, 1) + self.v_h(h + R, 1)
        dh = self.v_h(h, 1) + self.v_h(h + R, 1)
        return np.column_stack([dR, dh])

    def energy(self, R, h, d3=True, ref="atoms"):
        """ref='atoms': relative to two far-apart H atoms; ref='gas': to gas-phase H2."""
        shape = np.shape(np.broadcast_arrays(R, h)[0])
        R, h = [np.ravel(a) for a in np.broadcast_arrays(R, h)]
        X = self._x(R, h)
        E = self.v0(R, h) + _taper(R) * _taper(h) * self.gp.predict(X)
        if d3:
            E = E + self.gp_d3.predict(X)
        if ref == "gas":
            E = E - self.gas_ref(d3)
        return E.reshape(shape)

    def gradient(self, R, h, d3=True):
        """(dV/dR, dV/dh) in eV/A, shape (..., 2)."""
        shape = np.shape(np.broadcast_arrays(R, h)[0])
        R, h = [np.ravel(a) for a in np.broadcast_arrays(R, h)]
        X, J = self._x(R, h), self._dxdq(R, h)
        I, gI = self.gp.predict(X), _gp_mean_grad(self.gp, X) * J
        tR, th = _taper(R), _taper(h)
        g = self.v0(R, h, grad=True) + (tR * th)[:, None] * gI + I[:, None] * np.column_stack(
            [_taper(R, 1) * th, tR * _taper(h, 1)])
        if d3:
            g = g + _gp_mean_grad(self.gp_d3, X) * J
        return g.reshape(shape + (2,))

    def sigma(self, R, h):
        """GPR standard deviation of the PBE 3-body part (eV)."""
        shape = np.shape(np.broadcast_arrays(R, h)[0])
        R, h = [np.ravel(a) for a in np.broadcast_arrays(R, h)]
        sd = self.gp.predict(self._x(R, h), return_std=True)[1]
        return (_taper(R) * _taper(h) * sd).reshape(shape)

    # ---------- other coordinates
    def energy_z(self, z_low, z_up, **kw):
        """Energy from the two H heights above Si 13 (z_up > z_low)."""
        return self.energy(np.asarray(z_up) - np.asarray(z_low), z_low, **kw)

    def energy_com(self, R, Z, **kw):
        """Energy from H-H distance R and the H2 centre-of-mass height Z."""
        return self.energy(R, np.asarray(Z) - np.asarray(R) / 2, **kw)
