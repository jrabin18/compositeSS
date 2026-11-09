/*
 * compositess.js: a JavaScript port of the design, planning and analysis
 * functions of the compositeSS R package (version 0.1.0.9000).
 *
 * Function names follow the R package (csDesign = cs_design, and so on), and
 * each function reproduces the R computation step by step. tests/validate.html
 * compares the results with values computed by the R package itself
 * (tests/make_reference.R). Monte Carlo results (assurance, simulation,
 * bootstrap) use a different random number generator from R, so they agree
 * with R only up to Monte Carlo error.
 *
 * MIT licence. Works in the browser (global `CS`) and in Node (module.exports).
 */
(function (root) {
  'use strict';

  // ------------------------------------------------------------------
  // Distribution functions
  // ------------------------------------------------------------------
  const M_1_SQRT_2PI = 0.398942280401432677939946059934;
  const LOG_SQRT_2PI = 0.918938533204672741780329736406;

  // Normal distribution function: W. J. Cody's algorithm, as in R's pnorm.
  // Returns [lower tail, upper tail], each with full relative accuracy.
  function pnormBoth(x) {
    if (Number.isNaN(x)) return [NaN, NaN];
    if (x === Infinity) return [1, 0];
    if (x === -Infinity) return [0, 1];
    const a = [2.2352520354606839287, 161.02823106855587881, 1067.6894854603709582,
      18154.981253343561249, 0.065682337918207449113];
    const b = [47.20258190468824187, 976.09855173777669322, 10260.932208618978205,
      45507.789335026729956];
    const c = [0.39894151208813466764, 8.8831497943883759412, 93.506656132177855979,
      597.27027639480026226, 2494.5375852903726711, 6848.1904505362823326,
      11602.651437647350124, 9842.7148383839780218, 1.0765576773720192317e-8];
    const d = [22.266688044328115691, 235.38790178262499861, 1519.377599407554805,
      6485.558298266760755, 18615.571640885098091, 34900.952721145977266,
      38912.003286093271411, 19685.429676859990727];
    const p = [0.21589853405795699, 0.1274011611602473639, 0.022235277870649807,
      0.001421619193227893466, 2.9112874951168792e-5, 0.02307344176494017303];
    const q = [1.28426009614491121, 0.468238212480865118, 0.0659881378689285515,
      0.00378239633202758244, 7.29751555083966205e-5];
    const y = Math.abs(x);
    let cum, ccum, xnum, xden, temp, xsq, del;
    if (y <= 0.67448975) {
      if (y > 1.1102230246251565e-16) {
        xsq = x * x; xnum = a[4] * xsq; xden = xsq;
        for (let i = 0; i < 3; i++) { xnum = (xnum + a[i]) * xsq; xden = (xden + b[i]) * xsq; }
      } else { xnum = 0; xden = 0; }
      temp = x * (xnum + a[3]) / (xden + b[3]);
      return [0.5 + temp, 0.5 - temp];
    }
    if (y <= 5.656854249492380195206754896838) {
      xnum = c[8] * y; xden = y;
      for (let i = 0; i < 7; i++) { xnum = (xnum + c[i]) * y; xden = (xden + d[i]) * y; }
      temp = (xnum + c[7]) / (xden + d[7]);
    } else {
      if (y > 38.4674) return x > 0 ? [1, 0] : [0, 1];
      xsq = 1 / (x * x); xnum = p[5] * xsq; xden = xsq;
      for (let i = 0; i < 4; i++) { xnum = (xnum + p[i]) * xsq; xden = (xden + q[i]) * xsq; }
      temp = xsq * (xnum + p[4]) / (xden + q[4]);
      temp = (M_1_SQRT_2PI - temp) / y;
    }
    xsq = Math.trunc(y * 16) / 16; del = (y - xsq) * (y + xsq);
    cum = Math.exp(-xsq * xsq * 0.5) * Math.exp(-del * 0.5) * temp;
    ccum = 1 - cum;
    if (x > 0) { const t = cum; cum = ccum; ccum = t; }
    return [cum, ccum];
  }
  function pnorm(x, mean = 0, sd = 1, lower = true) {
    const r = pnormBoth((x - mean) / sd);
    return lower ? r[0] : r[1];
  }
  function dnorm(x, mean = 0, sd = 1) {
    const z = (x - mean) / sd;
    return Math.exp(-0.5 * z * z - LOG_SQRT_2PI) / sd;
  }
  // Normal quantile: Wichura's algorithm AS 241, as in R's qnorm.
  function poly(cf, r) { let s = cf[cf.length - 1]; for (let i = cf.length - 2; i >= 0; i--) s = s * r + cf[i]; return s; }
  const QA = [3.387132872796366608, 133.14166789178437745, 1971.5909503065514427,
    13731.693765509461125, 45921.953931549871457, 67265.770927008700853,
    33430.575583588128105, 2509.0809287301226727];
  const QB = [1, 42.313330701600911252, 687.1870074920579083, 5394.1960214247511077,
    21213.794301586595867, 39307.89580009271061, 28729.085735721942674,
    5226.495278852545925];
  const QC = [1.42343711074968357734, 4.6303378461565452959, 5.7694972214606914055,
    3.64784832476320460504, 1.27045825245236838258, 0.24178072517745061177,
    0.0227238449892691845833, 7.7454501427834140764e-4];
  const QD = [1, 2.05319162663775882187, 1.6763848301838038494, 0.68976733498510000455,
    0.14810397642748007459, 0.0151986665636164571966, 5.475938084995344946e-4,
    1.05075007164441684324e-9];
  const QE = [6.6579046435011037772, 5.4637849111641143699, 1.7848265399172913358,
    0.29656057182850489123, 0.026532189526576123093, 0.0012426609473880784386,
    2.71155556874348757815e-5, 2.01033439929228813265e-7];
  const QF = [1, 0.59983220655588793769, 0.13692988092273580531,
    0.0148753612908506148525, 7.868691311456132591e-4, 1.8463183175100546818e-5,
    1.4215117583164458887e-7, 2.04426310338993978564e-15];
  function qnorm(pr, mean = 0, sd = 1, lower = true) {
    if (Number.isNaN(pr) || pr < 0 || pr > 1) return NaN;
    // the quantile of the upper tail is minus the quantile of the lower tail
    if (!lower) return mean - sd * qnorm(pr, 0, 1, true);
    if (pr === 0) return -Infinity;
    if (pr === 1) return Infinity;
    const qq = pr - 0.5;
    let val;
    if (Math.abs(qq) <= 0.425) {
      const r = 0.180625 - qq * qq;
      val = qq * poly(QA, r) / poly(QB, r);
    } else {
      let r = qq < 0 ? pr : 1 - pr;
      r = Math.sqrt(-Math.log(r));
      if (r <= 5) { r -= 1.6; val = poly(QC, r) / poly(QD, r); }
      else { r -= 5; val = poly(QE, r) / poly(QF, r); }
      if (qq < 0) val = -val;
    }
    return mean + sd * val;
  }
  // log gamma: Lanczos approximation (g = 7, 9 terms)
  const LZ = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
    771.32342877765313, -176.61502916214059, 12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  function lgamma(x) {
    if (x < 0.5) return Math.log(Math.PI / Math.abs(Math.sin(Math.PI * x))) - lgamma(1 - x);
    x -= 1;
    let s = LZ[0];
    const t = x + 7.5;
    for (let i = 1; i < 9; i++) s += LZ[i] / (x + i);
    return LOG_SQRT_2PI + (x + 0.5) * Math.log(t) - t + Math.log(s);
  }
  function betacf(a, b, x) {
    const FPMIN = 1e-300, EPS = 1e-16;
    const qab = a + b, qap = a + 1, qam = a - 1;
    let c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    d = 1 / d;
    let h = d;
    for (let m = 1; m <= 10000; m++) {
      const m2 = 2 * m;
      let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d;
      const del = d * c;
      h *= del;
      if (Math.abs(del - 1) < EPS) break;
    }
    return h;
  }
  // regularized incomplete beta I_x(a, b); y = 1 - x supplied to keep accuracy
  function pbeta(x, a, b, y) {
    if (y === undefined) y = 1 - x;
    if (x <= 0) return 0;
    if (y <= 0) return 1;
    const lbt = lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log(y);
    if (x < (a + 1) / (a + b + 2)) return Math.exp(lbt) * betacf(a, b, x) / a;
    return 1 - Math.exp(lbt) * betacf(b, a, y) / b;
  }
  // Student t distribution function
  function pt(t, df, lower = true) {
    if (Number.isNaN(t)) return NaN;
    if (t === Infinity) return lower ? 1 : 0;
    if (t === -Infinity) return lower ? 0 : 1;
    if (df > 4e5) return pnorm(t, 0, 1, lower);
    const t2 = t * t;
    // tail = P(T > |t|)
    let tail;
    if (t2 < df) {
      // x = df/(df+t^2) near 1: use I_{t^2/(df+t^2)}(1/2, df/2) for accuracy
      const xx = t2 / (df + t2);
      tail = 0.5 * (1 - pbeta(xx, 0.5, df / 2, df / (df + t2)));
    } else {
      tail = 0.5 * pbeta(df / (df + t2), df / 2, 0.5, t2 / (df + t2));
    }
    const upper = t > 0 ? tail : 1 - tail;
    return lower ? 1 - upper : upper;
  }
  // P(T <= t) on the normal scale, qnorm(pt(t, df)), computed from the
  // smaller tail so that large statistics keep their accuracy
  function tToZ(t, df) {
    return t >= 0 ? -qnorm(pt(t, df, false)) : qnorm(pt(t, df, true));
  }

  // ------------------------------------------------------------------
  // Numerical tools
  // ------------------------------------------------------------------
  const XGK = [0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
    0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
    0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
    0.207784955007898467600689403773245, 0];
  const WGK = [0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
    0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
    0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
    0.204432940075298892414161999234649, 0.209482141084727828012999174891714];
  const WG = [0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
    0.381830050505118944950369775488975, 0.417959183673469387755102040816327];
  function gk15(f, a, b) {
    const c = (a + b) / 2, h = (b - a) / 2;
    const fc = f(c);
    let k = fc * WGK[7], g = fc * WG[3];
    for (let j = 0; j < 7; j++) {
      const dx = h * XGK[j];
      const s = f(c - dx) + f(c + dx);
      k += WGK[j] * s;
      if (j % 2 === 1) g += WG[(j - 1) / 2] * s;
    }
    return { k: k * h, err: Math.abs((k - g) * h) };
  }
  // adaptive Gauss-Kronrod quadrature (global subdivision of the worst interval)
  function integrate(f, a, b, relTol = 1e-10, absTol = 0, maxSub = 500) {
    if (!(b > a)) return 0;
    const first = gk15(f, a, b);
    const parts = [{ a, b, k: first.k, err: first.err }];
    let total = first.k, err = first.err;
    for (let it = 0; it < maxSub; it++) {
      if (err <= Math.max(absTol, relTol * Math.abs(total)) || err < 1e-300) break;
      let wi = 0;
      for (let i = 1; i < parts.length; i++) if (parts[i].err > parts[wi].err) wi = i;
      const w = parts[wi], m = (w.a + w.b) / 2;
      const L = gk15(f, w.a, m), R = gk15(f, m, w.b);
      parts.splice(wi, 1, { a: w.a, b: m, k: L.k, err: L.err }, { a: m, b: w.b, k: R.k, err: R.err });
      total = 0; err = 0;
      for (const p of parts) { total += p.k; err += p.err; }
    }
    return total;
  }
  // Brent's root finder (as R's uniroot / zeroin)
  function uniroot(f, lower, upper, tol = 1e-10, maxit = 1000) {
    let a = lower, b = upper, fa = f(a), fb = f(b);
    if (fa === 0) return a;
    if (fb === 0) return b;
    if (fa * fb > 0) throw new Error('f() values at end points not of opposite sign');
    let c = a, fc = fa, d = b - a, e = d;
    for (let it = 0; it < maxit; it++) {
      if (Math.abs(fc) < Math.abs(fb)) { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
      const tol1 = 2 * 2.220446e-16 * Math.abs(b) + 0.5 * tol;
      const xm = 0.5 * (c - b);
      if (Math.abs(xm) <= tol1 || fb === 0) return b;
      if (Math.abs(e) >= tol1 && Math.abs(fa) > Math.abs(fb)) {
        let p, q, r;
        const s = fb / fa;
        if (a === c) { p = 2 * xm * s; q = 1 - s; }
        else {
          q = fa / fc; r = fb / fc;
          p = s * (2 * xm * q * (q - r) - (b - a) * (r - 1));
          q = (q - 1) * (r - 1) * (s - 1);
        }
        if (p > 0) q = -q; else p = -p;
        if (2 * p < Math.min(3 * xm * q - Math.abs(tol1 * q), Math.abs(e * q))) { e = d; d = p / q; }
        else { d = xm; e = d; }
      } else { d = xm; e = d; }
      a = b; fa = fb;
      b += Math.abs(d) > tol1 ? d : (xm > 0 ? tol1 : -tol1);
      fb = f(b);
      if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) { c = a; fc = fa; d = b - a; e = d; }
    }
    return b;
  }
  // upper quantile of the chi-square distribution with 4 df:
  // P(X > x) = exp(-x/2) (1 + x/2)
  function qchisq4Upper(alpha) {
    const u = uniroot(v => Math.log1p(v) - v - Math.log(alpha), 0, 2000, 1e-14);
    return 2 * u;
  }

  // ------------------------------------------------------------------
  // Random numbers (seeded): xoshiro128** with SplitMix32 seeding
  // ------------------------------------------------------------------
  function makeRng(seed) {
    let x = (Math.floor(Math.abs(seed)) >>> 0) || 0x9e3779b9;
    const sm = () => {
      x = (x + 0x9e3779b9) >>> 0;
      let z = x;
      z = Math.imul(z ^ (z >>> 16), 0x85ebca6b);
      z = Math.imul(z ^ (z >>> 13), 0xc2b2ae35);
      return (z ^ (z >>> 16)) >>> 0;
    };
    const s = [sm(), sm(), sm(), sm()];
    const rotl = (v, k) => (v << k) | (v >>> (32 - k));
    const next = () => {
      const r = Math.imul(rotl(Math.imul(s[1], 5), 7), 9) >>> 0;
      const t = s[1] << 9;
      s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 11);
      return r;
    };
    let spare = null;
    const rng = {
      unif() { return ((next() >>> 5) * 67108864 + (next() >>> 6) + 0.5) / 9007199254740992; },
      norm(mu = 0, sd = 1) {
        if (spare !== null) { const z = spare; spare = null; return mu + sd * z; }
        let u, v, r2;
        do { u = 2 * rng.unif() - 1; v = 2 * rng.unif() - 1; r2 = u * u + v * v; } while (r2 >= 1 || r2 === 0);
        const f = Math.sqrt(-2 * Math.log(r2) / r2);
        spare = v * f;
        return mu + sd * u * f;
      },
      gamma(a) {
        if (a < 1) return rng.gamma(a + 1) * Math.pow(rng.unif(), 1 / a);
        const d = a - 1 / 3, c = 1 / Math.sqrt(9 * d);
        for (;;) {
          let z, v;
          do { z = rng.norm(); v = 1 + c * z; } while (v <= 0);
          v = v * v * v;
          const u = rng.unif();
          if (u < 1 - 0.0331 * z * z * z * z) return d * v;
          if (Math.log(u) < 0.5 * z * z + d * (1 - v + Math.log(v))) return d * v;
        }
      },
      beta(a, b) { const g1 = rng.gamma(a), g2 = rng.gamma(b); return g1 / (g1 + g2); },
      int(n) { return Math.floor(rng.unif() * n); }
    };
    return rng;
  }

  // ------------------------------------------------------------------
  // Small helpers
  // ------------------------------------------------------------------
  const sum = v => v.reduce((s, x) => s + x, 0);
  const mean = v => sum(v) / v.length;
  function variance(v) { const m = mean(v); let s = 0; for (const x of v) s += (x - m) * (x - m); return s / (v.length - 1); }
  function checkProb(p, name) {
    if (typeof p !== 'number' || !Number.isFinite(p) || p <= 0 || p >= 1)
      throw new Error(`'${name}' must be a single number strictly between 0 and 1.`);
  }
  function quantile7(v, probs) {
    const x = v.filter(Number.isFinite).sort((a, b) => a - b);
    const n = x.length;
    return probs.map(p => {
      if (!n) return NaN;
      const h = (n - 1) * p, lo = Math.floor(h), hi = Math.min(lo + 1, n - 1);
      return x[lo] + (h - lo) * (x[hi] - x[lo]);
    });
  }
  const upd = (d, o) => Object.assign({}, d, o);

  // ------------------------------------------------------------------
  // Design (design.R)
  // ------------------------------------------------------------------
  function alignedPlanningWeights(pA, pB, respA, respB, r) {
    const pw = (pA + r * pB) / (1 + r);
    const seD = Math.sqrt(pw * (1 - pw) * (1 + 1 / r));
    const qbar = (pA * respA + r * pB * respB) / (pA + r * pB);
    const seE = Math.sqrt(qbar * (1 - qbar) * (1 / pA + 1 / (r * pB)));
    return [(respA + respB) / 2 * seD, (pA + pB) / 2 * seE];
  }

  function csDesign(o) {
    const pA = o.pA, pB = o.pB, r = o.r === undefined ? 1 : o.r;
    const method = o.method || 'wz', efficacy = o.efficacy || 'continuous';
    const contingent = !!o.contingent, calibration = o.calibration || 'recalibrated';
    const alpha = o.alpha === undefined ? 0.025 : o.alpha;
    const gamma = o.gamma === undefined ? 'compute' : o.gamma;
    const margin_d = o.margin_d || 0, margin_e = o.margin_e || 0;
    const gate_ties = o.gate_ties || 'stop';
    let weights = o.weights === undefined ? [0.5, 0.5] : o.weights;
    let c = o.c === undefined ? 0 : o.c;
    checkProb(pA, 'pA'); checkProb(pB, 'pB');
    if (efficacy === 'continuous') {
      if (o.delta == null || o.sigma == null || !Number.isFinite(o.delta) || !Number.isFinite(o.sigma))
        throw new Error("Supply 'delta' and 'sigma'.");
      if (o.sigma <= 0) throw new Error("'sigma' must be positive.");
    } else {
      if (o.respA == null || o.respB == null)
        throw new Error("Supply 'respA' and 'respB' for a binary efficacy outcome.");
      checkProb(o.respA, 'respA'); checkProb(o.respB, 'respB');
    }
    const aligned = weights === 'aligned';
    if (aligned) {
      if (efficacy !== 'binary')
        throw new Error('Estimand-aligned weights need a binary (response) efficacy outcome: with a continuous outcome there is no response probability to align to. Use the expected benefit per randomized patient (theta) for the direct test instead.');
      if (method !== 'wz') throw new Error('Estimand-aligned weights apply to the weighted Z.');
      if (contingent) throw new Error('Estimand-aligned weights are for the unconditional weighted Z; the contingent test uses fixed weights.');
      if (margin_d !== 0 || margin_e !== 0) throw new Error('Estimand-aligned weights are implemented for superiority only.');
      weights = alignedPlanningWeights(pA, pB, o.respA, o.respB, r);
    }
    if (!Array.isArray(weights) || weights.length !== 2 || weights.some(w => !(w >= 0)) || sum(weights) === 0)
      throw new Error("'weights' must be two non-negative numbers, not both zero, or \"aligned\".");
    const nrm = Math.sqrt(weights[0] ** 2 + weights[1] ** 2);
    weights = [weights[0] / nrm, weights[1] / nrm];
    if (method === 'fisher' && Math.abs(weights[0] - weights[1]) > 1e-8)
      throw new Error("Fisher's method is implemented with equal weights only; use the weighted Z.");
    if (method === 'fisher' && (margin_d !== 0 || margin_e !== 0))
      throw new Error('Non-inferiority margins are implemented for the weighted Z only.');
    if (!(gamma === 'compute' || gamma === 'zero' || (typeof gamma === 'number' && Number.isFinite(gamma))))
      throw new Error("'gamma' must be 'compute', 'zero' or a number.");
    if (margin_d < 0 || margin_e < 0) throw new Error('Margins must be non-negative.');
    if (!(r > 0)) throw new Error("'r' must be positive.");
    if (!(alpha > 0 && alpha < 0.5)) throw new Error("'alpha' must be in (0, 0.5).");
    if (!contingent) c = -Infinity;
    return {
      type: 'cs', pA, pB, delta: o.delta, sigma: o.sigma, r, weights, method, contingent,
      c, calibration, alpha, gamma, efficacy, respA: o.respA, respB: o.respB,
      margin_d, margin_e, gate_ties, aligned
    };
  }

  function csAlignedWeights(pA, pB, muA, muB, seD, seE, s = 0) {
    checkProb(pA, 'pA'); checkProb(pB, 'pB');
    if (!(seD > 0) || !(seE > 0)) throw new Error("'se_d' and 'se_e' must be positive.");
    const mBar = (muA + muB) / 2 - s, pBar = (pA + pB) / 2;
    const w = [mBar * seD, pBar * seE], nrm = Math.sqrt(w[0] ** 2 + w[1] ** 2);
    return { weights: [w[0] / nrm, w[1] / nrm], c_sq: pBar / mBar, m_bar: mBar, p_bar: pBar,
      warning: mBar <= 0 ? "The completers' mean outcome is not above the dropout score s, so keeping patients in the trial lowers theta; the completion weight is negative." : null };
  }

  function gateOpen(zd, c, ties = 'stop') {
    if (!Number.isFinite(c)) return true;
    return ties === 'pass' ? zd >= c - 1e-12 : zd > c + 1e-12;
  }

  // ------------------------------------------------------------------
  // Components (components.R)
  // ------------------------------------------------------------------
  function stdEfficacy(d) {
    if (d.efficacy === 'continuous') return (d.delta + d.margin_e) / d.sigma;
    const qbar = (d.pA * d.respA + d.r * d.pB * d.respB) / (d.pA + d.r * d.pB);
    return (d.respA - d.respB + d.margin_e) / Math.sqrt(qbar * (1 - qbar));
  }
  function gammaDelta(dStd, mA, mB, pA, pB, se0) {
    const S = 1 / mA + 1 / mB;
    return 0.5 * dStd * Math.pow(S, -1.5) * (pA * (1 - pA) / (mA * mA) - pB * (1 - pB) / (mB * mB)) / se0;
  }
  function csComponents(d, n) {
    const nA = n, nB = d.r * n, pA = d.pA, pB = d.pB;
    const pbar = (nA * pA + nB * pB) / (nA + nB);
    let se0, eta_d, v_d;
    if (d.margin_d > 0) {
      se0 = Math.sqrt(pA * (1 - pA) / nA + pB * (1 - pB) / nB);
      eta_d = (pA - pB + d.margin_d) / se0; v_d = 1;
    } else {
      se0 = Math.sqrt(pbar * (1 - pbar) * (1 / nA + 1 / nB));
      eta_d = (pA - pB) / se0;
      v_d = (pA * (1 - pA) / nA + pB * (1 - pB) / nB) / (se0 * se0);
    }
    const d_std = stdEfficacy(d);
    const mA = nA * pA, mB = nB * pB;
    const eta_e = d_std / Math.sqrt(1 / mA + 1 / mB);
    const g = d.gamma === 'compute' ? gammaDelta(d_std, mA, mB, pA, pB, se0) : d.gamma === 'zero' ? 0 : d.gamma;
    return { eta_d, v_d, eta_e, gamma: g, nA, nB, se0, d_std };
  }
  function csGamma(d) { return csComponents(upd(d, { gamma: 'compute' }), 100).gamma; }

  // ------------------------------------------------------------------
  // Critical values (critical.R)
  // ------------------------------------------------------------------
  function boundary(z, k, d) {
    if (d.method === 'wz') return (k - d.weights[0] * z) / d.weights[1];
    const pd = pnorm(z, 0, 1, false), q = k / pd;
    return q < 1 ? qnorm(q, 0, 1, false) : -Infinity;
  }
  function nullSize(d, k) {
    const c0 = d.c;
    if (d.method === 'wz' && d.weights[1] === 0) return pnorm(Math.max(k, c0), 0, 1, false);
    if (d.method === 'wz' && !Number.isFinite(c0)) return pnorm(k, 0, 1, false);
    const lo = Number.isFinite(c0) ? c0 : -12;
    return integratePieces(z => dnorm(z) * pnorm(boundary(z, k, d), 0, 1, false), lo, 12, kinks(d, k), 1e-10);
  }
  // For Fisher's method the boundary jumps to -Inf where P(d) = k, a kink in
  // the integrand; integrating on each side of it keeps full accuracy.
  function kinks(d, k) { return d.method === 'fisher' && k > 0 && k < 1 ? [qnorm(k, 0, 1, false)] : []; }
  function integratePieces(f, lo, hi, cuts, tol) {
    const pts = [lo].concat(cuts.filter(x => x > lo && x < hi).sort((a, b) => a - b), [hi]);
    let s = 0;
    for (let i = 0; i < pts.length - 1; i++) s += integrate(f, pts[i], pts[i + 1], tol);
    return s;
  }
  function csCritical(d) {
    const alpha = d.alpha;
    if (d.method === 'wz') {
      const k0 = qnorm(1 - alpha);
      let k;
      if (!d.contingent || d.calibration === 'usual' || d.weights[1] === 0) k = k0;
      else k = uniroot(x => nullSize(d, x) - alpha, k0 - 8, k0, 1e-10);
      return { k, scale: 'Z', size: nullSize(d, k) };
    }
    const t0 = Math.exp(-qchisq4Upper(alpha) / 2);
    let t;
    if (!d.contingent || d.calibration === 'usual') t = t0;
    else t = uniroot(x => nullSize(d, x) - alpha, t0, Math.min(0.5, 50 * t0), 1e-12);
    return { k: t, scale: 'product', size: nullSize(d, t) };
  }

  // ------------------------------------------------------------------
  // Power and sample size (power.R)
  // ------------------------------------------------------------------
  function powerOne(d, cmp, k) {
    const wd = d.weights[0], we = d.weights[1];
    const sv = Math.sqrt(cmp.v_d), g = cmp.gamma;
    if (d.method === 'wz' && we === 0) {
      const thr = Math.max(k, d.c);
      return pnorm((thr - cmp.eta_d) / sv, 0, 1, false);
    }
    if (d.method === 'wz' && !Number.isFinite(d.c)) {
      const mu = wd * cmp.eta_d + we * cmp.eta_e;
      const tau = Math.sqrt(wd * wd * cmp.v_d + we * we + 2 * wd * we * g);
      return pnorm((k - mu) / tau, 0, 1, false);
    }
    const rho2 = g * g / cmp.v_d;
    if (rho2 >= 1) throw new Error('Implausible gamma: |gamma| must be < sqrt(v(d)).');
    const sr = Math.sqrt(1 - rho2);
    const f = z => {
      const me = cmp.eta_e + g * (z - cmp.eta_d) / cmp.v_d;
      return dnorm(z, cmp.eta_d, sv) * pnorm((boundary(z, k, d) - me) / sr, 0, 1, false);
    };
    let lo = cmp.eta_d - 12 * sv;
    if (Number.isFinite(d.c)) lo = Math.max(lo, d.c);
    const hi = cmp.eta_d + 12 * sv;
    if (lo >= hi) return 0;
    return integratePieces(f, lo, hi, kinks(d, k), 1e-9);
  }
  function csPower(d, ns) {
    const crit = csCritical(d), za = qnorm(1 - d.alpha);
    const rows = (Array.isArray(ns) ? ns : [ns]).map(n => {
      const cmp = csComponents(d, n);
      return {
        nA: n, nB: d.r * n, power: powerOne(d, cmp, crit.k),
        power_completion: pnorm((za - cmp.eta_d) / Math.sqrt(cmp.v_d), 0, 1, false),
        power_completers: pnorm(za - cmp.eta_e, 0, 1, false)
      };
    });
    rows.critical = crit;
    return rows;
  }
  function searchN(f, target, nMax) {
    let lo = 2;
    if (f(lo) >= target) return lo;
    let hi = 4;
    while (f(hi) < target) {
      lo = hi; hi *= 2;
      if (hi > nMax) throw new Error('Target power not reached below n_max.');
    }
    while (hi - lo > 1) {
      const mid = Math.floor((lo + hi) / 2);
      if (f(mid) >= target) hi = mid; else lo = mid;
    }
    return hi;
  }
  function csSampleSize(d, power = 0.80, nMax = 1e5) {
    const crit = csCritical(d);
    const pw = n => powerOne(d, csComponents(d, n), crit.k);
    const n = searchN(pw, power, nMax);
    const nB = Math.ceil(d.r * n);
    return { nA: n, nB, total: n + nB, power: pw(n), target: power, critical: crit, design: d };
  }
  function csOptWeights(d, n) {
    const cmp = csComponents(d, n);
    let w = [Math.max(cmp.eta_d, 0), Math.max(cmp.eta_e, 0)];
    if (sum(w) === 0) throw new Error('Neither component favours A under these planning values.');
    w = [w[0] / sum(w), w[1] / sum(w)];
    const nrm = Math.sqrt(w[0] ** 2 + w[1] ** 2);
    const dOpt = upd(d, { weights: [w[0] / nrm, w[1] / nrm] });
    return { weights: w, power_design: csPower(d, n)[0].power, power_optimal: csPower(dOpt, n)[0].power };
  }

  // ------------------------------------------------------------------
  // psi (psi.R)
  // ------------------------------------------------------------------
  function psiFun(pA, pB, delta, sigma, weights, r = 1, type = 'information') {
    const pbar = (pA + r * pB) / (1 + r);
    const h = (pA - pB) / Math.sqrt(pbar * (1 - pbar));
    const dd = delta / sigma;
    const nrm = Math.sqrt(weights[0] ** 2 + weights[1] ** 2);
    const w = [weights[0] / nrm, weights[1] / nrm];
    const pw = (1 + 1 / r) / (1 / pA + 1 / (r * pB));
    const f = type === 'information' ? Math.sqrt(pw) : 1;
    return { psi: w[0] * h + w[1] * dd * f, h, d: dd, pbar, pi_w: pw };
  }
  function csPsi(d, type = 'information') {
    if (d.efficacy !== 'continuous') throw new Error('psi is defined for a continuous efficacy outcome.');
    return psiFun(d.pA, d.pB, d.delta, d.sigma, d.weights, d.r, type);
  }

  // ------------------------------------------------------------------
  // theta (theta.R)
  // ------------------------------------------------------------------
  function thetaDesign(o) {
    const r = o.r === undefined ? 1 : o.r, alpha = o.alpha === undefined ? 0.025 : o.alpha;
    checkProb(o.pA, 'pA'); checkProb(o.pB, 'pB');
    if (!(o.sigma > 0)) throw new Error("'sigma' must be positive.");
    const margin = o.margin || 0, s = o.s || 0;
    if (margin < 0) throw new Error("'margin' must be non-negative.");
    return { type: 'theta', pA: o.pA, pB: o.pB, muA: o.muA, muB: o.muB, sigma: o.sigma, s, r, alpha, margin };
  }
  const thetaValue = d => (d.pA * d.muA + (1 - d.pA) * d.s) - (d.pB * d.muB + (1 - d.pB) * d.s);
  function thetaDecompose(d) {
    const comp = (d.pA - d.pB) * ((d.muA + d.muB) / 2 - d.s);
    const eff = ((d.pA + d.pB) / 2) * (d.muA - d.muB);
    return { completion: comp, completers: eff, theta: comp + eff };
  }
  const thetaV = (p, mu, sigma, s) => p * sigma * sigma + p * (1 - p) * (mu - s) * (mu - s);
  function thetaPower(d, n) {
    const VA = thetaV(d.pA, d.muA, d.sigma, d.s), VB = thetaV(d.pB, d.muB, d.sigma, d.s);
    const se = Math.sqrt(VA / n + VB / (d.r * n));
    return pnorm(qnorm(1 - d.alpha) - (thetaValue(d) + d.margin) / se, 0, 1, false);
  }
  function thetaSampleSize(d, power = 0.80) {
    const th = thetaValue(d) + d.margin;
    if (th <= 0) throw new Error('theta (+ margin) must be positive under the alternative.');
    const VA = thetaV(d.pA, d.muA, d.sigma, d.s), VB = thetaV(d.pB, d.muB, d.sigma, d.s);
    const nA = Math.ceil((qnorm(1 - d.alpha) + qnorm(power)) ** 2 * (VA + VB / d.r) / (th * th));
    const nB = Math.ceil(d.r * nA);
    return { nA, nB, total: nA + nB, theta: thetaValue(d), VA, VB, power: thetaPower(d, nA) };
  }
  function completionSampleSize(pA, pB, r = 1, alpha = 0.025, power = 0.80) {
    const pbar = (pA + r * pB) / (1 + r);
    const num = qnorm(1 - alpha) * Math.sqrt(pbar * (1 - pbar) * (1 + 1 / r)) +
      qnorm(power) * Math.sqrt(pA * (1 - pA) + pB * (1 - pB) / r);
    const nA = Math.ceil(num * num / ((pA - pB) * (pA - pB)));
    const nB = Math.ceil(r * nA);
    return { nA, nB, total: nA + nB };
  }

  // ------------------------------------------------------------------
  // MMRM comparator (mmrm_planning.R)
  // ------------------------------------------------------------------
  const NEWMEDS_CORR = [
    [1.00, 0.72, 0.59, 0.54, 0.49, 0.51],
    [0.72, 1.00, 0.79, 0.69, 0.66, 0.64],
    [0.59, 0.79, 1.00, 0.83, 0.76, 0.73],
    [0.54, 0.69, 0.83, 1.00, 0.86, 0.83],
    [0.49, 0.66, 0.76, 0.86, 1.00, 0.90],
    [0.51, 0.64, 0.73, 0.83, 0.90, 1.00]];
  function mmrmCorr(visits, structure, rho) {
    if (structure === 'newmeds') {
      if (visits !== 6) throw new Error('The NEWMEDS correlation matrix is for 6 weekly visits; use an AR(1) or compound-symmetry structure for another number of visits.');
      return NEWMEDS_CORR;
    }
    if (rho == null || !(rho >= 0) || !(rho < 1)) throw new Error("'rho' must be in [0, 1).");
    const R = [];
    for (let i = 0; i < visits; i++) {
      R.push([]);
      for (let j = 0; j < visits; j++)
        R[i].push(i === j ? 1 : structure === 'ar1' ? Math.pow(rho, Math.abs(i - j)) : rho);
    }
    return R;
  }
  function mmrmRetention(p, visits, pattern) {
    const out = [];
    for (let j = 1; j <= visits; j++) {
      let f = j / visits;
      if (pattern === 'early') f = 1 - (1 - f) ** 2; else if (pattern === 'late') f = f * f;
      out.push(1 - (1 - p) * f);
    }
    return out;
  }
  function solveSym(A, b) {
    const n = b.length, M = A.map((row, i) => row.slice(0, n).concat([b[i]]));
    for (let c = 0; c < n; c++) {
      let p = c;
      for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
      const t = M[c]; M[c] = M[p]; M[p] = t;
      if (Math.abs(M[c][c]) < 1e-14) throw new Error('The correlation matrix is singular.');
      for (let r = c + 1; r < n; r++) {
        const f = M[r][c] / M[c][c];
        for (let k = c; k <= n; k++) M[r][k] -= f * M[c][k];
      }
    }
    const x = new Array(n).fill(0);
    for (let r = n - 1; r >= 0; r--) {
      let s = M[r][n];
      for (let k = r + 1; k < n; k++) s -= M[r][k] * x[k];
      x[r] = s / M[r][r];
    }
    return x;
  }
  function mmrmPhi(R, ret) {
    const K = R.length;
    let prev = 0, phi = 0;
    const taus = [];
    for (let j = 1; j <= K; j++) {
      const s = R[K - 1].slice(0, j);
      const x = solveSym(R.slice(0, j).map(row => row.slice(0, j)), s);
      const V = s.reduce((a, si, i) => a + si * x[i], 0);
      taus.push(V - prev);
      phi += (V - prev) / ret[j - 1];
      prev = V;
    }
    return phi;
  }
  function csMmrmDesign(o) {
    const visits = o.visits === undefined ? 6 : o.visits;
    const structure = o.structure || 'newmeds', pattern = o.pattern || 'even';
    const r = o.r === undefined ? 1 : o.r, alpha = o.alpha === undefined ? 0.025 : o.alpha;
    checkProb(o.pA, 'pA'); checkProb(o.pB, 'pB');
    if (!Number.isFinite(o.d) || o.d === 0) throw new Error("'d' must be a non-zero number.");
    if (!(visits >= 1) || visits !== Math.round(visits)) throw new Error("'visits' must be a positive integer.");
    if (!(r > 0)) throw new Error("'r' must be positive.");
    if (!(alpha > 0 && alpha < 0.5)) throw new Error("'alpha' must be in (0, 0.5).");
    const R = mmrmCorr(visits, structure, o.rho);
    const phiA = mmrmPhi(R, mmrmRetention(o.pA, visits, pattern));
    const phiB = mmrmPhi(R, mmrmRetention(o.pB, visits, pattern));
    return { type: 'mmrm', d: o.d, pA: o.pA, pB: o.pB, visits, structure, rho: o.rho, pattern, r, alpha, phiA, phiB, R };
  }
  function csMmrmPower(d, nA) {
    return pnorm(Math.abs(d.d) * Math.sqrt(nA) / Math.sqrt(d.phiA + d.phiB / d.r) - qnorm(1 - d.alpha));
  }
  function csMmrmSampleSize(d, power = 0.90) {
    const nA = Math.ceil((qnorm(1 - d.alpha) + qnorm(power)) ** 2 * (d.phiA + d.phiB / d.r) / (d.d * d.d));
    const nB = Math.ceil(d.r * nA);
    return { nA, nB, total: nA + nB, power: csMmrmPower(d, nA), phiA: d.phiA, phiB: d.phiB };
  }

  // ------------------------------------------------------------------
  // Planning tools (planning.R)
  // ------------------------------------------------------------------
  function csPowerGrid(d, n, pAs, pBs, deltas) {
    pAs = pAs || [d.pA]; pBs = pBs || [d.pB]; deltas = deltas || [d.delta];
    const crit = csCritical(d), out = [];
    for (const delta of deltas) for (const pB of pBs) for (const pA of pAs) {
      const d2 = upd(d, { pA, pB });
      if (delta != null && Number.isFinite(delta)) d2.delta = delta;
      out.push({ pA, pB, delta, power: powerOne(d2, csComponents(d2, n), crit.k) });
    }
    return out;
  }
  function anyPower(d, n) {
    return d.type === 'theta' ? thetaPower(d, n) : csPower(d, n)[0].power;
  }
  function csAssurance(d, n, o = {}) {
    const nsim = o.nsim || 2000, rng = makeRng(o.seed === undefined ? 1 : o.seed);
    const draw = (pr, fixed) => {
      if (!pr) return new Array(nsim).fill(fixed);
      if (pr.length !== 2 || pr.some(v => !(v > 0))) throw new Error('Beta prior must be c(a, b) with a, b > 0.');
      return Array.from({ length: nsim }, () => rng.beta(pr[0], pr[1]));
    };
    const pA = draw(o.prior_pA, d.pA), pB = draw(o.prior_pB, d.pB);
    let eff, pw, p0;
    if (d.type === 'theta') {
      const eff0 = d.muA - d.muB;
      eff = o.prior_delta ? Array.from({ length: nsim }, () => rng.norm(o.prior_delta[0], o.prior_delta[1])) : new Array(nsim).fill(eff0);
      pw = eff.map((e, i) => thetaPower(upd(d, { pA: pA[i], pB: pB[i], muA: d.muB + e }), n));
      p0 = thetaPower(d, n);
    } else {
      eff = o.prior_delta ? Array.from({ length: nsim }, () => rng.norm(o.prior_delta[0], o.prior_delta[1])) : new Array(nsim).fill(d.delta);
      const crit = csCritical(d);
      pw = eff.map((e, i) => {
        const d2 = upd(d, { pA: pA[i], pB: pB[i], delta: e });
        return powerOne(d2, csComponents(d2, n), crit.k);
      });
      p0 = csPower(d, n)[0].power;
    }
    const a = mean(pw);
    return { assurance: a, mc_se: Math.sqrt(variance(pw) / nsim), power_at_planning_values: p0,
      draws: pA.map((x, i) => ({ pA: x, pB: pB[i], effect: eff[i], power: pw[i] })) };
  }
  function csBetaPrior(m, sd) {
    checkProb(m, 'mean');
    const k = m * (1 - m) / (sd * sd) - 1;
    if (!(sd > 0) || !(k > 0)) throw new Error("'sd' must be positive and smaller than sqrt(mean * (1 - mean)).");
    return [m * k, (1 - m) * k];
  }

  function csEstimandStatement(d, o = {}) {
    const population = o.population || "patients as defined by the trial's eligibility criteria";
    const treatment = o.treatment || 'the active treatment (A) versus the comparator (B), each as randomized';
    const endpoint = o.endpoint || 'improvement in the primary efficacy score at the final visit';
    const ice = 'Treatment discontinuation (including switching treatment or adding rescue medication not defined as part of the regimen) is handled by the composite strategy: it is part of the outcome, and no outcomes after discontinuation are modelled.';
    let variable, summary;
    if (d.type === 'theta') {
      variable = `Benefit per randomized patient: ${endpoint} for patients who complete the trial on assigned treatment; a score of ${fmtG(d.s)} for patients who discontinue.`;
      summary = `Difference in mean benefit per randomized patient (theta, in clinical units)${d.margin > 0 ? `; non-inferiority margin ${fmtG(d.margin)}` : ''}.`;
    } else {
      const ws = d.weights[0] + d.weights[1];
      const w = [d.weights[0] / ws, d.weights[1] / ws];
      variable = `Two components: completion of the trial on assigned treatment, and ${d.efficacy === 'continuous' ? endpoint : 'response among completers'} among completers.`;
      summary = `Composite of the difference in completion and the difference among completers, combined by ${d.method === 'wz' ? 'the weighted Z' : "Fisher's method"} with prespecified relative weights ${w[0].toFixed(2)} (completion) and ${w[1].toFixed(2)} (efficacy)${d.contingent ? `; efficacy counts only if A is not worse on completion (gate Z(d) >= ${fmtG(d.c)}, ${d.calibration} critical value)` : ''}; summarized by the standardized composite effect size psi with its confidence interval, reported with both components.`;
    }
    return `Population: ${population}.\nTreatment: ${treatment}.\nVariable: ${variable}\nIntercurrent events: ${ice}\nPopulation-level summary: ${summary}`;
  }
  function fmtG(x) { return String(+x.toPrecision(6)); }

  function csCompareDesigns(o) {
    const { pA, pB, delta, sigma, muB } = o;
    const s = o.s || 0, alpha = o.alpha === undefined ? 0.025 : o.alpha;
    const power = o.power === undefined ? 0.80 : o.power, gamma = o.gamma || 'zero';
    const c = o.c === undefined ? 0 : o.c;
    let opts = [
      ['Weighted Z, equal weights (0.5 / 0.5)', { w: [0.5, 0.5] }],
      ['Weighted Z, 0.7 / 0.3', { w: [0.7, 0.3] }],
      ['Weighted Z, 0.3 / 0.7', { w: [0.3, 0.7] }],
      ['Weighted Z, 0.3 / 0.7, contingent, recalibrated', { w: [0.3, 0.7], cont: true, cal: 'recalibrated' }],
      ['Weighted Z, 0.3 / 0.7, contingent, usual critical value', { w: [0.3, 0.7], cont: true, cal: 'usual' }],
      ['Fisher (Shih-Quan), equal weights', { w: [0.5, 0.5], m: 'fisher' }],
      ['Fisher, equal weights, contingent, recalibrated', { w: [0.5, 0.5], m: 'fisher', cont: true, cal: 'recalibrated' }],
      ['Completion rates only', { special: 'completion' }],
      ['Completers analysis only', { w: [0, 1] }],
      ['Expected benefit per randomized patient', { special: 'theta' }]];
    if (o.equal_contingent) {
      opts.splice(5, 0,
        ['Weighted Z, equal weights, contingent, recalibrated', { w: [0.5, 0.5], cont: true, cal: 'recalibrated' }],
        ['Weighted Z, equal weights, contingent, usual critical value', { w: [0.5, 0.5], cont: true, cal: 'usual' }]);
    }
    return opts.map(([name, op]) => {
      try {
        if (op.special === 'completion') {
          const r = completionSampleSize(pA, pB, 1, alpha, power);
          return { design: name, n_per_arm: r.nA, power: NaN };
        }
        if (op.special === 'theta') {
          const r = thetaSampleSize(thetaDesign({ pA, pB, muA: muB + delta, muB, sigma, s, alpha }), power);
          return { design: name, n_per_arm: r.nA, power: r.power };
        }
        const d = csDesign({ pA, pB, delta, sigma, weights: op.w, method: op.m || 'wz', contingent: !!op.cont,
          c, calibration: op.cal || 'recalibrated', alpha, gamma });
        const r = csSampleSize(d, power);
        return { design: name, n_per_arm: r.nA, power: r.power };
      } catch (e) {
        return { design: name, n_per_arm: NaN, power: NaN, error: e.message };
      }
    });
  }
  function csTable6(gamma = 'zero', equalContingent = false) {
    const p = csCompareDesigns({ pA: 0.65, pB: 0.50, delta: 6, sigma: 18, muB: 15, gamma, equal_contingent: equalContingent });
    const a = csCompareDesigns({ pA: 0.72, pB: 0.64, delta: 3, sigma: 18, muB: 18, gamma, equal_contingent: equalContingent });
    return p.map((row, i) => ({ design: row.design, placebo_controlled: row.n_per_arm, active_comparator: a[i].n_per_arm }));
  }

  // ------------------------------------------------------------------
  // Linear models for the analysis (analyze.R)
  // ------------------------------------------------------------------
  // Ordinary least squares with model-based and HC2 variances.
  function ols(X, y) {
    const n = X.length, p = X[0].length;
    const XtX = Array.from({ length: p }, () => new Array(p).fill(0));
    const Xty = new Array(p).fill(0);
    for (let i = 0; i < n; i++) {
      const xi = X[i];
      for (let a = 0; a < p; a++) {
        Xty[a] += xi[a] * y[i];
        for (let b = a; b < p; b++) XtX[a][b] += xi[a] * xi[b];
      }
    }
    for (let a = 0; a < p; a++) for (let b = 0; b < a; b++) XtX[a][b] = XtX[b][a];
    const XtXi = invertSPD(XtX);
    const beta = XtXi.map(row => row.reduce((s, v, j) => s + v * Xty[j], 0));
    const resid = y.map((yi, i) => yi - X[i].reduce((s, v, j) => s + v * beta[j], 0));
    const dfr = n - p;
    const rss = resid.reduce((s, e) => s + e * e, 0);
    const sigma = Math.sqrt(rss / dfr);
    const hat = X.map(xi => {
      let h = 0;
      for (let a = 0; a < p; a++) { let t = 0; for (let b = 0; b < p; b++) t += XtXi[a][b] * xi[b]; h += xi[a] * t; }
      return h;
    });
    return { beta, XtXi, resid, df: dfr, sigma, hat, X };
  }
  function invertSPD(A) {
    const n = A.length, M = A.map((r, i) => r.concat(Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))));
    for (let c = 0; c < n; c++) {
      let p = c;
      for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
      const t = M[c]; M[c] = M[p]; M[p] = t;
      const piv = M[c][c];
      const scale = Math.max(1, Math.abs(A[c][c]));
      if (Math.abs(piv) < 1e-10 * scale) throw new Error('The model matrix is not of full rank (a covariate or stratum is collinear with treatment, or a level has no patients).');
      for (let k = 0; k < 2 * n; k++) M[c][k] /= piv;
      for (let r = 0; r < n; r++) if (r !== c) {
        const f = M[r][c];
        if (f !== 0) for (let k = 0; k < 2 * n; k++) M[r][k] -= f * M[c][k];
      }
    }
    return M.map(r => r.slice(n));
  }
  function hc2Se(fit, j) {
    const p = fit.beta.length, meat = Array.from({ length: p }, () => new Array(p).fill(0));
    fit.X.forEach((xi, i) => {
      const u = fit.resid[i] / Math.sqrt(1 - fit.hat[i]);
      for (let a = 0; a < p; a++) for (let b = 0; b < p; b++) meat[a][b] += u * u * xi[a] * xi[b];
    });
    const B = fit.XtXi;
    let v = 0;
    for (let a = 0; a < p; a++) for (let b = 0; b < p; b++) v += B[j][a] * meat[a][b] * B[b][j];
    return Math.sqrt(v);
  }
  // rows: [{A, comp, y, cov: [..], strata}]; builds y ~ A + covariates (+ strata)
  function designMatrix(rows, nCov, useStrata) {
    let levels = [];
    if (useStrata) {
      levels = [...new Set(rows.map(r => r.strata))].sort(compareLevels);
    }
    const X = rows.map(r => {
      const x = [1, r.A].concat(r.cov.slice(0, nCov));
      for (let l = 1; l < levels.length; l++) x.push(r.strata === levels[l] ? 1 : 0);
      return x;
    });
    return { X, levels };
  }
  function compareLevels(a, b) {
    const na = Number(a), nb = Number(b);
    if (a !== '' && b !== '' && Number.isFinite(na) && Number.isFinite(nb)) return na - nb;
    return a < b ? -1 : a > b ? 1 : 0;
  }
  function cmhZ(rows, levels) {
    let num = 0, v = 0;
    for (const lv of levels) {
      const d = rows.filter(r => r.strata === lv);
      const N = d.length;
      if (N < 2) continue;
      const n1 = d.filter(r => r.A === 1).length, n0 = N - n1;
      const x = d.filter(r => r.comp).length, x1 = d.filter(r => r.comp && r.A === 1).length;
      num += x1 - n1 * x / N;
      v += n1 * n0 * x * (N - x) / (N * N * (N - 1));
    }
    return num / Math.sqrt(v);
  }
  function csEstimates(df, nCov, stratified, s, tests = true) {
    const A1 = df.filter(r => r.A === 1), A0 = df.filter(r => r.A === 0);
    const nA = A1.length, nB = A0.length;
    const pA = mean(A1.map(r => (r.comp ? 1 : 0))), pB = mean(A0.map(r => (r.comp ? 1 : 0)));
    const covOk = r => r.cov.slice(0, nCov).every(Number.isFinite);
    // completers' efficacy
    const cdf = df.filter(r => r.comp && Number.isFinite(r.y) && covOk(r));
    const cLev = stratified ? [...new Set(cdf.map(r => r.strata))] : [];
    const dm = designMatrix(cdf, nCov, stratified && cLev.length >= 2);
    const fit = ols(dm.X, cdf.map(r => r.y));
    const Delta = fit.beta[1], seDelta = Math.sqrt(fit.XtXi[1][1]) * fit.sigma;
    // theta: dropouts scored s
    const tdf = df.filter(r => !(r.comp && !Number.isFinite(r.y)) && covOk(r));
    const tLev = stratified ? [...new Set(tdf.map(r => r.strata))] : [];
    const dt = designMatrix(tdf, nCov, stratified && tLev.length >= 2);
    const tfit = ols(dt.X, tdf.map(r => (r.comp ? r.y : s)));
    const out = {
      nA, nB, r: nB / nA, pA, pB,
      muA: mean(cdf.filter(r => r.A === 1).map(r => r.y)), muB: mean(cdf.filter(r => r.A === 0).map(r => r.y)),
      Delta, seDelta, sigma: fit.sigma, df: fit.df, theta: tfit.beta[1], se_theta: hc2Se(tfit, 1)
    };
    if (tests) {
      if (stratified) {
        const lv = [...new Set(df.map(r => r.strata))].sort(compareLevels);
        out.Zd = cmhZ(df, lv);
      } else {
        const pbar = mean(df.map(r => (r.comp ? 1 : 0)));
        out.Zd = (pA - pB) / Math.sqrt(pbar * (1 - pbar) * (1 / nA + 1 / nB));
      }
      out.Ze = tToZ(Delta / seDelta, fit.df);
    }
    return out;
  }
  function psiSe(est, w, type) {
    const th = [est.pA, est.pB, est.Delta, est.sigma];
    const V = [est.pA * (1 - est.pA) / est.nA, est.pB * (1 - est.pB) / est.nB,
      est.seDelta * est.seDelta, est.sigma * est.sigma / (2 * est.df)];
    const f = x => psiFun(x[0], x[1], x[2], x[3], w, est.r, type).psi;
    let s = 0;
    for (let j = 0; j < 4; j++) {
      const h = 1e-6 * Math.max(1, Math.abs(th[j]));
      const up = th.slice(), dn = th.slice();
      up[j] += h; dn[j] -= h;
      const g = (f(up) - f(dn)) / (2 * h);
      s += g * g * V[j];
    }
    return Math.sqrt(s);
  }
  function parseLogical(v) {
    if (v === true || v === false) return v;
    if (v == null) return null;
    const t = String(v).trim();
    if (t === '' || t.toUpperCase() === 'NA') return null;
    if (['TRUE', 'true', 'True', 'T'].includes(t)) return true;
    if (['FALSE', 'false', 'False', 'F'].includes(t)) return false;
    const x = Number(t);
    if (Number.isFinite(x)) return x !== 0;
    return null;
  }
  function parseNum(v) {
    if (v == null) return NaN;
    const t = String(v).trim();
    if (t === '' || t.toUpperCase() === 'NA') return NaN;
    const x = Number(t);
    return Number.isFinite(x) ? x : NaN;
  }
  // data: array of row objects (values as read from CSV)
  function csAnalyze(data, o) {
    const method = o.method || 'wz', calibration = o.calibration || 'recalibrated';
    const psiType = o.psi_type || 'information', alpha = o.alpha === undefined ? 0.025 : o.alpha;
    const conf = o.conf_level || 0.95, s = o.s || 0, contingent = !!o.contingent;
    const c = o.c === undefined ? 0 : o.c, gateTies = o.gate_ties || 'stop';
    const covariates = o.covariates || [], strata = o.strata || null;
    const cols = data.length ? Object.keys(data[0]) : [];
    const miss = [o.arm, o.completed, o.outcome, ...covariates, ...(strata ? [strata] : [])].filter(v => !cols.includes(v));
    if (miss.length) throw new Error('Columns not found: ' + miss.join(', '));
    const df = data.map(row => ({
      A: String(row[o.arm]).trim() === String(o.active).trim() ? 1 : 0,
      comp: parseLogical(row[o.completed]),
      y: parseNum(row[o.outcome]),
      cov: covariates.map(v => parseNum(row[v])),
      strata: strata ? String(row[strata]).trim() : null
    }));
    if (df.some(r => r.comp === null)) throw new Error('Completion status must be known for every patient.');
    if (new Set(df.map(r => r.A)).size !== 2) throw new Error("'arm' must have exactly two groups.");
    const nMiss = df.filter(r => r.comp && !Number.isFinite(r.y)).length;
    const aligned = o.weights === 'aligned';
    if (aligned) {
      if (method !== 'wz') throw new Error('Estimand-aligned weights apply to the weighted Z.');
      if (contingent) throw new Error('Estimand-aligned weights are for the unconditional weighted Z; the contingent test uses fixed weights.');
    }
    const des = csDesign({ pA: 0.5, pB: 0.5, delta: 0, sigma: 1, weights: aligned ? [0.5, 0.5] : (o.weights || [0.5, 0.5]),
      method, contingent, c, calibration, alpha });
    const crit = csCritical(des);
    const nCov = covariates.length, stratified = !!strata;
    const est = csEstimates(df, nCov, stratified, s);
    const Zd = est.Zd, Ze = est.Ze;
    const pd = pnorm(Zd, 0, 1, false), pe = pnorm(Ze, 0, 1, false);
    let w = des.weights, alignedInfo = null;
    if (aligned) {
      const pbar = mean(df.map(r => (r.comp ? 1 : 0)));
      const seD = Math.sqrt(pbar * (1 - pbar) * (1 / est.nA + 1 / est.nB));
      alignedInfo = csAlignedWeights(est.pA, est.pB, est.muA, est.muB, seD, est.seDelta, s);
      w = alignedInfo.weights;
    }
    const open = !contingent || gateOpen(Zd, c, gateTies);
    let stat, pUncond, exceeds;
    if (method === 'wz') {
      stat = w[0] * Zd + w[1] * Ze;
      pUncond = pnorm(stat, 0, 1, false);
      exceeds = stat > crit.k;
    } else {
      stat = pd * pe;
      pUncond = stat * (1 - Math.log(stat));
      exceeds = stat < crit.k;
    }
    const pValue = !contingent ? pUncond : !open ? 1 : nullSize(des, stat);
    const reject = open && exceeds;
    const zq = qnorm(1 - (1 - conf) / 2);
    const psi = psiFun(est.pA, est.pB, est.Delta, est.sigma, w, est.r, psiType).psi;
    const pse = psiSe(est, w, psiType);
    const out = {
      n: { A: est.nA, B: est.nB }, completion: { A: est.pA, B: est.pB },
      completers_effect: { estimate: est.Delta, se: est.seDelta, sigma: est.sigma, df: est.df },
      components: { Z_d: Zd, Z_e: Ze, P_d: pd, P_e: pe },
      test: { method, weights: w, aligned, statistic: stat, critical: crit.k, scale: crit.scale,
        contingent, c, gate_ties: gateTies, gate_open: open, calibration, p_unconditional: pUncond,
        p_value: pValue, reject, alpha },
      psi: { estimate: psi, se: pse, lower: psi - zq * pse, upper: psi + zq * pse }, psi_type: psiType,
      theta: { estimate: est.theta, se: est.se_theta, lower: est.theta - zq * est.se_theta,
        upper: est.theta + zq * est.se_theta, p_one_sided: pnorm(est.theta / est.se_theta, 0, 1, false) },
      s, conf_level: conf, n_completers_missing_outcome: nMiss, aligned_info: alignedInfo
    };
    const B = o.bootstrap || 0;
    if (B > 0) {
      const rng = makeRng(o.seed === undefined ? 1 : o.seed);
      const iA = [], iB = [];
      df.forEach((r, i) => (r.A === 1 ? iA : iB).push(i));
      const bs = [];
      for (let b = 0; b < B; b++) {
        const rows = [];
        for (let i = 0; i < iA.length; i++) rows.push(df[iA[rng.int(iA.length)]]);
        for (let i = 0; i < iB.length; i++) rows.push(df[iB[rng.int(iB.length)]]);
        try {
          const e = csEstimates(rows, nCov, stratified, s, false);
          bs.push([psiFun(e.pA, e.pB, e.Delta, e.sigma, w, e.r, psiType).psi, e.theta]);
        } catch (err) { bs.push([NaN, NaN]); }
      }
      const pr = [(1 - conf) / 2, 1 - (1 - conf) / 2];
      out.bootstrap = { B, psi: quantile7(bs.map(x => x[0]), pr), theta: quantile7(bs.map(x => x[1]), pr) };
    }
    return out;
  }

  // ------------------------------------------------------------------
  // Simulation with dropout unrelated to outcome (simulate.R)
  // ------------------------------------------------------------------
  function csSimulation(d, n, o = {}) {
    if (d.efficacy !== 'continuous' || d.margin_d > 0 || d.margin_e > 0)
      throw new Error('Simulation currently supports continuous efficacy, superiority form.');
    const nsim = o.nsim || 1000, isNull = !!o.null, s = o.s || 0, conf = o.conf_level || 0.95;
    const muB = o.muB === undefined ? null : o.muB;
    const rng = makeRng(o.seed === undefined ? 2026 : o.seed);
    const pA = d.pA, pB = isNull ? d.pA : d.pB, delta = isNull ? 0 : d.delta, sigma = d.sigma;
    const mB0 = muB == null ? 0 : muB, muA0 = mB0 + delta;
    const nA = n, nB = Math.ceil(d.r * n);
    const crit = csCritical(d), w = d.weights, zq = qnorm(1 - (1 - conf) / 2), za = qnorm(1 - d.alpha);
    const psiTrue = psiFun(pA, pB, delta, sigma, w, d.r).psi;
    const thetaTrue = (pA * muA0 + (1 - pA) * s) - (pB * mB0 + (1 - pB) * s);
    const res = [];
    let done = 0;
    function arm(m, p, mu) {
      const ys = [];
      for (let i = 0; i < m; i++) if (rng.unif() < p) ys.push(rng.norm(mu, sigma));
      return ys;
    }
    function one() {
      const yA = arm(nA, pA, muA0), yB = arm(nB, pB, mB0);
      const mA = yA.length, mB = yB.length;
      if (mA < 2 || mB < 2) return;
      const pAh = mA / nA, pBh = mB / nB, pbar = (mA + mB) / (nA + nB);
      const Zd = (pAh - pBh) / Math.sqrt(pbar * (1 - pbar) * (1 / nA + 1 / nB));
      const Dh = mean(yA) - mean(yB), dfe = mA + mB - 2;
      const sp = Math.sqrt(((mA - 1) * variance(yA) + (mB - 1) * variance(yB)) / dfe);
      const seD = sp * Math.sqrt(1 / mA + 1 / mB);
      const Ze = tToZ(Dh / seD, dfe);
      const gate = !d.contingent || gateOpen(Zd, d.c, d.gate_ties || 'stop');
      const ex = d.method === 'wz' ? (w[0] * Zd + w[1] * Ze) > crit.k
        : pnorm(Zd, 0, 1, false) * pnorm(Ze, 0, 1, false) < crit.k;
      const est = { pA: pAh, pB: pBh, Delta: Dh, seDelta: seD, sigma: sp, df: dfe, nA, nB, r: nB / nA };
      const ps = psiFun(pAh, pBh, Dh, sp, w, nB / nA).psi;
      const pse = psiSe(est, w, 'information');
      let thCov = NaN, thRej = NaN;
      if (muB != null) {
        const vA = yA.concat(new Array(nA - mA).fill(s)), vB = yB.concat(new Array(nB - mB).fill(s));
        const th = mean(vA) - mean(vB);
        const sth = Math.sqrt(variance(vA) / nA + variance(vB) / nB);
        thCov = Math.abs(th - thetaTrue) <= zq * sth ? 1 : 0;
        thRej = th / sth > za ? 1 : 0;
      }
      res.push([Zd, Ze, gate && ex ? 1 : 0, Zd > za ? 1 : 0, Ze > za ? 1 : 0,
        Math.abs(ps - psiTrue) <= zq * pse ? 1 : 0, thCov, thRej]);
    }
    return {
      total: nsim,
      get done() { return done; },
      step(k) { const end = Math.min(nsim, done + k); while (done < end) { one(); done++; } return done >= nsim; },
      result() {
        const col = j => res.map(r => r[j]);
        const rj = mean(col(2));
        const zd = col(0), ze = col(1), mzd = mean(zd), mze = mean(ze);
        let sxy = 0, sxx = 0, syy = 0;
        for (let i = 0; i < res.length; i++) { const a = zd[i] - mzd, b = ze[i] - mze; sxy += a * b; sxx += a * a; syy += b * b; }
        return {
          nsim: res.length, n: { A: nA, B: nB }, null: isNull, dropout: 'random',
          reject: [rj, Math.sqrt(rj * (1 - rj) / res.length)],
          reject_completion: mean(col(3)), reject_completers: mean(col(4)),
          cor_Zd_Ze: sxy / Math.sqrt(sxx * syy), psi_true: psiTrue, psi_coverage: mean(col(5)),
          theta_true: muB != null ? thetaTrue : NaN,
          theta_coverage: muB != null ? mean(col(6)) : NaN,
          theta_reject: muB != null ? mean(col(7)) : NaN,
          calculated_power: isNull ? NaN : csPower(d, n)[0].power
        };
      }
    };
  }

  // ------------------------------------------------------------------
  // CSV parsing (RFC 4180, comma or semicolon separated)
  // ------------------------------------------------------------------
  function parseCSV(text) {
    text = text.replace(/^﻿/, '');
    const first = text.split(/\r?\n/, 1)[0] || '';
    const sep = (first.match(/;/g) || []).length > (first.match(/,/g) || []).length ? ';' : ',';
    const rows = [];
    let row = [], field = '', inQ = false;
    for (let i = 0; i < text.length; i++) {
      const ch = text[i];
      if (inQ) {
        if (ch === '"') { if (text[i + 1] === '"') { field += '"'; i++; } else inQ = false; }
        else field += ch;
      } else if (ch === '"') inQ = true;
      else if (ch === sep) { row.push(field); field = ''; }
      else if (ch === '\n' || ch === '\r') {
        if (ch === '\r' && text[i + 1] === '\n') i++;
        row.push(field); field = '';
        if (row.length > 1 || row[0] !== '') rows.push(row);
        row = [];
      } else field += ch;
    }
    if (field !== '' || row.length) { row.push(field); if (row.length > 1 || row[0] !== '') rows.push(row); }
    if (!rows.length) throw new Error('The file is empty.');
    const header = rows[0].map(h => h.trim());
    return rows.slice(1).map(r => Object.fromEntries(header.map((h, j) => [h, r[j] === undefined ? '' : r[j]])));
  }

  const api = {
    version: '0.1.0.9000',
    pnorm, dnorm, qnorm, pt, pbeta, lgamma, tToZ, integrate, uniroot, qchisq4Upper, makeRng, quantile7,
    csDesign, csAlignedWeights, csComponents, csGamma, csCritical, csPower, csSampleSize, csOptWeights,
    csPsi, psiFun, thetaDesign, thetaValue, thetaDecompose, thetaPower, thetaSampleSize, completionSampleSize,
    csMmrmDesign, csMmrmPower, csMmrmSampleSize, csPowerGrid, anyPower, csAssurance, csBetaPrior,
    csEstimandStatement, csCompareDesigns, csTable6, csAnalyze, csSimulation, parseCSV, nullSize, gateOpen
  };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.CS = api;
})(typeof self !== 'undefined' ? self : this);
