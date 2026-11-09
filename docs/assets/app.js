/* compositeSS calculator: user interface. Calculations are in compositess.js. */
(function () {
  'use strict';
  const $ = id => document.getElementById(id);
  const num = id => { const v = $(id).value.trim(); return v === '' ? NaN : Number(v); };
  const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const f3 = x => Number.isFinite(x) ? x.toFixed(3) : '–';
  const f4 = x => Number.isFinite(x) ? x.toFixed(4) : '–';
  function sig(x, d = 4) {
    if (!Number.isFinite(x)) return x === Infinity ? '∞' : x === -Infinity ? '−∞' : '–';
    if (x === 0) return '0';
    if (Math.abs(x) < 1e-4) return x.toExponential(d - 1);
    if (Math.abs(x) >= Math.pow(10, d)) return String(Math.round(x));
    return String(+x.toPrecision(d));
  }
  const intf = n => Number.isFinite(n) ? Math.round(n).toLocaleString('en-US') : '–';
  const css = name => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  function showErr(id, msg) { const e = $(id); if (msg) { e.textContent = msg; e.hidden = false; } else e.hidden = true; }

  // ------------------------------------------------------------------ tooltip
  const tip = $('tip');
  function showTip(html, ev) {
    tip.innerHTML = html; tip.hidden = false;
    const w = tip.offsetWidth, h = tip.offsetHeight;
    let x = ev.clientX + 14, y = ev.clientY + 14;
    if (x + w > window.innerWidth - 8) x = ev.clientX - w - 14;
    if (y + h > window.innerHeight - 8) y = ev.clientY - h - 14;
    tip.style.left = Math.max(8, x) + 'px'; tip.style.top = Math.max(8, y) + 'px';
  }
  const hideTip = () => { tip.hidden = true; };

  // ------------------------------------------------------------------ charts
  function ticks(lo, hi, n = 5) {
    const span = hi - lo || 1, step0 = span / n, mag = Math.pow(10, Math.floor(Math.log10(step0)));
    const step = [1, 2, 2.5, 5, 10].map(m => m * mag).find(s => s >= step0);
    const out = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(+v.toFixed(10));
    return out;
  }
  function frame(el, H) {
    const W = Math.max(280, el.clientWidth || 600);
    return { W, H, m: { l: 48, r: 16, t: 12, b: 40 } };
  }
  // line chart with crosshair tooltip
  function lineChart(el, o) {
    const { W, H, m } = frame(el, o.height || 300);
    const x0 = o.xDomain[0], x1 = o.xDomain[1], y0 = o.yDomain[0], y1 = o.yDomain[1];
    const sx = v => m.l + (v - x0) / (x1 - x0) * (W - m.l - m.r);
    const sy = v => H - m.b - (v - y0) / (y1 - y0) * (H - m.t - m.b);
    let s = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(o.label)}"><g class="grid">`;
    for (const t of ticks(y0, y1, 5)) s += `<line x1="${m.l}" x2="${W - m.r}" y1="${sy(t)}" y2="${sy(t)}"/><text x="${m.l - 8}" y="${sy(t) + 4}" text-anchor="end">${o.yFmt ? o.yFmt(t) : t}</text>`;
    s += '</g>';
    for (const t of ticks(x0, x1, Math.max(3, Math.floor((W - m.l) / 90)))) s += `<text x="${sx(t)}" y="${H - m.b + 16}" text-anchor="middle">${o.xFmt ? o.xFmt(t) : t}</text>`;
    s += `<line class="axis" x1="${m.l}" x2="${W - m.r}" y1="${H - m.b}" y2="${H - m.b}"/>`;
    s += `<text class="title-x" x="${(m.l + W - m.r) / 2}" y="${H - 6}" text-anchor="middle">${esc(o.xLabel)}</text>`;
    s += `<text class="title-y" transform="translate(12 ${(m.t + H - m.b) / 2}) rotate(-90)" text-anchor="middle">${esc(o.yLabel)}</text>`;
    if (o.target != null) s += `<line class="target" x1="${m.l}" x2="${W - m.r}" y1="${sy(o.target)}" y2="${sy(o.target)}"/><text class="lbl" x="${W - m.r}" y="${sy(o.target) - 5}" text-anchor="end">target ${o.target}</text>`;
    if (o.marker != null && o.marker >= x0 && o.marker <= x1) s += `<line class="marker" x1="${sx(o.marker)}" x2="${sx(o.marker)}" y1="${m.t}" y2="${H - m.b}"/><text class="lbl" x="${sx(o.marker) + 5}" y="${m.t + 10}">${esc(o.markerLabel || '')}</text>`;
    o.series.forEach(se => {
      const pts = o.x.map((xv, i) => Number.isFinite(se.values[i]) ? `${sx(xv).toFixed(1)},${sy(se.values[i]).toFixed(1)}` : null).filter(Boolean);
      s += `<polyline class="series" style="stroke:var(${se.color})" points="${pts.join(' ')}"/>`;
    });
    s += `<g class="hover" visibility="hidden"><line class="cross" y1="${m.t}" y2="${H - m.b}"/>${o.series.map(se => `<circle class="dot" r="4" style="fill:var(${se.color})"/>`).join('')}</g>`;
    s += `<rect class="hit" x="${m.l}" y="${m.t}" width="${W - m.l - m.r}" height="${H - m.t - m.b}"/></svg>`;
    el.innerHTML = s;
    const svg = el.querySelector('svg'), hov = svg.querySelector('.hover'), hit = svg.querySelector('.hit');
    hit.addEventListener('pointermove', ev => {
      const r = svg.getBoundingClientRect(), px = (ev.clientX - r.left) * W / r.width;
      const xv = x0 + (px - m.l) / (W - m.l - m.r) * (x1 - x0);
      let i = 0; for (let j = 1; j < o.x.length; j++) if (Math.abs(o.x[j] - xv) < Math.abs(o.x[i] - xv)) i = j;
      hov.setAttribute('visibility', 'visible');
      const cx = sx(o.x[i]);
      hov.querySelector('line').setAttribute('x1', cx); hov.querySelector('line').setAttribute('x2', cx);
      hov.querySelectorAll('circle').forEach((c, k) => {
        const v = o.series[k].values[i];
        c.setAttribute('cx', cx); c.setAttribute('cy', Number.isFinite(v) ? sy(v) : -99);
      });
      showTip(o.tip(i), ev);
    });
    hit.addEventListener('pointerleave', () => { hov.setAttribute('visibility', 'hidden'); hideTip(); });
  }
  // horizontal bars, label above each bar
  function barChart(el, rows, o) {
    const W = Math.max(280, el.clientWidth || 600), rowH = 40, top = 6, H = top + rows.length * rowH + 22;
    const max = Math.max(...rows.map(r => r.value).filter(Number.isFinite), 1);
    const sx = v => 4 + v / max * (W - 70);
    let s = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(o.label)}">`;
    rows.forEach((r, i) => {
      const y = top + i * rowH;
      s += `<text class="lbl" x="4" y="${y + 12}">${esc(r.label)}</text>`;
      if (Number.isFinite(r.value)) {
        const w = Math.max(2, sx(r.value) - 4);
        s += `<path class="bar${r.ref ? ' ref' : ''}" d="M4 ${y + 18}h${w - 4}a4 4 0 0 1 4 4v6a4 4 0 0 1 -4 4h${-(w - 4)}z"/>`;
        s += `<text class="lbl" x="${4 + w + 6}" y="${y + 29}">${intf(r.value)}</text>`;
      } else s += `<text class="lbl" x="4" y="${y + 29}">not reached</text>`;
      s += `<rect class="hit" x="0" y="${y}" width="${W}" height="${rowH}" data-i="${i}"/>`;
    });
    s += `<text class="title-x" x="4" y="${H - 4}">${esc(o.xLabel)}</text></svg>`;
    el.innerHTML = s;
    el.querySelectorAll('.hit').forEach(h => {
      h.addEventListener('pointermove', ev => showTip(o.tip(+h.dataset.i), ev));
      h.addEventListener('pointerleave', hideTip);
    });
  }
  // sequential colour scale from the --seq-* tokens
  function seqColor(v) {
    const stops = [0, 1, 2, 3, 4, 5, 6].map(i => css('--seq-' + i));
    const rgb = h => [1, 3, 5].map(i => parseInt(h.slice(i, i + 2), 16));
    const t = Math.min(1, Math.max(0, v)) * 6, i = Math.min(5, Math.floor(t)), f = t - i;
    const a = rgb(stops[i]), b = rgb(stops[i + 1]);
    return `rgb(${a.map((x, k) => Math.round(x + (b[k] - x) * f)).join(',')})`;
  }
  // marching squares: line segments of the level set z = lev on a grid z[i][j] (i along x, j along y)
  function contourSegments(xs, ys, z, lev) {
    const segs = [];
    const interp = (p1, p2, v1, v2) => [p1[0] + (lev - v1) / (v2 - v1) * (p2[0] - p1[0]), p1[1] + (lev - v1) / (v2 - v1) * (p2[1] - p1[1])];
    for (let i = 0; i < xs.length - 1; i++) for (let j = 0; j < ys.length - 1; j++) {
      const P = [[xs[i], ys[j]], [xs[i + 1], ys[j]], [xs[i + 1], ys[j + 1]], [xs[i], ys[j + 1]]];
      const V = [z[i][j], z[i + 1][j], z[i + 1][j + 1], z[i][j + 1]];
      const pts = [];
      for (let k = 0; k < 4; k++) {
        const a = V[k], b = V[(k + 1) % 4];
        if ((a < lev) !== (b < lev)) pts.push(interp(P[k], P[(k + 1) % 4], a, b));
      }
      if (pts.length === 2) segs.push(pts);
      else if (pts.length === 4) { segs.push([pts[0], pts[1]]); segs.push([pts[2], pts[3]]); }
    }
    return segs;
  }
  function heatmap(el, xs, ys, z, o) {
    const W = Math.max(280, Math.min(el.clientWidth || 520, 560)), m = { l: 52, r: 12, t: 10, b: 42 };
    const H = Math.round((W - m.l - m.r) * 0.85) + m.t + m.b;
    const dx = (xs[1] - xs[0]) / 2, dy = (ys[1] - ys[0]) / 2;
    const X0 = xs[0] - dx, X1 = xs[xs.length - 1] + dx, Y0 = ys[0] - dy, Y1 = ys[ys.length - 1] + dy;
    const sx = v => m.l + (v - X0) / (X1 - X0) * (W - m.l - m.r);
    const sy = v => H - m.b - (v - Y0) / (Y1 - Y0) * (H - m.t - m.b);
    let s = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(o.label)}">`;
    const cw = sx(xs[1]) - sx(xs[0]), ch = sy(ys[0]) - sy(ys[1]);
    xs.forEach((x, i) => ys.forEach((y, j) => {
      s += `<rect x="${(sx(x) - cw / 2).toFixed(2)}" y="${(sy(y) - ch / 2).toFixed(2)}" width="${(cw + 0.6).toFixed(2)}" height="${(ch + 0.6).toFixed(2)}" fill="${seqColor(z[i][j])}" data-i="${i}" data-j="${j}"/>`;
    }));
    for (const lev of [0.5, 0.7, 0.8, 0.9]) {
      const segs = contourSegments(xs, ys, z, lev);
      if (!segs.length) continue;
      s += `<path class="contour" d="${segs.map(g => `M${sx(g[0][0]).toFixed(1)} ${sy(g[0][1]).toFixed(1)}L${sx(g[1][0]).toFixed(1)} ${sy(g[1][1]).toFixed(1)}`).join('')}"/>`;
      // label each level near the upper end of its line, clear of the planning-value dot
      const bySum = segs.slice().sort((a, b) => (a[0][0] + a[0][1]) - (b[0][0] + b[0][1]));
      const g = bySum[Math.floor(bySum.length * 0.82)];
      s += `<text class="clabel" x="${sx((g[0][0] + g[1][0]) / 2).toFixed(1)}" y="${sy((g[0][1] + g[1][1]) / 2).toFixed(1)}" text-anchor="middle" dy="3">${lev}</text>`;
    }
    if (o.point) s += `<circle class="plan" r="5" cx="${sx(o.point[0])}" cy="${sy(o.point[1])}"><title>Planning values</title></circle>`;
    for (const t of ticks(xs[0], xs[xs.length - 1], 5)) s += `<text x="${sx(t)}" y="${H - m.b + 15}" text-anchor="middle">${t.toFixed(2)}</text>`;
    for (const t of ticks(ys[0], ys[ys.length - 1], 5)) s += `<text x="${m.l - 6}" y="${sy(t) + 4}" text-anchor="end">${t.toFixed(2)}</text>`;
    s += `<text class="title-x" x="${(m.l + W - m.r) / 2}" y="${H - 6}" text-anchor="middle">${esc(o.xLabel)}</text>`;
    s += `<text class="title-y" transform="translate(12 ${(m.t + H - m.b) / 2}) rotate(-90)" text-anchor="middle">${esc(o.yLabel)}</text></svg>`;
    el.innerHTML = s;
    const svg = el.querySelector('svg');
    svg.addEventListener('pointermove', ev => {
      const t = ev.target;
      if (t.dataset && t.dataset.i !== undefined) {
        const i = +t.dataset.i, j = +t.dataset.j;
        showTip(o.tip(xs[i], ys[j], z[i][j]), ev);
      } else hideTip();
    });
    svg.addEventListener('pointerleave', hideTip);
  }
  function histogram(el, values, o) {
    const W = Math.max(280, el.clientWidth || 600), H = 200, m = { l: 44, r: 12, t: 10, b: 38 }, nb = 30;
    const counts = new Array(nb).fill(0);
    values.forEach(v => { counts[Math.min(nb - 1, Math.max(0, Math.floor(v * nb)))]++; });
    const max = Math.max(...counts, 1);
    const sx = v => m.l + v * (W - m.l - m.r), sy = c => H - m.b - c / max * (H - m.t - m.b);
    let s = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(o.label)}"><g class="grid">`;
    for (const t of ticks(0, max, 4)) s += `<line x1="${m.l}" x2="${W - m.r}" y1="${sy(t)}" y2="${sy(t)}"/><text x="${m.l - 6}" y="${sy(t) + 4}" text-anchor="end">${t}</text>`;
    s += '</g>';
    const bw = (W - m.l - m.r) / nb;
    counts.forEach((c, i) => {
      if (c > 0) {
        const x = sx(i / nb) + 1, w = bw - 2, y = sy(c), h = H - m.b - y, r = Math.min(4, w / 2, h);
        s += `<path class="bar" d="M${x} ${H - m.b}V${y + r}a${r} ${r} 0 0 1 ${r} ${-r}h${w - 2 * r}a${r} ${r} 0 0 1 ${r} ${r}V${H - m.b}z"/>`;
      }
      s += `<rect class="hit" x="${sx(i / nb)}" y="${m.t}" width="${bw}" height="${H - m.t - m.b}" data-i="${i}"/>`;
    });
    s += `<line class="axis" x1="${m.l}" x2="${W - m.r}" y1="${H - m.b}" y2="${H - m.b}"/>`;
    for (const t of [0, 0.2, 0.4, 0.6, 0.8, 1]) s += `<text x="${sx(t)}" y="${H - m.b + 15}" text-anchor="middle">${t}</text>`;
    if (o.marker != null) s += `<line class="marker" x1="${sx(o.marker)}" x2="${sx(o.marker)}" y1="${m.t}" y2="${H - m.b}"/><text class="lbl" x="${sx(o.marker) + 5}" y="${m.t + 10}">assurance ${o.marker.toFixed(3)}</text>`;
    s += `<text class="title-x" x="${(m.l + W - m.r) / 2}" y="${H - 6}" text-anchor="middle">${esc(o.xLabel)}</text>`;
    s += `<text class="title-y" transform="translate(12 ${(m.t + H - m.b) / 2}) rotate(-90)" text-anchor="middle">Prior draws</text></svg>`;
    el.innerHTML = s;
    el.querySelectorAll('.hit').forEach(h => {
      const i = +h.dataset.i;
      h.addEventListener('pointermove', ev => showTip(`Power ${(i / nb).toFixed(2)}–${((i + 1) / nb).toFixed(2)}<br>${counts[i]} of ${values.length} draws`, ev));
      h.addEventListener('pointerleave', hideTip);
    });
  }

  // ------------------------------------------------------------------ design form
  function syncForm() {
    const b = document.body.classList;
    const est = document.querySelector('input[name=estimand]:checked').value;
    const binary = $('efficacy').value === 'binary', wz = $('method').value === 'wz';
    const aligned = est === 'cs' && binary && wz && $('aligned').checked;
    b.toggle('est-cs', est === 'cs'); b.toggle('est-theta', est === 'theta');
    b.toggle('eff-binary', binary); b.toggle('method-wz', wz);
    b.toggle('aligned-ok', binary && wz); b.toggle('is-aligned', aligned);
    b.toggle('fixedw', wz && !aligned);
    b.toggle('contingent', $('contingent').checked && !aligned);
    $('wdOut').textContent = Number($('wd').value).toFixed(2);
    b.toggle('mm-param', $('mmStruct').value !== 'newmeds');
    const anWz = $('anMethod').value === 'wz', anAl = anWz && $('anAligned').checked;
    b.toggle('an-wz', anWz); b.toggle('an-aligned', anAl); b.toggle('an-fixed', anWz && !anAl);
    b.toggle('an-gate', $('anCont').checked && !anAl);
    $('anWdOut').textContent = Number($('anWd').value).toFixed(2);
  }
  function readDesign() {
    const est = document.querySelector('input[name=estimand]:checked').value;
    const g = { pA: num('pA'), pB: num('pB'), muA: num('muA'), muB: num('muB'), sigma: num('sigma'),
      alpha: num('alpha'), power: num('power'), r: num('r') };
    if (!(g.power > 0 && g.power < 1)) throw new Error('Power must be between 0 and 1.');
    if (est === 'theta') {
      for (const k of ['muA', 'muB']) if (!Number.isFinite(g[k])) throw new Error('Enter the completers’ mean improvement in each arm.');
      if (!(g.r > 0)) throw new Error("'r' must be positive.");
      if (!(g.alpha > 0 && g.alpha < 0.5)) throw new Error("'alpha' must be in (0, 0.5).");
      const d = CS.thetaDesign({ pA: g.pA, pB: g.pB, muA: g.muA, muB: g.muB, sigma: g.sigma, s: num('s') || 0,
        r: g.r, alpha: g.alpha, margin: num('margin_t') || 0 });
      return { d, g, est };
    }
    const binary = $('efficacy').value === 'binary', method = $('method').value;
    const aligned = binary && method === 'wz' && $('aligned').checked;
    const wd = method === 'fisher' ? 0.5 : Number($('wd').value);
    const d = CS.csDesign({
      pA: g.pA, pB: g.pB, delta: binary ? null : g.muA - g.muB, sigma: binary ? null : g.sigma, r: g.r,
      weights: aligned ? 'aligned' : [wd, 1 - wd], method,
      contingent: $('contingent').checked && !aligned, c: num('gate') || 0, calibration: $('cal').value,
      alpha: g.alpha, gamma: $('gamma').value, efficacy: binary ? 'binary' : 'continuous',
      respA: num('respA'), respB: num('respB'),
      margin_d: method === 'wz' ? (num('margin_d') || 0) : 0, margin_e: method === 'wz' ? (num('margin_e') || 0) : 0 });
    return { d, g, est };
  }
  function sampleSize(D) {
    return D.d.type === 'theta' ? CS.thetaSampleSize(D.d, D.g.power) : CS.csSampleSize(D.d, D.g.power);
  }

  // ------------------------------------------------------------------ cache of the current design
  let cur = null;
  function compute() {
    try {
      const D = readDesign();
      const ss = sampleSize(D);
      cur = { D, ss, err: null };
    } catch (e) {
      let D = null;
      try { D = readDesign(); } catch (e2) { /* invalid design */ }
      cur = { D, ss: null, err: e.message };
    }
    return cur;
  }

  // ------------------------------------------------------------------ sample size view
  function renderSampleSize() {
    const c = cur;
    showErr('ssErr', c.err);
    $('ssBody').hidden = !!c.err;
    if (c.err) return;
    const { D, ss } = c, d = D.d;
    $('ssTotal').textContent = intf(ss.total); $('ssA').textContent = intf(ss.nA);
    $('ssB').textContent = intf(ss.nB); $('ssPow').textContent = f3(ss.power);
    const facts = [];
    if (d.type === 'theta') {
      const dc = CS.thetaDecompose(d);
      facts.push(['θ', `${f3(dc.theta)} = completion part ${f3(dc.completion)} + completers part ${f3(dc.completers)}`]);
      facts.push(['Variance per patient', `A ${f3(ss.VA)}, B ${f3(ss.VB)}`]);
    } else {
      if (d.aligned) facts.push(['Estimand-aligned weights', `${f3(d.weights[0])} / ${f3(d.weights[1])}; the test targets the probability of completing and responding`]);
      else if (d.method === 'wz') {
        try {
          const ow = CS.csOptWeights(d, ss.nA);
          facts.push(['Power-maximizing weights', `${ow.weights[0].toFixed(2)} / ${ow.weights[1].toFixed(2)} (power ${f3(ow.power_optimal)} vs ${f3(ow.power_design)} with the chosen weights)`]);
        } catch (e) { facts.push(['Power-maximizing weights', e.message]); }
      }
      facts.push(['Critical value', `${ss.critical.k.toFixed(ss.critical.scale === 'Z' ? 4 : 6)} (${ss.critical.scale === 'Z' ? 'Z scale' : 'product of P values'})`]);
      facts.push(['γ', CS.csComponents(d, ss.nA).gamma.toFixed(4)]);
      if (d.efficacy === 'continuous') facts.push(['ψ', `${f3(CS.csPsi(d).psi)} (information weights), ${f3(CS.csPsi(d, 'clinical').psi)} (clinical weights)`]);
    }
    facts.push(['Dropout', 'No inflation for dropout is needed: dropout is part of the endpoint.']);
    $('ssFacts').innerHTML = facts.map(([k, v]) => `<dt>${k}</dt><dd>${esc(v)}</dd>`).join('');

    const nHi = Math.max(50, Math.ceil(2 * ss.nA));
    const ns = [...new Set(Array.from({ length: 80 }, (_, i) => Math.round(5 + i * (nHi - 5) / 79)))];
    const pw = d.type === 'theta' ? ns.map(n => CS.thetaPower(d, n)) : CS.csPower(d, ns).map(r => r.power);
    lineChart($('curve'), {
      label: 'Power by number of patients on A', x: ns, series: [{ values: pw, color: '--series-1' }],
      xDomain: [ns[0], ns[ns.length - 1]], yDomain: [0, 1], xLabel: 'Patients on A', yLabel: 'Power',
      target: D.g.power, marker: ss.nA, markerLabel: `n = ${ss.nA}`,
      tip: i => `n<sub>A</sub> = ${ns[i]} (n<sub>B</sub> = ${Math.ceil(d.r * ns[i])})<br>Power ${f3(pw[i])}`
    });
    if (d.type === 'theta') {
      $('pwTab').innerHTML = `<thead><tr><th>n<sub>A</sub></th><th>Power</th></tr></thead><tbody><tr><td>${ss.nA}</td><td>${f4(ss.power)}</td></tr></tbody>`;
    } else {
      const p = CS.csPower(d, ss.nA)[0];
      $('pwTab').innerHTML = `<thead><tr><th>n<sub>A</sub></th><th>n<sub>B</sub></th><th>Composite test</th><th>Completion alone</th><th>Completers alone</th></tr></thead>
        <tbody><tr><td>${p.nA}</td><td>${sig(p.nB, 6)}</td><td>${f4(p.power)}</td><td>${f4(p.power_completion)}</td><td>${f4(p.power_completers)}</td></tr></tbody>`;
    }
    $('estimandTxt').textContent = CS.csEstimandStatement(d);
  }

  // ------------------------------------------------------------------ compare view
  function renderCompare() {
    const g = cur.D ? cur.D.g : null;
    try {
      if (!g) throw new Error(cur.err);
      if (!Number.isFinite(g.muA) || !Number.isFinite(g.muB) || !(g.sigma > 0)) throw new Error('Enter the completers’ means and SD in the Design panel.');
      const rows = CS.csCompareDesigns({ pA: g.pA, pB: g.pB, delta: g.muA - g.muB, sigma: g.sigma, muB: g.muB,
        s: num('s') || 0, alpha: g.alpha, power: g.power, gamma: $('cmpGamma').value });
      showErr('cmpErr', null);
      $('cmpTab').innerHTML = `<thead><tr><th>Design</th><th>Patients per arm</th><th>Power</th></tr></thead><tbody>${rows.map(r =>
        `<tr><td>${esc(r.design)}</td><td>${Number.isFinite(r.n_per_arm) ? intf(r.n_per_arm) : '<span class="na">not reached</span>'}</td><td class="${Number.isFinite(r.power) ? '' : 'na'}">${Number.isFinite(r.power) ? f3(r.power) : '–'}</td></tr>`).join('')}</tbody>`;
      barChart($('cmpChart'), rows.map(r => ({ label: r.design, value: r.n_per_arm, ref: / only$/.test(r.design) })), {
        label: 'Patients per arm by design', xLabel: 'Patients per arm (1:1 allocation); grey bars are reference designs',
        tip: i => `${esc(rows[i].design)}<br>${Number.isFinite(rows[i].n_per_arm) ? intf(rows[i].n_per_arm) + ' per arm' : 'not reached'}${Number.isFinite(rows[i].power) ? '<br>power ' + f3(rows[i].power) : ''}`
      });
    } catch (e) {
      showErr('cmpErr', e.message); $('cmpTab').innerHTML = ''; $('cmpChart').innerHTML = '';
    }
    if (!$('t6Tab').dataset.done) {
      const t6 = CS.csTable6();
      $('t6Tab').innerHTML = `<thead><tr><th>Design</th><th>Placebo-controlled</th><th>Active comparator</th></tr></thead><tbody>${t6.map(r =>
        `<tr><td>${esc(r.design)}</td><td>${intf(r.placebo_controlled)}</td><td>${intf(r.active_comparator)}</td></tr>`).join('')}</tbody>`;
      $('t6Tab').dataset.done = '1';
    }
  }

  // ------------------------------------------------------------------ sensitivity view
  function senN() {
    const n = num('senN');
    if (Number.isFinite(n)) { if (!(n >= 2)) throw new Error('Patients on A must be at least 2.'); return Math.round(n); }
    if (!cur.ss) throw new Error(cur.err);
    return cur.ss.nA;
  }
  function renderSensitivity() {
    try {
      if (!cur.D) throw new Error(cur.err);
      const d = cur.D.d, n = senN();
      const rng = (a, b) => { if (!(a > 0 && b < 1 && a < b)) throw new Error('Completion ranges must satisfy 0 < from < to < 1.'); return Array.from({ length: 21 }, (_, i) => a + (b - a) * i / 20); };
      const pa = rng(num('senA1'), num('senA2')), pb = rng(num('senB1'), num('senB2'));
      let z;
      if (d.type === 'theta') z = pa.map(a => pb.map(b => CS.thetaPower(Object.assign({}, d, { pA: a, pB: b }), n)));
      else {
        const grid = CS.csPowerGrid(d, n, pa, pb, [d.delta]);
        z = pa.map((_, i) => pb.map((_, j) => grid[j * pa.length + i].power));
      }
      showErr('senErr', null);
      $('heatTitle').textContent = `Power at n_A = ${n} by true completion rates`;
      heatmap($('heat'), pa, pb, z, {
        label: 'Power by true completion rates', xLabel: 'True completion, A', yLabel: 'True completion, B', point: [d.pA, d.pB],
        tip: (a, b, v) => `Completion A ${a.toFixed(3)}, B ${b.toFixed(3)}<br>Power ${f3(v)}`
      });
      $('heatLegend').innerHTML = `<span>Power 0</span><span class="ramp" style="background:linear-gradient(90deg,${[0, 1, 2, 3, 4, 5, 6].map(i => css('--seq-' + i)).join(',')})"></span><span>1</span><span>· contours at 0.5, 0.7, 0.8, 0.9 · dot: planning values</span>`;
      // completers' effect
      const g = cur.D.g, eff0 = g.muA - g.muB;
      if (d.type === 'cs' && d.efficacy === 'binary') {
        $('deltaCurve').innerHTML = '<p class="hint">Not available for a response outcome: the planning value is a pair of response rates, not a difference in means.</p>';
      } else {
        const eff = Array.from({ length: 50 }, (_, i) => i * 2 * Math.max(Math.abs(eff0), 1) / 49);
        const pw = eff.map(e => d.type === 'theta' ? CS.thetaPower(Object.assign({}, d, { muA: d.muB + e }), n)
          : CS.csPowerGrid(d, n, [d.pA], [d.pB], [e])[0].power);
        lineChart($('deltaCurve'), {
          label: 'Power by true completers’ difference', x: eff, series: [{ values: pw, color: '--series-1' }],
          xDomain: [0, eff[eff.length - 1]], yDomain: [0, 1], xLabel: 'True completers’ difference in mean improvement',
          yLabel: 'Power', target: g.power, marker: eff0, markerLabel: 'planning value', height: 300, xFmt: t => +t.toFixed(2),
          tip: i => `Difference ${eff[i].toFixed(2)}<br>Power ${f3(pw[i])} at n<sub>A</sub> = ${n}`
        });
      }
    } catch (e) {
      showErr('senErr', e.message); $('heat').innerHTML = ''; $('deltaCurve').innerHTML = ''; $('heatLegend').innerHTML = '';
    }
  }
  function runAssurance() {
    try {
      if (!cur.D) throw new Error(cur.err);
      const d = cur.D.d, n = senN();
      const bp = (m, k) => { if (!(m > 0 && m < 1) || !(k > 0)) throw new Error('Prior means must be between 0 and 1 and effective n positive.'); return [m * k, (1 - m) * k]; };
      const sd = num('prDsd');
      const a = CS.csAssurance(d, n, { prior_pA: bp(num('prAm'), num('prAn')), prior_pB: bp(num('prBm'), num('prBn')),
        prior_delta: sd > 0 ? [num('prDm'), sd] : null, nsim: Math.max(100, Math.round(num('prN') || 2000)), seed: num('prSeed') || 1 });
      $('assurOut').innerHTML = `<dl class="facts"><dt>Assurance at n<sub>A</sub> = ${n}</dt><dd><strong>${f3(a.assurance)}</strong> (Monte Carlo SE ${f4(a.mc_se)})</dd><dt>Power at the planning values</dt><dd>${f3(a.power_at_planning_values)}</dd></dl>`;
      histogram($('assurHist'), a.draws.map(x => x.power), { label: 'Power over prior draws', xLabel: 'Power over prior draws', marker: a.assurance });
    } catch (e) { $('assurOut').innerHTML = `<div class="error">${esc(e.message)}</div>`; $('assurHist').innerHTML = ''; }
  }

  // ------------------------------------------------------------------ simulation view
  let simJob = null;
  function runSim() {
    showErr('simErr', null); $('simOut').innerHTML = '';
    try {
      if (!cur.D) throw new Error(cur.err);
      const D = cur.D, d = D.d;
      if (d.type !== 'cs' || d.efficacy !== 'continuous' || d.margin_d !== 0 || d.margin_e !== 0)
        throw new Error('Simulation needs the composite test with a continuous outcome (superiority).');
      let n = num('simN');
      if (!Number.isFinite(n)) { if (!cur.ss) throw new Error(cur.err); n = cur.ss.nA; }
      if (!(n >= 4)) throw new Error('Patients on A must be at least 4.');
      const nsim = Math.max(100, Math.min(100000, Math.round(num('simNsim') || 2000)));
      const job = CS.csSimulation(d, Math.round(n), { nsim, null: $('simNull').checked, muB: D.g.muB, s: num('s') || 0, seed: num('simSeed') || 2026 });
      simJob = job;
      $('goSim').disabled = true; $('stopSim').hidden = false; $('simProg').hidden = false;
      let chunk = 50;
      const tick = () => {
        if (simJob !== job) return;
        const t0 = performance.now();
        const fin = job.step(chunk);
        const dt = performance.now() - t0;
        chunk = Math.max(10, Math.round(chunk * Math.min(4, Math.max(0.25, 40 / Math.max(dt, 1)))));
        $('simBar').style.width = (70 * job.done / job.total) + '%';
        $('simPct').textContent = `${job.done.toLocaleString('en-US')} of ${job.total.toLocaleString('en-US')} trials`;
        if (fin) finishSim(job); else setTimeout(tick, 0);
      };
      setTimeout(tick, 0);
    } catch (e) { showErr('simErr', e.message); }
  }
  function finishSim(job) {
    simJob = null; $('goSim').disabled = false; $('stopSim').hidden = true;
    const r = job.result();
    if (!r.nsim) { showErr('simErr', 'No usable trials: too few completers.'); return; }
    const lab = r.null ? 'Type I error' : 'Power';
    const rows = [
      ['Trials', `${r.nsim.toLocaleString('en-US')} (n<sub>A</sub> = ${r.n.A}, n<sub>B</sub> = ${r.n.B}), dropout unrelated to outcome, ${r.null ? 'under H<sub>0</sub>' : 'under H<sub>1</sub>'}`],
      [`${lab}, composite test`, `<strong>${f4(r.reject[0])}</strong> (Monte Carlo SE ${f4(r.reject[1])})${r.null ? '' : `; calculated ${f4(r.calculated_power)}`}`],
      ['Rejection by each component alone', `completion ${f4(r.reject_completion)}, completers ${f4(r.reject_completers)}`],
      ['Correlation of Z(d) and Z(e)', f4(r.cor_Zd_Ze)],
      ['ψ', `true ${f4(r.psi_true)}, CI coverage ${f4(r.psi_coverage)}`],
      ['θ (dropouts scored s)', `true ${f4(r.theta_true)}, CI coverage ${f4(r.theta_coverage)}, rejection ${f4(r.theta_reject)}`]];
    $('simOut').innerHTML = `<dl class="facts">${rows.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('')}</dl>`;
  }

  // ------------------------------------------------------------------ MMRM view
  function renderMmrm() {
    try {
      if (!cur.D) throw new Error(cur.err);
      const g = cur.D.g, structure = $('mmStruct').value;
      const md = CS.csMmrmDesign({ d: num('mmD'), pA: g.pA, pB: g.pB, visits: structure === 'newmeds' ? 6 : Math.round(num('mmK')),
        structure, rho: num('mmRho'), pattern: $('mmPat').value, r: g.r, alpha: g.alpha });
      if (!(g.power > 0 && g.power < 1)) throw new Error('Power must be between 0 and 1.');
      const ss = CS.csMmrmSampleSize(md, g.power);
      showErr('mmErr', null); $('mmBody').hidden = false;
      $('mmTotal').textContent = intf(ss.total); $('mmA').textContent = intf(ss.nA); $('mmB').textContent = intf(ss.nB); $('mmPow').textContent = f3(ss.power);
      const facts = [['Variance inflation over complete data', `A ${f3(ss.phiA)}, B ${f3(ss.phiB)}`],
        ['Model', `${md.visits} post-baseline visits; ${structure === 'newmeds' ? 'NEWMEDS correlations' : structure === 'ar1' ? `AR(1), ρ = ${md.rho}` : `equal correlation ρ = ${md.rho}`}; dropout ${md.pattern === 'even' ? 'evenly over visits' : md.pattern === 'early' ? 'mostly early' : 'mostly late'}`]];
      if (cur.ss) facts.push(['Composite design in the Design panel', `${intf(cur.ss.total)} patients in total. The two answer different questions, so choose between them by the estimand, not by the number.`]);
      $('mmFacts').innerHTML = facts.map(([k, v]) => `<dt>${k}</dt><dd>${esc(v)}</dd>`).join('');
    } catch (e) { showErr('mmErr', e.message); $('mmBody').hidden = true; }
  }

  // ------------------------------------------------------------------ analyze view
  let data = null;
  function loadData(rows, label) {
    if (!rows.length) throw new Error('The file has no data rows.');
    data = rows;
    const cols = Object.keys(rows[0]);
    const guess = (re, k) => cols.find(c => re.test(c)) || cols[Math.min(k, cols.length - 1)];
    const fill = (id, opts, sel) => { $(id).innerHTML = opts.map(o => `<option${o === sel ? ' selected' : ''}>${esc(o)}</option>`).join(''); };
    fill('cArm', cols, guess(/^(arm|treat|group|trt)/i, 0));
    fill('cComp', cols, guess(/compl|finish/i, 1));
    fill('cY', cols, guess(/improv|change|outcome|^y$/i, 2));
    fill('cStrata', ['(none)'].concat(cols), '(none)');
    $('cCov').innerHTML = cols.map(c => `<label class="check"><input type="checkbox" value="${esc(c)}"> ${esc(c)}</label>`).join('');
    fillActive();
    $('colUi').hidden = false; $('goAn').disabled = false;
    $('fileInfo').textContent = `${label}: ${rows.length.toLocaleString('en-US')} rows, ${cols.length} columns.`;
    $('anHeadBox').hidden = false;
    $('anHead').innerHTML = `<thead><tr>${cols.map(c => `<th>${esc(c)}</th>`).join('')}</tr></thead><tbody>${rows.slice(0, 8).map(r => `<tr>${cols.map(c => `<td>${esc(r[c])}</td>`).join('')}</tr>`).join('')}</tbody>`;
    $('anOut').innerHTML = ''; showErr('anErr', null);
  }
  function fillActive() {
    const col = $('cArm').value;
    const lv = [...new Set(data.map(r => String(r[col]).trim()))].sort();
    $('cActive').innerHTML = lv.map(v => `<option>${esc(v)}</option>`).join('');
  }
  function runAnalysis() {
    showErr('anErr', null); $('anOut').innerHTML = '';
    try {
      const method = $('anMethod').value, aligned = method === 'wz' && $('anAligned').checked;
      const wd = method === 'fisher' ? 0.5 : Number($('anWd').value);
      const cov = [...$('cCov').querySelectorAll('input:checked')].map(i => i.value);
      const st = $('cStrata').value;
      const a = CS.csAnalyze(data, {
        arm: $('cArm').value, active: $('cActive').value, completed: $('cComp').value, outcome: $('cY').value,
        covariates: cov, strata: st === '(none)' ? null : st, weights: aligned ? 'aligned' : [wd, 1 - wd], method,
        contingent: $('anCont').checked && !aligned, c: num('anC') || 0, alpha: num('anAlpha'), s: num('anS') || 0,
        psi_type: $('anPsi').value, bootstrap: Math.max(0, Math.round(num('anBoot') || 0)), seed: 1 });
      const t = a.test, cl = Math.round(100 * a.conf_level);
      const testName = t.method !== 'wz' ? 'Fisher' : t.aligned ? 'weighted Z, estimand-aligned' : 'weighted Z';
      const rows = [
        ['Randomized', `A ${a.n.A}, B ${a.n.B}`],
        ['Completion', `A ${sig(a.completion.A)}, B ${sig(a.completion.B)}; Z(d) = ${sig(a.components.Z_d)} (P = ${sig(a.components.P_d)})`],
        ['Completers’ effect', `${sig(a.completers_effect.estimate)} (SE ${sig(a.completers_effect.se)}); Z(e) = ${sig(a.components.Z_e)} (P = ${sig(a.components.P_e)})`],
        ['Composite test', `${testName}, weights ${sig(t.weights[0])} / ${sig(t.weights[1])}${t.contingent ? `, contingent c = ${t.c}, ${t.calibration}` : ''}`],
        ['Statistic', `${sig(t.statistic)} (critical value ${sig(t.critical)}); P = ${sig(t.p_value)}`],
        ['Decision', `<span class="decision ${t.reject ? 'yes' : ''}">${t.reject ? 'Reject H<sub>0</sub>' : 'Do not reject H<sub>0</sub>'}</span> at one-sided α = ${t.alpha}`]];
      if (t.contingent && !t.gate_open) rows.push(['Gate', 'Closed: A observed worse on completion; testing stopped.']);
      rows.push([`ψ (${a.psi_type})`, `${sig(a.psi.estimate)}, ${cl}% CI ${sig(a.psi.lower)} to ${sig(a.psi.upper)}`]);
      rows.push([`θ (s = ${a.s})`, `${sig(a.theta.estimate)}, ${cl}% CI ${sig(a.theta.lower)} to ${sig(a.theta.upper)}; one-sided P = ${sig(a.theta.p_one_sided)}`]);
      if (a.bootstrap) rows.push([`Bootstrap (B = ${a.bootstrap.B})`, `ψ ${sig(a.bootstrap.psi[0])} to ${sig(a.bootstrap.psi[1])}; θ ${sig(a.bootstrap.theta[0])} to ${sig(a.bootstrap.theta[1])}`]);
      let notes = '';
      if (a.n_completers_missing_outcome > 0) notes += `<p class="warn">${a.n_completers_missing_outcome} completer(s) with a missing final outcome were excluded.</p>`;
      if (a.aligned_info && a.aligned_info.warning) notes += `<p class="warn">${esc(a.aligned_info.warning)}</p>`;
      $('anOut').innerHTML = `<div class="panel"><h3>Composite-strategy analysis</h3><dl class="facts">${rows.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('')}</dl></div>${notes}`;
    } catch (e) { showErr('anErr', e.message); }
  }

  // ------------------------------------------------------------------ tabs and wiring
  const DESIGN_TABS = ['sample-size', 'compare', 'sensitivity', 'simulation', 'mmrm'];
  const RENDER = { 'sample-size': renderSampleSize, compare: renderCompare, sensitivity: renderSensitivity, mmrm: renderMmrm };
  let tab = 'sample-size';
  const dirty = new Set(DESIGN_TABS);
  function showTab(name) {
    if (!document.querySelector(`[data-view="${name}"]`)) name = 'sample-size';
    tab = name;
    document.querySelectorAll('[data-view]').forEach(v => { v.hidden = v.dataset.view !== name; });
    document.querySelectorAll('[data-tab]').forEach(a => a.setAttribute('aria-selected', a.dataset.tab === name ? 'true' : 'false'));
    $('layout').classList.toggle('nodesign', !DESIGN_TABS.includes(name));
    if (dirty.has(name) && RENDER[name]) { RENDER[name](); dirty.delete(name); }
  }
  let timer = null;
  function designChanged() {
    syncForm();
    clearTimeout(timer);
    timer = setTimeout(() => {
      compute();
      DESIGN_TABS.forEach(t => dirty.add(t));
      if (RENDER[tab]) { RENDER[tab](); dirty.delete(tab); }
    }, 120);
  }
  $('designPanel').addEventListener('input', designChanged);
  $('designPanel').addEventListener('change', designChanged);
  $('cmpGamma').addEventListener('change', renderCompare);
  $('loadT6').addEventListener('click', () => {
    Object.entries({ pA: 0.65, pB: 0.5, muA: 21, muB: 15, sigma: 18, alpha: 0.025, power: 0.8 }).forEach(([k, v]) => { $(k).value = v; });
    designChanged();
  });
  ['senN', 'senA1', 'senA2', 'senB1', 'senB2'].forEach(id => $(id).addEventListener('change', renderSensitivity));
  $('goAssur').addEventListener('click', runAssurance);
  $('goSim').addEventListener('click', runSim);
  $('stopSim').addEventListener('click', () => { if (simJob) { const j = simJob; finishSim(j); $('simPct').textContent += ' (stopped)'; } });
  ['mmD', 'mmStruct', 'mmRho', 'mmK', 'mmPat'].forEach(id => $(id).addEventListener('input', () => { syncForm(); renderMmrm(); }));
  document.querySelector('#view-analyze').addEventListener('input', syncForm);
  document.querySelector('#view-analyze').addEventListener('change', syncForm);
  $('file').addEventListener('change', () => {
    const f = $('file').files[0]; if (!f) return;
    const rd = new FileReader();
    rd.onload = () => { try { loadData(CS.parseCSV(String(rd.result)), f.name); } catch (e) { showErr('anErr', e.message); } };
    rd.onerror = () => showErr('anErr', 'The file could not be read.');
    rd.readAsText(f);
  });
  $('loadEx').addEventListener('click', async () => {
    try {
      const res = await fetch('example_trial.csv');
      if (!res.ok) throw new Error('The example data set could not be loaded.');
      loadData(CS.parseCSV(await res.text()), 'Example (simulated) trial');
    } catch (e) { showErr('anErr', e.message); }
  });
  $('cArm').addEventListener('change', fillActive);
  $('goAn').addEventListener('click', runAnalysis);
  document.querySelectorAll('[data-copy]').forEach(b => b.addEventListener('click', async () => {
    const el = $(b.dataset.copy), txt = el.textContent;
    try { await navigator.clipboard.writeText(txt); b.textContent = 'Copied'; }
    catch (e) { const r = document.createRange(); r.selectNodeContents(el); const s = getSelection(); s.removeAllRanges(); s.addRange(r); b.textContent = 'Selected: press Ctrl+C'; }
    setTimeout(() => { b.textContent = 'Copy'; }, 2000);
  }));
  window.addEventListener('hashchange', () => showTab(location.hash.slice(1)));
  let rz = null;
  window.addEventListener('resize', () => {
    clearTimeout(rz);
    rz = setTimeout(() => { DESIGN_TABS.forEach(t => dirty.add(t)); if (RENDER[tab]) { RENDER[tab](); dirty.delete(tab); } }, 200);
  });
  window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => { DESIGN_TABS.forEach(t => dirty.add(t)); if (RENDER[tab]) { RENDER[tab](); dirty.delete(tab); } });

  // On GitHub Pages (https://<user>.github.io/<repo>/) the repository is known
  // from the address, so the install line and repository link fill themselves in.
  (function repoInfo() {
    const m = location.hostname.match(/^([^.]+)\.github\.io$/i);
    const repo = location.pathname.split('/').filter(Boolean)[0];
    if (!m || !repo) return;
    const slug = `${m[1]}/${repo}`;
    $('installCode').textContent = `remotes::install_github("${slug}")\ncompositeSS::run_calculator()`;
    $('repoLink').innerHTML = `, <a href="https://github.com/${esc(slug)}">github.com/${esc(slug)}</a>`;
  })();

  syncForm();
  compute();
  showTab(location.hash.slice(1) || 'sample-size');
})();
