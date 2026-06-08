// Phase 2 analysis - foundations claimed by the lab atlas v1.0:
//   |F| triangles = 256
//   Betti numbers beta = (1, 0, 14, 0) of the primal simplicial complex
//   B_op (vertex-edge boundary, O_h-equivariant) rank = 12
//   dim Hom_G(V_edge240, V_dual144 (x) chi) = 744
//   cycle space H_1^dual irrep decomposition under O_h
//   shell-ratio phi-proximity diagnostic (5% bound)
//
// All routines are pure functions of the geometry; no gauge or convention
// from the corpus is invented.  Each routine is honest about its definition.

import { build } from "./geometry.mjs";
import { analyze } from "./analysis.mjs";
import {
  OH_CLASS_NAMES, OH_CLASS_SIZES, OH_GROUP_ORDER,
  OH_IRREP_NAMES, OH_IRREP_DIMS,
  OH_CHARACTER_TABLE, classifyOh, decomposeIrrep, formatDecomposition,
} from "./oh-character-table.mjs";

const k6 = (p) => p.map((x) => (Object.is(x, -0) ? 0 : +x).toFixed(6)).join(",");
const applyG = (g, v) => [g.s[0]*v[g.pi[0]], g.s[1]*v[g.pi[1]], g.s[2]*v[g.pi[2]]];

// ---- the 48 elements of O_h --------------------------------------------
function buildOh() {
  const perms = [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]];
  const sgn = (p) => { let s = 1; for (let i = 0; i < 3; i++) for (let j = i+1; j < 3; j++) if (p[i] > p[j]) s = -s; return s; };
  const G = [];
  for (const pi of perms) for (const sx of [-1,1]) for (const sy of [-1,1]) for (const sz of [-1,1])
    G.push({ pi, s: [sx,sy,sz], det: sgn(pi)*sx*sy*sz });
  return G;
}

// ---- Gaussian elimination, returns rank over R ------------------------
function rankReal(rows, cols, mat, tol = 1e-9) {
  const A = mat.map((r) => r.slice());
  let rank = 0, row = 0;
  for (let col = 0; col < cols && row < rows; col++) {
    let piv = -1, max = tol;
    for (let i = row; i < rows; i++) if (Math.abs(A[i][col]) > max) { max = Math.abs(A[i][col]); piv = i; }
    if (piv < 0) continue;
    if (piv !== row) [A[row], A[piv]] = [A[piv], A[row]];
    for (let i = 0; i < rows; i++) {
      if (i !== row && Math.abs(A[i][col]) > tol) {
        const f = A[i][col] / A[row][col];
        for (let k = col; k < cols; k++) A[i][k] -= f * A[row][k];
      }
    }
    rank++; row++;
  }
  return rank;
}

// ========================================================================
//  TRIANGLES, BETTI NUMBERS, BOUNDARY MAPS
// ========================================================================
export function trianglesAndBetti(M, A) {
  // map each vertex position to a stable integer id
  const vmap = new Map(A.vertices.map((v, i) => [k6(v.pos), i]));
  // ---- triangles (2-simplices) = the 4 face-triples of every tetrahedron -
  const triKey = (ids) => ids.slice().sort((a,b)=>a-b).join(",");
  const triMap = new Map(); // key -> ids sorted
  for (const t of M.tetrahedra) {
    const ids = t.vertices.map((v) => vmap.get(k6(v)));
    for (let i=0;i<4;i++) for (let j=i+1;j<4;j++) for (let k=j+1;k<4;k++) {
      const sorted = [ids[i], ids[j], ids[k]].sort((a,b)=>a-b);
      triMap.set(triKey(sorted), sorted);
    }
  }
  const triangles = [...triMap.values()];
  // also count 3-cliques in the bare 1-skeleton (sanity reference) ---------
  const adj = new Map(A.vertices.map((v) => [v.id, new Set()]));
  for (const e of A.edges) {
    const a = vmap.get(k6(e.a)), b = vmap.get(k6(e.b));
    adj.get(a).add(b); adj.get(b).add(a);
  }
  let cliques = 0;
  for (const v of A.vertices) {
    const nbrs = [...adj.get(v.id)].filter((u) => u > v.id);
    for (let i = 0; i < nbrs.length; i++) for (let j = i+1; j < nbrs.length; j++) {
      if (adj.get(nbrs[i]).has(nbrs[j])) cliques++;
    }
  }
  // ---- boundary matrices ----------------------------------------------
  const V = A.vertices.length;          // 63
  const E = A.edges.length;             // 240
  const F = triangles.length;           // 256 expected
  const T = M.tetrahedra.length;        // 64
  // edge id map (sorted vertex-id pair -> edge id)
  const edgeId = new Map();
  A.edges.forEach((e, i) => {
    const a = vmap.get(k6(e.a)), b = vmap.get(k6(e.b));
    const lo = Math.min(a,b), hi = Math.max(a,b);
    edgeId.set(`${lo},${hi}`, i);
  });
  // d_1 : edges -> vertices  (graph boundary, signed)
  const D1 = Array.from({ length: V }, () => new Array(E).fill(0));
  A.edges.forEach((e, i) => {
    const a = vmap.get(k6(e.a)), b = vmap.get(k6(e.b));
    const lo = Math.min(a,b), hi = Math.max(a,b);
    D1[hi][i] += 1; D1[lo][i] -= 1;
  });
  // d_2 : triangles -> edges
  //   for tri [a,b,c] sorted (a<b<c): d2 = [b,c] - [a,c] + [a,b]
  const triIndex = new Map(triangles.map((tri, i) => [triKey(tri), i]));
  const D2 = Array.from({ length: E }, () => new Array(F).fill(0));
  triangles.forEach((tri, ti) => {
    const [a,b,c] = tri;
    D2[edgeId.get(`${b},${c}`)][ti] += 1;
    D2[edgeId.get(`${a},${c}`)][ti] -= 1;
    D2[edgeId.get(`${a},${b}`)][ti] += 1;
  });
  // d_3 : tetrahedra -> triangles
  //   for tet [a,b,c,d] sorted (a<b<c<d):
  //   d3 = [b,c,d] - [a,c,d] + [a,b,d] - [a,b,c]
  const D3 = Array.from({ length: F }, () => new Array(T).fill(0));
  M.tetrahedra.forEach((t, ti) => {
    const ids = t.vertices.map((v) => vmap.get(k6(v))).sort((a,b)=>a-b);
    const [a,b,c,d] = ids;
    D3[triIndex.get(triKey([b,c,d]))][ti] += 1;
    D3[triIndex.get(triKey([a,c,d]))][ti] -= 1;
    D3[triIndex.get(triKey([a,b,d]))][ti] += 1;
    D3[triIndex.get(triKey([a,b,c]))][ti] -= 1;
  });
  // d_n ranks
  const r1 = rankReal(V, E, D1);
  const r2 = rankReal(E, F, D2);
  const r3 = rankReal(F, T, D3);
  // Betti numbers
  // beta_0 = V - r1
  // beta_1 = (E - r1) - r2
  // beta_2 = (F - r2) - r3
  // beta_3 = T - r3
  const beta = [V - r1, E - r1 - r2, F - r2 - r3, T - r3];
  return { triangles, cliques, V, E, F, T, ranks: { r1, r2, r3 }, beta };
}

// ========================================================================
//  B_op  (vertex-edge boundary restricted to Shell 1, rank claim = 12)
// ========================================================================
export function bopRank(M, A) {
  // identify Shell 1 vertex ids (radius sqrt 2, i.e. r2 == 2)
  const shell1 = A.shells.find((s) => Math.abs(s.r2 - 2) < 1e-9);
  if (!shell1) return { rank: 0, rows: 0, cols: 0 };
  // analysis.mjs stores ids under `vertexIds`
  const shell1Ids = new Set(shell1.vertexIds);
  const vmap = new Map(A.vertices.map((v, i) => [k6(v.pos), i]));
  // project d_1 codomain onto Shell 1 only (zero out non-Shell-1 rows)
  const rows = shell1Ids.size;
  const cols = A.edges.length;
  const idxOf = new Map([...shell1Ids].map((id, i) => [id, i]));
  const Bop = Array.from({ length: rows }, () => new Array(cols).fill(0));
  A.edges.forEach((e, ei) => {
    const a = vmap.get(k6(e.a)), b = vmap.get(k6(e.b));
    const lo = Math.min(a,b), hi = Math.max(a,b);
    if (shell1Ids.has(hi)) Bop[idxOf.get(hi)][ei] += 1;
    if (shell1Ids.has(lo)) Bop[idxOf.get(lo)][ei] -= 1;
  });
  const rank = rankReal(rows, cols, Bop);
  return { rank, rows, cols, definition: "d_1 with codomain projected onto Shell 1 (12 vertices)" };
}

// ========================================================================
//  Hom_G(V_edge240, V_dual144 (x) chi_chiral)   character-formula
//
//  Computed under several natural rep conventions because the corpus's
//  exact orientation/signing of V_edge and V_dual is not given.  We report
//  every variant honestly so the matching one (if any) to 744 is visible.
// ========================================================================
export function homGdim(M, A) {
  const G = buildOh();
  const edgeKey = (a, b) => {
    const ka = k6(a), kb = k6(b);
    return ka < kb ? ka + "|" + kb : kb + "|" + ka;
  };
  const tetraKey = (t) => t.vertices.map(k6).sort().join("/");
  const tetraKeyToId = new Map(M.tetrahedra.map((t, i) => [tetraKey(t), i]));
  const dualPairs = A.dualGraphs.edge.pairs;

  const variants = {}; // name -> sum over G
  const init = (k) => { variants[k] = 0; };
  for (const k of ["edge_uns_dual_uns_det", "edge_sig_dual_uns_det", "edge_uns_dual_sig_det",
                   "edge_sig_dual_sig_det", "edge_uns_dual_uns_NoChi", "edge_sig_dual_sig_NoChi",
                   "edge_uns_dual_uns_A2g", "edge_uns_dual_uns_A2u"]) init(k);

  // helper: character of A2g (= sign on 6C2' and 6C4 and 6sigma_d and 6S4)
  // A2g classes signs: E +, 8C3 +, 6C2' -, 6C4 -, 3C2 +, i +, 8S6 +, 6sd -, 6S4 -, 3sh +
  function chiA2g(g) { const c = classifyOh(g); return OH_CHARACTER_TABLE[1][c]; }
  function chiA2u(g) { const c = classifyOh(g); return OH_CHARACTER_TABLE[6][c]; }

  for (const g of G) {
    // signed-vs-unsigned edge character
    let chiE_uns = 0, chiE_sig = 0;
    for (const e of A.edges) {
      const ka = k6(e.a), kb = k6(e.b);
      const lo = ka < kb ? ka : kb, hi = ka < kb ? kb : ka;
      const ga = k6(applyG(g, e.a)), gb = k6(applyG(g, e.b));
      const gLo = ga < gb ? ga : gb, gHi = ga < gb ? gb : ga;
      if (lo === gLo && hi === gHi) {
        chiE_uns += 1;
        // signed: +1 if orientation preserved, -1 if reversed
        //   orientation: from low canonical key to high; preserved iff g maps low->low
        const origLowKey = ka < kb ? ka : kb;
        const newLowKey  = ga < gb ? ga : gb;
        chiE_sig += (origLowKey === newLowKey) ? 1 : -1;
      }
    }
    // image of every tetra under g
    const tetraImage = M.tetrahedra.map((t) =>
      tetraKeyToId.get(t.vertices.map((v) => applyG(g, v)).map(k6).sort().join("/"))
    );
    // signed-vs-unsigned dual character
    let chiD_uns = 0, chiD_sig = 0;
    for (const [i, j] of dualPairs) {
      const ii = tetraImage[i], jj = tetraImage[j];
      if (ii === i && jj === j) { chiD_uns++; chiD_sig += 1; }
      else if (ii === j && jj === i) { chiD_uns++; chiD_sig += -1; }
    }
    const chiChir_det = g.det;
    const chiChir_A2g = chiA2g(g);
    const chiChir_A2u = chiA2u(g);
    variants.edge_uns_dual_uns_det   += chiE_uns * chiD_uns * chiChir_det;
    variants.edge_sig_dual_uns_det   += chiE_sig * chiD_uns * chiChir_det;
    variants.edge_uns_dual_sig_det   += chiE_uns * chiD_sig * chiChir_det;
    variants.edge_sig_dual_sig_det   += chiE_sig * chiD_sig * chiChir_det;
    variants.edge_uns_dual_uns_NoChi += chiE_uns * chiD_uns;
    variants.edge_sig_dual_sig_NoChi += chiE_sig * chiD_sig;
    variants.edge_uns_dual_uns_A2g   += chiE_uns * chiD_uns * chiChir_A2g;
    variants.edge_uns_dual_uns_A2u   += chiE_uns * chiD_uns * chiChir_A2u;
  }
  const dims = Object.fromEntries(Object.entries(variants).map(([k, v]) => [k, v / OH_GROUP_ORDER]));
  return { dims, claim: 744, match: Object.entries(dims).find(([, v]) => Math.round(v) === 744)?.[0] || null };
}

// ========================================================================
//  Irrep decomposition of H_1(dual graph) of dimension 81 under O_h
// ========================================================================
export function dualCycleSpaceIrreps(M, A) {
  const G = buildOh();
  const tetraKey = (t) => t.vertices.map(k6).sort().join("/");
  const tetraKeyToId = new Map(M.tetrahedra.map((t, i) => [tetraKey(t), i]));
  // for each g, compute chi_C0(g) = # tetra fixed and chi_C1(g) = signed # dual edges fixed.
  // chi_H1(g) = chi_C1(g) - chi_C0(g) + 1   (since beta_0 = 1 for connected dual)
  const dualPairs = A.dualGraphs.edge.pairs;
  const classChars = new Array(10).fill(0);
  const classCounts = new Array(10).fill(0);
  for (const g of G) {
    const tetraImage = M.tetrahedra.map((t) => {
      const k = t.vertices.map((v) => applyG(g, v)).map(k6).sort().join("/");
      return tetraKeyToId.get(k);
    });
    let chiC0 = 0;
    for (let i = 0; i < M.tetrahedra.length; i++) if (tetraImage[i] === i) chiC0++;
    let chiC1 = 0;
    for (const [i, j] of dualPairs) {
      const ii = tetraImage[i], jj = tetraImage[j];
      if (ii === i && jj === j) chiC1 += 1;
      else if (ii === j && jj === i) chiC1 -= 1; // orientation reversed
    }
    const chiH1 = chiC1 - chiC0 + 1;
    const c = classifyOh(g);
    classChars[c] += chiH1;
    classCounts[c]++;
  }
  // collapse to per-class character (divide by class count, but they should be equal for one element)
  const chars = classChars.map((sum, c) => classCounts[c] > 0 ? sum / classCounts[c] : 0);
  const mults = decomposeIrrep(chars);
  return { chars, mults, formatted: formatDecomposition(mults) };
}

// ========================================================================
//  Forman-Ricci curvature per edge (graph version)
//
//  Simple combinatorial form for an unweighted graph:
//     F(e=(u,v)) = 4 - deg(u) - deg(v)
//  cf. Sreejith-Forman 2016, the "augmented" Forman curvature for graphs.
//  This is a graph-level approximation of the full Forman 2003 cell-complex
//  formula; sufficient for a structural "gravity-readout" diagnostic.
// ========================================================================
export function formanRicci(M, A) {
  const vmap = new Map(A.vertices.map((v, i) => [k6(v.pos), i]));
  const degree = new Array(A.vertices.length).fill(0);
  for (const e of A.edges) {
    degree[vmap.get(k6(e.a))]++;
    degree[vmap.get(k6(e.b))]++;
  }
  const F = A.edges.map((e, idx) => {
    const a = vmap.get(k6(e.a)), b = vmap.get(k6(e.b));
    return { edgeId: idx, deg_u: degree[a], deg_v: degree[b], F: 4 - degree[a] - degree[b] };
  });
  const hist = new Map();
  for (const f of F) hist.set(f.F, (hist.get(f.F) || 0) + 1);
  const histArr = [...hist.entries()].sort((a, b) => a[0] - b[0]);
  const mean = F.reduce((s, f) => s + f.F, 0) / F.length;
  return { F, degree, histogram: histArr, mean };
}

// ========================================================================
//  Closed walks  tr(A^k)  for primal 1-skeleton and dual graph
// ========================================================================
function matmul(A, B, n) {
  const C = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let i = 0; i < n; i++) for (let k = 0; k < n; k++) {
    if (A[i][k] === 0) continue;
    const aik = A[i][k];
    for (let j = 0; j < n; j++) C[i][j] += aik * B[k][j];
  }
  return C;
}
function trace(M) { let t = 0; for (let i = 0; i < M.length; i++) t += M[i][i]; return t; }

export function closedWalks(M, A, kMax = 8) {
  // ----- primal adjacency (63 x 63) -----
  const vmap = new Map(A.vertices.map((v, i) => [k6(v.pos), i]));
  const V = A.vertices.length;
  const Ap = Array.from({ length: V }, () => new Array(V).fill(0));
  for (const e of A.edges) {
    const a = vmap.get(k6(e.a)), b = vmap.get(k6(e.b));
    Ap[a][b] = 1; Ap[b][a] = 1;
  }
  // ----- dual adjacency (64 x 64) -----
  const T = M.tetrahedra.length;
  const Ad = Array.from({ length: T }, () => new Array(T).fill(0));
  for (const [i, j] of A.dualGraphs.edge.pairs) { Ad[i][j] = 1; Ad[j][i] = 1; }

  function spectrumTrace(Adj, n, kMax) {
    const out = [];
    let Ak = Adj.map((r) => r.slice());
    out.push({ k: 1, trace: trace(Ak) });
    for (let k = 2; k <= kMax; k++) {
      Ak = matmul(Ak, Adj, n);
      out.push({ k, trace: trace(Ak) });
    }
    return out;
  }
  return {
    primal: spectrumTrace(Ap, V, kMax),
    dual:   spectrumTrace(Ad, T, kMax),
  };
}

// ========================================================================
//  O_h orbits on dual edges (each orbit = a "type" of dual edge)
// ========================================================================
export function dualEdgeOrbits(M, A) {
  const G = buildOh();
  const dualPairs = A.dualGraphs.edge.pairs;
  const tetraKey = (t) => t.vertices.map(k6).sort().join("/");
  const tetraKeyToId = new Map(M.tetrahedra.map((t, i) => [tetraKey(t), i]));
  const dualKey = (i, j) => i < j ? `${i}|${j}` : `${j}|${i}`;
  const dualKeyToIdx = new Map(dualPairs.map((p, idx) => [dualKey(p[0], p[1]), idx]));
  // precompute image of each tetra under each g for speed
  const tetraImage = G.map((g) => M.tetrahedra.map((t) =>
    tetraKeyToId.get(t.vertices.map((v) => applyG(g, v)).map(k6).sort().join("/"))
  ));
  const visited = new Array(dualPairs.length).fill(false);
  const orbits = [];
  for (let idx = 0; idx < dualPairs.length; idx++) {
    if (visited[idx]) continue;
    const orbit = [idx];
    visited[idx] = true;
    for (let gi = 0; gi < G.length; gi++) {
      const [i, j] = dualPairs[idx];
      const ii = tetraImage[gi][i], jj = tetraImage[gi][j];
      const k = dualKeyToIdx.get(dualKey(ii, jj));
      if (k !== undefined && !visited[k]) { visited[k] = true; orbit.push(k); }
    }
    orbits.push(orbit);
  }
  return orbits.map((o, i) => ({ orbit: i, size: o.length, dualEdgeIds: o }));
}

// ========================================================================
//  Shell-ratio diagnostic vs phi
// ========================================================================
export function shellRatiosPhi(A) {
  const phi = (1 + Math.sqrt(5)) / 2;
  const tol = 0.05; // 5%
  const radii = A.shells.map((s) => ({ name: s.name, r: s.r }));
  // all (i,j) ratios with j > i and r_i > 0
  const out = [];
  for (let i = 0; i < radii.length; i++) for (let j = i + 1; j < radii.length; j++) {
    if (radii[i].r < 1e-9) continue;
    const ratio = radii[j].r / radii[i].r;
    const dev = Math.abs(ratio - phi) / phi;
    out.push({ from: radii[i].name, to: radii[j].name, ratio, deviation_from_phi: dev, within_5pct: dev < tol });
  }
  const within = out.filter((r) => r.within_5pct);
  return { phi, ratios: out, count_within_5pct: within.length, within };
}
