// Structural analysis of the 64-tetrahedron grid:
//   - chirality bipartition (32+ / 32- via signed volume)
//   - radial shells S0..S5 of the 63 unique vertices
//   - bulk edges (length sqrt(2)) and their tetrahedron multiplicities
//   - several dual-graph candidates (face / edge / vertex sharing,
//     and the tetrahedron <-> central-octahedron bipartite incidence)
//   - antipodal pairs under spatial inversion P
//   - candidate boundary matrices and their R^240 -> R^64 kernel dimensions
//   - Oh and O (rotation-only) orbits on tetrahedra and on vertices
//
// All numbers come from the geometry; nothing is asserted physically.

import { build } from "./geometry.mjs";

const eps = 1e-9;
const k6 = (p) => p.map((x) => (Object.is(x, -0) ? 0 : +x).toFixed(6)).join(",");
const add = (a, b) => [a[0]+b[0], a[1]+b[1], a[2]+b[2]];
const sub = (a, b) => [a[0]-b[0], a[1]-b[1], a[2]-b[2]];
const dot = (a, b) => a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
const dist2 = (a, b) => { const d=sub(a,b); return dot(d,d); };
const det3 = (a, b, c) =>
  a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0]);

// --- the 48 elements of O_h as signed permutations of (x,y,z) ----------
function buildOh() {
  const perms = [
    [0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0],
  ];
  const sgn = (p) => { // sign of permutation
    let s = 1; for (let i = 0; i < 3; i++) for (let j = i+1; j < 3; j++) if (p[i] > p[j]) s = -s; return s;
  };
  const G = []; // each element: {pi, s, det = sgn(pi)*prod(s)}
  for (const pi of perms) for (const sx of [-1,1]) for (const sy of [-1,1]) for (const sz of [-1,1])
    G.push({ pi, s: [sx,sy,sz], det: sgn(pi)*sx*sy*sz });
  return G;
}
function applyG(g, v) {
  return [g.s[0]*v[g.pi[0]], g.s[1]*v[g.pi[1]], g.s[2]*v[g.pi[2]]];
}

// ---- Gaussian elimination over R (partial pivoting) -------------------
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

// ---- main analysis ----------------------------------------------------
export function analyze(M = build()) {
  // ===== 1. CHIRALITY (32+ / 32-) =====
  // chirality = sign of det(b1-apex, b2-apex, b3-apex);
  // closed form: chir(s) = -prod(s_x,s_y,s_z), independent of stella centre P
  const chirality = M.tetrahedra.map((t) => {
    const [a, b, c, d] = t.vertices;
    const sv = det3(sub(b, a), sub(c, a), sub(d, a));
    return { id: t.id, stella: t.stella, sign: t.sign, signedVolume: sv, chirality: Math.sign(sv) };
  });
  const chiralityCounts = { plus: 0, minus: 0 };
  for (const c of chirality) { if (c.chirality > 0) chiralityCounts.plus++; else chiralityCounts.minus++; }

  // ===== 2. UNIQUE VERTICES & RADIAL SHELLS =====
  const vmap = new Map();
  for (const t of M.tetrahedra) for (const v of t.vertices) {
    const k = k6(v);
    if (!vmap.has(k)) vmap.set(k, { pos: v, tetra: [] });
    vmap.get(k).tetra.push(t.id);
  }
  const vertices = [...vmap.values()].map((o, i) => ({ id: i, pos: o.pos, tetra: o.tetra, r2: dot(o.pos, o.pos) }));
  const shellMap = new Map();
  for (const v of vertices) {
    const k = v.r2.toFixed(6);
    if (!shellMap.has(k)) shellMap.set(k, []);
    shellMap.get(k).push(v.id);
  }
  const shells = [...shellMap.entries()]
    .map(([r2, ids]) => ({ r2: +r2, r: Math.sqrt(+r2), count: ids.length, vertexIds: ids }))
    .sort((a, b) => a.r - b.r);
  shells.forEach((s, i) => (s.name = "S" + i));

  // ===== 3. BULK EDGES (length sqrt(2)) =====
  const eMap = new Map();
  for (const t of M.tetrahedra) for (let i = 0; i < 4; i++) for (let j = i + 1; j < 4; j++) {
    const a = k6(t.vertices[i]), b = k6(t.vertices[j]);
    const key = a < b ? a + "|" + b : b + "|" + a;
    if (!eMap.has(key)) eMap.set(key, { a: t.vertices[i], b: t.vertices[j], tetra: [] });
    eMap.get(key).tetra.push(t.id);
  }
  const edges = [...eMap.values()].map((e, id) => ({ id, ...e, mult: e.tetra.length }));
  const edgeMultHist = {};
  for (const e of edges) edgeMultHist[e.mult] = (edgeMultHist[e.mult] || 0) + 1;

  // ===== 4. DUAL-GRAPH CANDIDATES =====
  // for each tetra: set of vertex keys
  const tVerts = M.tetrahedra.map((t) => new Set(t.vertices.map(k6)));
  // for each tetra: set of edge keys
  const tEdgeKeys = M.tetrahedra.map((t) => {
    const s = new Set();
    for (let i=0;i<4;i++) for (let j=i+1;j<4;j++) {
      const a=k6(t.vertices[i]),b=k6(t.vertices[j]); s.add(a<b?a+"|"+b:b+"|"+a);
    }
    return s;
  });
  // for each tetra: set of face keys
  const tFaceKeys = M.tetrahedra.map((t) => {
    const s = new Set();
    for (let i=0;i<4;i++) for (let j=i+1;j<4;j++) for (let k=j+1;k<4;k++) {
      const ks=[k6(t.vertices[i]),k6(t.vertices[j]),k6(t.vertices[k])].sort().join("|");
      s.add(ks);
    }
    return s;
  });
  const intersect = (A, B) => { let n = 0; for (const x of A) if (B.has(x)) n++; return n; };
  const candidate = (label, predicate) => {
    const pairs = [];
    for (let i = 0; i < 64; i++) for (let j = i + 1; j < 64; j++) if (predicate(i, j)) pairs.push([i, j]);
    return { label, count: pairs.length, pairs };
  };
  const dualFace   = candidate("face-sharing",   (i, j) => intersect(tFaceKeys[i],   tFaceKeys[j])   >= 1);
  const dualEdge   = candidate("edge-sharing",   (i, j) => intersect(tEdgeKeys[i],   tEdgeKeys[j])   >= 1);
  const dualVertex = candidate("vertex-sharing", (i, j) => intersect(tVerts[i],      tVerts[j])      >= 1);
  // tetra <-> octahedron bipartite incidence: each tetra's base face equals an octahedron face
  const octFaceSets = M.octahedra.map((o) => {
    const verts = o.vertices.map(k6);
    const s = new Set();
    for (let i=0;i<6;i++) for (let j=i+1;j<6;j++) for (let k=j+1;k<6;k++)
      s.add([verts[i],verts[j],verts[k]].sort().join("|"));
    return s;
  });
  let tetraOctaEdges = 0;
  for (let i=0;i<64;i++) for (let oc=0;oc<8;oc++)
    if (intersect(tFaceKeys[i], octFaceSets[oc]) >= 1) tetraOctaEdges++;
  const tetraOcta = { label: "tetra<->octahedron bipartite (faces)", count: tetraOctaEdges };

  // independent cycles per candidate (assuming connected): E - V + 1
  const cyclesOf = (E, V = 64) => E - V + 1;
  for (const d of [dualFace, dualEdge, dualVertex]) d.cycles = cyclesOf(d.count);

  // ===== 5. ANTIPODAL PAIRS (inversion P) =====
  const vkey = new Map(vertices.map((v) => [k6(v.pos), v.id]));
  const antipodal = [];
  for (const v of vertices) {
    const j = vkey.get(k6(v.pos.map((x) => -x)));
    if (j !== undefined && j > v.id) antipodal.push([v.id, j]);
  }
  const selfAntipodal = vertices.filter((v) => vkey.get(k6(v.pos.map((x) => -x))) === v.id).length;
  // tetra antipodal pairs: t mapped to inversion of all its vertices
  const tetraKey = M.tetrahedra.map((t) => t.vertices.map(k6).sort().join("/"));
  const tetraInvKey = M.tetrahedra.map((t) => t.vertices.map((v)=>k6(v.map((x)=>-x))).sort().join("/"));
  const keyToId = new Map(tetraKey.map((k, i) => [k, i]));
  const antipodalTetra = [];
  for (let i = 0; i < 64; i++) {
    const j = keyToId.get(tetraInvKey[i]);
    if (j !== undefined && j > i) antipodalTetra.push([i, j]);
  }

  // ===== 6. BOUNDARY MATRICES & ker(B^T) =====
  // candidate B1: unsigned edges x tetra incidence (240 rows, 64 cols)
  const E = edges.length, T = M.tetrahedra.length;
  const edgeKeyToId = new Map();
  [...eMap.keys()].forEach((k, i) => edgeKeyToId.set(k, i));
  const B1 = Array.from({ length: E }, () => new Array(T).fill(0));
  M.tetrahedra.forEach((t, ti) => {
    for (let i=0;i<4;i++) for (let j=i+1;j<4;j++) {
      const a=k6(t.vertices[i]),b=k6(t.vertices[j]);
      const key = a<b?a+"|"+b:b+"|"+a;
      const ei = edgeKeyToId.get(key);
      B1[ei][ti] = 1;
    }
  });
  // candidate B2: same with column sign = chirality of tetra
  const B2 = B1.map((row) => row.map((x, ti) => x * chirality[ti].chirality));
  const rankB1 = rankReal(E, T, B1);
  const rankB2 = rankReal(E, T, B2);
  const kerB1T = E - rankB1, kerB2T = E - rankB2;
  // candidate B3: standard graph boundary  d1 : edges -> vertices,  +1 at head, -1 at tail.
  // For a connected graph with |V| vertices, rank(d1) = |V| - 1, so ker(d1^T) = |E| - |V| + 1
  // = the cycle space (1st homology of the bulk 1-skeleton).
  const V = vertices.length;
  const vkeyMap = new Map(vertices.map((v) => [k6(v.pos), v.id]));
  const B3 = Array.from({ length: V }, () => new Array(E).fill(0));
  edges.forEach((e, ei) => {
    const ka = k6(e.a), kb = k6(e.b);
    const va = vkeyMap.get(ka), vb = vkeyMap.get(kb);
    // canonical orientation: lower key -> higher key
    const [tail, head] = ka < kb ? [va, vb] : [vb, va];
    B3[head][ei] += 1;
    B3[tail][ei] -= 1;
  });
  const rankB3 = rankReal(V, E, B3);
  const kerB3T = E - rankB3; // dimension of cycle space

  // ===== 7. O_h and O ORBITS =====
  const G = buildOh();          // 48 elements
  const Orot = G.filter((g) => g.det > 0);   // 24 rotations
  const orbits = (group, items, equal) => {
    const seen = new Array(items.length).fill(false);
    const out = [];
    for (let i = 0; i < items.length; i++) {
      if (seen[i]) continue;
      const orb = [i];
      seen[i] = true;
      for (const g of group) {
        const gx = g.act(items[i]);
        for (let j = 0; j < items.length; j++) if (!seen[j] && equal(items[j], gx)) { seen[j] = true; orb.push(j); }
      }
      out.push(orb);
    }
    return out;
  };
  const actVertex = (v) => (g) => applyG(g, v);
  const vertexOrbitsOh = orbits(G.map((g) => ({ act: (v) => applyG(g, v) })), vertices.map((v) => v.pos), (a, b) => k6(a) === k6(b));
  const vertexOrbitsO  = orbits(Orot.map((g) => ({ act: (v) => applyG(g, v) })), vertices.map((v) => v.pos), (a, b) => k6(a) === k6(b));
  const tetraKeySet = (verts) => verts.map(k6).sort().join("/");
  const tetraOrbitsOh = orbits(G.map((g) => ({ act: (verts) => verts.map((v) => applyG(g, v)) })),
                               M.tetrahedra.map((t) => t.vertices), (a, b) => tetraKeySet(a) === tetraKeySet(b));
  const tetraOrbitsO  = orbits(Orot.map((g) => ({ act: (verts) => verts.map((v) => applyG(g, v)) })),
                               M.tetrahedra.map((t) => t.vertices), (a, b) => tetraKeySet(a) === tetraKeySet(b));

  return {
    M,
    chirality, chiralityCounts,
    vertices, shells,
    edges, edgeMultHist,
    dualGraphs: { face: dualFace, edge: dualEdge, vertex: dualVertex, tetraOcta },
    antipodal: { vertexPairs: antipodal, vertexSelf: selfAntipodal, tetraPairs: antipodalTetra },
    boundary: {
      candidates: [
        { name: "B1: unsigned edges x tetra",     rows: E, cols: T, rank: rankB1, kerColT: E - rankB1 },
        { name: "B2: signed-by-chirality",         rows: E, cols: T, rank: rankB2, kerColT: E - rankB2 },
        { name: "B3: graph boundary d1 (V x E)",  rows: V, cols: E, rank: rankB3, kerColT: E - rankB3 },
      ],
    },
    orbits: {
      vertices_Oh: vertexOrbitsOh.map((o) => o.length),
      vertices_O:  vertexOrbitsO.map((o) => o.length),
      tetra_Oh:    tetraOrbitsOh.map((o) => o.length),
      tetra_O:     tetraOrbitsO.map((o) => o.length),
    },
  };
}
