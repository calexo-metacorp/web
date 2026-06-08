// 64-tetrahedron grid with a central cuboctahedral "vector equilibrium" nucleus.
//
// Construction (isotropic vector matrix / FCC lattice, integer coordinates):
//
//   * The 8 star tetrahedra (stellae octangulae) are centered on the octant
//     points  P = (sx, sy, sz),  s in {-1,+1}^3.
//   * Each stella is a central regular octahedron (vertices P +- e_i) wearing a
//     regular tetrahedron on every one of its 8 faces.  Those 8 spike tetrahedra
//     give  8 stellae x 8 = 64 tetrahedra,  all regular with edge sqrt(2).
//   * Every stella sends one spike straight at the origin.  The bases of those
//     spikes are the octahedron vertices nearest the origin, and across all 8
//     stellae they land on the 12 points  perm(0, +-1, +-1).  Those 12 points
//     are a cuboctahedron whose circumradius (sqrt 2) equals its edge (sqrt 2):
//     Buckminster Fuller's vector equilibrium.  It is the empty nuclear space
//     that balances the array (sum of its 12 vertex vectors = 0).
//
// This module is dependency-free and is shared, verbatim in spirit, with the
// in-browser viewer so the picture and the data are the same object.

export const EDGE = Math.SQRT2; // FCC nearest-neighbour distance, integer coords

const E = [
  [1, 0, 0],
  [0, 1, 0],
  [0, 0, 1],
];

const SIGNS = []; // the 8 octant sign-vectors
for (const sx of [-1, 1]) for (const sy of [-1, 1]) for (const sz of [-1, 1]) SIGNS.push([sx, sy, sz]);

const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const scale = (a, k) => [a[0] * k, a[1] * k, a[2] * k];
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dist2 = (a, b) => { const d = sub(a, b); return d[0] * d[0] + d[1] * d[1] + d[2] * d[2]; };
const len = (a) => Math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
const centroid = (pts) => scale(pts.reduce(add, [0, 0, 0]), 1 / pts.length);
const key = (p) => p.map((x) => (Object.is(x, -0) ? 0 : x).toFixed(6)).join(",");

export function build() {
  const stellaCenters = SIGNS.map((s) => [...s]);

  const tetrahedra = [];
  const octahedra = [];

  stellaCenters.forEach((P, si) => {
    // central octahedron of this stella: P +- e_i
    const octaVerts = [];
    for (const e of E) { octaVerts.push(add(P, e)); octaVerts.push(sub(P, e)); }
    octahedra.push({ stella: si, center: [...P], vertices: octaVerts });

    // one regular tetrahedron per octahedron face <-> per octant sign s
    SIGNS.forEach((s) => {
      const apex = add(P, s);
      const base = [add(P, scale(E[0], s[0])), add(P, scale(E[1], s[1])), add(P, scale(E[2], s[2]))];
      const vertices = [apex, ...base];
      tetrahedra.push({
        id: tetrahedra.length,
        stella: si,
        sign: [...s],
        vertices,
        centroid: centroid(vertices),
      });
    });
  });

  // --- the nuclear cuboctahedron (vector equilibrium): perm(0, +-1, +-1) ---
  const veVerts = [];
  for (const a of [-1, 1]) for (const b of [-1, 1]) {
    veVerts.push([0, a, b]);
    veVerts.push([a, 0, b]);
    veVerts.push([a, b, 0]);
  }
  // dedupe (each is generated once here, but be safe)
  const seen = new Set();
  const ve = [];
  for (const v of veVerts) { const k = key(v); if (!seen.has(k)) { seen.add(k); ve.push(v); } }

  // cuboctahedron edges: vertex pairs at distance sqrt(2)
  const edges = [];
  for (let i = 0; i < ve.length; i++) for (let j = i + 1; j < ve.length; j++) {
    if (Math.abs(dist2(ve[i], ve[j]) - EDGE * EDGE) < 1e-9) edges.push([i, j]);
  }
  // 8 triangular faces (one per octant) + 6 square faces (one per +-axis)
  const triFaces = SIGNS.map((s) => [
    ve.findIndex((v) => key(v) === key([0, s[1], s[2]])),
    ve.findIndex((v) => key(v) === key([s[0], 0, s[2]])),
    ve.findIndex((v) => key(v) === key([s[0], s[1], 0])),
  ]);
  const sqFaces = [];
  for (let axis = 0; axis < 3; axis++) for (const sgn of [-1, 1]) {
    const idx = ve.map((v, i) => [v, i]).filter(([v]) => v[axis] === sgn).map(([, i]) => i);
    // order the 4 indices into a proper square (sort by angle in the face plane)
    const c = centroid(idx.map((i) => ve[i]));
    const ref = sub(ve[idx[0]], c);
    const nrm = E[axis];
    const ang = (i) => {
      const w = sub(ve[i], c);
      const x = w[0] * ref[0] + w[1] * ref[1] + w[2] * ref[2];
      const cr = [ref[1] * w[2] - ref[2] * w[1], ref[2] * w[0] - ref[0] * w[2], ref[0] * w[1] - ref[1] * w[0]];
      const y = cr[0] * nrm[0] + cr[1] * nrm[1] + cr[2] * nrm[2];
      return Math.atan2(y, x);
    };
    idx.sort((p, q) => ang(p) - ang(q));
    sqFaces.push(idx);
  }

  return { stellaCenters, tetrahedra, octahedra, vectorEquilibrium: { vertices: ve, edges, triFaces, sqFaces } };
}

// Verification of every property the artifact claims.
export function verify(model) {
  const { tetrahedra, vectorEquilibrium: ve, stellaCenters } = model;
  const checks = {};
  const eq = (a, b, t = 1e-9) => Math.abs(a - b) < t;
  const vEq = (a, b, t = 1e-9) => eq(a[0], b[0], t) && eq(a[1], b[1], t) && eq(a[2], b[2], t);

  checks.tetra_count_is_64 = tetrahedra.length === 64;

  // every tetra is regular with edge sqrt(2)
  let allRegular = true;
  for (const t of tetrahedra) {
    const v = t.vertices;
    for (let i = 0; i < 4; i++) for (let j = i + 1; j < 4; j++) {
      if (!eq(dist2(v[i], v[j]), EDGE * EDGE)) allRegular = false;
    }
  }
  checks.all_tetrahedra_regular_edge_sqrt2 = allRegular;

  // vector equilibrium: 12 vertices, radius == edge
  checks.ve_vertex_count_is_12 = ve.vertices.length === 12;
  checks.ve_radius_equals_edge =
    ve.vertices.every((v) => eq(len(v), EDGE)) && eq(EDGE, EDGE);
  checks.ve_edge_count_is_24 = ve.edges.length === 24;
  checks.ve_faces_8_tri_6_sq = ve.triFaces.length === 8 && ve.sqFaces.length === 6;

  // balance: vectors sum to zero
  checks.ve_vertices_sum_to_zero = vEq(ve.vertices.reduce(add, [0, 0, 0]), [0, 0, 0]);
  checks.stella_centers_sum_to_zero = vEq(stellaCenters.reduce(add, [0, 0, 0]), [0, 0, 0]);
  checks.array_centroid_is_origin = vEq(centroid(tetrahedra.map((t) => t.centroid)), [0, 0, 0], 1e-9);

  // each VE vertex is shared by exactly two stellae (the binding of the array)
  const shareCounts = ve.vertices.map((v) =>
    stellaCenters.filter((P) => eq(dist2(P, v), 1)).length
  );
  checks.each_ve_vertex_shared_by_two_stellae = shareCounts.every((c) => c === 2);

  return checks;
}

export const metrics = (model) => {
  const tetVol = (Math.pow(EDGE, 3)) / (6 * Math.SQRT2); // regular tetra edge sqrt2 = 1/3
  const veVol = (5 / 3) * Math.SQRT2 * Math.pow(EDGE, 3); // cuboctahedron edge sqrt2 = 20/3
  return {
    edge: EDGE,
    tetrahedron_volume: tetVol,
    total_tetrahedra_volume: tetVol * model.tetrahedra.length,
    nucleus_cuboctahedron_volume: veVol,
    symmetry_group: "Oh (full octahedral)",
    symmetry_order: 48,
  };
};
