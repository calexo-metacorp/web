#!/usr/bin/env node
// Builds the 64-tetrahedron / cuboctahedral-vector-equilibrium artifact,
// verifies all of its geometric claims, and writes tetra64.json.
//
//   node tetra64/generate.mjs

import { writeFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";
import { build, verify, metrics, EDGE } from "./geometry.mjs";

const here = dirname(fileURLToPath(import.meta.url));
const model = build();
const checks = verify(model);
const m = metrics(model);

const allPass = Object.values(checks).every(Boolean);

const artifact = {
  name: "64-tetrahedron grid with a cuboctahedral vector-equilibrium nucleus",
  description:
    "Eight star tetrahedra (stellae octangulae) on the octant points (+-1,+-1,+-1) " +
    "of the isotropic vector matrix. Each contributes 8 regular tetrahedra (64 total). " +
    "Their innermost vertices coincide on a cuboctahedron whose circumradius equals its " +
    "edge - Fuller's vector equilibrium - the balanced nuclear space of the array.",
  lattice: "FCC / isotropic vector matrix, integer coordinates, edge = sqrt(2)",
  edge_length: EDGE,
  metrics: m,
  stella_centers: model.stellaCenters,
  tetrahedra: model.tetrahedra,
  central_octahedra: model.octahedra,
  vector_equilibrium: {
    vertices: model.vectorEquilibrium.vertices,
    edges: model.vectorEquilibrium.edges,
    triangular_faces: model.vectorEquilibrium.triFaces,
    square_faces: model.vectorEquilibrium.sqFaces,
    circumradius: EDGE,
    edge: EDGE,
    radius_equals_edge: true,
  },
  verification: checks,
  verification_passed: allPass,
};

writeFileSync(join(here, "tetra64.json"), JSON.stringify(artifact, null, 2));

// ---- human-readable report ----
const ok = (b) => (b ? "PASS" : "FAIL");
console.log("64-tetrahedron grid / cuboctahedral vector equilibrium");
console.log("=".repeat(56));
console.log(`edge length .................. ${EDGE.toFixed(6)} (= sqrt 2)`);
console.log(`tetrahedra ................... ${model.tetrahedra.length}`);
console.log(`star tetrahedra (stellae) .... ${model.stellaCenters.length}  (8 x 8 = 64)`);
console.log(`central octahedra ............ ${model.octahedra.length}`);
console.log(`nucleus ...................... cuboctahedron, ${model.vectorEquilibrium.vertices.length} vertices`);
console.log(`  circumradius / edge ........ ${EDGE.toFixed(6)} / ${EDGE.toFixed(6)}  -> equilibrium`);
console.log(`  tetra volume (each) ........ ${m.tetrahedron_volume.toFixed(6)}  (= 1/3)`);
console.log(`  64-tetra volume ............ ${m.total_tetrahedra_volume.toFixed(6)}  (= 64/3)`);
console.log(`  nucleus volume ............. ${m.nucleus_cuboctahedron_volume.toFixed(6)}  (= 20/3)`);
console.log(`  symmetry ................... ${m.symmetry_group}, order ${m.symmetry_order}`);
console.log("-".repeat(56));
for (const [k, v] of Object.entries(checks)) console.log(`  [${ok(v)}] ${k}`);
console.log("-".repeat(56));
console.log(`ALL CHECKS: ${ok(allPass)}`);
console.log(`wrote ${join("tetra64", "tetra64.json")}`);

process.exit(allPass ? 0 : 1);
