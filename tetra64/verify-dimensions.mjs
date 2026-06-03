#!/usr/bin/env node
// Cross-checks the structural numbers in the rector-intuition text against
// what the geometry actually produces.  Emits a human report and a JSON file.
//
//   node tetra64/verify-dimensions.mjs

import { writeFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";
import { analyze } from "./analysis.mjs";

const here = dirname(fileURLToPath(import.meta.url));
const A = analyze();

const CLAIMS = {
  "R^64_tetra (tetrahedra)":            { claim: 64,  got: A.M.tetrahedra.length },
  "R^63_vertex (unique tetra vertices)":{ claim: 63,  got: A.vertices.length },
  "chirality split 32+ / 32-":          { claim: "32 / 32", got: `${A.chiralityCounts.plus} / ${A.chiralityCounts.minus}` },
  "Shell 1 (cuboctahedron, VE) size":   { claim: 12,  got: A.shells.find((s) => Math.abs(s.r2 - 2) < 1e-9)?.count ?? 0 },
  "shells S0..S5 count":                { claim: 6,   got: A.shells.length },
  "R^240_edge (bulk edges, |e|=sqrt 2)":{ claim: 240, got: A.edges.length },
  "R^144_dual_edge (dual graph edges)": { claim: 144, gotByDef: {
    "face-sharing":   A.dualGraphs.face.count,
    "edge-sharing":   A.dualGraphs.edge.count,
    "vertex-sharing": A.dualGraphs.vertex.count,
    "tetra<->octa":   A.dualGraphs.tetraOcta.count,
  } },
  "R^81_dual_cycle (independent cycles)": { claim: 81, gotByDef: {
    "from face-sharing":   A.dualGraphs.face.cycles,
    "from edge-sharing":   A.dualGraphs.edge.cycles,
    "from vertex-sharing": A.dualGraphs.vertex.cycles,
  } },
  "ker(B^T)_64TG dim = 178":            { claim: 178, gotByDef: Object.fromEntries(
    A.boundary.candidates.map((c) => [c.name, c.kerColT])
  ) },
  "O_h orbits on 63 vertices = 6 shells": { claim: A.shells.length, got: A.orbits.vertices_Oh.length },
  "tetra orbits under O_h":               { claim: "(no explicit claim)", got: `${A.orbits.tetra_Oh.length} orbits, sizes ${JSON.stringify(A.orbits.tetra_Oh)}` },
  "tetra orbits under O":                 { claim: "(no explicit claim)", got: `${A.orbits.tetra_O.length} orbits, sizes ${JSON.stringify(A.orbits.tetra_O)}` },
};

const report = {
  generated_at: new Date().toISOString(),
  shells: A.shells.map((s) => ({ name: s.name, r: +s.r.toFixed(6), r2: s.r2, count: s.count })),
  chirality: A.chiralityCounts,
  edges: { count: A.edges.length, multiplicity_histogram: A.edgeMultHist },
  dual_graphs: {
    face_sharing:    { edges: A.dualGraphs.face.count,    cycles: A.dualGraphs.face.cycles },
    edge_sharing:    { edges: A.dualGraphs.edge.count,    cycles: A.dualGraphs.edge.cycles },
    vertex_sharing:  { edges: A.dualGraphs.vertex.count,  cycles: A.dualGraphs.vertex.cycles },
    tetra_octa:      { incidences: A.dualGraphs.tetraOcta.count },
  },
  antipodal: { vertex_pairs: A.antipodal.vertexPairs.length, self: A.antipodal.vertexSelf, tetra_pairs: A.antipodal.tetraPairs.length },
  boundary_matrices: A.boundary.candidates,
  orbits: {
    vertices_Oh: A.orbits.vertices_Oh,
    vertices_O:  A.orbits.vertices_O,
    tetra_Oh:    A.orbits.tetra_Oh,
    tetra_O:     A.orbits.tetra_O,
  },
  claims: CLAIMS,
};
writeFileSync(join(here, "tetra64-dimensions.json"), JSON.stringify(report, null, 2));

// ---- pretty print ----
const ok = (b) => (b ? "PASS" : "fail");
const line = (s = "") => console.log(s);
line("Verificacion de dimensiones del 64-TG");
line("=".repeat(60));
line(`vertices unicos:          ${A.vertices.length}`);
line(`tetraedros:               ${A.M.tetrahedra.length}`);
line(`bipartition quiral:       ${A.chiralityCounts.plus}+ / ${A.chiralityCounts.minus}-`);
line(`aristas (|e|=sqrt2):      ${A.edges.length}    histograma de multiplicidad: ${JSON.stringify(A.edgeMultHist)}`);
line(`cascarones radiales:`);
for (const s of A.shells) line(`  ${s.name}  r=${s.r.toFixed(4)}  count=${s.count}`);
line();
line("Grafos duales (candidatos):");
line(`  face-sharing:    aristas=${A.dualGraphs.face.count.toString().padStart(4)}  ciclos=${A.dualGraphs.face.cycles}`);
line(`  edge-sharing:    aristas=${A.dualGraphs.edge.count.toString().padStart(4)}  ciclos=${A.dualGraphs.edge.cycles}`);
line(`  vertex-sharing:  aristas=${A.dualGraphs.vertex.count.toString().padStart(4)}  ciclos=${A.dualGraphs.vertex.cycles}`);
line(`  tetra<->octa:    incidencias=${A.dualGraphs.tetraOcta.count}`);
line();
line("Matrices de borde candidatas (240 x 64):");
for (const c of A.boundary.candidates)
  line(`  ${c.name.padEnd(34)} rank=${c.rank}  ker(col^T)=${c.kerColT}`);
line();
line("Orbitas:");
line(`  O_h sobre vertices:   ${JSON.stringify(A.orbits.vertices_Oh)}  -> ${A.orbits.vertices_Oh.length} orbitas`);
line(`  O   sobre vertices:   ${JSON.stringify(A.orbits.vertices_O)}  -> ${A.orbits.vertices_O.length} orbitas`);
line(`  O_h sobre 64 tetra:   ${JSON.stringify(A.orbits.tetra_Oh)}  -> ${A.orbits.tetra_Oh.length} orbitas`);
line(`  O   sobre 64 tetra:   ${JSON.stringify(A.orbits.tetra_O)}  -> ${A.orbits.tetra_O.length} orbitas`);
line(`  antipodal vertex pairs: ${A.antipodal.vertexPairs.length}  (self-antipodal: ${A.antipodal.vertexSelf})`);
line(`  antipodal tetra pairs:  ${A.antipodal.tetraPairs.length}`);
line();
line("Comparacion con el texto rector (claims, no claims):");
line("-".repeat(60));
for (const [name, info] of Object.entries(CLAIMS)) {
  if (info.gotByDef) {
    line(`  ${name}`);
    line(`     claim: ${info.claim}`);
    for (const [def, val] of Object.entries(info.gotByDef)) {
      const match = String(val) === String(info.claim);
      line(`     ${match ? "[MATCH]" : "[    ]"} ${def}: ${val}`);
    }
  } else {
    const match = String(info.got) === String(info.claim);
    line(`  [${match ? "MATCH" : "diff "}] ${name}  claim=${info.claim}  got=${info.got}`);
  }
}
line();
line("Nota: las metricas no marcadas MATCH son honestas: depende de que");
line("definicion (gauge, signo, grafo dual exacto) se use; varias de las");
line("cifras del corpus rector son construcciones especificas (B con un signo");
line("particular, dual graph con cierta convencion) que no se reproducen aqui");
line("salvo que coincidan estructuralmente.");
line();
line(`escrito  tetra64/tetra64-dimensions.json`);
