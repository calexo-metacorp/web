#!/usr/bin/env node
// Phase 2 verification - cross-check the lab atlas v1.0 foundations table
// against what the geometry produces from first principles.
//
//   node tetra64/verify-foundations.mjs

import { writeFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";
import { build } from "./geometry.mjs";
import { analyze } from "./analysis.mjs";
import { trianglesAndBetti, bopRank, homGdim, dualCycleSpaceIrreps, shellRatiosPhi,
         formanRicci, closedWalks, dualEdgeOrbits } from "./analysis-phase2.mjs";
import { OH_IRREP_NAMES, OH_IRREP_DIMS } from "./oh-character-table.mjs";

const here = dirname(fileURLToPath(import.meta.url));
const M = build();
const A = analyze(M);

console.log("Phase 2 foundations - lab atlas v1.0 cross-check");
console.log("=".repeat(64));

// --- baseline preserved -----------------------------------------------
const baseline = {
  "tetrahedra |T|":       { claim: 64,  got: M.tetrahedra.length },
  "vertices |V|":         { claim: 63,  got: A.vertices.length },
  "edges |E|":            { claim: 240, got: A.edges.length },
  "dual edges (E-sharing)": { claim: 144, got: A.dualGraphs.edge.count },
  "dual cycles":          { claim: 81,  got: A.dualGraphs.edge.cycles },
  "chirality 32/32":      { claim: "32/32", got: `${A.chiralityCounts.plus}/${A.chiralityCounts.minus}` },
  "shells":               { claim: 6,   got: A.shells.length },
};
const okIcon = (b) => b ? "PASS" : "FAIL";
console.log("\n-- baseline (Phase 1 unchanged) --");
for (const [name, v] of Object.entries(baseline)) {
  const eq = String(v.got) === String(v.claim);
  console.log(`  [${okIcon(eq)}] ${name.padEnd(28)} claim=${v.claim} got=${v.got}`);
}

// --- triangles + Betti ------------------------------------------------
console.log("\n-- triangles + Betti numbers --");
const TB = trianglesAndBetti(M, A);
console.log(`  triangles |F| (= 4 faces per tetra, 64*4):  ${TB.F}    claim 256  ${okIcon(TB.F === 256)}`);
console.log(`  3-cliques in 1-skeleton (sanity reference): ${TB.cliques}`);
console.log(`  rank d_1 (V x E):  ${TB.ranks.r1}    rank d_2 (E x F):  ${TB.ranks.r2}    rank d_3 (F x T):  ${TB.ranks.r3}`);
console.log(`  beta = (${TB.beta.join(", ")})    claim (1, 0, 14, 0)  ${okIcon(TB.beta[0]===1 && TB.beta[1]===0 && TB.beta[2]===14 && TB.beta[3]===0)}`);

// --- B_op rank --------------------------------------------------------
console.log("\n-- B_op rank (O_h-equivariant vertex-edge boundary) --");
const BO = bopRank(M, A);
console.log(`  ${BO.definition}`);
console.log(`  rank = ${BO.rank}   (matrix ${BO.rows} x ${BO.cols})   claim 12  ${okIcon(BO.rank === 12)}`);

// --- Hom_G dimension --------------------------------------------------
console.log("\n-- dim Hom_G(V_edge240, V_dual144 (x) chi) via character formula --");
const HG = homGdim(M, A);
console.log(`  variants tried (each = (1/|G|) sum_g chi_edge(g) * chi_dual(g) * chi(g)):`);
for (const [k, v] of Object.entries(HG.dims)) {
  const hit = Math.round(v) === HG.claim;
  console.log(`     ${hit ? "[MATCH 744]" : "[         ]"} ${k.padEnd(28)} = ${v}`);
}
if (HG.match) console.log(`  -> match found: ${HG.match}`);
else console.log(`  -> no variant matches 744; corpus's V_edge/V_dual rep signing differs from these conventions.`);

// --- cycle space irrep decomposition ----------------------------------
console.log("\n-- H_1(dual graph) ~ R^81 irrep decomposition under O_h --");
const CYC = dualCycleSpaceIrreps(M, A);
const total = CYC.mults.reduce((acc, m, i) => acc + m * OH_IRREP_DIMS[i], 0);
console.log(`  characters per class: ${CYC.chars.map((c) => c.toFixed(2)).join("  ")}`);
console.log(`  decomposition: ${CYC.formatted}`);
console.log(`  sum a_i * dim_i = ${total.toFixed(3)}    must equal 81  ${okIcon(Math.abs(total - 81) < 1e-6)}`);
console.log(`  multiplicities (rounded):`);
for (let i = 0; i < OH_IRREP_NAMES.length; i++) {
  const m = CYC.mults[i];
  if (Math.abs(m) < 1e-9) continue;
  console.log(`     ${OH_IRREP_NAMES[i].padEnd(4)} (dim ${OH_IRREP_DIMS[i]}):  ${m.toFixed(4)}`);
}

// --- shell ratios vs phi ----------------------------------------------
console.log("\n-- shell-ratio diagnostic vs phi (5% bound) --");
const SR = shellRatiosPhi(A);
console.log(`  phi = ${SR.phi.toFixed(6)}`);
console.log(`  ratios within 5% of phi: ${SR.count_within_5pct}    atlas diagnostic claim 5`);
if (SR.count_within_5pct > 0) {
  for (const r of SR.within) console.log(`     ${r.from} -> ${r.to}: ${r.ratio.toFixed(4)}   |dev|=${(r.deviation_from_phi*100).toFixed(2)}%`);
} else {
  console.log("     (none of the shell-radius ratios fall within 5% of phi under the geometric shell definition used here)");
  console.log("     atlas note: 'NO claim derivacion canonical de phi'; Q24 SACRED anti-heterodoxy preserved.");
}

// --- bipartite chirality structure of the dual graph -------------------
console.log("\n-- bipartite chirality structure of the dual graph --");
const chir = new Map(A.chirality.map((c) => [c.id, c.chirality]));
let edgesPP = 0, edgesMM = 0, edgesPM = 0;
for (const [i, j] of A.dualGraphs.edge.pairs) {
  const a = chir.get(i), b = chir.get(j);
  if (a > 0 && b > 0) edgesPP++;
  else if (a < 0 && b < 0) edgesMM++;
  else edgesPM++;
}
const bipartiteOk = edgesPP === 0 && edgesMM === 0;
console.log(`  +/+ dual edges: ${edgesPP}    -/- : ${edgesMM}    +/- : ${edgesPM}`);
console.log(`  bipartite (every dual edge crosses chirality)?  ${okIcon(bipartiteOk)}    atlas [OBS-BIPARTITE-STRUCTURE]`);

// --- Forman-Ricci curvature per edge (graph version) ------------------
console.log("\n-- Forman-Ricci curvature per edge (graph version) --");
const FR = formanRicci(M, A);
console.log(`  formula: F(e=(u,v)) = 4 - deg(u) - deg(v)`);
console.log(`  mean F = ${FR.mean.toFixed(4)}    histogram of values:`);
for (const [val, count] of FR.histogram) {
  const bar = "#".repeat(Math.min(60, Math.round(count/2)));
  console.log(`     F = ${String(val).padStart(4)}  count=${String(count).padStart(3)}  ${bar}`);
}
console.log(`  degree distribution per shell:`);
for (const s of A.shells) {
  const degs = s.vertexIds.map((id) => FR.degree[id]);
  const minD = Math.min(...degs), maxD = Math.max(...degs);
  console.log(`     ${s.name}  r=${s.r.toFixed(3)}  count=${s.count}   deg in [${minD}, ${maxD}]`);
}

// --- closed walks  tr(A^k)  -------------------------------------------
console.log("\n-- closed walks  tr(A^k)  on primal (63x63) and dual (64x64) --");
const CW = closedWalks(M, A, 8);
const fmt = (n) => String(n).padStart(11);
console.log("        k= " + CW.primal.map((r) => String(r.k).padStart(11)).join(" "));
console.log("  primal:  " + CW.primal.map((r) => fmt(r.trace)).join(" "));
console.log("    dual:  " + CW.dual.map((r) => fmt(r.trace)).join(" "));
const tr2 = CW.primal[1].trace, tr3 = CW.primal[2].trace, atlasE = 2 * A.edges.length;
const triangles_count = tr3 / 6;
console.log(`  tr(A^2) primal = 2|E| ?      ${tr2} vs ${atlasE}            ${okIcon(tr2 === atlasE)}`);
console.log(`  tr(A^3) primal / 6 = #triangles? ${tr3}/6 = ${triangles_count} vs 256 ${okIcon(triangles_count === 256)}`);
const oddDual = CW.dual.filter((r) => r.k % 2 === 1).every((r) => r.trace === 0);
console.log(`  dual graph bipartite signature (tr(A^k)=0 for all odd k)?  ${okIcon(oddDual)}`);

// --- O_h orbits on dual edges -----------------------------------------
console.log("\n-- O_h orbits on the 144 dual edges --");
const DO = dualEdgeOrbits(M, A);
console.log(`  ${DO.length} orbit(s) of sizes ${JSON.stringify(DO.map((o) => o.size))}    (total ${DO.reduce((s, o) => s + o.size, 0)})`);

// --- bipartite check (already above) ----------------------------------

// ---- export JSON -----------------------------------------------------
const out = {
  generated_at: new Date().toISOString(),
  baseline: Object.fromEntries(Object.entries(baseline).map(([k, v]) => [k, { ...v, pass: String(v.got) === String(v.claim) }])),
  triangles_and_betti: TB,
  bop: BO,
  hom_g: HG,
  cycle_space_irreps_H1_dual: { chars: CYC.chars, multiplicities: Object.fromEntries(OH_IRREP_NAMES.map((n, i) => [n, CYC.mults[i]])), formatted: CYC.formatted },
  shell_ratios_phi: SR,
  forman_ricci: { mean: FR.mean, histogram: FR.histogram },
  closed_walks: CW,
  dual_edge_oh_orbits: { count: DO.length, sizes: DO.map((o) => o.size) },
  bipartite_dual: { pp: edgesPP, mm: edgesMM, pm: edgesPM, is_bipartite: bipartiteOk },
};
writeFileSync(join(here, "tetra64-foundations.json"), JSON.stringify(out, null, 2));
console.log("\nwrote tetra64/tetra64-foundations.json");
