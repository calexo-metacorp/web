#!/usr/bin/env node
// Curvature snapshot: 240 bulk edges coloured by Forman-Ricci F = 4-deg(u)-deg(v),
// 144 dual edges coloured by O_h orbit, nucleus + radials kept for orientation.
//   node tetra64/render-curvature-svg.mjs

import { writeFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";
import { build } from "./geometry.mjs";
import { analyze } from "./analysis.mjs";
import { formanRicci, dualEdgeOrbits } from "./analysis-phase2.mjs";

const here = dirname(fileURLToPath(import.meta.url));
const M = build();
const A = analyze(M);
const FR = formanRicci(M, A);
const DO = dualEdgeOrbits(M, A);
const formanVals = [...new Set(FR.F.map((f) => f.F))].sort((a, b) => a - b);
const formanIdx = new Map(formanVals.map((v, i) => [v, i]));
const dualOrbitOf = new Array(A.dualGraphs.edge.pairs.length).fill(-1);
DO.forEach((o, i) => o.dualEdgeIds.forEach((e) => (dualOrbitOf[e] = i)));

const W = 1180, H = 800, PADX = 50, TOP = 80, BOT = 80;
const CX = W / 2, CY = TOP + (H - TOP - BOT) / 2;
const yaw = 0.62, pitch = -0.46, F = 7;
const sub = (a, b) => [a[0]-b[0], a[1]-b[1], a[2]-b[2]];
const cross = (a, b) => [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
const llen = (a) => Math.hypot(a[0], a[1], a[2]);
const ddot = (a, b) => a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
const cen = (p) => { const s=p.reduce((a,b)=>[a[0]+b[0],a[1]+b[1],a[2]+b[2]],[0,0,0]); return [s[0]/p.length,s[1]/p.length,s[2]/p.length]; };
function rot(p){
  const cy=Math.cos(yaw), sy=Math.sin(yaw), cx=Math.cos(pitch), sx=Math.sin(pitch);
  const x=p[0]*cy+p[2]*sy, z=-p[0]*sy+p[2]*cy, y=p[1];
  return [x, y*cx-z*sx, y*sx+z*cx];
}
const ALL = [...M.tetrahedra.flatMap(t=>t.vertices), ...M.vectorEquilibrium.vertices];
const SCALE = (() => {
  let mu=0, mv=0;
  for(const p of ALL){ const r=rot(p), k=F/(F-r[2]); mu=Math.max(mu,Math.abs(r[0]*k)); mv=Math.max(mv,Math.abs(r[1]*k)); }
  return Math.min((W/2-PADX)/mu, ((H-TOP-BOT)/2)/mv);
})();
function proj(p){
  const r=rot(p), persp=F/(F-r[2]), s=SCALE*persp;
  return { x:CX+r[0]*s, y:CY-r[1]*s, z:r[2] };
}

const F_PAL = ["#3d4ea6","#4a87c8","#46c6c6","#8fd07e","#f0b450","#ff7a59"];
const O_PAL = ["#6fdc8c","#c894ff","#ff7a59","#ffd86b","#5ec0ff","#ff6fb5"];

const items = [];

// 240 bulk edges coloured by Forman-Ricci
A.edges.forEach((e, ei) => {
  const pa = proj(e.a), pb = proj(e.b), z = (pa.z + pb.z) / 2 - 0.01;
  const col = F_PAL[formanIdx.get(FR.F[ei].F)];
  items.push({ z, s: `<line x1="${pa.x.toFixed(1)}" y1="${pa.y.toFixed(1)}" x2="${pb.x.toFixed(1)}" y2="${pb.y.toFixed(1)}" stroke="${col}" stroke-width="1.1" stroke-opacity="0.85"/>` });
});

// 144 dual edges coloured by O_h orbit
A.dualGraphs.edge.pairs.forEach(([i, j], ei) => {
  const a = M.tetrahedra[i].centroid, b = M.tetrahedra[j].centroid;
  const pa = proj(a), pb = proj(b), z = (pa.z + pb.z) / 2;
  const col = O_PAL[dualOrbitOf[ei] % O_PAL.length];
  items.push({ z, s: `<line x1="${pa.x.toFixed(1)}" y1="${pa.y.toFixed(1)}" x2="${pb.x.toFixed(1)}" y2="${pb.y.toFixed(1)}" stroke="${col}" stroke-width="0.9" stroke-opacity="0.7"/>` });
});

// nucleus VE
const VEv = M.vectorEquilibrium.vertices, Pve = VEv.map(proj), Rve = VEv.map(rot);
function shade(rgb, v3, emis) {
  const n = cross(sub(v3[1], v3[0]), sub(v3[2], v3[0])), m = llen(n) || 1;
  const lit = Math.abs(ddot([n[0]/m, n[1]/m, n[2]/m], [0.35, 0.55, 0.9]) / llen([0.35, 0.55, 0.9])), k = emis + (1 - emis) * lit;
  return `rgb(${Math.round(rgb[0]*k)},${Math.round(rgb[1]*k)},${Math.round(rgb[2]*k)})`;
}
const faces = [...M.vectorEquilibrium.triFaces.map(f => ({ f, sq: false })), ...M.vectorEquilibrium.sqFaces.map(f => ({ f, sq: true }))];
for (const { f, sq } of faces) {
  const v3 = f.map(i => Rve[i]), z = cen(v3)[2] + 0.002;
  const pts = f.map(i => `${Pve[i].x.toFixed(1)},${Pve[i].y.toFixed(1)}`).join(" ");
  items.push({ z, s: `<polygon points="${pts}" fill="${shade([245,196,81], v3, 0.42)}" fill-opacity="${sq?0.35:0.5}" stroke="#f5c451" stroke-width="1.2"/>` });
}
const o = proj([0, 0, 0]);
items.push({ z: 99, s: `<circle cx="${o.x.toFixed(1)}" cy="${o.y.toFixed(1)}" r="3.4" fill="#ffffff"/>` });

items.sort((a, b) => a.z - b.z);

// legends
const Fhist = new Map();
for (const f of FR.F) Fhist.set(f.F, (Fhist.get(f.F) || 0) + 1);
const Fleg = formanVals.map((v, i) => {
  const y = TOP + 4 + i * 18;
  return `<rect x="${W-260}" y="${y-6}" width="14" height="3" fill="${F_PAL[i]}"/><text x="${W-240}" y="${y}" fill="#cdd3e0" font-size="11">F=${v}  (${Fhist.get(v)} aristas)</text>`;
}).join("\n");
const Oleg = DO.map((o, i) => {
  const y = TOP + 4 + (i + 7) * 18;
  return `<rect x="${W-260}" y="${y-6}" width="14" height="3" fill="${O_PAL[i % O_PAL.length]}"/><text x="${W-240}" y="${y}" fill="#cdd3e0" font-size="11">orbit ${i+1}  (size ${o.size})</text>`;
}).join("\n");

const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" font-family="ui-monospace,monospace">
<rect width="${W}" height="${H}" fill="#0a0c12"/>
<text x="24" y="38" fill="#e8ecf4" font-size="20">64-TG  ·  Forman-Ricci  +  O_h orbits duales</text>
<text x="24" y="58" fill="#8b93a7" font-size="12">240 aristas coloreadas por F = 4 - deg(u) - deg(v)  ·  144 dual edges coloreadas por orbita O_h</text>
<text x="24" y="${H-58}" fill="#cdd3e0" font-size="11">cross-checks (atlas v1.0):  256 = tr(A^3)/6  ·  beta = (1, 0, 14, 0)  ·  B_op rank = 12  ·  bipartite +/- = 0,0,144</text>
<text x="24" y="${H-40}" fill="#cdd3e0" font-size="11">H_1(dual) bajo O_h:  A2g + 2Eg + 7T1g + 5T2g + 2A1u + 2A2u + 6Eu + 4T1u + 4T2u  =  81</text>
<text x="24" y="${H-22}" fill="#8b93a7" font-size="10.5">dual tr(A^k) odd k = 0 (bipartite chirality)  ·  6 valores distintos de F-Ricci  ·  6 orbitas O_h sobre 144 dual edges</text>
${Fleg}
${Oleg}
${items.map(i=>i.s).join("\n")}
</svg>`;
writeFileSync(join(here, "tetra64-curvature.svg"), svg);
console.log("wrote tetra64/tetra64-curvature.svg");
