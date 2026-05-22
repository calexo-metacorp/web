#!/usr/bin/env node
// Renders a static SVG snapshot of the 64-tetrahedron / vector-equilibrium
// artifact using the shared geometry (no dependencies, painter's algorithm).
//   node tetra64/render-svg.mjs

import { writeFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";
import { build } from "./geometry.mjs";

const here = dirname(fileURLToPath(import.meta.url));
const M = build();

const W = 960, H = 680;
const PADX = 40, TOP = 80, BOT = 44;            // margins reserved for labels
const CX = W / 2, CY = TOP + (H - TOP - BOT) / 2; // drawing centre
const yaw = 0.62, pitch = -0.46;
const F = 7;
const sub = (a, b) => [a[0]-b[0], a[1]-b[1], a[2]-b[2]];
const dot = (a, b) => a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
const cross = (a, b) => [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
const llen = (a) => Math.hypot(a[0], a[1], a[2]);
const cen = (p) => { const s=p.reduce((a,b)=>[a[0]+b[0],a[1]+b[1],a[2]+b[2]],[0,0,0]); return [s[0]/p.length,s[1]/p.length,s[2]/p.length]; };

function rot(p){
  const cy=Math.cos(yaw), sy=Math.sin(yaw), cx=Math.cos(pitch), sx=Math.sin(pitch);
  const x=p[0]*cy+p[2]*sy, z=-p[0]*sy+p[2]*cy, y=p[1];
  return [x, y*cx-z*sx, y*sx+z*cx];
}
// auto-fit: measure pre-scale extent, then pick a scale that fits the frame
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
const L=(()=>{const l=[0.35,0.55,0.9],m=llen(l);return [l[0]/m,l[1]/m,l[2]/m];})();
function shade(rgb, v3, emis){
  const n=cross(sub(v3[1],v3[0]),sub(v3[2],v3[0])), m=llen(n)||1;
  const lit=Math.abs(dot([n[0]/m,n[1]/m,n[2]/m],L)), k=emis+(1-emis)*lit;
  return `rgb(${Math.round(rgb[0]*k)},${Math.round(rgb[1]*k)},${Math.round(rgb[2]*k)})`;
}
const items=[];

// 64 tetrahedra (wireframe)
const TE=[[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]];
for(const t of M.tetrahedra){
  const P=t.vertices.map(proj), R=t.vertices.map(rot), z=cen(R)[2];
  const d=TE.map(([a,b])=>`M${P[a].x.toFixed(1)} ${P[a].y.toFixed(1)}L${P[b].x.toFixed(1)} ${P[b].y.toFixed(1)}`).join("");
  items.push({z, s:`<path d="${d}" stroke="#7aa2ff" stroke-opacity="0.5" stroke-width="0.7" fill="none"/>`});
}

// nucleus cuboctahedron
const V=M.vectorEquilibrium.vertices, P=V.map(proj), R=V.map(rot);
const faces=[...M.vectorEquilibrium.triFaces.map(f=>({f,sq:false})),...M.vectorEquilibrium.sqFaces.map(f=>({f,sq:true}))];
for(const {f,sq} of faces){
  const v3=f.map(i=>R[i]), z=cen(v3)[2]+0.002;
  const pts=f.map(i=>`${P[i].x.toFixed(1)},${P[i].y.toFixed(1)}`).join(" ");
  items.push({z, s:`<polygon points="${pts}" fill="${shade([245,196,81],v3,0.42)}" fill-opacity="${sq?0.55:0.7}" stroke="#f5c451" stroke-width="1.6"/>`});
}
// 12 radial equilibrium vectors
const o=proj([0,0,0]);
for(let i=0;i<V.length;i++){
  const p=P[i], z=rot(V[i])[2]/2;
  items.push({z, s:`<line x1="${o.x.toFixed(1)}" y1="${o.y.toFixed(1)}" x2="${p.x.toFixed(1)}" y2="${p.y.toFixed(1)}" stroke="#4fd0e0" stroke-opacity="0.85" stroke-width="1.1"/>`});
}
for(let i=0;i<V.length;i++){const p=P[i],z=rot(V[i])[2]+0.01;
  items.push({z, s:`<circle cx="${p.x.toFixed(1)}" cy="${p.y.toFixed(1)}" r="3.2" fill="#ffe49a"/>`});}
items.push({z:99, s:`<circle cx="${o.x.toFixed(1)}" cy="${o.y.toFixed(1)}" r="3.6" fill="#ffffff"/>`});

items.sort((a,b)=>a.z-b.z);

const svg=`<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" font-family="ui-monospace,monospace">
<rect width="${W}" height="${H}" fill="#0a0c12"/>
<text x="24" y="40" fill="#e8ecf4" font-size="20">64-Tetrahedron Grid</text>
<text x="24" y="62" fill="#8b93a7" font-size="13">8 star tetrahedra (8x8) around a cuboctahedral vector-equilibrium nucleus  ·  radius = edge = sqrt2</text>
<text x="24" y="${H-24}" fill="#8b93a7" font-size="12">Sigma(vertices)=0  ·  64 regular tetrahedra  ·  nucleus: 12 vertices, 24 edges, 8 tri + 6 squares  ·  symmetry Oh (48)</text>
${items.map(i=>i.s).join("\n")}
</svg>`;
writeFileSync(join(here,"tetra64.svg"), svg);
console.log("wrote tetra64/tetra64.svg");
