#!/usr/bin/env node
// Extended snapshot: chirality 32+/32-, shells S0..S5 dots, 8 stella centres,
// 240 bulk edges, 144 dual edges, nucleus + radials. Painter's algorithm SVG.
//   node tetra64/render-extended-svg.mjs

import { writeFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";
import { build } from "./geometry.mjs";
import { analyze } from "./analysis.mjs";

const here = dirname(fileURLToPath(import.meta.url));
const M = build();
const Araw = analyze(M);
// adapt analysis.mjs's shape to the field names this renderer expects
const A = {
  verts: Araw.vertices,
  shells: Araw.shells,
  edges: Araw.edges,
  dualEdges: Araw.dualGraphs.edge.pairs,
};

const W = 1100, H = 760, PADX = 40, TOP = 80, BOT = 60;
const CX = W / 2, CY = TOP + (H - TOP - BOT) / 2;
const yaw = 0.62, pitch = -0.46, F = 7;

const sub = (a, b) => [a[0]-b[0], a[1]-b[1], a[2]-b[2]];
const ddot = (a, b) => a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
const cross = (a, b) => [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
const llen = (a) => Math.hypot(a[0], a[1], a[2]);
const cen = (p) => { const s=p.reduce((a,b)=>[a[0]+b[0],a[1]+b[1],a[2]+b[2]],[0,0,0]); return [s[0]/p.length,s[1]/p.length,s[2]/p.length]; };
function rot(p){
  const cy=Math.cos(yaw), sy=Math.sin(yaw), cx=Math.cos(pitch), sx=Math.sin(pitch);
  const x=p[0]*cy+p[2]*sy, z=-p[0]*sy+p[2]*cy, y=p[1];
  return [x, y*cx-z*sx, y*sx+z*cx];
}
const ALL = [...M.tetrahedra.flatMap(t=>t.vertices), ...M.vectorEquilibrium.vertices, ...M.stellaCenters];
const SCALE = (() => {
  let mu=0, mv=0;
  for(const p of ALL){ const r=rot(p), k=F/(F-r[2]); mu=Math.max(mu,Math.abs(r[0]*k)); mv=Math.max(mv,Math.abs(r[1]*k)); }
  return Math.min((W/2-PADX)/mu, ((H-TOP-BOT)/2)/mv);
})();
function proj(p){
  const r=rot(p), persp=F/(F-r[2]), s=SCALE*persp;
  return { x:CX+r[0]*s, y:CY-r[1]*s, z:r[2] };
}

const SHELL = ["#ffffff","#f5c451","#4fd0e0","#ff7ad9","#7be094","#ffae5a"];
const CHIRP="#ff8c5a", CHIRM="#5ec0ff";
const items = [];

// 240 bulk edges (faint)
for(const e of A.edges){
  const pa=proj(e.a), pb=proj(e.b), z=(pa.z+pb.z)/2 - 0.01;
  items.push({z, s:`<line x1="${pa.x.toFixed(1)}" y1="${pa.y.toFixed(1)}" x2="${pb.x.toFixed(1)}" y2="${pb.y.toFixed(1)}" stroke="#3f475a" stroke-opacity="0.55" stroke-width="0.5"/>`});
}

// 144 dual edges (between tetra centroids)
for(const [i,j] of A.dualEdges){
  const a=M.tetrahedra[i].centroid, b=M.tetrahedra[j].centroid;
  const pa=proj(a), pb=proj(b), z=(pa.z+pb.z)/2 - 0.005;
  items.push({z, s:`<line x1="${pa.x.toFixed(1)}" y1="${pa.y.toFixed(1)}" x2="${pb.x.toFixed(1)}" y2="${pb.y.toFixed(1)}" stroke="#6fdc8c" stroke-opacity="0.45" stroke-width="0.6"/>`});
}

// 64 tetra wireframe coloured by chirality
const TE=[[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]];
for(const t of M.tetrahedra){
  const v=t.vertices;
  const sv = (b,c,d)=>{
    const v1=sub(b,v[0]),v2=sub(c,v[0]),v3=sub(d,v[0]);
    return v1[0]*(v2[1]*v3[2]-v2[2]*v3[1])-v1[1]*(v2[0]*v3[2]-v2[2]*v3[0])+v1[2]*(v2[0]*v3[1]-v2[1]*v3[0]);
  };
  const chir = Math.sign(sv(v[1],v[2],v[3]));
  const col = chir>0?CHIRP:CHIRM;
  const P=v.map(proj), R=v.map(rot), z=cen(R)[2];
  const d=TE.map(([a,b])=>`M${P[a].x.toFixed(1)} ${P[a].y.toFixed(1)}L${P[b].x.toFixed(1)} ${P[b].y.toFixed(1)}`).join("");
  items.push({z, s:`<path d="${d}" stroke="${col}" stroke-opacity="0.7" stroke-width="0.8" fill="none"/>`});
}

// nucleus cuboctahedron (gold faces + edges)
const VEv=M.vectorEquilibrium.vertices, Pve=VEv.map(proj), Rve=VEv.map(rot);
function shade(rgb,v3,emis){
  const n=cross(sub(v3[1],v3[0]),sub(v3[2],v3[0])), m=llen(n)||1;
  const lit=Math.abs(ddot([n[0]/m,n[1]/m,n[2]/m],[0.35,0.55,0.9])/llen([0.35,0.55,0.9])), k=emis+(1-emis)*lit;
  return `rgb(${Math.round(rgb[0]*k)},${Math.round(rgb[1]*k)},${Math.round(rgb[2]*k)})`;
}
const faces=[...M.vectorEquilibrium.triFaces.map(f=>({f,sq:false})),...M.vectorEquilibrium.sqFaces.map(f=>({f,sq:true}))];
for(const{f,sq} of faces){
  const v3=f.map(i=>Rve[i]), z=cen(v3)[2]+0.002;
  const pts=f.map(i=>`${Pve[i].x.toFixed(1)},${Pve[i].y.toFixed(1)}`).join(" ");
  items.push({z, s:`<polygon points="${pts}" fill="${shade([245,196,81],v3,0.42)}" fill-opacity="${sq?0.5:0.65}" stroke="#f5c451" stroke-width="1.4"/>`});
}

// 12 radials
const o=proj([0,0,0]);
for(let i=0;i<VEv.length;i++){
  const p=Pve[i], z=Rve[i][2]/2;
  items.push({z, s:`<line x1="${o.x.toFixed(1)}" y1="${o.y.toFixed(1)}" x2="${p.x.toFixed(1)}" y2="${p.y.toFixed(1)}" stroke="#4fd0e0" stroke-opacity="0.7" stroke-width="0.9"/>`});
}

// 63 vertex dots coloured by shell
for(const v of A.verts){
  const p=proj(v.pos), z=rot(v.pos)[2]+0.005;
  const r= v.shell===0?5: (v.shell===1?3.2:2.4);
  items.push({z, s:`<circle cx="${p.x.toFixed(1)}" cy="${p.y.toFixed(1)}" r="${r}" fill="${SHELL[v.shell]}"/>`});
}

// 8 stella centres
for(const c of M.stellaCenters){
  const p=proj(c), z=rot(c)[2]+0.006;
  items.push({z, s:`<circle cx="${p.x.toFixed(1)}" cy="${p.y.toFixed(1)}" r="3.6" fill="#ffffff" stroke="#000" stroke-width="0.4"/>`});
}
// stella cube edges
for(let i=0;i<8;i++) for(let j=i+1;j<8;j++){
  const d=sub(M.stellaCenters[i],M.stellaCenters[j]);
  if(Math.abs(ddot(d,d)-4)<1e-9){
    const pa=proj(M.stellaCenters[i]), pb=proj(M.stellaCenters[j]);
    const z=(pa.z+pb.z)/2 - 0.001;
    items.push({z, s:`<line x1="${pa.x.toFixed(1)}" y1="${pa.y.toFixed(1)}" x2="${pb.x.toFixed(1)}" y2="${pb.y.toFixed(1)}" stroke="#ffffff" stroke-opacity="0.5" stroke-width="0.6" stroke-dasharray="3 3"/>`});
  }
}
items.push({z:99, s:`<circle cx="${o.x.toFixed(1)}" cy="${o.y.toFixed(1)}" r="3.4" fill="#ffffff"/>`});

items.sort((a,b)=>a.z-b.z);

// legend block
const LEG = [
  ["S0  origen",              SHELL[0], 1],
  ["S1  VE (12)",             SHELL[1], 12],
  ["S2  axis +-2 (6)",        SHELL[2], 6],
  ["S3  perm(+-2,+-1,+-1) (24)", SHELL[3], 24],
  ["S4  perm(+-2,+-2,0) (12)",   SHELL[4], 12],
  ["S5  (+-2,+-2,+-2) (8)",      SHELL[5], 8],
];
let legendY = TOP + 4;
const legendItems = LEG.map((row, i) => {
  const y = legendY + i*18;
  return `<circle cx="${W-180}" cy="${y}" r="5" fill="${row[1]}"/><text x="${W-168}" y="${y+4}" fill="#cdd3e0" font-size="11">${row[0]}</text>`;
}).join("\n");

const svg=`<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" font-family="ui-monospace,monospace">
<rect width="${W}" height="${H}" fill="#0a0c12"/>
<text x="24" y="38" fill="#e8ecf4" font-size="20">64-TG  +  capas extendidas</text>
<text x="24" y="58" fill="#8b93a7" font-size="12">quiralidad 32+/32-  ·  cascarones S0..S5  ·  240 aristas  ·  144 dual edges  ·  8 centros estrella</text>
<text x="24" y="${H-36}" fill="#cdd3e0" font-size="11">verificacion:  64 tetra  ·  63 vertices  ·  32+/32-  ·  240 edges (R240)  ·  144 dual edges (R144)  ·  81 dual cycles (R81)  ·  cycle-space dim = E-V+1 = 178 (R178)</text>
<text x="24" y="${H-18}" fill="#8b93a7" font-size="10.5">naranja = quiralidad +    azul = quiralidad -    dorado = nucleo VE    cian = radiales    verde = grafo dual</text>
${legendItems}
${items.map(i=>i.s).join("\n")}
</svg>`;
writeFileSync(join(here,"tetra64-extended.svg"), svg);
console.log("wrote tetra64/tetra64-extended.svg");
