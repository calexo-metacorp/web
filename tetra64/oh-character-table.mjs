// Character table for O_h (octahedral group with inversion), 48 elements,
// 10 conjugacy classes, 10 irreducible representations.
//
// Class order (index 0..9):
//   0: E         (1  element)
//   1: 8 C_3     (8  elements: 3-cycles in coords, det=+1)
//   2: 6 C_2'    (6  elements: edge-axis rotations, transp perm, det=+1, trace=-1)
//   3: 6 C_4     (6  elements: face-axis 90 deg rotations, transp perm, det=+1, trace=+1)
//   4: 3 C_2     (3  elements: face-axis 180 deg rotations, id perm, two -1 signs)
//   5: i         (1  element : inversion, det=-1, s=(-1,-1,-1))
//   6: 8 S_6     (8  elements: i*C_3, det=-1, trace=0)
//   7: 6 sigma_d (6  elements: dihedral mirrors, transp perm, det=-1, trace=+1)
//   8: 6 S_4     (6  elements: rotoreflections, transp perm, det=-1, trace=-1)
//   9: 3 sigma_h (3  elements: face mirrors, id perm, one -1 sign)
//
// Irrep order: A1g A2g Eg T1g T2g  A1u A2u Eu T1u T2u  (dims 1 1 2 3 3 1 1 2 3 3)

export const OH_CLASS_NAMES = ["E","8C3","6C2'","6C4","3C2","i","8S6","6sigma_d","6S4","3sigma_h"];
export const OH_CLASS_SIZES = [1, 8, 6, 6, 3, 1, 8, 6, 6, 3];
export const OH_GROUP_ORDER = 48;

export const OH_IRREP_NAMES = ["A1g","A2g","Eg","T1g","T2g","A1u","A2u","Eu","T1u","T2u"];
export const OH_IRREP_DIMS  = [1, 1, 2, 3, 3, 1, 1, 2, 3, 3];

// rows = irreps, cols = classes in order above
export const OH_CHARACTER_TABLE = [
  // E   8C3  6C2'  6C4  3C2   i   8S6  6sd  6S4  3sh
  [  1,   1,   1,   1,   1,   1,   1,   1,   1,   1 ], // A1g
  [  1,   1,  -1,  -1,   1,   1,   1,  -1,  -1,   1 ], // A2g
  [  2,  -1,   0,   0,   2,   2,  -1,   0,   0,   2 ], // Eg
  [  3,   0,  -1,   1,  -1,   3,   0,  -1,   1,  -1 ], // T1g
  [  3,   0,   1,  -1,  -1,   3,   0,   1,  -1,  -1 ], // T2g
  [  1,   1,   1,   1,   1,  -1,  -1,  -1,  -1,  -1 ], // A1u
  [  1,   1,  -1,  -1,   1,  -1,  -1,   1,   1,  -1 ], // A2u
  [  2,  -1,   0,   0,   2,  -2,   1,   0,   0,  -2 ], // Eu
  [  3,   0,  -1,   1,  -1,  -3,   0,   1,  -1,   1 ], // T1u
  [  3,   0,   1,  -1,  -1,  -3,   0,  -1,   1,   1 ], // T2u
];

// classify a signed-permutation O_h element by its conjugacy class
// g = { pi: [a,b,c], s: [+-1,+-1,+-1], det: +-1 }
export function classifyOh(g) {
  const fixed = [0,1,2].filter((i) => g.pi[i] === i);
  const trace = fixed.reduce((acc, i) => acc + g.s[i], 0);
  const fixedCount = fixed.length;
  if (g.det === +1) {
    if (fixedCount === 3) return trace === 3 ? 0 : 4;      // E / 3C2
    if (fixedCount === 1) return trace === 1 ? 3 : 2;      // 6C4 / 6C2'
    /* fixedCount === 0 */    return 1;                     // 8C3
  } else {
    if (fixedCount === 3) return trace === -3 ? 5 : 9;     // i / 3sigma_h
    if (fixedCount === 1) return trace === 1 ? 7 : 8;      // 6sigma_d / 6S4
    /* fixedCount === 0 */    return 6;                     // 8S6
  }
}

// Inner-product based irrep decomposition.
//   classChars: array of length 10 giving the character of the representation per class
// Returns array of multiplicities a_i for each irrep.
export function decomposeIrrep(classChars) {
  const out = [];
  for (let i = 0; i < OH_IRREP_NAMES.length; i++) {
    let sum = 0;
    for (let c = 0; c < 10; c++) sum += OH_CLASS_SIZES[c] * OH_CHARACTER_TABLE[i][c] * classChars[c];
    out.push(sum / OH_GROUP_ORDER);
  }
  return out;
}

export function formatDecomposition(mults) {
  const parts = [];
  for (let i = 0; i < mults.length; i++) {
    const m = mults[i];
    if (Math.abs(m) < 1e-9) continue;
    const r = Math.round(m);
    parts.push(`${r === 1 ? "" : (Math.abs(r - m) < 1e-6 ? r : m.toFixed(3))}${OH_IRREP_NAMES[i]}`);
  }
  return parts.join(" + ");
}
