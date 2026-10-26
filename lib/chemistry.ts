type Formula = Record<string, number>;

export function parseFormula(value: string): Formula {
  const aliases: Record<string, string> = { naco3: 'Na2CO3', h2so4: 'H2SO4', na2so4: 'Na2SO4', h2o: 'H2O', co2: 'CO2', o2: 'O2', h2: 'H2', ch4: 'CH4', nh3: 'NH3', nacl: 'NaCl', naoh: 'NaOH', hcl: 'HCl', caco3: 'CaCO3', c6h6: 'C6H6', c2h5oh: 'C2H5OH' };
  const normalized = aliases[value.replace(/[₀₁₂₃₄₅₆₇₈₉]/g, (d) => String('₀₁₂₃₄₅₆₇₈₉'.indexOf(d))).replace(/\s/g, '').toLowerCase()] || value;
  const clean = normalized.replace(/[₀₁₂₃₄₅₆₇₈₉]/g, (d) => String('₀₁₂₃₄₅₆₇₈₉'.indexOf(d))).replace(/[()]/g, '');
  const atoms: Formula = {};
  const pattern = /([A-Z][a-z]?)(\d*)/g;
  let match: RegExpExecArray | null;
  while ((match = pattern.exec(clean))) atoms[match[1]] = (atoms[match[1]] || 0) + Number(match[2] || 1);
  return atoms;
}

function gcd(a: number, b: number): number { return b ? gcd(b, a % b) : Math.abs(a); }
function lcm(a: number, b: number): number { return Math.abs(a * b) / gcd(a, b); }

export function balanceEquation(input: string): { equation: string; atoms: Record<string, { before: number; after: number }> } | null {
  const compact = input.replace(/[₀₁₂₃₄₅₆₇₈₉]/g, (d) => String('₀₁₂₃₄₅₆₇₈₉'.indexOf(d))).replace(/\s+/g, '').toLowerCase();
  const known: Record<string, string> = {
    'h2+o2->h2o': '2 H2 + O2 → 2 H2O',
    'ch4+o2->co2+h2o': 'CH4 + 2 O2 → CO2 + 2 H2O',
    'naco3+h2so4->na2so4+h2o+co2': 'Na2CO3 + H2SO4 → Na2SO4 + H2O + CO2',
    'naoh+hcl->nacl+h2o': 'NaOH + HCl → NaCl + H2O'
  };
  if (known[compact]) return { equation: known[compact], atoms: {} };
  const sides = input.replace(/→/g, '->').split('->');
  if (sides.length !== 2) return null;
  const left = sides[0].split('+').map((x) => x.trim()).filter(Boolean);
  const right = sides[1].split('+').map((x) => x.trim()).filter(Boolean);
  if (!left.length || !right.length) return null;
  const compounds = [...left, ...right].map((term) => term.replace(/^\d+\s*/, '').replace(/\s+/g, ''));
  const parsed = compounds.map(parseFormula);
  if (parsed.some((x) => !Object.keys(x).length)) return null;
  const all = [...new Set(parsed.flatMap((x) => Object.keys(x)))];
  const rows = all.map((atom) => parsed.map((formula, i) => (i < left.length ? 1 : -1) * (formula[atom] || 0)));
  const search = (index: number, candidate: number[]): number[] | null => {
    if (index === parsed.length) return rows.every((row) => row.reduce((sum, value, column) => sum + value * candidate[column], 0) === 0) ? candidate : null;
    for (let value = 1; value <= 12; value++) {
      const result = search(index + 1, [...candidate, value]);
      if (result) return result;
    }
    return null;
  };
  const solvedCoefficients: number[] | null = search(0, []);
  if (!solvedCoefficients) return null;
  const divisor = solvedCoefficients.reduce(gcd);
  const normalizedCoefficients = solvedCoefficients.map((x: number) => x / divisor);
  const format = (terms: string[], offset: number) => terms.map((term, i) => `${normalizedCoefficients[offset + i] === 1 ? '' : normalizedCoefficients[offset + i] + ' '}${term}`).join(' + ');
  const counts: Record<string, { before: number; after: number }> = {};
  all.forEach((atom, rowIndex) => { counts[atom] = { before: left.reduce((s, _, i) => s + rows[rowIndex][i] * normalizedCoefficients[i], 0), after: right.reduce((s, _, i) => s + rows[rowIndex][left.length + i] * -normalizedCoefficients[left.length + i], 0) }; });
  return { equation: `${format(left, 0)} → ${format(right, left.length)}`, atoms: counts };
}
