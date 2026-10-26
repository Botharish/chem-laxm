"""ChemSphere scientific service.

The Next.js app remains the primary UI. This small service provides validated,
reproducible scientific operations where Python's chemistry ecosystem is useful.
"""
from __future__ import annotations

import re
from typing import Any

import numpy as np
import sympy as sp
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from pint import UnitRegistry
from scipy.integrate import solve_ivp

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
except ImportError:  # keeps the API startup error explicit if RDKit is unavailable
    Chem = None

ureg = UnitRegistry()
app = FastAPI(title="ChemSphere Scientific API", version="1.0.0")
app.add_middleware(CORSMiddleware, allow_origins=["http://localhost:3000"], allow_methods=["*"], allow_headers=["*"])

class SmilesRequest(BaseModel):
    smiles: str = Field(min_length=1, max_length=500)

class BalanceRequest(BaseModel):
    equation: str = Field(min_length=3, max_length=1000)

class QuantityRequest(BaseModel):
    value: float
    from_unit: str
    to_unit: str

class KineticsRequest(BaseModel):
    initial_a: float = Field(ge=0)
    initial_b: float = Field(ge=0)
    rate_constant: float = Field(gt=0)
    duration: float = Field(gt=0, le=10000)
    points: int = Field(default=100, ge=10, le=2000)

@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok", "rdkit": "available" if Chem else "missing"}

@app.post("/api/structure")
def structure(payload: SmilesRequest) -> dict[str, Any]:
    if Chem is None:
        raise HTTPException(503, "RDKit is not installed in the API environment")
    molecule = Chem.MolFromSmiles(payload.smiles)
    if molecule is None:
        raise HTTPException(422, "Invalid SMILES: RDKit could not parse this structure")
    molecule = Chem.AddHs(molecule)
    if AllChem.EmbedMolecule(molecule, randomSeed=0xC0FFEE) != 0:
        AllChem.Compute2DCoords(molecule)
    else:
        AllChem.UFFOptimizeMolecule(molecule)
    conformer = molecule.GetConformer()
    coordinates = [{"element": atom.GetSymbol(), "x": conformer.GetAtomPosition(atom.GetIdx()).x, "y": conformer.GetAtomPosition(atom.GetIdx()).y, "z": 0.0} for atom in molecule.GetAtoms()]
    return {"valid": True, "canonical_smiles": Chem.MolToSmiles(Chem.RemoveHs(molecule)), "formula": rdMolDescriptors.CalcMolFormula(molecule), "molar_mass": round(Descriptors.MolWt(molecule), 4), "atoms": coordinates, "molfile": Chem.MolToMolBlock(molecule)}

@app.post("/api/balance")
def balance(payload: BalanceRequest) -> dict[str, Any]:
    sides = re.split(r"(?:->|→|=)", payload.equation)
    if len(sides) != 2:
        raise HTTPException(422, "Use a complete equation with reactants and products separated by → or ->")
    left = [x.strip() for x in sides[0].split("+") if x.strip()]
    right = [x.strip() for x in sides[1].split("+") if x.strip()]
    formulas = [re.sub(r"^\d+\s*", "", x).replace(" ", "") for x in left + right]
    atom_names = sorted(set(re.findall(r"[A-Z][a-z]?", "".join(formulas))))
    matrix = sp.Matrix([[sp.Integer((1 if i < len(left) else -1) * _count(formula, atom)) for i, formula in enumerate(formulas)] for atom in atom_names])
    null = matrix.nullspace()
    if not null:
        raise HTTPException(422, "No atom-conserving integer solution was found")
    vector = null[0]
    multiplier = sp.ilcm(*[term.q for term in vector])
    coefficients = [int(term * multiplier) for term in vector]
    divisor = int(sp.igcd(*coefficients)); coefficients = [x // divisor for x in coefficients]
    formatted = " + ".join(f"{'' if c == 1 else str(c) + ' '}{f}" for c, f in zip(coefficients, formulas[:len(left)])) + " → " + " + ".join(f"{'' if c == 1 else str(c) + ' '}{f}" for c, f in zip(coefficients[len(left):], formulas[len(left):]))
    return {"equation": formatted, "coefficients": coefficients, "atoms": {atom: {"reactants": sum(coefficients[i] * _count(formulas[i], atom) for i in range(len(left))), "products": sum(coefficients[len(left) + i] * _count(formulas[len(left) + i], atom) for i in range(len(right)))} for atom in atom_names}}

def _count(formula: str, atom: str) -> int:
    match = re.search(rf"{atom}(\d*)", formula)
    return int(match.group(1) or 1) if match else 0

@app.post("/api/units")
def convert_units(payload: QuantityRequest) -> dict[str, Any]:
    try:
        quantity = payload.value * ureg(payload.from_unit)
        converted = quantity.to(payload.to_unit)
    except Exception as exc:
        raise HTTPException(422, f"Unsupported or incompatible units: {exc}") from exc
    return {"value": converted.magnitude, "unit": str(converted.units)}

@app.post("/api/kinetics/first-order")
def kinetics(payload: KineticsRequest) -> dict[str, list[float]]:
    def derivative(_time: float, state: np.ndarray) -> list[float]:
        rate = payload.rate_constant * state[0] * state[1]
        return [-rate, -rate, rate]
    times = np.linspace(0, payload.duration, payload.points)
    solution = solve_ivp(derivative, (0, payload.duration), [payload.initial_a, payload.initial_b, 0.0], t_eval=times)
    if not solution.success:
        raise HTTPException(500, "Kinetics solver failed")
    return {"time": solution.t.tolist(), "A": solution.y[0].tolist(), "B": solution.y[1].tolist(), "product": solution.y[2].tolist()}
