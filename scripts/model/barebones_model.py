import json
import re

from pathlib import Path
from cobra import Model, Reaction, Metabolite

# --- JSON parsing functions
def load_json(path):
    """
    Load a JSON file from disk.
    """
    with Path(path).open("r", encoding="utf-8") as f:
        return json.load(f)
    
# Parsed metabolite information should be added into a dictionary, 
# id and name comes from DEFINITON and EQUATION
# Formula comes from KEGG compund id request(FORMULA)

# Parsed reaction information should be used directly
# id part of it --> Reaction('string') comes from KEGG reaction ID
# name comes from genbank parsing --> "Protein": part
# lower and upper bound limits come from EQUATITION value 
#   <=> for -1000, 1000 
#   <= for -1000, 0
#   => for 0, 1000 <-- Default
# when adding metabolites, use EQUATITION to get the minus/plus values and get the metabolite from a dictionary
# gene rules comes from "Gene": part

def parse_reversibility(equation_str):
    if "<=>" in equation_str:
        return -1000, 1000
    elif "<=" in equation_str:
        return -1000, 0
    elif "=>" in equation_str:
        return 0, 1000
    else:
        raise ValueError(f"Unknown reaction direction: {equation_str}")

def split_equation(equation_str):
    if "<=>" in equation_str:
        left, right = equation_str.split("<=>")
    elif "=>" in equation_str:
        left, right = equation_str.split("=>")
    elif "<=" in equation_str:
        left, right = equation_str.split("<=")
    else:
        raise ValueError("No valid reaction arrow found")

    return left.strip(), right.strip()

def parse_equation_ids(equation_str):
    """
    Returns:
      left:  list of (compound_id, coeff)
      right: list of (compound_id, coeff)
    """
    left, right = split_equation(equation_str)

    def parse_side(side):
        parts = []
        for token in side.split("+"):
            token = token.strip()
            match = re.match(r"^(\d+)?\s*(C\d{5})$", token)
            if not match:
                raise ValueError(f"Bad KEGG token: {token}")

            coeff = int(match.group(1)) if match.group(1) else 1
            cid = match.group(2)
            parts.append((cid, coeff))
        return parts

    return parse_side(left), parse_side(right)


# --- Metabolite Dictionary
def get_or_create_metabolite(kegg_id, name, compound_metadata, metabolite_registry, compartment="c"):
    if kegg_id in metabolite_registry:
        return metabolite_registry[kegg_id]

    formula = compound_metadata.get(kegg_id, {}).get("formula")

    met = Metabolite(
        id=f"{kegg_id}_{compartment}",
        name=name,
        formula=formula,
        compartment=compartment
    )

    metabolite_registry[kegg_id] = met
    return met

# --- Plus/minus values for the reaction
def build_stoichiometry(equation_str, compound_metadata, metabolite_registry, compartment="c"):
    left, right = parse_equation_ids(equation_str)

    stoich = {}

    for cid, coeff in left:
        meta = compound_metadata.get(cid, {})
        name = meta.get("name", cid)

        met = get_or_create_metabolite(
            cid,
            name,
            compound_metadata,
            metabolite_registry,
            compartment
        )
        stoich[met] = -coeff

    for cid, coeff in right:
        meta = compound_metadata.get(cid, {})
        name = meta.get("name", cid)

        met = get_or_create_metabolite(
            cid,
            name,
            compound_metadata,
            metabolite_registry,
            compartment
        )
        stoich[met] = coeff

    return stoich


# --- Reaction Creation
def parse_reaction(reaction_id, equation, cmpd_ids, genes, compound_metadata, metabolite_registry, subsystem=""):
    lb, ub = parse_reversibility(equation)

    reaction = Reaction(
        id=reaction_id,
        name=equation,
        lower_bound=lb,
        upper_bound=ub,
        subsystem=subsystem
    )

    reaction.gene_reaction_rule = genes

    stoich = build_stoichiometry(
        equation,
        cmpd_ids,
        compound_metadata,
        metabolite_registry
    )

    reaction.add_metabolites(stoich)
    return reaction
