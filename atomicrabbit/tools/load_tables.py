from astroquery.nist import Nist
import astropy.units as u
import yaml
from atomicrabbit.physics.states import FineState
from atomicrabbit.physics.transitions import Transition
import pandas as pd

'''
This code snippet uses astroquery to import data from nist and turns the data for
relevant transitions into a .yaml file.
It's fine as long as it works so it's written by Claude

'''
import re
import yaml
import numpy as np
from astroquery.nist import Nist
import astropy.units as u

L_LETTERS = {'S':0, 'P':1, 'D':2, 'F':3, 'G':4, 'H':5}

def parse_term_symbol(term_str):
    """Parse e.g. '1P*_1' or '3D_2' into (S, L, J, parity)"""
    
    if not isinstance(term_str, str):
        return None
    parity = 'odd' if '*' in term_str else 'even'
    match = re.match(r'(\d+)([SPDFGH])\*?_([\d.]+)', term_str.strip())
    
    if match:
        multiplicity = int(match.group(1))
        S = (multiplicity - 1) / 2
        L = L_LETTERS[match.group(2)]
        J = float(match.group(3))
        return S, L, J, parity
    return None

def parse_state(level_str, energy):
    """
    Parse a level string like '6s.6p  | 1P*  | 1' into a State.
    """
    if not isinstance(level_str, str):
        return None
    
    parts = [p.strip() for p in level_str.split('|')]
    if len(parts) != 3:
        return None
    
    config, term, J_str = parts
    parity = 'odd' if '*' in term else 'even'
    
    match = re.match(r'(\d+)([SPDFGH])\*?', term.strip())
    if match is None:
        return None
    
    multiplicity = int(match.group(1))
    S = (multiplicity - 1) / 2
    L = L_LETTERS[match.group(2)]
    J = float(J_str)
    
    return FineState(
        configuration=config.strip(),
        S=S, L=L, J=J,
        term=term,
        energy_cm=float(energy),
        parity=parity
    )

def load_transitions_from_table(table):
    """
    Parse an astroquery NIST table into a list of Transition objects,
    and return a YAML-friendly list in parallel.
    """
    transitions = []
    yaml_data = []
    ei_ek_col = [c for c in table.colnames if 'Ei' in c][0]
    
    
    for row in table:
        try:
            ei, ek = row[ei_ek_col].split('-')
            ei = float(ei.strip())
            ek = float(ek.strip())
            
            lower = parse_state(row['Lower level'], energy=ei)
            upper = parse_state(row['Upper level'], energy=ek)
            
            if lower is None or upper is None:
                continue

            A = float(row['Aki']) if row['Aki'] else None
            t = Transition(lower=lower, upper=upper, A=A)
            transitions.append(t)

            yaml_data.append({
                'lower': lower.spectroscopic_name(),
                'upper': upper.spectroscopic_name(),
                'lower_energy_cm-1': lower.energy,
                'upper_energy_cm-1': upper.energy,
                'wavelength_nm': float(row['Observed']),  # adjust colname
                'Aki': A,
            })

        except Exception as e:
            print(f"Skipping row due to error: {e}")
            continue

    return transitions, yaml_data

# --- main ---
table = Nist.query(500 * u.nm, 600 * u.nm, linename="Ba I")
table = table[pd.to_numeric(table['Rel.'], errors='coerce') > 100]
print(table)
print(table['Lower level'][:3])
print(table['Upper level'][:3])

transitions, yaml_data = load_transitions_from_table(table)

with open('ba_transitions.yaml', 'w') as f:
    yaml.dump(yaml_data, f, allow_unicode=True, sort_keys=False)

print(f"Loaded {len(transitions)} transitions")

