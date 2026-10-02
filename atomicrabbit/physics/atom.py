from dataclasses import dataclass
from typing import Optional
import numpy as np
from qutip import *
from scipy.constants import h, hbar, epsilon_0, c, k, physical_constants
from atomicrabbit.physics import states
from atomicrabbit.physics import lasers
from atomicrabbit.physics import transitions
import networkx as nx

@dataclass(frozen=True)
class Atom:
    species: str
    isotope: str
    transition_list: list[transitions.FineTransition]
    laser_list: list[lasers.Laser]
    I: Optional[float] = 0
    mu_I: Optional[float] = 0
    B: Optional[float] = 0
    
    # def load_transition(self, table):

    @classmethod
    def make_atom(cls, element, mass_number, db_path, laser_list, B=0):
        """
        Construct an Atom from spectroscopy.db, replacing hand-coded FineState/Transition setup.
        Usage: Yb174 = Atom.make_atom('Yb', 174, db_path, laser_list, B=1e-4)

        Loads is_recommended energy levels for the isotope (energies and hyperfine
        constants from paper sources) and matches them to NIST A coefficients by
        energy proximity (isotope shifts are <0.01 cm⁻¹, so a 1 cm⁻¹ tolerance is safe).

        Parameters
        ----------
        element : str           Element symbol, e.g. 'Yb'.
        mass_number : int       Isotope mass number, e.g. 174.
        db_path : str           Path to spectroscopy.db.
        laser_list : list       List of Laser objects for the simulation.
        B : float               Magnetic field in Tesla.

        Returns
        -------
        Atom
        """
        import sqlite3 as _sqlite3
        from atomicrabbit.db import queries as _q

        conn = _sqlite3.connect(db_path)
        conn.row_factory = _sqlite3.Row
        try:
            atom_props = _q.get_atom_properties(conn, element, mass_number)
            state_rows = _q.get_fine_states(conn, element, mass_number)
            nist_rows  = _q.get_nist_transitions(conn, element)
        finally:
            conn.close()

        # Build FineState objects keyed by (J, energy_cm1).
        # Multiple papers may contribute the same level; keep only the first seen.
        # For jj-coupled states (L or S is None in DB), fall back to L=S=0 — gJ
        # will be approximate but those states have no NIST A coefficient anyway.
        seen = set()
        state_map = {}
        for r in state_rows:
            key = (float(r["J"]), float(r["energy_cm1"]))
            if key in seen:
                continue
            seen.add(key)
            L = float(r["L"]) if r["L"] is not None else 0.0
            S = float(r["S"]) if r["S"] is not None else 0.0
            J = float(r["J"])
            # J=0 states have no electronic angular momentum so A=B=0 always.
            # hyperfine_split raises if both are None, so supply 0.0 explicitly.
            A_hf = r["hyperfine_A_MHz"] if J > 0 else 0.0
            B_hf = r["hyperfine_B_MHz"] if J > 0 else 0.0
            fs = states.FineState(
                configuration=r["configuration"],
                L=L,
                S=S,
                J=J,
                energy_cm=float(r["energy_cm1"]),
                parity=r["parity"],
                A_hyperfine=A_hf,
                B_hyperfine=B_hf,
            )
            state_map[key] = fs

        def _find_state(J, E, tol=1.0):
            J, E = float(J), float(E)
            for (sJ, sE), fs in state_map.items():
                if sJ == J and abs(sE - E) < tol:
                    return fs
            return None

        transition_list = []
        for t in nist_rows:
            lower = _find_state(t["lower_J"], t["lower_E"])
            upper = _find_state(t["upper_J"], t["upper_E"])
            if lower is None or upper is None:
                continue
            transition_list.append(transitions.FineTransition(
                lower=lower,
                upper=upper,
                A=float(t["einstein_A_s1"]),
            ))

        I    = float(atom_props.get("nuclear_spin") or 0)
        mu_I = float(atom_props.get("mu_I") or 0)

        return cls(
            species=f"{element} I",
            isotope=str(mass_number),
            transition_list=transition_list,
            laser_list=laser_list,
            I=I,
            mu_I=mu_I,
            B=B,
        )

    def add_fine_states(self):
        fine_states = set()
        for transition in self.transition_list:
            fine_states.add(transition.lower)
            fine_states.add(transition.upper)
        return list(fine_states)
    
    def add_zeeman_states(self):
        fine_states = self.add_fine_states()
        
        if self.I == 0:
            return states.batch_zeeman_split(fine_states, self.B)
        else:
            hyperfine_states = states.batch_hyperfine_split(fine_states, I=self.I, mu_I=self.mu_I)
            return states.batch_zeeman_split(hyperfine_states, self.B)
        
    def basis(self):
        states = self.add_zeeman_states()
        index = {s: i for i, s in enumerate(states)}
        return states, index
    
    def magnetic_sublevels(self, state):
        if self.I == 0:
            return states.zeeman_split(state, self.B)
    
        return states.batch_zeeman_split(
            states.hyperfine_split(
                state,
                I=self.I,
                mu_I=self.mu_I,
            ),
            self.B,
        )

    def zeeman_transitions(self):
        zeeman_transitions = []
    
        for transition in self.transition_list:
    
            lower_states = self.magnetic_sublevels(transition.lower)
            upper_states = self.magnetic_sublevels(transition.upper)
    
            for lower in lower_states:
                for upper in upper_states:
                    zeeman_transitions.append(
                        transitions.Transition(
                            lower=lower,
                            upper=upper,
                            A=transition.A,
                            linewidth=transition.linewidth,
                            lifetime=transition.lifetime,
                        )
                    )

        return zeeman_transitions
    
    def branching_ratios(self):
        """
        Returns normalized branching ratios for every Zeeman transition.
        """
    
        ratios = {}
    
        for transition in self.transition_list:
    
            lower_states = self.magnetic_sublevels(transition.lower)
            upper_states = self.magnetic_sublevels(transition.upper)
    
            # Collect all allowed decay amplitudes
            amplitudes = {}
    
            for upper in upper_states:
                total = 0.0
    
                for lower in lower_states:
    
                    decay = transitions.Transition(
                        lower=lower,
                        upper=upper,
                        A=transition.A,
                        linewidth=transition.linewidth,
                        lifetime=transition.lifetime,
                    )
    
                    amp2 = abs(decay.branching_amplitude)**2
    
                    if amp2 > 0:
                        amplitudes[decay] = amp2
                        total += amp2
    
                # Normalize excited state's decay channels
                for decay, amp2 in amplitudes.items():
                    if decay.upper == upper:
                        ratios[decay] = amp2 / total
    
        return ratios

    def driven_transitions(self):
        # Add situations where the transitions are not zeeman states later
        driven_transition = []
        zeeman_transitions = self.zeeman_transitions()
        for zeeman_transition in zeeman_transitions:
            for laser in self.laser_list:
                dipole = zeeman_transition.transition_dipole(laser)
                detuning = zeeman_transition.detuning(laser)
                
                # Dropping the dipole-forbidden terms
                if dipole != 0:
                    
                    # Dropping the far-off-resonant transitions
                    if abs(detuning) <= 1e13:
                        driven_transition.append(transitions.DrivenTransition(
                            transition = zeeman_transition,
                            laser = laser,
                            dipole = dipole
                            ))
                        
                    print(f'Driven transition {zeeman_transition.lower.spectroscopic_name} - {zeeman_transition.upper.spectroscopic_name}')
                    print(f'Transition dipole = {dipole}')
                    print(f'Laser = {laser}')                
                        
        return driven_transition
    
    def driven_transition_debug(self):
        for driven in self.driven_transitions():
            print(
                f"{driven.lower.spectroscopic_name}"
                f" --> "
                f"{driven.upper.spectroscopic_name}"
                f"   q={driven.laser.polarization}"
                f"   Ω={driven.rabi_frequency:.3e}"
                f"   Δ={driven.detuning:.3e}"
            )
        
    @property
    def system_size(self):
        # add non-zeeman states later too
        # this thing calls zeeman_transitions twice and might slow the calculations down!
        return len(self.zeeman_transitions())
    
    def atom_graph(self):
        """
        Building the Grotrian diagram for rotating frame energy assignment.
        The edges are transitions driven by the lasers.
        """
        G = nx.Graph()
        G.add_nodes_from(self.add_zeeman_states())
        for driven in self.driven_transitions():
            G.add_edge(driven.transition.lower, driven.transition.upper,
                       detuning=driven.laser.frequency - driven.transition.frequency)
        return G
    
    def assign_rotating_frame_energies(self, G=None):
        """
        Assigning the rotating frame energies.
        The rotating frame is determined by a tree diagram, where the nodes are separated by the laser frequencies.
        The off-diagonal terms are the cumulative detunings to reach the state from ground state.
        """
        if G is None:
            G = self.atom_graph()
        energies = {}
        for component in nx.connected_components(G):
            ref = min(component, key=lambda s: s.energy_omega)
            sub_energies = {ref: ref.energy_omega}

            for u, v in nx.bfs_edges(G, ref):
                detuning = G[u][v]['detuning']
                sub_energies[v] = sub_energies[u] - 2 * np.pi * detuning
            energies.update(sub_energies)

        return energies

    def build_hamiltonian(self):
        states, index = self.basis()
        n = len(states)
        H = np.zeros((n, n), dtype=complex)

        # Build graph once; shared by diagonal and edge check
        G = self.atom_graph()

        # Diagonal: driven states (have edges) get rotating-frame energies.
        # Isolated nodes are undriven — their raw energy_omega is optical-scale; skip them.
        energies = self.assign_rotating_frame_energies(G)
        for state, E in energies.items():
            if G.degree(state) > 0:
                H[index[state], index[state]] = E
                print(f'diagonal detuning = {E} added')

        # Off-diagonal: Rabi couplings
        for driven in self.driven_transitions():
            i, j = index[driven.transition.lower], index[driven.transition.upper]
            H[i, j] += driven.rabi_frequency / 2
            print(f'rabi frequency = {driven.rabi_frequency}')
            H[j, i] = np.conj(H[i, j])

        return H, index
                
    def build_collapse_operators(self):

        states, index = self.basis()
    
        branching = self.branching_ratios()
    
        c_ops = []
    
        for transition, ratio in branching.items():
    
            gamma = transition.A * ratio
    
            if gamma == 0:
                continue
    
            lower = index[transition.lower]
            upper = index[transition.upper]
    
            c = (
                np.sqrt(gamma)
                * basis(len(states), lower)
                * basis(len(states), upper).dag()
            )
    
            c_ops.append(c)
    
        return c_ops
    
