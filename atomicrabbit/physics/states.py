from dataclasses import dataclass
from typing import Optional
import numpy as np
from scipy.constants import h, hbar, epsilon_0, c, k, physical_constants

'''
example use:
    
d52 = FineState(...)

F3 = HyperfineState(
    parent=d52,
    F=3
)

mF2 = ZeemanState(
    parent=F3,
    m=2
)

'''
mu_B = physical_constants["Bohr magneton"][0]

def joule_to_cm(energy_joule):
    energy_cm = energy_joule / (h * c * 100)
    return energy_cm
    
def cm_to_joule(energy_cm):
    energy_joule = energy_cm * h * c * 100
    return energy_joule
    
@dataclass(frozen=True)
class FineState:
    configuration: str
    L: float
    S: float
    J: float
    energy_cm: float
    parity: str
    gJ: Optional[float] = None
    A_hyperfine: Optional[float] = None
    B_hyperfine: Optional[float] = None

    
    @property
    def g(self):
       return self.gJ
   
    @property
    def degeneracy(self):
        return 2 * self.J +1
    
    @property
    def term_symbol(self):
        L_letters = {
            0: "S",
            1: "P",
            2: "D",
            3: "F",
            4: "G",
        }

        multiplicity = int(2 * self.S + 1)
        L_letter = L_letters.get(self.L, "?")

        return f"{multiplicity}{L_letter}_{self.J}"
    
    @property
    def spectroscopic_name(self):
        return f"{self.configuration} {self.term_symbol}"


@dataclass(frozen=True)
class HyperfineState:
    parent: FineState
    F: float
    I: Optional[float] = None

    @property
    def J(self):
        return self.parent.J
    
    @property
    def A_hyperfine(self):
        return self.parent.A_hyperfine
    
    @property
    def B_hyperfine(self):
        return self.parent.B_hyperfine
    
    @property
    def g(self):
        return self.gF
    
    @property
    def hyperfine_shift(self):
        '''
        calculate hyperfine shift
        '''
        J = self.J
        I = self.I
        F = self.F
        A = self.A_hyperfine
        B = self.B_hyperfine

        C = F * (F + 1) - I * (I + 1) - J * (J + 1)
    
        # Magnetic dipole term
        delta_E = 0.5 * A * C
    
        # Electric quadrupole term
        if B != None and I > 0.5 and J > 0.5:
            numerator = (
                1.5 * C * (C + 1)
                - 2 * I * (I + 1) * J * (J + 1)
            )
    
            denominator = (
                I * (2 * I - 1)
                * J * (2 * J - 1)
            )
    
            delta_E += (B / 4) * numerator / denominator
    
        return delta_E
    
    @property
    def energy_cm(self):
        return self.parent.energy_cm + self.hyperfine_shift

    
@dataclass(frozen=True)
class ZeemanState:
    parent: FineState | HyperfineState
    B : float   # Magnetic field strength
    mF: Optional[float] = None
    gF: Optional[float] = None

    @property
    def zeeman_shift(self):
        delta_E = mu_B * self.gF * self.mF * self.B
        return joule_to_cm(delta_E)
    
    @property
    def energy_cm(self):
        return self.parent.energy_cm + self.zeeman_shift

def allowed_F(I, J):
    # Finding all allowed F values(F between |I-J| and I+J)
    F_min2 = round(2 * abs(I - J))
    F_max2 = round(2 * (I + J))
    
    return [f2 / 2 for f2 in range(F_min2, F_max2 + 1, 2)]

def hyperfine_split(state, I=None):
    '''
    I: The nuclear spin(induces hyperfine splitting)
    '''
    state_list = []
    
    if I == None:
        print(f'No hyperfine structure for {state.spectroscopic_name}')
        return
    
    F_list = allowed_F(I, state.J)
    
    for F in F_list:     
        if state.A_hyperfine is None and state.B_hyperfine is None:
            print('Hyperfine constants A and B missing when creating the fine structure')
        
        state_list.append(HyperfineState(
                parent = state,
                F = F,
                I = I
                ))
            
    return state_list

def allowed_mF(F):
    # Finding all allowed mF values between -F and F
    F2 = round(2 * F)
    return [mF2 / 2 for mF2 in range(-F2, F2 + 1, 2)]

def zeeman_split(state, B=None, gF=None):
    state_list = []
    
    if B == None:
        print('No Zeeman splitting due to zero magnetic field')
        return
    
    if isinstance(state, FineState):
        mF_list = allowed_mF(state.J)
        
    elif isinstance(state, HyperfineState):
        mF_list = allowed_mF(state.F)
        
    else:
        print('Only FineState and HyperfineState instances can have zeeman splitting')
        
    for mF in mF_list:
        if gF is None:
            print('Please enter gF to calculated Zeeman splitting')
        
        state_list.append(ZeemanState(
            parent = state,
            B = B,
            mF = mF,
            gF = gF))
        
    return state_list


'''
Functions to split a list of states in one go
'''
def flatten(lst):
    result = []
    for item in lst:
        if isinstance(item, list):
            result.extend(flatten(item))
        else:
            result.append(item)
    return result

def batch_hyperfine_split(state_list, I=None):
    state_list = flatten(state_list)
    hyperfine_list = []
    for state in state_list:
        hyperfine_state = hyperfine_split(state, I)
        hyperfine_list.append(hyperfine_state)
    
    hyperfine_list = flatten(hyperfine_list)
    return hyperfine_list

def batch_zeeman_split(state_list, B=None, gF=None):
    state_list = flatten(state_list)
    zeeman_list = []
    for state in state_list:
        zeeman_state = zeeman_split(state, B, gF)
        zeeman_list.append(zeeman_state)
        
    zeeman_list = flatten(zeeman_list)
    return zeeman_list
    
        
        
            
            
            
            