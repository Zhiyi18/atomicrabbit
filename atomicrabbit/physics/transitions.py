from dataclasses import dataclass
import numpy as np
from scipy.constants import h, hbar, epsilon_0, c, k, physical_constants

@dataclass(frozen=True)
class Transition:
    lower: str
    upper: str
    A: float