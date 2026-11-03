from dataclasses import dataclass
from typing import Optional
from typing import Literal
import warnings
import numpy as np
from scipy.constants import h, hbar, epsilon_0, c, k, physical_constants

def unit(v, name="vector"):
    v = np.asarray(v, dtype=float)
    if v.shape != (3,):
        raise ValueError(f"{name} must have 3 components.")
    n = np.linalg.norm(v)
    if n == 0:
        raise ValueError(f"{name} must be non-zero.")
    return v / n

def transverse_basis(z):
    # Right-handed orthonormal (e1, e2) with e1 x e2 = z.
    ref = np.array([1.0, 0, 0]) if abs(z[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = ref - ref @ z * z
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(z, e1)
    return e1, e2


def cartesian_to_spherical(E):
    # Cartesian components (Ex, Ey, Ez) -> (c_-1, c_0, c_+1).
    Ex, Ey, Ez = E
    s = np.sqrt(2.0)
    return np.array([(Ex + 1j * Ey) / s, Ez, -(Ex - 1j * Ey) / s])


@dataclass(frozen=True, kw_only=True)
class Laser:
    """
    Properties
    ---------
    frequency: Hz
    intensity: W * m^-2   The peak intensity of the laser.
    linewidth: Hz
    ---------
    All axes are given in Cartesian coordinates.
    If not specified, the default propagation direction is in z axis direction.
    -----------------------------
    Polarization is given in Spherical basis or Jones basis.
    Spherical basis: +/- 1 for sigma +/- polarization, and 0 for pi polarization.
                     The components are in order (sigma-, pi, sigma+).
    Jones basis: Polarization represented in Jones vector.
    -----------------------------


    """

    frequency: float
    intensity: float
    linewidth: float
    polarization_basis: Literal["spherical", "jones"]
    polarization: int | tuple[complex, ...]
    propagation_direction: tuple[float, float, float] = (0.0, 0.0, 1.0)
    quantization_axis: tuple[float, float, float] | None = None

    def __post_init__(self):
        unit(self.propagation_direction, "propagation_direction")

        if self.polarization_basis == "spherical":
            p = self.polarization
            if isinstance(p, int) and not isinstance(p, bool):
                if p not in (-1, 0, 1):
                    raise ValueError("Integer spherical polarization must be -1, 0 or +1.")
                
            elif isinstance(p, tuple) and len(p) == 3:
                if np.linalg.norm(np.asarray(p, dtype=complex)) == 0:
                    raise ValueError("Spherical polarization vector must be non-zero.")
                
            else:
                raise ValueError("Spherical polarization must be q = -1/0/+1 or a 3-component tuple.")
            
            if self.quantization_axis is not None:
                warnings.warn("quantization_axis is ignored for spherical-basis polarization.")

        elif self.polarization_basis == "jones":
            p = self.polarization
            if not (isinstance(p, tuple) and len(p) == 2):
                raise ValueError("Jones polarization must be a 2-component tuple.")
            
            if np.linalg.norm(np.asarray(p, dtype=complex)) == 0:
                raise ValueError("Jones vector must be non-zero.")
            
            if self.quantization_axis is None:
                raise ValueError("Jones polarization requires a quantization_axis.")
            
        else:
            raise ValueError("polarization_basis must be 'spherical' or 'jones'.")

    def electric_field_amplitude(self):
        return np.sqrt(2 * self.intensity / (c * epsilon_0))

    def spherical_components(self, quantization_axis=None):
        #Normalized components (sigma-, pi, sigma+) w.r.t. the quantization axis.
        if self.polarization_basis == "spherical":
            if isinstance(self.polarization, int):
                out = np.zeros(3, dtype=complex)
                out[self.polarization + 1] = 1.0   # -1 -> idx 0, 0 -> 1, +1 -> 2
                return tuple(out)

            v = np.asarray(self.polarization, dtype=complex)
            return tuple(v / np.linalg.norm(v))

        # Jones -> spherical
        if self.polarization_basis == "jones":
            axis = self.quantization_axis

            J = np.asarray(self.polarization, dtype=complex)
            J = J / np.linalg.norm(J)

            k = unit(self.propagation_direction, "propagation_direction")
            e1, e2 = transverse_basis(k)
            E_lab = J[0] * e1 + J[1] * e2

            zq = unit(axis, "quantization_axis")
            xq, yq = transverse_basis(zq)
            E_q = np.array([E_lab @ xq, E_lab @ yq, E_lab @ zq])
            
            return tuple(cartesian_to_spherical(E_q))
            

@dataclass(frozen=True, kw_only=True)  
class LaserPulse:


    
    