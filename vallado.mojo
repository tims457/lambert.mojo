from tensor import Tensor
from math import sqrt, pow, cos, abs
from testing import assert_true

from stumpff import c2, c3
from Vector import Vector, norm, mul
from angles import get_transfer_angle

from memory import memset_zero
from random import rand, random_float64

alias PI = 3.141592653589793
alias maxiter = 100


fn lambert(
    mu: Float64,
    # r1: Tensor[DType.float64],
    # r2: Tensor[DType.float64],
    # r1: SIMD[DType.float64, 3],
    # r2: SIMD[DType.float64, 3],
    r1: Vector,
    r2: Vector,
    tof: Float64,
    n: Int64,
    prograde: Bool = True,
) raises -> Tuple[Vector, Vector]:
    # let maxiter: Int64 = 100
    let rtol: Float64 = 1e-7
    var y: Float64

    let r1_norm: Float64 
    r1_norm = norm(r1)
    
    let r2_norm: Float64 
    r2_norm = norm(r2)
    
    let c_norm: Float64 = norm(r2 - r1)

    let dtheta: Float64 = get_transfer_angle(r1, r2, prograde)

    # Compute Vallado's transfer angle parameter
    let A: Float64 
    A = _get_A(r1_norm, r2_norm, dtheta)
    if A == 0.0:
        raise Error("Cannot compute orbit, phase angle is 180 degrees")

    # The initial guess and limits for the bisection method
    var psi: Float64 = 0.0
    var psi_low: Float64 = -4.0 * PI**2
    var psi_up: Float64 = 4.0 * PI**2

    var tof_new: Float64 = 0.0
    var X: Float64 = 0.0

    for i in range( maxiter + 1):
        y = _y_at_psi(psi, r1_norm, r2_norm, A)

        if A > 0.0: 
            while y < 0.0:
                psi_low = psi
                psi = (
                    0.8
                    * (1.0 / c3(psi))
                    * (1.0 - (r1_norm * r2_norm) * sqrt(c2(psi)) / A)
                )
                y = _y_at_psi(psi, r1_norm, r2_norm, A)

        X = _X_at_psi(psi, y)
        tof_new = _tof_vallado(mu, psi, X, A, y)

        # Convergence check
        if abs((tof_new - tof) / tof) < rtol:
            break

        # Bisection check
        let condition: Bool = tof_new <= tof
        psi_low = psi_low + (psi - psi_low) * condition
        psi_up = psi_up + (psi - psi_up) * (not condition)

        psi = (psi_up + psi_low) / 2
    else:
        raise Error("Exceeded maximum number of iterations!")

    let f: Float64 = 1.0 - y / r1_norm
    let g: Float64 = A * sqrt(y / mu)

    let gdot: Float64 = 1.0 - y / r2_norm

    let v1: Vector = (r2 - mul( r1, f)) / g
    let v2: Vector = (mul(r2,  gdot) - r1) / g

    # v1.print()
    # v2.print()
    
    # return [v1, v2 ]
    return (v1, v2 )
    


fn _get_A(r1_norm: Float64, r2_norm: Float64, dtheta: Float64) -> Float64:
    """Computes the value of the A constant.

    Parameters
    ----------
    r1_norm: float
        Initial position vector norm.
    r2_norm: float
        Final position vector norm.
    dtheta: float
        The transfer angle in radians.

    Returns
    -------
    A: float
        The transfer angle parameter.

    """
    let t_m: Float64 = 1.0 if dtheta < PI else -1.0
    let A: Float64 = t_m * sqrt(r1_norm * r2_norm * (1 + cos(dtheta)))
    return A


def _y_at_psi(psi: Float64, r1_norm: Float64, r2_norm: Float64, A: Float64) -> Float64:
    """Evaluates the value of y at given psi.

    Parameters
    ----------
    psi: float
        The free-parameter or independent variable.
    r1_norm: float
        Initial position vector norm.
    r2_norm: float
        Final position vector norm.
    A: float
        The transfer angle parameter.

    Returns
    -------
    y: float
        Auxiliary variable.

    Notes
    -----
    This is equation (7-59) simplified, similarly as made in [1].

    """
    return (r1_norm + r2_norm) + A * (psi * c3(psi) - 1) / sqrt(c2(psi))


fn _tof_vallado(
    mu: Float64, psi: Float64, X: Float64, A: Float64, y: Float64
) -> Float64:
    """Evaluates universal Kepler's equation.

    Parameters
    ----------
    mu: float
        The gravitational parameter.
    psi: float
        The free-parameter or independent variable.
    X: float
        Auxiliary variable.
    A: float
        The transfer angle parameter.
    y: float
        Auxiliary variable.

    Returns
    -------
    tof: float
        The computed time of flight.

    """
    return (X**3 * c3(psi) + A * sqrt(y)) / sqrt(mu)


fn _X_at_psi(psi: Float64, y: Float64) -> Float64:
    """Computes the value of X at given psi.

    Parameters
    ----------
    psi: float
        The free-parameter or independent variable.
    y: float
        Auxiliary variable.

    Returns
    -------
    X: float
        Auxiliary variable.

    """
    return sqrt(y / c2(psi))
