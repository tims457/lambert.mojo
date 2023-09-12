from math import acos
from Vector import Vector, cross, dot, norm, sign

alias PI = 3.141592653589793


fn get_transfer_angle(r1:Vector, r2:Vector, prograge:Bool) raises -> Float64:

    # Solves for the transfer angle being known the sense of rotation.

    # Parameters
    # ----------
    # r1: np.array
    #     Initial position vector.
    # r2: np.array
    #     Final position vector.
    # prograde: bool
    #     If True, it assumes prograde motion, otherwise assumes retrograde.

    # Returns
    # -------
    # dtheta: float
    #     Transfer angle in radians.


    # Check if both position vectors are collinear. If so, check if the transfer
    # angle is 0 or pi.


    if cross(r1, r2) == 0.0:
        if sign(r1) == sign(r2):
            return 0.0
        else:
            return PI

    # Solve for a unitary vector normal to the vector plane. Its direction and
    # sense the one given by the cross product (right-hand) from r1 to r2.
    let h: Vector = cross(r1, r2) / norm(cross(r1, r2))

    # Compute the projection of the normal vector onto the reference plane.
    let zhat = Vector(3)
    zhat[0] = 0; zhat[1] = 0; zhat[2] = 1
    
    let alpha: Float64 = dot(zhat, h)

    let r1_norm: Float64 = norm(r1)
    let r2_norm: Float64 = norm(r2)

    let theta0: Float64 = acos(dot(r1, r2) / (r1_norm * r2_norm))

    if prograge:
        if alpha > 0:
            return theta0
        else:
            return 2 * PI - theta0
    else:
        if alpha < 0:
            return theta0
        else:
            return 2 * PI - theta0
        
