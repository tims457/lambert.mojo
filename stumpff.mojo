from math import sin, cos, sqrt, sinh, cosh, tgamma


fn c2(psi: Float64) -> Float64:
    # Second Stumpff function.
    # For positive arguments:
    # .. math::
    #     c_2(\psi) = \frac{1 - \cos{\sqrt{\psi}}}{\psi}

    var res: Float64
    let eps: Float64 = 1.0

    if psi > eps:
        res = (1 - cos(sqrt(psi))) / psi
    elif psi < -eps:
        res = (cosh(sqrt(-psi)) - 1) / (-psi)
    else:
        res = 1.0 / 2.0
        var delta: Float64 = (-psi) / tgamma(SIMD[DType.float64, 1](5))
        var k: Float64 = 1
        while res + delta != res:
            res = res + delta
            k += 1
            delta = (-psi) ** k / tgamma(SIMD[DType.float64, 1](2 * k + 2 + 1))
    return res


fn c3(psi: Float64) -> Float64:
    # Third Stumpff function.
    # For positive arguments:
    # .. math::
    #     c_3(\psi) = \frac{\sqrt{\psi} - \sin{\sqrt{\psi}}}{\sqrt{\psi^3}}

    var res: Float64
    let eps: Float64 = 1.0

    if psi > eps:
        res = (sqrt(psi) - sin(sqrt(psi))) / (psi * sqrt(psi))
    elif psi < -eps:
        res = (sinh(sqrt(-psi)) - sqrt(-psi)) / (-psi * sqrt(-psi))
    else:
        res = 1.0 / 6.0
        var delta: Float64 = (-psi) / tgamma(SIMD[DType.float64, 1](6))
        var k: Float64 = 1
        while res + delta != res:
            res = res + delta
            k += 1
            delta = (-psi) ** k / tgamma(SIMD[DType.float64, 1](2 * k + 3 + 1))
    return res