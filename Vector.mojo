from math import sqrt, pow
from testing import assert_true

from stumpff import c2, c3


from memory import memset_zero
from random import rand, random_float64


struct Vector:
    var data: DTypePointer[DType.float64]
    var rows: Int

    fn __init__(inout self, rows: Int):
        self.data = DTypePointer[DType.float64].alloc(rows)
        rand(self.data, rows)
        self.rows = rows

    fn __del__(owned self):
        self.data.free()

    fn __eq__(self, rhs: Vector) raises -> Bool:
        if self.rows != rhs.rows:
            return False
        for i in range(self.rows):
            if self[i] != rhs[i]:
                return False
        return True
    
    fn __eq__(self, rhs: Float64) raises -> Bool:
        for i in range(self.rows):
            if self[i] != rhs:
                return False
        return True


    fn zero(inout self):
        memset_zero(self.data, self.rows)

    @always_inline
    fn __getitem__(self, i: Int) raises -> Float64:
        if i < 0 or i >= self.rows:
            raise Error("Out of bounds")
        return self.load[1](i)

    @always_inline
    fn __setitem__(self, i: Int, val: Float64) raises:
        if i < 0 or i >= self.rows:
            raise Error("Out of bounds")
        return self.store[1](i, val)

    fn __copyinit__(inout self, other: Self) -> None:
            self.__init__(other.rows)
            memcpy[DType.float64](self.data, other.data, self.rows)


    fn __sub__(self, rhs: Vector) raises -> Vector:
        let z = Vector(self.rows)
        for i in range(self.rows):
            z[i] = self[i] - rhs[i]
    
        return z

    fn __truediv__(self, rhs: Float64) raises -> Vector:
        let z = Vector(self.rows)
        for i in range(self.rows):
            z[i] = self[i] / rhs
    
        return z

    fn __mul__(self, rhs: Float64) raises -> Vector:
        let z = Vector(self.rows)
        for i in range(self.rows):
            z[i] = self[i] * rhs
    
        return z

    

    @always_inline
    fn load[nelts: Int](self, i: Int) -> SIMD[DType.float64, nelts]:
        return self.data.simd_load[nelts](i)

    @always_inline
    fn store[nelts: Int](self, i: Int, val: SIMD[DType.float64, nelts]):
        return self.data.simd_store[nelts](i, val)

    fn print(self) raises ->None:
        for i in range(self.rows):
            print(self[i])


fn norm(x: Vector) raises -> Float64:
    var s: Float64 = 0.0
    for i in range(x.rows):
        s += pow(x[i],2)
    return sqrt(s)


fn cross(a:Vector, b:Vector) raises -> Vector:

    # Computes the cross product between two vectors.

    # Parameters
    # ----------
    # r1: Vector
    #     First vector.
    # r2: Vector
    #     Second vector.

    # Returns
    # -------
    # r3: Vector
    #     Cross product between r1 and r2.
    let c = Vector(3)

    c[0] = a[1]*b[2] - a[2]*b[1]
    c[1] = a[2]*b[0] - a[0]*b[2]
    c[2] = a[0]*b[1] - a[1]*b[0]
    return c

fn dot(a:Vector, b:Vector) raises -> Float64:
    # Computes the dot product between two vectors.

    # Parameters
    # ----------
    # r1: Vector
    #     First vector.
    # r2: Vector
    #     Second vector.

    # Returns
    # -------
    # r3: Float64
    #     Dot product between r1 and r2.
    var c: Float64 = 0.0
    for i in range(a.rows):
        c += a[i]*b[i]
    return c


fn sign(x:Vector) raises -> Vector:
    let y = Vector(x.rows)
    for i in range(x.rows):
        if x[i] > 0.0:
            y[i] = 1.0
        elif x[i] < 0.0:
            y[i] = -1.0
        else:
            y[i] = 0.0
    return y

fn mul(x:Vector, y:Float64) raises -> Vector:
    let z = Vector(x.rows)
    for i in range(x.rows):
        z[i] = x[i] * y
    return z

