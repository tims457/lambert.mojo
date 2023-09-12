
import time
from tensor import Tensor
from math import tgamma, lgamma, exp, log, sqrt

from python import Python


from stumpff import c2, c3
from vallado import lambert
from Vector import Vector


def main():
    let np = Python.import_module("numpy")
    
    let r1p = np.array([7250.4622962932, 3380.946093925595, 0])
    let r2p = np.array([0, 8500, 0])
    let r1 = Vector(3)
    let r2 = Vector(3)

    let mu = 398600.4415
    let tof = 36000

    for i in range(3):
        r1[i] = r1p[i].to_float64()
        r2[i] = r2p[i].to_float64()

    let start_time = time.now()

    lambert(mu, r1, r2, tof, 0, True)

    print("Elapsed time", (time.now() - start_time)/1e9, "sec")



