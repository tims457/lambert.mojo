import time
import numpy as np
# from lamberthub import izzo2015 as solver
from vallado import vallado2013 as solver


r1 = np.array([7250.4622962932, 3380.946093925595, 0])
r2 = np.array([0, 8500, 0])
tof = 10*3600

start_time = time.time()
v1, v2 = solver(398600.4415, r1, r2, tof, prograde=True, low_path=True,
                     maxiter=100, atol=1e-5, rtol=1e-7, full_output=False)
end_time = time.time()
# print(v1, v2)
print(f"Elapsed time: {(end_time - start_time):.2e} seconds.")