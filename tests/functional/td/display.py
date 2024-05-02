from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import get_times_from_h5
from pyphare.pharesee.hierarchy import flat_finest_field

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np    
import os

mpl.use("Agg")

run_path = "./noRefinement"
run = Run(run_path)

time = get_times_from_h5(os.path.join(run_path, "EM_B.h5"))

plot_time = time[int(len(time) / 2)]


B = run.GetB(plot_time)
by, xby = flat_finest_field(B, "By")

velocity_ions = run.GetVi(plot_time)
velocity_electrons = run.GetVe(plot_time)
print(velocity_ions)
print(velocity_electrons)

ion_density = run.GetNi(plot_time)
electron_density = run.GetNe(plot_time)
print(ion_density)
print(electron_density)

Vix, xVix = flat_finest_field(velocity_ions, "Vx")
Vex, xVex = flat_finest_field(velocity_electrons, "Vx")
Ni, xNi = flat_finest_field(ion_density, "rho")
Ne, xNe = flat_finest_field(electron_density, "rho")


print(velocity_electrons)

fig, axarr = plt.subplots(nrows=4, figsize=(8, 10))

def S(x, x0, l):
    return 0.5 * (1 + np.tanh((x - x0) / l))

def BY(x):
    L = 100
    v1 = -1
    v2 = 1
    return v1 + (v2 - v1) * (S(x, L * 0.25, 1) - S(x, L * 0.75, 1))


ax0, ax1, ax2, ax3 = axarr

ax0.plot(xby, by, color="k", ls="-")
ax0.plot(xby, BY(xby), color="darkorange", ls="--")

ax1.plot(xNi, Ni, color="k", ls="-")
ax1.plot(xNe, Ne, color="r", ls="-")


ax2.plot(xVix, Vix, color="k", ls="-")
ax2.plot(xVex, Vex, color="r", ls="--")

ax3.plot(xVex, Vex, color="k", ls="-")



fig.savefig("td1d_pic_info.png")
