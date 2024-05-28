from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import get_times_from_h5
from pyphare.pharesee.hierarchy import flat_finest_field
from pyphare.pharesee.hierarchy import hierarchy_from

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np    
import os

mpl.use("Agg")

run_path = "./noFlowPic"
run = Run(run_path)

time = get_times_from_h5(os.path.join(run_path, "EM_B.h5"))

plot_time = time[int(len(time)-1)]
first_time = time[1]
half_time = time[int(len(time)/2)]


B = run.GetB(plot_time)
by, xby = flat_finest_field(B, "By")

B1 = run.GetB(half_time)
by1, xby1 = flat_finest_field(B1, "By")


velocity_ions = run.GetVi(plot_time)
velocity_electrons = run.GetVe(plot_time)

flux_ions = run.GetFlux(plot_time, "protons")
flux_electrons = run.GetElectronFlux(plot_time, "electrons")

ion_density = run.GetNi(plot_time)
electron_density = run.GetNe(plot_time)

ion_particles= run.GetParticles(plot_time, "protons")
electron_particles = run.GetElectronParticles(plot_time, "electrons")

Vix, xVix = flat_finest_field(velocity_ions, "Vx")
Vex, xVex = flat_finest_field(velocity_electrons, "Vx")
Ni, xNi = flat_finest_field(ion_density, "rho")
Ne, xNe = flat_finest_field(electron_density, "rho")
Fix, xFix = flat_finest_field(flux_ions, "protons_Fx")
Fex, xFex = flat_finest_field(flux_electrons, "pop_Fx")

hier = hierarchy_from(h5_filename="./noFlowPic/EM_B.h5")
"""
patches = hier.level(0).patches
for patch in patches:
    print(patch.patch_datas.items())
    print(patch.patch_datas.keys())
"""

fig, axarr = plt.subplots(nrows=5, figsize=(8, 10))

def S(x, x0, l):
    return 0.5 * (1 + np.tanh((x - x0) / l))

def BY(x):
    L = 100
    v1 = -1
    v2 = 1
    return v1 + (v2 - v1) * (S(x, L * 0.25, 1) - S(x, L * 0.75, 1))


ax0, ax1, ax2, ax3, ax4 = axarr

ax0.plot(xby, BY(xby), color="darkorange", ls="--", label="t=0")
ax0.plot(xby, by, color="k", ls="-")
ax0.legend()

ax4.set_xlabel("x")
ax0.set_ylabel(r"$B_y$")

ax1.plot(xNi, Ni, color="k", ls="-", label="Ions")
ax1.plot(xNe, Ne, color="r", ls="-", label="Electrons")
ax1.set_ylabel("Density")

ax2.plot(xVix, Vix, color="k", ls="-", label="Ions")
ax2.plot(xVex, Vex, color="r", ls="-", label="Electrons", zorder=1)
ax2.set_ylabel("Bulk Velocity")

ion_particles.dist_plot(
                axis=("x", "Vx"),
                ax=ax3,
            )
ax3.set_ylabel(r"Ion $V_x$")
electron_particles.dist_plot(
                axis=("x", "Vx"),
                ax=ax4,
            )
ax4.set_ylabel(r"Electron $V_x$")
"""
ax5.plot(xFix, Fix, color="k", ls="-")
ax5.plot(xFex, Fex, color="r", ls="--")
ax5.set_ylabel(r"$F_x$")
"""
ax1.legend(loc="upper left")
ax2.legend(loc="upper left")
#fig.tight_layout()

fig.savefig("td1d_pic_noflow.png")

E = run.GetE(plot_time)
Ex, xEx = flat_finest_field(E, "Ex")
Ey, xEy = flat_finest_field(E, "Ey")
Ez, xEz = flat_finest_field(E, "Ez")

E0 = run.GetE(first_time)
Ex0, xEx0 = flat_finest_field(E0, "Ex")
Ey0, xEy0 = flat_finest_field(E0, "Ey")
Ez0, xEz0 = flat_finest_field(E0, "Ez")

E1 = run.GetE(half_time)
Ex1, xEx1 = flat_finest_field(E1, "Ex")
Ey1, xEy1 = flat_finest_field(E1, "Ey")
Ez1, xEz1 = flat_finest_field(E1, "Ez")

fig, axarr = plt.subplots(nrows=3, figsize=(8, 10))
ax0, ax1, ax2 = axarr

ax0.plot(xEx, Ex, color="k", ls="-")
ax0.plot(xEx0, Ex0, color="r", ls="--")
ax0.plot(xEx1, Ex1, color="b", ls="--")
ax1.plot(xEy, Ey, color="k", ls="-")
ax1.plot(xEy0, Ey0, color="r", ls="--")
ax1.plot(xEy1, Ey1, color="b", ls="--")
ax2.plot(xEz, Ez, color="k", ls="-")
ax2.plot(xEz0, Ez0, color="r", ls="--")
ax2.plot(xEz1, Ez1, color="b", ls="--")

fig.savefig("td1d_pic_E_noflow.png")
