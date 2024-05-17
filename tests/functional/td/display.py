from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import get_times_from_h5
from pyphare.pharesee.hierarchy import flat_finest_field
from pyphare.pharesee.hierarchy import hierarchy_from

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np    
import os

mpl.use("Agg")

run_path = "./noRefinement"
run = Run(run_path)

time = get_times_from_h5(os.path.join(run_path, "EM_B.h5"))

plot_time = time[int(len(time)-1)]
first_time = time[1]


B = run.GetB(plot_time)
by, xby = flat_finest_field(B, "By")

E = run.GetE(plot_time)
Ex, xEx = flat_finest_field(E, "Ex")
Ey, xEy = flat_finest_field(E, "Ey")
Ez, xEz = flat_finest_field(E, "Ez")

E0 = run.GetE(first_time)
Ex0, xEx0 = flat_finest_field(E0, "Ex")
Ey0, xEy0 = flat_finest_field(E0, "Ey")
Ez0, xEz0 = flat_finest_field(E0, "Ez")

velocity_ions = run.GetVi(plot_time)
velocity_electrons = run.GetVe(plot_time)

flux_ions = run.GetFlux(plot_time, "protons")
flux_electrons = run.GetElectronFlux(plot_time, "electrons")

ion_density = run.GetNi(plot_time)
ion_mass_density = run.GetMassDensity(plot_time)
electron_density = run.GetNe(plot_time)

ion_particles= run.GetParticles(plot_time, "protons")
electron_particles = run.GetElectronParticles(plot_time, "electrons")

Vix, xVix = flat_finest_field(velocity_ions, "Vx")
Vex, xVex = flat_finest_field(velocity_electrons, "Vx")
Ni, xNi = flat_finest_field(ion_density, "rho")
Ne, xNe = flat_finest_field(electron_density, "rho")
Fix, xFix = flat_finest_field(flux_ions, "protons_Fx")
Fex, xFex = flat_finest_field(flux_electrons, "pop_Fx")

hier = hierarchy_from(h5_filename="./noRefinement/EM_B.h5")
patches = hier.level(0).patches
for patch in patches:
    print(patch.patch_datas.items())
    print(patch.patch_datas.keys())


fig, axarr = plt.subplots(nrows=6, figsize=(8, 10))

def S(x, x0, l):
    return 0.5 * (1 + np.tanh((x - x0) / l))

def BY(x):
    L = 100
    v1 = -1
    v2 = 1
    return v1 + (v2 - v1) * (S(x, L * 0.25, 1) - S(x, L * 0.75, 1))


ax0, ax1, ax2, ax3, ax4, ax5 = axarr

ax0.plot(xby, by, color="k", ls="-")
ax0.plot(xby, BY(xby), color="darkorange", ls="--")

ax1.plot(xNi, Ni, color="k", ls="-")
ax1.plot(xNe, Ne, color="r", ls="-")


ax2.plot(xVix, Vix, color="k", ls="-")
ax2.plot(xVex, Vex, color="r", ls="--")

ion_particles.dist_plot(
                axis=("x", "Vx"),
                ax=ax3,
            )

electron_particles.dist_plot(
                axis=("x", "Vx"),
                ax=ax4,
            )

ax5.plot(xFix, Fix, color="k", ls="-")
ax5.plot(xFex, Fex, color="r", ls="--")

fig.savefig("td1d_pic_info.png")

fig, axarr = plt.subplots(nrows=3, figsize=(8, 10))
ax0, ax1, ax2 = axarr

ax0.plot(xEx, Ex, color="k", ls="-")
ax0.plot(xEx0, Ex0, color="r", ls="--")
ax1.plot(xEy, Ey, color="k", ls="-")
ax1.plot(xEy0, Ey0, color="r", ls="--")
ax2.plot(xEz, Ez, color="k", ls="-")
ax2.plot(xEz0, Ez0, color="r", ls="--")

fig.savefig("td1d_pic_E.png")
