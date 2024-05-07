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


B = run.GetB(plot_time)
by, xby = flat_finest_field(B, "By")

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
Vex, xVex = flat_finest_field(velocity_electrons, "Vex")
Ni, xNi = flat_finest_field(ion_density, "rho")
Nmi, xNmi = flat_finest_field(ion_mass_density, "rho")
Ne, xNe = flat_finest_field(electron_density, "rho")
Fix, xFix = flat_finest_field(flux_ions, "protons_Fx")
Fex, xFex = flat_finest_field(flux_electrons, "pop_Fx")

hier = hierarchy_from(h5_filename="./noRefinement/EM_B.h5")
patches = hier.level(0).patches
for patch in patches:
    print(patch.patch_datas.items())
    print(patch.patch_datas.keys())


fig, axarr = plt.subplots(nrows=5, figsize=(8, 10))

def S(x, x0, l):
    return 0.5 * (1 + np.tanh((x - x0) / l))

def BY(x):
    L = 100
    v1 = -1
    v2 = 1
    return v1 + (v2 - v1) * (S(x, L * 0.25, 1) - S(x, L * 0.75, 1))


ax0, ax1, ax2, ax3, ax4 = axarr

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

ax4.plot(xFix, Fix, color="k", ls="-")
ax4.plot(xFex, Fex, color="r", ls="--")

fig.savefig("td1d_pic_info.png")
