from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import get_times_from_h5
from pyphare.pharesee.hierarchy import flat_finest_field
from pyphare.pharesee.hierarchy import hierarchy_from

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np    
import os

mpl.use("Agg")

run_path = "./wavePropagation"
run = Run(run_path)

time = get_times_from_h5(os.path.join(run_path, "EM_B.h5"))

plot_time = time[int(len(time)-1)]
first_time = time[1]
half_time = time[int(len(time)/2)]


B = run.GetB(plot_time)
by, xby = flat_finest_field(B, "By")

B1 = run.GetB(half_time)
by1, xby1 = flat_finest_field(B1, "By")


E0 = run.GetE(first_time)
Ez0, xEz0 = flat_finest_field(E0, "Ez")

E1 = run.GetE(plot_time)
Ez1, xEz1 = flat_finest_field(E1, "Ez")


fig, axarr = plt.subplots(nrows=2, figsize=(8, 3))



ax0, ax1 = axarr

ax0.plot(xEz0, Ez0, color="k", ls="-")
ax0.plot(xby, by, color="r", ls="--")

ax1.plot(xby1, by1, color="r", ls="--")
ax1.plot(xEz1, Ez1, color="k", ls="-")



fig.savefig("td1d_pic_waveprop.png")
