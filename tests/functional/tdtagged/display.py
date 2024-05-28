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
half_time = time[int(len(time)/2)]



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

fig.savefig("td1d_E.png")
