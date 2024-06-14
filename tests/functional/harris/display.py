from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import get_times_from_h5
from pyphare.pharesee.hierarchy import flat_finest_field
from pyphare.pharesee.hierarchy import hierarchy_from

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np    
import os

mpl.use("Agg")

run_path = "./phare_outputs/test/harris/2d"
run = Run(run_path)

time = get_times_from_h5(os.path.join(run_path, "EM_B.h5"))

plot_time = time[int(len(time)-1)]
first_time = time[0]
half_time = time[int(len(time)/2)]


B = run.GetB(plot_time)
by, xby = flat_finest_field(B, "By")

B1 = run.GetB(half_time)
by1, xby1 = flat_finest_field(B1, "By")


def S(x, x0, l):
    return 0.5 * (1 + np.tanh((x - x0) / l))

def BY(x):
    L = 100
    v1 = -1
    v2 = 1
    return v1 + (v2 - v1) * (S(x, L * 0.25, 1) - S(x, L * 0.75, 1))




from pyphare.pharesee.plotting import finest_field_plot


finest_field_plot(run_path, "Bx", filename="Bx.png", time=plot_time)
finest_field_plot(run_path, "By", filename="By.png", time=plot_time)
finest_field_plot(run_path, "Bz", filename="Bz.png", time=plot_time, title="plot time = {}".format(plot_time))  

finest_field_plot(run_path, "Ex", filename="Ex.png", time=plot_time)
finest_field_plot(run_path, "Ey", filename="Ey.png", time=plot_time)
finest_field_plot(run_path, "Ez", filename="Ez.png", time=plot_time)

finest_field_plot(run_path, "Bx", filename="Bx0.png", time=first_time)
finest_field_plot(run_path, "By", filename="By0.png", time=first_time)
finest_field_plot(run_path, "Bz", filename="Bz0.png", time=first_time)
finest_field_plot(run_path, "Ex", filename="Ex0.png", time=first_time)
finest_field_plot(run_path, "Ey", filename="Ey0.png", time=first_time)
finest_field_plot(run_path, "Ez", filename="Ez0.png", time=first_time)
