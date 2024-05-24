#!/usr/bin/env python3


import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.run import Run


import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.use("Agg")


from tests.diagnostic import all_timestamps


def density(x):
    return 1.0

def ion_density(x):
    return 1.


def S(x, x0, l):
    return 0.5 * (1 + np.tanh((x - x0) / l))


def bx(x):
    return 0.0


def by(x):
    from pyphare.pharein.global_vars import sim

    L = sim.simulation_domain()[0]
    L1 = L * 0.25
    return np.sin(2 * np.pi * x / L1)


def bz(x):
    return 0.0




def simulation_params(diagdir, **extra):
    params = {
        "interp_order": 1,
        "time_step_nbr": 200,
        "time_step": 0.05,
        "boundary_types": "periodic",
        "cells": 100,
        "dl": 1.0,
        "normalized_c": 1.0,
        "diag_options": {
            "format": "phareh5",
            "options": {"dir": diagdir, "mode": "overwrite"},
        },
    }
    params.update(**extra)
    return params


def config(**options):
    sim = ph.Simulation(**options)
    ph.MaxwellianFluidModel(
        bx=bx, by=by, bz=bz
    )


    timestamps = all_timestamps(sim)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )
    

    return sim



def noRefinement(diagdir):
    return config(**simulation_params(diagdir))



def main():
    Simulator(noRefinement(diagdir="wavePropagation")).run()
    ph.global_vars.sim = None

def fig():
    from pyphare.pharesee.run import Run
    from pyphare.pharesee.hierarchy import get_times_from_h5
    from pyphare.pharesee.hierarchy import flat_finest_field
  
    import os

    mpl.use("Agg")

    run_path = "./wavePropagation"
    run = Run(run_path)

    time = get_times_from_h5(os.path.join(run_path, "EM_B.h5"))

    plot_time = time[int(len(time)-1)]
    first_time = time[1]
    half_time = time[int(len(time)/2)]


    B = run.GetB(first_time)
    by, xby = flat_finest_field(B, "By")

    B1 = run.GetB(half_time)
    by1, xby1 = flat_finest_field(B1, "By")

    B2 = run.GetB(plot_time)
    by2, xby2 = flat_finest_field(B2, "By")


    E0 = run.GetE(first_time)
    Ez0, xEz0 = flat_finest_field(E0, "Ez")

    E1 = run.GetE(half_time)
    Ez1, xEz1 = flat_finest_field(E1, "Ez")

    E2 = run.GetE(plot_time)
    Ez2, xEz2 = flat_finest_field(E2, "Ez")


    fig, axarr = plt.subplots(nrows=3, figsize=(8, 3))



    ax0, ax1, ax2 = axarr

    ax0.plot(xEz0, Ez0, color="k", ls="-")
    ax0.plot(xby, by, color="r", ls="--")

    ax1.plot(xby1, by1, color="r", ls="--")
    ax1.plot(xEz1, Ez1, color="k", ls="-")

    ax2.plot(xby2, by2, color="r", ls="--")
    ax2.plot(xEz2, Ez2, color="k", ls="-")



    fig.savefig("td1d_pic_waveprop.png")


if __name__ == "__main__":
    main()
    fig()
