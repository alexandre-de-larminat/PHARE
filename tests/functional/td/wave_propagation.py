#!/usr/bin/env python3


import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.run import Run

from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import get_times_from_h5
from pyphare.pharesee.hierarchy import flat_finest_field
  
import os

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import matplotlib.animation as animation

mpl.use("Agg")

from tests.diagnostic import all_timestamps


def bx(x):
    return 0.0

def by(x):
    from pyphare.pharein.global_vars import sim

    L = sim.simulation_domain()[0]
    L1 = L #* 0.5
    return np.sin(2 * np.pi * x / L1)

def bz(x):
    return 0.


def simulation_params(diagdir, **extra):
    params = {
        "interp_order": 1,
        "time_step_nbr": 100,
        "time_step": 1.,
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


def wave(x, a0, k, phi):
    return a0 * np.cos(k * x + phi)


def phase_speed(run_path, ampl, xmax):
    from scipy.signal import medfilt
    from scipy.optimize import curve_fit
    import os

    time = get_times_from_h5(os.path.join(run_path, "EM_B.h5"))
    r = Run(run_path)
    phase = np.zeros_like(time)
    amplitude = np.zeros_like(time)
    wave_vec = np.zeros_like(time)

    for it, t in enumerate(time):
        B = r.GetB(t)
        by, xby = flat_finest_field(B, "By")
        a, k, phi = curve_fit(wave, xby, by, p0=(ampl, 2 * np.pi / xmax, 0))[0]
        phase[it] = phi
        amplitude[it] = a
        wave_vec[it] = k

    vphi = medfilt(np.gradient(phase, time) / wave_vec, kernel_size=7)
    return vphi, time, phase, amplitude, wave_vec


def fig():

    run_path = "./wavePropagation"
    run = Run(run_path)

    time = get_times_from_h5(os.path.join(run_path, "EM_B.h5"))

    B = run.GetB(time[0])
    by, xby = flat_finest_field(B, "By")

    E = run.GetE(time[0])
    Ez, xEz = flat_finest_field(E, "Ez")


    fig, axarr = plt.subplots(nrows=2, figsize=(6, 4))

    ax0, ax1 = axarr


    Bline, = ax0.plot(xby, by, color="r", ls="-")
    ax0.set_ylabel(r"$B_y$")
    Eline, = ax1.plot(xEz, Ez, color="k", ls="-")
    ax1.set_ylabel(r"$E_z$")
    ax1.set_xlabel("x")

    vphi, t, phi, a, k = phase_speed(run_path, 1, 100)

    ax0.set_title(r"$V_\phi = {:6.4f}$".format(vphi.mean()))

    def update(frame):
        Bline.set_ydata(flat_finest_field(run.GetB(time[frame]), "By")[0])
        Eline.set_ydata(flat_finest_field(run.GetE(time[frame]), "Ez")[0])
        return Bline, Eline

    anim = animation.FuncAnimation(fig, update, frames=len(time), interval=100, blit=True)
    anim.save("td1d_pic_waveprop.gif", writer="imagemagick", fps=10)

    ax0.plot(xby, flat_finest_field(run.GetB(time[int(len(time)/4)]), "By")[0], color="r", ls="-")
    ax1.plot(xEz, flat_finest_field(run.GetE(time[int(len(time)/4)]), "Ez")[0], color="k", ls="-")
    fig.savefig("td1d_pic_waveprop.png")


if __name__ == "__main__":
    main()
    fig()
