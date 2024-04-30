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


def S(x, x0, l):
    return 0.5 * (1 + np.tanh((x - x0) / l))


def bx(x):
    return 0.0


def by(x):
    from pyphare.pharein.global_vars import sim

    L = sim.simulation_domain()[0]
    v1 = -1
    v2 = 1.0
    return v1 + (v2 - v1) * (S(x, L * 0.25, 1) - S(x, L * 0.75, 1))


def bz(x):
    return 0.5


def b2(x):
    return bx(x) ** 2 + by(x) ** 2 + bz(x) ** 2


def T(x):
    K = 1
    return 1 / density(x) * (K - b2(x) * 0.5)


def vx(x):
    return 2.0


def vy(x):
    return 0.0


def vz(x):
    return 0.0


def vthx(x):
    return T(x)


def vthy(x):
    return T(x)


def vthz(x):
    return T(x)

mass_electron = 1./10

def vthxe(x):
    return T(x)/mass_electron

def vthye(x):
    return T(x)/mass_electron

def vthze(x):
    return T(x)/mass_electron


vvv = {
    "vbulkx": vx,
    "vbulky": vy,
    "vbulkz": vz,
    "vthx": vthx,
    "vthy": vthy,
    "vthz": vthz,
}

vvv_electrons = {
    "vbulkx": vx,
    "vbulky": vy,
    "vbulkz": vz,
    "vthx": vthxe,
    "vthy": vthye,
    "vthz": vthze,
}

# used to only test on the early particle diagnostic files
particle_diagnostics = {"count": 10, "idx": 0}


def simulation_params(diagdir, **extra):
    params = {
        "interp_order": 1,
        "time_step_nbr": 500,
        "time_step": 0.04,
        "boundary_types": "periodic",
        "cells": 200,
        "hyper_resistivity": 0.01,
        "dl": 1.0,
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
        bx=bx, by=by, bz=bz, 
        protons={"charge": 1., "density": density, **vvv}, 
        electrons={"charge": -1., "mass":mass_electron, "density": density, **vvv_electrons}
    )
    #ph.ElectronModel(closure="isothermal", Te=0.12)

    timestamps = all_timestamps(sim)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )
        
    for quantity in ["density", "bulkVelocity"]:
        ph.FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )
        """
    for pop in sim.model.populations:
        for quantity in ["domain"]:
            ph.ParticleDiagnostics(
                quantity=quantity,
                compute_timestamps=timestamps[: particle_diagnostics["count"] + 1],
                write_timestamps=timestamps[: particle_diagnostics["count"] + 1],
                population_name=pop,
            )
"""

    return sim



def noRefinement(diagdir):
    return config(**simulation_params(diagdir))


def make_figure():
    from scipy.optimize import curve_fit

    rNoRef = Run("./noRefinement")

    plot_time = 11
    v = 2

    BNoRef = rNoRef.GetB(plot_time, merged=True, interp="linear")
    JNoRef = rNoRef.GetJ(plot_time, merged=True, interp="linear")

    xbyNoRef = BNoRef["By"][1][0]
    byNoRef = BNoRef["By"][0](xbyNoRef)
    xjzNoRef = JNoRef["Jz"][1][0]
    jzNoRef = JNoRef["Jz"][0](xjzNoRef)

    fig, axarr = plt.subplots(nrows=3, figsize=(8, 8))

    def S(x, x0, l):
        return 0.5 * (1 + np.tanh((x - x0) / l))

    def by(x):
        L = 200
        v1 = -1
        v2 = 1
        return v1 + (v2 - v1) * (S(x, L * 0.25, 1) - S(x, L * 0.75, 1))

    wT0 = 150.0

    ax0, ax1, ax2 = axarr

    ax0.plot(xbyNoRef, byNoRef, color="k", ls="-")
    ax0.plot(xbyNoRef, by(xbyNoRef), color="darkorange", ls="--")

    ax1.set_xlim((wT0, 195))
    ax1.set_ylim((-1.5, 2))
    ax1.plot(xbyNoRef, byNoRef, color="k", ls="-")

    ax2.plot(xjzNoRef, jzNoRef, color="k")
    ax2.set_xlim((wT0, 195))
    ax2.set_ylim((-1.5, 0.5))


    from pyphare.pharesee.plotting import zoom_effect

    zoom_effect(ax0, ax1, wT0, 195)

    for ax in (ax0, ax1, ax2):
        ax.axvline(wT0 + plot_time * v, color="r")

    fig.savefig("td1d_pic.png")

    # select data around the rightward TD
    idx = np.where((xbyNoRef > 150) & (xbyNoRef < 190))
    xx = xbyNoRef[idx]
    bby = byNoRef[idx]

    # now we will fit by_fit to the data
    # and we expect to find x0=172 and L=1
    # or close enough
    def by_fit(x, x0, L):
        v1 = 1
        v2 = -1
        return v1 + (v2 - v1) * S(x, x0, L)

    popt, pcov = curve_fit(by_fit, xx, bby, p0=(150, 1))
    x0, L = popt

    #if np.abs(L - 1) > 0.5:
    #    raise RuntimeError(f"L (={L}) too far from 1.O")
    #if np.abs(x0 - (150 + plot_time * v)) > 0.5:
    #    raise RuntimeError(f"x0 (={x0}) too far from 172")



def main():
    Simulator(noRefinement(diagdir="noRefinement")).run()
    ph.global_vars.sim = None
    
    make_figure()


if __name__ == "__main__":
    main()
