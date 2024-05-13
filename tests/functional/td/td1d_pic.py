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
    return 0.1


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

def vthe(x):
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
    "vthx": vthe,
    "vthy": vthe,
    "vthz": vthe,
}

steps = 250
# used to only test on the early particle diagnostic files
particle_diagnostics = {"count": steps, "idx": 0}


def simulation_params(diagdir, **extra):
    params = {
        "interp_order": 1,
        "time_step_nbr": steps,
        "time_step": 0.04,
        "boundary_types": "periodic",
        "cells": 100,
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
        protons={"charge": 1., "mass": 1., "density": ion_density, **vvv}, 
        electrons={"charge": -1., "mass":mass_electron, "density": density, **vvv_electrons}
    )


    timestamps = all_timestamps(sim)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )
        
    for quantity in ["density", "mass_density", "bulkVelocity"]:
        ph.FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )
        
    for quantity in ["density", "electron_bulkVelocity"]:
        ph.FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
            population_name="all_electrons",
        )


    for pop in sim.model.populations:
        for quantity in ["density", "flux"]:
            ph.FluidDiagnostics(
                quantity=quantity,
                write_timestamps=timestamps,
                compute_timestamps=timestamps,
                population_name=pop,
            )
    
    for pop in sim.model.populations:
        for quantity in ["domain"]:
            ph.ParticleDiagnostics(
                quantity=quantity,
                compute_timestamps=timestamps[: particle_diagnostics["count"] + 1],
                write_timestamps=timestamps[: particle_diagnostics["count"] + 1],
                population_name=pop,
            )
    

    return sim



def noRefinement(diagdir):
    return config(**simulation_params(diagdir))



def main():
    Simulator(noRefinement(diagdir="noRefinement")).run()
    ph.global_vars.sim = None


if __name__ == "__main__":
    main()
