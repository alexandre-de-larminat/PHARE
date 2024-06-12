#!/usr/bin/env python3

import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.hierarchy import get_times_from_h5
from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import flat_finest_field

from tests.diagnostic import all_timestamps

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.use("Agg")


####################################################################
#
#                     Simulation configuration
#
####################################################################
v = .00
b = .01
cell_nbr = 100
dl = .1
x_length = cell_nbr * dl

def config():
    # configure the simulation

    sim = ph.Simulation(
        time_step_nbr=100,  # number of time steps (not specified if time_step and final_time provided)
        final_time=10,  # simulation final time (not specified if time_step and time_step_nbr provided)
        boundary_types="periodic",  # boundary condition, string or tuple, length == len(cell) == len(dl)
        cells=cell_nbr,  # integer or tuple length == dimension
        dl=dl,  # mesh size of the root level, float or tuple
        normalized_c=1.,
        refinement_boxes={},
        diag_options={
            "format": "phareh5",
            "options": {"dir": "alt", "mode": "overwrite"},
        },
    )

    def density(x):
        return 0.1

    def by(x):
        L = sim.simulation_domain()
        return b * np.cos(2 * np.pi * x / L[0])

    def bz(x):
        L = sim.simulation_domain()
        return b * np.sin(2 * np.pi * x / L[0])

    def bx(x):
        return 0.01

    def vx(x):
        return 0.0


    def vy(x):
        L = sim.simulation_domain()
        return 0#v * np.cos(2 * np.pi * x / L[0])

    def vz(x):
        L = sim.simulation_domain()
        return 0#v * np.sin(2 * np.pi * x / L[0])
    
    def vye(x):
        L = sim.simulation_domain()
        return -vy(x) 

    def vze(x):
        L = sim.simulation_domain()
        return -vz(x)

    def vth(x):
        return .0
    
    def vthe(x):
        return vth(x)#/electron_mass

    electron_mass = 1

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vth,
        "vthy": vth,
        "vthz": vth,
    }

    vvv_e = {
        "vbulkx": vx,
        "vbulky": vye,
        "vbulkz": vze,
        "vthx": vthe,
        "vthy": vthe,
        "vthz": vthe,
    }

    ph.MaxwellianFluidModel(
        bx=bx, by=by, bz=bz, protons={"charge": 1, "density": density, **vvv},
        electrons={"charge": -1, "mass":electron_mass, "density": density, **vvv_e}
    )

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

    return sim


####################################################################
#                      post processing
####################################################################


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


def main():
    from pyphare.cpp import cpp_lib

    cpp = cpp_lib()

    from pyphare.pharesee.run import Run
    from pyphare.pharesee.hierarchy import flat_finest_field

    sim = config()
    Simulator(sim).run()

def figure():
        vphi, t, phi, a, k = phase_speed(".", b, x_length)

        r = Run(".")
        t = get_times_from_h5("EM_B.h5")
        fig, ax = plt.subplots(figsize=(9, 5), nrows=1)

        B = r.GetB(t[int(len(t) / 2)])
        by, xby = flat_finest_field(B, "By")
        ax.plot(xby, by, label="t = 500", alpha=0.6)

 

        B = r.GetB(t[-1])
        by, xby = flat_finest_field(B, "By")
        ax.plot(xby, by, label="t = 1000", alpha=0.6)
        ax.plot(
            xby,
            wave(xby, b, 2 * np.pi / x_length, 2 * np.pi / x_length * x_length/2),
            color="k",
            ls="--",
            label="T=500 (theory)",
        )

        B = r.GetB(t[0])
        by, xby = flat_finest_field(B, "By")
        ax.plot(xby, by, label="t = 0", color="k")

        ax.set_xlabel("x")
        ax.set_ylabel(r"$B_y$")
        ax.legend(ncol=4, loc="upper center")
        ax.set_ylim((-1.2*b, 1.3*b))
        ax.set_title(r"$V_\phi = {:6.4f}$".format(vphi.mean()))

        fig.tight_layout()

        fig.savefig("alfven_wave.png", dpi=200)

        assert np.mean(np.abs(vphi - 1) < 5e-2)


if __name__ == "__main__":
    main()
    #figure()
