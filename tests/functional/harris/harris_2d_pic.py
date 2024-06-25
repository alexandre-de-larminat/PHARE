#!/usr/bin/env python3

import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator, startMPI

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use("Agg")

from pyphare.cpp import cpp_lib

cpp = cpp_lib()
startMPI()

diag_outputs = "phare_outputs/test/harris/2d"
from datetime import datetime


def config():
    sim = ph.Simulation(
        time_step_nbr=4000,
        time_step=0.001,
        # boundary_types="periodic",
        cells=(200, 400),
        dl=(0.2, 0.2),
        refinement_boxes={},
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_outputs, "mode": "overwrite"},
        },
        strict=False,
    )

    def density(x, y):
        L = sim.simulation_domain()[1]
        return (
            0.2
            + 1.0 / np.cosh((y - L * 0.3) / 0.5) ** 2
            + 1.0 / np.cosh((y - L * 0.7) / 0.5) ** 2
        )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def by(x, y):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        w1 = 0.2
        w2 = 1.0
        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly
        w3 = np.exp(-(x0 * x0 + y1 * y1) / (w2 * w2))
        w4 = np.exp(-(x0 * x0 + y2 * y2) / (w2 * w2))
        w5 = 2.0 * w1 / w2
        return (w5 * x0 * w3) + (-w5 * x0 * w4)

    def bx(x, y):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        w1 = 0.2
        w2 = 1.0
        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly
        w3 = np.exp(-(x0 * x0 + y1 * y1) / (w2 * w2))
        w4 = np.exp(-(x0 * x0 + y2 * y2) / (w2 * w2))
        w5 = 2.0 * w1 / w2
        v1 = -1
        v2 = 1.0
        return (
            v1
            + (v2 - v1) * (S(y, Ly * 0.3, 0.5) - S(y, Ly * 0.7, 0.5))
            + (-w5 * y1 * w3)
            + (+w5 * y2 * w4)
        )

    def bz(x, y):
        return 0.0

    def b2(x, y):
        return bx(x, y) ** 2 + by(x, y) ** 2 + bz(x, y) ** 2

    def Ti(x, y):
        K = 1
        temp = 1.0 / density(x, y) * (K - b2(x, y) * 0.5)
        assert np.all(temp > 0)
        return temp
    
    def Te(x, y):
        K = 1
        temp = Ti(x, y) * 5
        assert np.all(temp > 0)
        return temp

    def vx(x, y):
        return 0.0

    def vy(x, y):
        return 0.0

    def vz(x, y):
        return 0.0
    
    def dB_y_dx(x, y):
        Lx = sim.simulation_domain()[0]  # Assuming sim is the simulation object
        Ly = sim.simulation_domain()[1]
        w1 = 0.2
        w2 = 1.0
        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly
        w3 = np.exp(-(x0 * x0 + y1 * y1) / (w2 * w2))
        w4 = np.exp(-(x0 * x0 + y2 * y2) / (w2 * w2))
        w5 = 2.0 * w1 / w2

        # Partial derivatives
        dx0_dx = 1.0
        dw3_dx0 = -2.0 * x0 / (w2 * w2) * w3
        dw4_dx0 = -2.0 * x0 / (w2 * w2) * w4

        # Chain rule
        dby_dx = w5 * (dx0_dx * w3 - x0 * dw3_dx0) + (-w5) * (dx0_dx * w4 - x0 * dw4_dx0)

        return dby_dx
    
    def dB_x_dy(x, y):
        Lx = sim.simulation_domain()[0]  # Assuming sim is the simulation object
        Ly = sim.simulation_domain()[1]
        w1 = 0.2
        w2 = 1.0
        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly
        v1 = -1
        v2 = 1.0
        w3 = np.exp(-(x0 * x0 + y1 * y1) / (w2 * w2))
        w4 = np.exp(-(x0 * x0 + y2 * y2) / (w2 * w2))
        w5 = 2.0 * w1 / w2

        # Partial derivatives
        dS_dy1 = 0.5 / 0.5 * np.cosh((y - Ly * 0.3) / 0.5) ** (-2)
        dS_dy2 = 0.5 / 0.5 * np.cosh((y - Ly * 0.7) / 0.5) ** (-2)
        dw3_dy = -2.0 * y1 / (w2 * w2) * w3
        dw4_dy = -2.0 * y2 / (w2 * w2) * w4

        # Calculate dBx_dy
        dBx_dy = (v2 - v1) * (dS_dy1 - dS_dy2) - w5 * y1 * dw3_dy + w5 * y2 * dw4_dy

        return dBx_dy
    
    def vez(x,y):
        J0 = dB_y_dx(x,y) - dB_x_dy(x, y)
        return -J0/density(x,y)

    def vth_ions(x, y):
        return np.sqrt(Ti(x, y))

    def vth_electrons(x, y):
        return np.sqrt(Te(x, y))
    
    mass_electron = .2

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vth_ions,
        "vthy": vth_ions,
        "vthz": vth_ions,
        "nbr_part_per_cell": 100,
    }

    vvv_electrons = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vez,
        "vthx": vth_electrons,
        "vthy": vth_electrons,
        "vthz": vth_electrons,
        "nbr_part_per_cell": 100,
    }

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={"charge": 1, "density": density, **vvv, "init": {"seed": 12334}},
        electrons={"charge": -1, "mass": mass_electron, "density": density, **vvv_electrons, "init": {"seed": 12334}},
    )

    dt = 10 * sim.time_step
    nt = sim.final_time / dt + 1
    timestamps = dt * np.arange(nt)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )

    return sim


def get_time(path, time, datahier=None):
    time = "{:.10f}".format(time)
    from pyphare.pharesee.hierarchy import hierarchy_from

    datahier = hierarchy_from(h5_filename=path + "/EM_E.h5", time=time, hier=datahier)
    datahier = hierarchy_from(h5_filename=path + "/EM_B.h5", time=time, hier=datahier)
    return datahier


def post_advance(new_time):
    if cpp.mpi_rank() == 0:
        print(f"running tests at time {new_time}")
        from tests.simulator.test_advance import AdvanceTestBase

        test = AdvanceTestBase()
        test.base_test_overlaped_fields_are_equal(
            get_time(diag_outputs, new_time), new_time
        )
        print(f"tests passed")


def main():
    s = Simulator(config())#, post_advance=post_advance)
    s.initialize()
    s.run()


if __name__ == "__main__":
    main()
