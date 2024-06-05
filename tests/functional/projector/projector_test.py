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
        time_step_nbr=1,
        time_step=0.1,
        # boundary_types="periodic",
        cells=(10, 10),
        dl=(1.0, 1.0),
        refinement_boxes={},
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_outputs, "mode": "overwrite"},
        },
        strict=False,
    )

    def density(x, y):
        return 1.0


    def by(x, y):
        return 0.0

    def bx(x, y):
        return 0.0

    def bz(x, y):
        return 0.0

    def b2(x, y):
        return bx(x, y) ** 2 + by(x, y) ** 2 + bz(x, y) ** 2

    def T(x, y):
        K = 1
        temp = 1.0 / density(x, y) * (K - b2(x, y) * 0.5)
        assert np.all(temp > 0)
        return temp

    def vx(x, y):
        return 1.0

    def vy(x, y):
        return 0.0

    def vz(x, y):
        return 0.0

    def vthx(x, y):
        return np.sqrt(T(x, y))

    def vthy(x, y):
        return np.sqrt(T(x, y))

    def vthz(x, y):
        return np.sqrt(T(x, y))
    

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthx,
        "vthy": vthy,
        "vthz": vthz,
        "nbr_part_per_cell": 10,
    }

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={"charge": 1, "density": density, **vvv},
        electrons={"charge": -1, "density": density, **vvv},
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




def main():
    s = Simulator(config())
    s.initialize()
    s.run()


if __name__ == "__main__":
    main()
