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

diag_outputs = "output"
from datetime import datetime
from tests.diagnostic import all_timestamps


def config():
    sim = ph.Simulation(
        time_step_nbr=10,
        time_step=0.1,
        # boundary_types="periodic",
        cells=(10,10),
        dl=(1.0,1.0),
        refinement_boxes={},
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_outputs, "mode": "overwrite"},
        },
        strict=False,
    )

    def density(x,y):
        return 1.0


    def by(x,y):
        return 0.0

    def bx(x,y):
        return 0.0

    def bz(x,y):
        return 0.0


    def T(x,y):
        return 1.0

    def vx(x,y):
        return 1.0

    def vy(x,y):
        return 0.0

    def vz(x,y):
        return 0.0

    def vthx(x,y):
        return np.sqrt(T(x,y))

    def vthy(x,y):
        return np.sqrt(T(x,y))

    def vthz(x,y):
        return np.sqrt(T(x,y))
    

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


    timestamps = all_timestamps(sim)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )

    for quantity in ["density"]:
        ph.FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )
        
    for quantity in ["density"]:
        ph.FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
            population_name="all_electrons",
        )


    return sim


def figure():
    from pyphare.pharesee.run import Run
    from pyphare.pharesee.hierarchy import get_times_from_h5
    from pyphare.pharesee.hierarchy import flat_finest_field
    from pyphare.pharesee.hierarchy import hierarchy_from
  
    import os


    run_path = "./output"
    run = Run(run_path)

    time = get_times_from_h5(os.path.join(run_path, "EM_B.h5"))

    plot_time = time[int(len(time)-1)]
    first_time = time[0]


    J = run.GetJ(plot_time, merged=True, interp="bilinear")

    print(J["Jx"]) 


    B = run.GetB(plot_time)

    print(B)
    print(B["Bx"]) 


    xyjx = J["Jz"][1][0]
    jx = J["Jz"][0](xyjx)

    print(jx)
    print(xyjx)

    #xjx = xyjx[:,0]
    #yjx = xyjx[:,1]

    ion_density = run.GetNi(plot_time)
    electron_density = run.GetNe(plot_time)

    Ni, xNi = flat_finest_field(ion_density, "rho")
    Ne, xNe = flat_finest_field(electron_density, "rho")


    fig, axarr = plt.subplots(nrows=2, figsize=(8, 10))




    ax0, ax1 = axarr

    ax0.plot(jx, color="k", ls="-")
    ax0.set_ylabel("Jx")

    charge_density = Ni
    for i in range(len(Ni)):
        charge_density[i] = Ni[i] - Ne[i]


    ax1.plot(xNi, charge_density, color="k", ls="-")
    ax1.set_ylabel("Charge density")

    #fig.tight_layout()
    fig.savefig("test.png")



def main():
    s = Simulator(config())
    s.initialize()
    s.run()


if __name__ == "__main__":
    main()
    figure()
