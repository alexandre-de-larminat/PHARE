#!/usr/bin/env python3

import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.hierarchy import get_times_from_h5
from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import flat_finest_field

import matplotlib.animation as animation

from tests.diagnostic import all_timestamps

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.use("Agg")



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

    from pyphare.pharesee.run import Run
    from pyphare.pharesee.hierarchy import flat_finest_field

    vphi, t, phi, a, k = phase_speed(".", 1., 100)

    r = Run(".")
    t = get_times_from_h5("EM_B.h5")
    fig, ax = plt.subplots(figsize=(9, 5), nrows=1)

    #ax, ax2 = ax
    B = r.GetB(t[int(len(t) / 2)])
    by, xby = flat_finest_field(B, "By")
    ax.plot(xby, by, label="t = 500", alpha=0.6)



    B = r.GetB(t[-1])
    by, xby = flat_finest_field(B, "By")
    ax.plot(xby, by, label="t = 1000", alpha=0.6)
    ax.plot(
            xby,
            wave(xby, 1, 2 * np.pi / 100.0, 2 * np.pi / 100 * 50),
            color="k",
            ls="--",
            label="T=500 (theory)",
        )

    B = r.GetB(t[0])
    by, xby = flat_finest_field(B, "By")
    #ax.plot(xby, by, label="t = 0", color="k")
    for j in range(10):
        i=j+1
        ax.plot(xby, flat_finest_field(r.GetB(t[int(len(t) * i*0.1)-1]), "By")[0], alpha=.1*i, c="k")

    ax.set_xlabel("x")
    ax.set_ylabel(r"$B_y$")
    ax.legend(ncol=4, loc="upper center")
    ax.set_ylim((-1.2, 1.3))
    ax.set_title(r"$V_\phi = {:6.4f}$".format(vphi.mean()))

    E = r.GetE(t[int(len(t) / 2)])
    ey, xey = flat_finest_field(E, "Ez")

    #ax2.plot(xey, ey, label="t = 500", alpha=0.6, color='r')
    #ax2.plot(xey, flat_finest_field(r.GetE(t[0]), "Ez")[0], alpha=0.6, color='k')
    #ax2.plot(xey, flat_finest_field(r.GetE(t[-1]), "Ez")[0], alpha=0.6, color='b')
    fig.tight_layout()

    fig.savefig("alfven_wave.png", dpi=200)


"""
    fig, ax = plt.subplots(nrows=1, figsize=(8, 4))

    Bline, = ax.plot(xby, by, color="r", ls="-")
    title = ax.set_title("")
    ax.set_ylabel(r"$B_y$")
    ax.set_xlabel("x")

    def update(frame):
        Bline.set_ydata(flat_finest_field(r.GetB(t[frame]), "By")[0])
        title.set_text(frame)
        return Bline,

    anim = animation.FuncAnimation(
        fig, update, frames=len(t), interval=10, blit=True
    )

    anim.save("alfven.gif", writer="imagemagick", fps=20)
"""
    #assert np.mean(np.abs(vphi - 1) < 5e-2)


if __name__ == "__main__":
    main()
