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

    #vphi, t, phi, a, k = phase_speed(".", .01, 10)

    r = Run("alt")
    t = get_times_from_h5("alt/EM_B.h5")
    fig, ax = plt.subplots(figsize=(9, 5), nrows=1)

    #ax, ax2 = ax

    B = r.GetB(t[0])
    by, xby = flat_finest_field(B, "By")
    ax.plot(xby, by, label="t = 0", alpha=0.6)
    B = r.GetB(t[int(len(t) * 0.1)])
    by, xby = flat_finest_field(B, "By")
    ax.plot(xby, by, label="t = 1", alpha=0.6)

    B = r.GetB(t[int(len(t) * 0.2)])
    by, xby = flat_finest_field(B, "By")
    ax.plot(xby, by, label="t = 2", alpha=0.6)


    B = r.GetB(t[int(len(t) * 0.3)])
    by, xby = flat_finest_field(B, "By")
    ax.plot(xby, by, label="t = 3", alpha=0.6)

    B = r.GetB(t[int(len(t) * 0.4)])
    by, xby = flat_finest_field(B, "By")
    ax.plot(xby, by, label="t = 4", alpha=0.6)

    B = r.GetB(t[int(len(t) * 0.5)])
    by, xby = flat_finest_field(B, "By")
    ax.plot(xby, by, label="t = 5", alpha=0.6)

    B = r.GetB(t[int(len(t) *.75)])
    by, xby = flat_finest_field(B, "By")
    ax.plot(xby, by, label="t = 7.5", alpha=0.6)

    B = r.GetB(t[int(len(t) -1)])
    by, xby = flat_finest_field(B, "By")
    ax.plot(xby, by, label="t = 10", alpha=0.6)


    #B = r.GetB(t[-1])
    #by, xby = flat_finest_field(B, "By")
    #ax.plot(xby, by, label="t = 100", alpha=0.6)

    ax.plot(
            xby,
            wave(xby, .01, 2 * np.pi / 10.0, 2 * np.pi / 10 * 5),
            color="k",
            ls="--",
            label="T=50 (theory)",
        )
    """
    B = r.GetB(t[0])
    by, xby = flat_finest_field(B, "By")
    #ax.plot(xby, by, label="t = 0", color="k")
    for j in range(2):
        i=j+1
        ax.plot(xby, flat_finest_field(r.GetB(t[int(len(t) * i*0.5)-1]), "By")[0], alpha=.5*i, c="k")
    """
    ax.set_xlabel("x")
    ax.set_ylabel(r"$B_y$")
    ax.legend(ncol=4, loc="upper center")
    ax.set_ylim((-1.2*0.01, 1.3*0.01))
    #ax.set_title(r"$V_\phi = {:6.4f}$".format(vphi.mean()))

    E = r.GetE(t[int(len(t) / 2)])
    ey, xey = flat_finest_field(E, "Ez")

    #ax2.plot(xey, ey, label="t = 500", alpha=0.6, color='r')
    #ax2.plot(xey, flat_finest_field(r.GetE(t[0]), "Ez")[0], alpha=0.6, color='k')
    #ax2.plot(xey, flat_finest_field(r.GetE(t[-1]), "Ez")[0], alpha=0.6, color='b')
    fig.tight_layout()

    fig.savefig("alfven_wave.png", dpi=200)


def anim():
    r = Run(".")
    t = get_times_from_h5("EM_B.h5")
    B = r.GetB(t[-1])
    by, xby = flat_finest_field(B, "By")

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
        fig, update, frames=int(len(t)/10), interval=10, blit=True
    )

    anim.save("alfven.gif", writer="imagemagick", fps=20)

    #assert np.mean(np.abs(vphi - 1) < 5e-2)


if __name__ == "__main__":
    main()
    #anim()
