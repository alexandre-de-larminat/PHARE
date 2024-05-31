import numpy as np
import matplotlib.pyplot as plt

traj1 = np.loadtxt("pusher_test_in.txt", delimiter = " ")
traj2 = np.loadtxt("pusher_alt_trajectory.txt", delimiter = " ")

reduced1 = traj1[::1]
reduced2 = traj2[::1]

tstart = 0.0
tend = 10.0
dt = 0.0001
nt = int((tend - tstart) / dt) + 1
t = np.arange(0, nt * dt, dt)
t = t[::1]

plt.plot(t, reduced1[:, 0], "-", label="non corrected")
plt.plot(t, reduced2[:, 0], "--", label="corrected")
#plt.plot(t, reduced1[:, 0] - reduced2[:, 0], "-", label="difference")
plt.xlabel("Time")
plt.ylabel("Position in x")
#plt.xlim(9.9975, 10.0)
#plt.ylim(30.017, 30.031)
#plt.legend()
#plt.show()
plt.tight_layout()
plt.savefig("pusher_test.pdf")