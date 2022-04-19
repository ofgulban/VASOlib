"""[!!!WIP!!!] VASO."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


# =============================================================================
def stage1(time, T1):
    """Stage 1: time < Ti."""
    return 1 - 2 * np.exp(-time / T1)


def stage2(time, T1, Ti):
    """Stage 2: time < Tr+Tr."""
    return 1 - np.exp((-time+Ti) / T1)


def stage3(time, T1, Ti):
    """Stage 3: time < Tr+Tr+Ti."""
    return 1 - np.exp((-time+Tr+Tr) / T1) * (1 + (1 - np.exp((-(Tr+Tr) + Ti) / T1)))


def stage4(time, T1, Ti):
    """Stage 4: time < Tr+Tr+Tr+Tr."""
    return 1 - np.exp((-time+Ti+Tr+Tr) / T1)


def stage5(time, T1, Ti):
    """Stage 5: time < Tr+Tr+Ti+Tr+Tr."""
    return 1 - np.exp((-time+Tr+Tr+Tr+Tr) / T1) * (1 + (1 - np.exp((-(Tr+Tr+Tr+Tr)+Ti+Tr+Tr) / T1)))


def stage6(time, T1, Ti):
    """Stage 6: time < Tr+Tr+Tr+Tr+Tr+Tr."""
    return 1 - np.exp((-x+Ti+Tr+Tr+Tr+Tr) / T1)


# =============================================================================
# Initial parameters
T1csf = 5.
T1gm = 1.9
T1b = 2.1
Tr = 2.
T1 = T1b  # für Tr= 3
Ti = (np.log(2) - np.log(1 + np.exp(-2 * Tr / T1))) * T1

x = 5 * Tr + 0.06
time = np.linspace(0, x, 1001)

# -----------------------------------------------------------------------------
signal1 = np.zeros(time.shape)
signal1[time < Tr+Tr+Tr+Tr+Tr+Tr] = stage6(time[time < Tr+Tr+Tr+Tr+Tr+Tr], T1b, Ti)
signal1[time < Tr+Tr+Ti+Tr+Tr] = stage5(time[time < Tr+Tr+Ti+Tr+Tr], T1b, Ti)
signal1[time < Tr+Tr+Tr+Tr] = stage4(time[time < Tr+Tr+Tr+Tr], T1b, Ti)
signal1[time < Tr+Tr+Ti] = stage3(time[time < Tr+Tr+Ti], T1b, Ti)
signal1[time < Tr+Tr] = stage2(time[time < Tr+Tr], T1b, Ti)
signal1[time < Ti] = stage1(time[time < Ti], T1b)

# -----------------------------------------------------------------------------
signal2 = np.zeros(time.shape)
signal2[time < Tr+Tr+Tr+Tr+Tr+Tr] = stage6(time[time < Tr+Tr+Tr+Tr+Tr+Tr], T1gm, Ti)
signal2[time < Tr+Tr+Ti+Tr+Tr] = stage5(time[time < Tr+Tr+Ti+Tr+Tr], T1gm, Ti)
signal2[time < Tr+Tr+Tr+Tr] = stage4(time[time < Tr+Tr+Tr+Tr], T1gm, Ti)
signal2[time < Tr+Tr+Ti] = stage3(time[time < Tr+Tr+Ti], T1gm, Ti)
signal2[time < Tr+Tr] = stage2(time[time < Tr+Tr], T1gm, Ti)
signal2[time < Ti] = stage1(time[time < Ti], T1gm)

plt.plot(time, signal1)
plt.plot(time, signal2)

plt.title("Original VASO")
plt.xlabel("Time [s]")
plt.ylabel(r"$M_z$")
plt.ylim([-1, 1])
plt.grid(True)

plt.show()
