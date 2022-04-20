"""Replicate Renzo's thesis Fig. 3.2 Panel A Page 47 from gnuplot source."""

import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# Gray matter (direct translation of gnuplot conditions)
# =============================================================================
def stage1_gm(time, T1):
    """Compute Mz for condition: time < Ti."""
    return 1 - 2 * np.exp(-time / T1)


def stage2_gm(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr."""
    return 1 - np.exp((-time+Ti) / T1)


def stage3_gm(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr+Ti."""
    return 1 - np.exp((-time+Tr+Tr) / T1) * (1 + (1 - np.exp((-(Tr+Tr)+Ti) / T1)))


def stage4_gm(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr+Tr+Tr."""
    return 1 - np.exp((-time+Ti+Tr+Tr) / T1)


def stage5_gm(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr+Ti+Tr+Tr."""
    return 1 - np.exp((-time+Tr+Tr+Tr+Tr) / T1) * (1 + (1 - np.exp((-(Tr+Tr+Tr+Tr)+Ti+Tr+Tr) / T1gm)))


def stage6_gm(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr+Tr+Tr+Tr+Tr."""
    return 1 - np.exp((-time+Ti+Tr+Tr+Tr+Tr) / T1)


# =============================================================================
# Blood (direct translation of gnuplot conditions)
# =============================================================================
# NOTE[Faruk]: Renzo seems to skip two stages here because 90 deg pulse has no
# effect on the longitudinal recovery as it happens on zero crossing for blood
def stage1_blood(time, T1):
    """Compute Mz for condition: time < Ti."""
    return 1 - 2 * np.exp(-time / T1)


def stage2_blood(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr."""
    return 1 - np.exp((-time+Ti) / T1)


def stage3_blood(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr+Tr+Tr."""
    return 1 - np.exp((-(time-2*Tr))/T1) * (1 + (1 - np.exp((-(Tr+Tr)+Ti) / T1)))

def stage4_blood(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr+Tr+Tr+Tr+Tr."""
    return 1 - np.exp((-(time-4*Tr)) / T1) * (1 + 1 - np.exp((-(Tr+Tr+Tr+Tr-2*Tr)) / T1) * (1 + (1 -2 * np.exp(-(Tr+Tr)/T1))))

# =============================================================================
# Initial parameters (direct translation of gnuplot conditions)
# =============================================================================
T1csf = 5.
T1gm = 1.9
T1b = 2.1  # steady state blood
Tr = 2.
T1 = T1b  # für Tr= 3
Ti = (np.log(2) - np.log(1 + np.exp(-2 * Tr / T1))) * T1

x = 5 * Tr + 0.06
time = np.linspace(0, x, 1001)

# -----------------------------------------------------------------------------
# Compute plot data
# -----------------------------------------------------------------------------
signal1 = np.zeros(time.shape)
idx = time < Tr+Tr+Tr+Tr+Tr+Tr
signal1[idx] = stage6_gm(time[idx], T1gm, Ti)
idx = time < Tr+Tr+Ti+Tr+Tr
signal1[idx] = stage5_gm(time[idx], T1gm, Ti)
idx = time < Tr+Tr+Tr+Tr
signal1[idx] = stage4_gm(time[idx], T1gm, Ti)
idx = time < Tr+Tr+Ti
signal1[idx] = stage3_gm(time[idx], T1gm, Ti)
idx = time < Tr+Tr
signal1[idx] = stage2_gm(time[idx], T1gm, Ti)
idx = time < Ti
signal1[idx] = stage1_gm(time[idx], T1gm)

signal2 = np.zeros(time.shape)
idx = time < Tr+Tr+Tr+Tr+Tr+Tr
signal2[idx] = stage4_blood(time[idx], T1b, Ti)
idx = time < Tr+Tr+Tr+Tr
signal2[idx] = stage3_blood(time[idx], T1b, Ti)
idx = time < Tr+Tr
signal2[idx] = stage2_blood(time[idx], T1b, Ti)
idx = time < Ti
signal2[idx] = stage1_blood(time[idx], T1b)

# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
plt.plot(time, signal1, linewidth=2, color="blue")
plt.plot(time, signal2, linewidth=2, color="red")

plt.title("Original VASO")
plt.xlabel("Time [s]")
plt.ylabel(r"$M_z$")
plt.ylim([-1, 1])
plt.grid(True)
plt.legend(["Gray Matter", "Steady State Blood"])

plt.show()
