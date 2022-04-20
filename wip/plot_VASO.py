"""Reproduce Renzo's thesis Fig. 3.2 Panel A Page 47."""

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
signal1[time < Tr+Tr+Tr+Tr+Tr+Tr] = stage6_gm(time[time < Tr+Tr+Tr+Tr+Tr+Tr], T1gm, Ti)
signal1[time < Tr+Tr+Ti+Tr+Tr] = stage5_gm(time[time < Tr+Tr+Ti+Tr+Tr], T1gm, Ti)
signal1[time < Tr+Tr+Tr+Tr] = stage4_gm(time[time < Tr+Tr+Tr+Tr], T1gm, Ti)
signal1[time < Tr+Tr+Ti] = stage3_gm(time[time < Tr+Tr+Ti], T1gm, Ti)
signal1[time < Tr+Tr] = stage2_gm(time[time < Tr+Tr], T1gm, Ti)
signal1[time < Ti] = stage1_gm(time[time < Ti], T1gm)

signal2 = np.zeros(time.shape)
signal2[time < Tr+Tr+Tr+Tr+Tr+Tr] = stage4_blood(time[time < Tr+Tr+Tr+Tr+Tr+Tr], T1b, Ti)
signal2[time < Tr+Tr+Tr+Tr] = stage3_blood(time[time < Tr+Tr+Tr+Tr], T1b, Ti)
signal2[time < Tr+Tr] = stage2_blood(time[time < Tr+Tr], T1b, Ti)
signal2[time < Ti] = stage1_blood(time[time < Ti], T1b)



# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Faruk's rework start from here
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def Mz(time, M0_equi, M0_init, FA_rad, T1):
    """Longitudinal magnetization."""
    return M0_equi - (M0_equi - M0_init * np.cos(FA_rad)) * np.exp(-time/T1)


# Prepare condition array
cond = np.full(time.shape, 2)
idx1 = time % (Tr*2) < Ti  # Stages after 180 deg pulse
cond[idx1] = 1

# =============================================================================
# Compute signal array
# =============================================================================
signal3 = np.zeros(time.shape)

# Prepare initial parameters
M0_equi = 1.  # This never changes
M0_init = 1.
for i, t in enumerate(time):
    t %= 2*Tr

    # -------------------------------------------------------------------------
    # Handle first signal separately
    # -------------------------------------------------------------------------
    if i == 0:
        signal3[i] = Mz(time=t, M0_equi=M0_equi, M0_init=M0_init,
                        FA_rad=np.deg2rad(180), T1=T1gm)

    # -------------------------------------------------------------------------
    # 180 degree pulse
    # -------------------------------------------------------------------------
    elif cond[i] == 1:
        if cond[i] != cond[i-1]:  # Update M0 upon condition switch
            M0_init = Mz(time=Tr+Tr-Ti, M0_equi=M0_equi, M0_init=M0_init,
                         FA_rad=np.deg2rad(90), T1=T1gm)
        signal3[i] = Mz(time=t, M0_equi=M0_equi, M0_init=M0_init,
                        FA_rad=np.deg2rad(180), T1=T1gm)

    # -------------------------------------------------------------------------
    # 90 degree pulse
    # -------------------------------------------------------------------------
    else:
        signal3[i] = Mz(time=t-Ti, M0_equi=M0_equi, M0_init=M0_init,
                        FA_rad=np.deg2rad(90), T1=T1gm)

# =============================================================================
# Compute signal array
# =============================================================================
signal4 = np.zeros(time.shape)

# Prepare initial parameters
M0_equi = 1.  # This never changes
M0_init = 1.
for i, t in enumerate(time):
    t %= 2*Tr

    # -------------------------------------------------------------------------
    # Handle first signal separately
    # -------------------------------------------------------------------------
    if i == 0:
        signal4[i] = Mz(time=t, M0_equi=M0_equi, M0_init=M0_init,
                        FA_rad=np.deg2rad(180), T1=T1b)

    # -------------------------------------------------------------------------
    # 180 degree pulse
    # -------------------------------------------------------------------------
    elif t < Ti:
        if cond[i] != cond[i-1]:  # Update M0 upon condition switch
            M0_init = Mz(time=Tr+Tr-Ti, M0_equi=M0_equi, M0_init=M0_init,
                         FA_rad=np.deg2rad(90), T1=T1b)
        signal4[i] = Mz(time=t, M0_equi=M0_equi, M0_init=M0_init,
                            FA_rad=np.deg2rad(180), T1=T1b)

    # -------------------------------------------------------------------------
    # 90 degree pulse
    # -------------------------------------------------------------------------
    else:
        signal4[i] = Mz(time=t-Ti, M0_equi=M0_equi, M0_init=M0_init,
                        FA_rad=np.deg2rad(90), T1=T1b)


# =============================================================================
# Plotting
# =============================================================================
plt.plot(time, signal1, linewidth=3, color="blue")
plt.plot(time, signal2, linewidth=3, color="red")

plt.plot(time, signal3, linewidth=2, color="cyan")
plt.plot(time, signal4, linewidth=2, color="orange")

plt.title("Original VASO")
plt.xlabel("Time [s]")
plt.ylabel(r"$M_z$")
plt.ylim([-1, 1])
plt.grid(True)
plt.legend(["R: Gray Matter", "R: Steady State Blood",
            "F: Gray Matter", "F: Steady State Blood"])

plt.show()
