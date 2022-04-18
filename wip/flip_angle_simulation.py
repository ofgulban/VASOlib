"""[!!!WIP!!!] Flip angle simulation from Renzo."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


# =============================================================================
def Mz(time, T1, FA_rad):
    """Longitudinal decay."""
    return (1 - np.exp(-time/T1)) / (1 - np.cos(FA_rad) * np.exp(-time/T1))


def Mxy(time, T1, FA_rad):
    """Transversal decay."""
    return np.sin(FA_rad) * Mz(time, T1, FA_rad)


def update(val):
    """Update plot data after each slider interaction."""
    time = stime.val

    signal1 = Mxy(time, T1_Blood, FA_rad)
    signal2 = Mxy(time, T1_GM, FA_rad)
    signal3 = Mxy(time, T1_CSF, FA_rad)

    signal4 = Mz(time, T1_Blood, FA_rad)
    signal5 = Mz(time, T1_GM, FA_rad)
    signal6 = Mz(time, T1_CSF, FA_rad)

    line1.set_ydata(signal1)
    line2.set_ydata(signal2)
    line3.set_ydata(signal3)
    line4.set_ydata(signal4)
    line5.set_ydata(signal5)
    line6.set_ydata(signal6)

    fig.canvas.draw_idle()


# =============================================================================
# Initial parameters
# NOTE: T1_Blood = 2122. T1_GM = 1700, T1_CSF = 5500.
FA_deg = np.linspace(0, 60, 101)  # ms
FA_rad = np.deg2rad(FA_deg)

T1_Blood = 2122
T1_GM = 1700
T1_CSF = 5500

time = 72  # ms

signal1 = Mxy(time, T1_Blood, FA_rad)
signal2 = Mxy(time, T1_GM, FA_rad)
signal3 = Mxy(time, T1_CSF, FA_rad)

signal4 = Mz(time, T1_Blood, FA_rad)
signal5 = Mz(time, T1_GM, FA_rad)
signal6 = Mz(time, T1_CSF, FA_rad)

# -----------------------------------------------------------------------------
# Prepare figure
fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.set_title(r"Transversal decay")
ax1.set_xlabel(r"FA [$\degree$]")
ax1.set_ylabel(r"$M_xy$")
# ax1.set_xlim([0, 100])
# ax1.set_ylim([0, 0.2])
line1, = ax1.plot(FA_deg, signal1, lw=2, color="red")
line2, = ax1.plot(FA_deg, signal2, lw=2, color="gray")
line3, = ax1.plot(FA_deg, signal3, lw=2, color="black")

ax2.set_title(r"Longitudinal decay")
ax2.set_xlabel(r"FA [$\degree$]")
ax2.set_ylabel(r"$M_z$")
# ax1.set_xlim([0, 100])
# ax1.set_ylim([0, 0.2])
line4, = ax2.plot(FA_deg, signal4, lw=2, color="red")
line5, = ax2.plot(FA_deg, signal5, lw=2, color="gray")
line6, = ax2.plot(FA_deg, signal6, lw=2, color="black")


plt.tight_layout()
# -----------------------------------------------------------------------------
# Sliders
fig2 = plt.figure()
axcolor = 'lightgoldenrodyellow'

# [left, bottom, width, height]
axtime = plt.axes([0.1, 0.9, 0.8, 0.03], facecolor=axcolor)

stime = Slider(axtime, r"$time$", 1, 100, valinit=time, valstep=1)

stime.on_changed(update)

plt.show()
