"""Interactive VASO longitudinal magnetization simulation.

Reference
---------
- Renzo Huber's PhD Thesis Fig. 3.2 Panel A Page 47.

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


# =============================================================================
# Faruk's reformulation of Mz equation from Renzo's gnuplot scripts
# =============================================================================
def Mz(time, M0_equi, M0_init, FA_rad, T1):
    """Longitudinal magnetization."""
    return M0_equi - (M0_equi - M0_init * np.cos(FA_rad)) * np.exp(-time/T1)


def compute_VASO_Mz_signal(time, T1, Tr, Ti):
    """Compute VASO Mz signal."""
    signal = np.zeros(time.shape)
    M0_equi = 1.  # This never changes
    M0_init = 1.

    # Prepare condition array
    cond = np.full(time.shape, 2)
    idx1 = time % (Tr*2) < Ti  # Stages after 180 deg pulse
    cond[idx1] = 1

    for i, t in enumerate(time):
        t %= 2*Tr
        # ---------------------------------------------------------------------
        # Handle first signal separately
        # ---------------------------------------------------------------------
        if i == 0:
            signal[i] = Mz(time=t, M0_equi=M0_equi, M0_init=M0_init,
                           FA_rad=np.deg2rad(180), T1=T1)
        # ---------------------------------------------------------------------
        # 180 degree pulse
        # ---------------------------------------------------------------------
        elif cond[i] == 1:
            if cond[i] != cond[i-1]:  # Update M0 upon condition switch
                M0_init = Mz(time=Tr+Tr-Ti, M0_equi=M0_equi, M0_init=M0_init,
                             FA_rad=np.deg2rad(90), T1=T1)
            signal[i] = Mz(time=t, M0_equi=M0_equi, M0_init=M0_init,
                           FA_rad=np.deg2rad(180), T1=T1)
        # ----------------------------------------------------------------------
        # 90 degree pulse
        # ----------------------------------------------------------------------
        else:
            signal[i] = Mz(time=t-Ti, M0_equi=M0_equi, M0_init=M0_init,
                           FA_rad=np.deg2rad(90), T1=T1)
    return signal


def update(val):
    """Update plot data after each slider interaction."""
    T1gm = sT1.val
    max_time = sTime.val
    Ti = sTi.val
    Tr = sTr.val

    time = np.linspace(0, max_time, 501)
    signal1 = compute_VASO_Mz_signal(time, T1gm, Tr, Ti)
    signal2 = compute_VASO_Mz_signal(time, T1b, Tr, Ti)

    # -------------------------------------------------------------------------
    # Clear axis and re-draw elements
    # -------------------------------------------------------------------------
    ax1.cla()
    ax1.set_title("Original VASO")
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel(r"$M_z$")
    ax1.set_xlim([0, max_time])
    ax1.set_ylim([-1, 1])
    ax1.plot(time, signal1, lw=2, color="blue")
    ax1.plot(time, signal2, lw=2, color="red")
    ax1.legend(['Tissue X', 'Blood'], loc="upper left")

    # -------------------------------------------------------------------------
    # Reference lines
    # -------------------------------------------------------------------------
    # Horizontal
    ax1.hlines([0], 0, max_time, linestyle='solid',
               color='lightgray', zorder=0)

    # Vertical lines
    event_180deg = np.arange(0, max_time, 2*Tr)
    ax1.vlines(event_180deg, -1, 1, linestyle=':', color='gray', zorder=0)
    trans = ax1.get_xaxis_transform()
    for x in event_180deg:
        ax1.text(x, 0.02, r"$180\degree$ pulse", rotation=90, transform=trans)

    event_90deg = np.arange(+Ti, max_time, 2*Tr)
    for x in event_90deg:
        ax1.text(x, 0.02, r"$90\degree$ pulse", rotation=90, transform=trans)
    ax1.vlines(event_90deg, -1, 1, linestyle=':', color='gray', zorder=0)

    fig1.canvas.draw_idle()


# =============================================================================
# Initial parameters
# =============================================================================
T1csf = 5.
T1gm = 1.9
T1b = 2.1  # steady state blood
Tr = 2.
T1 = T1b  # für Tr= 3
Ti = (np.log(2) - np.log(1 + np.exp(-2 * Tr / T1))) * T1

max_time = 5 * Tr
time = np.linspace(0, max_time, 1001)

signal1 = compute_VASO_Mz_signal(time, T1gm, Tr, Ti)
signal2 = compute_VASO_Mz_signal(time, T1b, Tr, Ti)

# =============================================================================
# Plotting
# =============================================================================
# Prepare figure
fig1, (ax1) = plt.subplots(1, 1)

ax1.set_title("Original VASO")
ax1.set_xlabel("Time [s]")
ax1.set_ylabel(r"$M_z$")
ax1.set_xlim([0, max_time])
ax1.set_ylim([-1, 1])
ax1.plot(time, signal1, lw=2, color="blue")
ax1.plot(time, signal2, lw=2, color="red")
ax1.legend(['Tissue X', 'Blood'], loc="upper left")

# -----------------------------------------------------------------------------
# Reference lines
# -----------------------------------------------------------------------------
# Horizontal
ax1.hlines([0], 0, max_time, linestyle='solid',
           color='lightgray', zorder=0)

# Vertical lines
event_180deg = np.arange(0, max_time, 2*Tr)
ax1.vlines(event_180deg, -1, 1, linestyle=':', color='gray', zorder=0)
trans = ax1.get_xaxis_transform()
for x in event_180deg:
    ax1.text(x, 0.02, r"$180\degree$ pulse", rotation=90, transform=trans)

event_90deg = np.arange(+Ti, max_time, 2*Tr)
for x in event_90deg:
    ax1.text(x, 0.02, r"$90\degree$ pulse", rotation=90, transform=trans)
ax1.vlines(event_90deg, -1, 1, linestyle=':', color='gray', zorder=0)

# -----------------------------------------------------------------------------
# Sliders
# -----------------------------------------------------------------------------
fig2 = plt.figure()
axcolor = 'lightgoldenrodyellow'

# [left, bottom, width, height]
axT1 = plt.axes([0.15, 0.9, 0.70, 0.03], facecolor=axcolor)
axTime = plt.axes([0.15, 0.85, 0.70, 0.03], facecolor=axcolor)
axTi = plt.axes([0.15, 0.80, 0.70, 0.03], facecolor=axcolor)
axTr = plt.axes([0.15, 0.75, 0.70, 0.03], facecolor=axcolor)

sT1 = Slider(axT1, r"$T_1$", 0, 6.0, valinit=T1gm, valstep=0.1)
sTime = Slider(axTime, r"$Max. Time$", 1, 30, valinit=max_time, valstep=1)
sTi = Slider(axTi, r"$Ti$", 0, 6.0, valinit=Ti, valstep=0.1)
sTr = Slider(axTr, r"$Tr$", 0, 6.0, valinit=Tr, valstep=0.1)

sT1.on_changed(update)
sTime.on_changed(update)
sTi.on_changed(update)
sTr.on_changed(update)

plt.show()
