"""Interactive VASO simulations."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


# =============================================================================
def signal_bold(time, S0=100, T2star=28):
    """Equation 7 from Akbari et al. (2021)."""
    return S0 * np.exp(-time/T2star)


def signal_nulled(time, S0=100, T2star=28, TI=10, T1=10, TR=10, epsilon=1):
    """Equation 11 from Akbari et al. (2021)."""
    inner_term = 1 - (1 + epsilon) * np.exp(-TI/T1) + epsilon * np.exp(-TR/T1)
    return S0 * inner_term * np.exp(-time/T2star)


def update(val):
    """Update plot data after each slider interaction."""
    S0 = sS0.val
    T2star = sT2star.val
    TI = sTI.val
    T1 = sT1.val
    TR = sTR.val
    Ep = sEp.val
    NS = sNS.val
    signal1 = signal_bold(time, S0=S0, T2star=T2star) + NS * noise1
    signal2 = signal_nulled(time, S0=S0, T2star=T2star, TI=TI, TR=TR, T1=T1,
                            epsilon=Ep) + NS * noise2
    line1.set_ydata(signal1)
    line2.set_ydata(signal2)
    line3.set_ydata(signal2 / signal1)
    fig.canvas.draw_idle()


# =============================================================================
time = np.linspace(0, 100, 101)  # ms

# Prepare noise
NS = 0
noise1 = np.random.normal(loc=0.0, scale=1.0, size=101)
noise2 = np.random.normal(loc=0.0, scale=1.0, size=101)

# Prepare signal
S0 = 100
T2star = 28  # ms
TI, T1 = 1500, 2000  # ms
TR = TI * 2 + 1000  # ms
Ep = 1  # Inversion efficiency
signal1 = signal_bold(time, S0) + NS * noise1
signal2 = signal_nulled(time, S0, T2star, TI, T1, TR, Ep) + NS * noise2

# -----------------------------------------------------------------------------
# Prepare figure
fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.set_title(r"BOLD (red) & Nulled")
ax1.set_xlabel("Time [ms]")
ax1.set_ylabel("MRI signal")
line1, = ax1.plot(time, signal1, lw=2, color="red")
line2, = ax1.plot(time, signal2, lw=2, color="black")

ax2.set_title(r"Division (BOCO)")
ax2.set_xlabel("Time [ms]")
# ax2.set_ylabel("Signal")
line3, = ax2.plot(time, signal2 / signal1, lw=2)
ax2.set_ylim([0, 1])

plt.tight_layout()
# -----------------------------------------------------------------------------
# Iso-T1 lines
for i in np.linspace(1000, 3000, 21):
    signal2 = signal_nulled(time, S0, T2star, TI, i, TR, Ep)
    ax2.plot(time, signal2 / signal1, lw=0.5, color="gray")

# -----------------------------------------------------------------------------
# Sliders
axcolor = 'lightgoldenrodyellow'

fig2 = plt.figure()
# [left, bottom, width, height]
axS0 = plt.axes([0.1, 0.9, 0.8, 0.03], facecolor=axcolor)
axT2star = plt.axes([0.1, 0.85, 0.8, 0.03], facecolor=axcolor)
axTI = plt.axes([0.1, 0.80, 0.8, 0.03], facecolor=axcolor)
axT1 = plt.axes([0.1, 0.75, 0.8, 0.03], facecolor=axcolor)
axTR = plt.axes([0.1, 0.70, 0.8, 0.03], facecolor=axcolor)
axEp = plt.axes([0.1, 0.65, 0.8, 0.03], facecolor=axcolor)

axNS = plt.axes([0.1, 0.55, 0.8, 0.03], facecolor=axcolor)

sS0 = Slider(axS0, "$S_0$", 0, 200, valinit=S0, valstep=10)
sT2star = Slider(axT2star, r"$T_{2}^*$", 1, 100, valinit=T2star, valstep=1)
sTI = Slider(axTI, r"$TI$", 0, 4000, valinit=TI, valstep=100)
sT1 = Slider(axT1, r"$T_1$", 0, 4000, valinit=T1, valstep=100)
sTR = Slider(axTR, r"$TR$", 0, 6000, valinit=TR, valstep=100)
sEp = Slider(axEp, r"$\epsilon$", 0, 1, valinit=Ep, valstep=0.01)

sNS = Slider(axNS, r"$Noise$", 0, 10, valinit=NS, valstep=0.1)

sS0.on_changed(update)
sT2star.on_changed(update)
sTI.on_changed(update)
sT1.on_changed(update)
sTR.on_changed(update)
sEp.on_changed(update)
sNS.on_changed(update)

plt.show()
