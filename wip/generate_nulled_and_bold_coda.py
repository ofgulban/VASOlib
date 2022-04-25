"""Interactive VASO simulations."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import compoda.core as coda


# =============================================================================
def signal_bold(time, S0=100, T2star=28):
    """Equation 7 from Akbari et al. (2021)."""
    return S0 * np.exp(-time/T2star)


def signal_nulled(time, S0=100, T2star=28, TI=10, T1=10, TR=10, epsilon=1):
    """Equation 11 from Akbari et al. (2021)."""
    inner_term = 1 - (1 + epsilon) * np.exp(-TI/T1) + epsilon * np.exp(-TR/T1)
    return S0 * inner_term * np.exp(-time/T2star)


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
fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4)

ax1.set_title(r"BOLD (red) & Nulled")
ax1.set_xlabel("Time [ms]")
ax1.set_ylabel("MRI signal")
ax1.plot(time, signal1, lw=2, color="red")
ax1.plot(time, signal2, lw=2, color="black")

ax2.set_title(r"BOCO ~ (1-CBV)")
ax2.set_xlabel("Time [ms]")
ax2.set_ylim([0, 1])

ax3.set_title(r"(1-BOCO) ~ CBV")
ax3.set_xlabel("Time [ms]")
ax3.set_ylim([0, 1])

ax4.set_title(r"Centered logratio")
ax4.set_xlabel("Time [ms]")

plt.tight_layout()
# -----------------------------------------------------------------------------
# Iso-T1 lines
for i in np.linspace(1000, 3000, 21):
    signal2 = signal_nulled(time, S0, T2star, TI, i, TR, Ep)
    relative_blood = signal2 / signal1
    relative_notblood = 1 - relative_blood

    ax2.plot(time, relative_blood, lw=0.5, color="gray")
    ax3.plot(time, relative_notblood, lw=0.5, color="gray")

    # -------------------------------------------------------------------------
    # Compositional perspective
    signal3 = signal1 - signal2
    bary = np.vstack([signal2, signal3]).T
    bary = coda.closure(bary)
    clr = coda.clr_transformation(bary)[:, 0]
    ax4.plot(time, clr, lw=0.5, color="gray")

plt.show()
