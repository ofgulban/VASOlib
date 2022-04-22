"""Emulate SS-SI-VASO images from a T1 input image."""

import numpy as np
import nibabel as nb


# =============================================================================
def Mz(time, M0_equi, M0_init, FA_rad, T1):
    """Longitudinal magnetization."""
    return M0_equi - (M0_equi - M0_init * np.cos(FA_rad)) * np.exp(-time/T1)


def compute_SI_SS_VASO_Mz_signal(time, T1, Tr, Ti1, Ti2, mode_nonblood=False):
    """Compute VASO Mz signal."""
    signal = np.zeros(time.shape)
    M0_equi = 1.  # This never changes
    M0_init = 1.

    # Prepare condition array
    cond = np.full(time.shape, 3)
    idx2 = time % (Tr*2) < (Ti2+Tr)  # Stages after 180 deg pulse
    idx1 = time % (Tr*2) < Ti1  # Stages after the first 90 deg pulse
    cond[idx2] = 2
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
        # After 180 degree pulse
        # ---------------------------------------------------------------------
        elif cond[i] == 1:
            if mode_nonblood:
                if cond[i] != cond[i-1]:  # Update M0 upon condition switch
                    M0_init = Mz(time=Tr-Ti2, M0_equi=M0_equi, M0_init=M0_init,
                                 FA_rad=np.deg2rad(90), T1=T1)
            signal[i] = Mz(time=t, M0_equi=M0_equi, M0_init=M0_init,
                           FA_rad=np.deg2rad(180), T1=T1)
        # ---------------------------------------------------------------------
        # After the first 90 degree pulse
        # ---------------------------------------------------------------------
        elif cond[i] == 2:
            signal[i] = Mz(time=t-Ti1, M0_equi=M0_equi, M0_init=M0_init,
                           FA_rad=np.deg2rad(90), T1=T1)
        # ---------------------------------------------------------------------
        # After the second 90 degree pulse
        # ---------------------------------------------------------------------
        else:
            signal[i] = Mz(time=t-Tr-Ti2, M0_equi=M0_equi, M0_init=M0_init,
                           FA_rad=np.deg2rad(90), T1=T1)
    return signal


# =============================================================================
# Initial parameters
# =============================================================================
Ti1 = 1.45561 * 1000
Ti2 = 1.7 * 1000
Tr = 2. * 1000
max_time = 6 * Tr

# Compute readout times
time_readout1 = np.arange(Ti1, max_time, Tr*2)
time_readout2 = np.arange(Tr+Ti2, max_time, Tr*2)

nr_timepoints = time_readout1.size + time_readout2.size
time_readouts = np.zeros(nr_timepoints)
time_readouts[0::2] = time_readout1
time_readouts[1::2] = time_readout2

compute_SI_SS_VASO_Mz_signal(time_readouts-1, T1=1900, Tr=Tr, Ti1=Ti1, Ti2=Ti2, mode_nonblood=True)

# Load nifti
FILE = "/home/faruk/gdrive/proj-BOCO+/data/sub-01_T1_avg.nii.gz"
nii = nb.load(FILE)

# Select a single slice
data_T1 = np.asarray(nii.dataobj[:, :, 72], dtype="int")

# Impute zeros
data_T1[data_T1 < 2] = 2
dims = data_T1.shape

data_T1 = data_T1.flatten()
data = np.zeros((data_T1.size, nr_timepoints))
for i, v in enumerate(data_T1):
    data[i, :] = compute_SI_SS_VASO_Mz_signal(
        time_readouts-1, T1=v, Tr=Tr, Ti1=Ti1, Ti2=Ti2, mode_nonblood=True)
data = np.abs(data)

OUTFILE = "/home/faruk/gdrive/proj-BOCO+/data/test.nii.gz"
data = data.reshape([dims[0], dims[1], 1, nr_timepoints])
img = nb.Nifti1Image(data, affine=np.eye(4))
nb.save(img, OUTFILE)

# BOLD correction
nulled = data[..., 0::2]
bold = data[..., 1::2]
data2 =  nulled / bold

# data2[data2>5] = 5
OUTFILE = "/home/faruk/gdrive/proj-BOCO+/data/test2.nii.gz"
img = nb.Nifti1Image(data2, affine=np.eye(4))
nb.save(img, OUTFILE)

print("Finished.")
