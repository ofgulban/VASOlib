"""Demonstrate scale of proportions using the rainy days example."""

import numpy as np
import compoda.core as coda

# =============================================================================
d_rainy = [0.1, 0.2]  # desert rainy days in a year [%]
m_rainy = [20, 20.1]  # mountain rainy days in a year [%]

# -----------------------------------------------------------------------------
# Let's compute the change in 'rainy days in a year'.
# -----------------------------------------------------------------------------
d_change = d_rainy[1] - d_rainy[0]
m_change = m_rainy[1] - m_rainy[0]

print("Change in rainy days in a year: [desert | mountain]")
print("    {:.1f}% | {:.1f}%".format(d_change, m_change))

# -----------------------------------------------------------------------------
# Let's compute the rainy days increment of ratios (percent change)
# -----------------------------------------------------------------------------
d_increment = (d_rainy[1] - d_rainy[0]) / d_rainy[0] * 100  # %
m_increment = (m_rainy[1] - m_rainy[0]) / m_rainy[0] * 100  # %

print("Increment in rainy days in a year: [desert | mountain]")
print("    {:.1f}% | {:.1f}% \n".format(d_increment, m_increment))

# =============================================================================
# Let's compute the increment of ratios like in MRI (percent change)
# =============================================================================
d_nonrainy = [100-i for i in d_rainy]  # %
m_nonrainy = [100-i for i in m_rainy]  # %

d_change = d_nonrainy[1] - d_nonrainy[0]
m_change = m_nonrainy[1] - m_nonrainy[0]

print("Change in non-rainy days in a year: [desert | mountain]")
print("    {:.1f}% | {:.1f}%".format(d_change, m_change))

# -----------------------------------------------------------------------------
# Let's compute the non-rainy days increment of ratios (percent change)
# -----------------------------------------------------------------------------
d_increment = (d_nonrainy[1] - d_nonrainy[0]) / d_nonrainy[0] * 100  # %
m_increment = (m_nonrainy[1] - m_nonrainy[0]) / m_nonrainy[0] * 100  # %

print("Increment in non-rainy days in a year: [desert | mountain]")
print("    {:.3f}% | {:.3f}% \n".format(d_increment, m_increment))

# =============================================================================
# Let's recognize the sampling space (as simplex) and treat it accordingly
# =============================================================================
comp1 = [d_rainy, d_nonrainy]
comp2 = [m_rainy, m_nonrainy]

d_clr = coda.clr_transformation(np.asarray(comp1).T)
m_clr = coda.clr_transformation(np.asarray(comp2).T)

# -----------------------------------------------------------------------------
# Let's compute the change in 'rainy days in a year'.
# -----------------------------------------------------------------------------
d_clr_change = d_clr[1, 0] - d_clr[0, 0]
m_clr_change = m_clr[1, 0] - m_clr[0, 0]

print("Change in rainy days in a year: [desert | mountain]")
print("    {:.3f}% | {:.3f}%".format(d_clr_change, m_clr_change))

# -----------------------------------------------------------------------------
# Let's compute the non-rainy days increment of ratios (percent change)
# -----------------------------------------------------------------------------
d_clr_increment = (d_nonrainy[1] - d_nonrainy[0]) / d_nonrainy[0] * 100  # %
m_clr_increment = (m_nonrainy[1] - m_nonrainy[0]) / m_nonrainy[0] * 100  # %

print("Increment in non-rainy days in a year: [desert | mountain]")
print("    {:.3f}% | {:.3f}% \n".format(d_increment, m_increment))
