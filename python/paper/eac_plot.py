# plotting electrical power against power angle delta combined with the steady turbine power. Thus showing acceleration and deceleration areas

from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
from matplotlib.patches import Polygon

# redefining plot save parameterss
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Charter"],
    "font.size": 17
})

# # update savefig options for latex export
# mpl.use("pgf")

# inputs:
P_m = 0.6
P_max = 2.1
delta_c = 68.01

# power functions
def P_e(x):
    # determing the P_e curve under input in degrees
    return P_max*np.sin(x*2*np.pi/360)

def P_t(x):
    return P_m*np.ones(np.size(x))

# Integral limits
delta_0 = np.arcsin(P_m/P_max)*360/(2*np.pi)
delta_max = 180-np.arcsin(P_m/P_max)*360/(2*np.pi)

# calculation of P_e and P_t
x = np.linspace(0, 180)
y1 = P_e(x)
y2 = P_t(x)

# setting up the plot and plotting P_e and P_t
fig, ax = plt.subplots()
ax.plot(x, y1, linewidth=2, label='$P_\mathrm{e}$')# of the SG')
ax.plot(x, y2, linewidth=2, label='$P_\mathrm{m}$')# of the turbine')
ax.set_ylim(bottom=0)

# Make the shaded region for area_acc
ix1 = np.linspace(delta_0, delta_c)
iy1 = P_t(ix1)
verts = [(delta_0, 0), *zip(ix1, iy1), (delta_c, 0)]
poly = Polygon(verts, facecolor='0.9', edgecolor='0.5')
ax.add_patch(poly)

# Make the shaded region for area_dec, https://matplotlib.org/stable/gallery/lines_bars_and_markers/fill_between_demo.html
ix2 = np.linspace(delta_c, delta_max)
iy2 = P_e(ix2)
ax.fill_between(ix2, iy2, P_m, facecolor='0.9', edgecolor='0.5')
# ax.vlines(x = delta_max, ymin = 0, ymax = P_m, colors = '0.5') # vertical line for marking delta_max

# adding area descriptions
ax.text(0.5*(delta_0 + delta_c), 0.5*(P_m), r"$A_\mathrm{acc}$", horizontalalignment='center', fontsize=17)
ax.text(0.47*(delta_max + delta_c), P_m+0.3*(P_max-P_m), r"$A_\mathrm{dec}$", horizontalalignment='center', fontsize=17)

# axis text and optical manipulation of the plot
# fig.text(0.9, 0.05, '$\delta$ in $°$')
# fig.text(0.1, 0.9, '$P$ in $pu$')

# ax.spines[['top', 'right']].set_visible(False)
ax.set_xticks([0,90,180,delta_0, delta_c, delta_max], labels=[0,90,180,'$\delta_\mathrm{0}$', '$\delta_\mathrm{c}$', '$\delta_\mathrm{max}$'])
ax.set_yticks(np.arange(0, P_max, 0.3))

plt.xlabel("power angle $\delta$ in deg")
plt.ylabel("power $P$ in pu")

# plt.grid()
plt.legend()
# plt.savefig("plots/eac.png")
plt.show()