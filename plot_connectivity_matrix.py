"""
=================================
            SYNESNET
=================================
Plot coreness results in .png and also save for complementary 3D oplot in matlab
"""
from scipy import io
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from tools import load_node_names
from config import DATA_DIR


# Node names
n_name, n_name_full = load_node_names()

# Subject types
subjects = {"synes": np.arange(1, 18), "ctr": np.arange(18, 35)}
labels = ["Left", "Right"]
x = [round(len(n_name) / 4), round(len(n_name) / 4 * 3)]

for sub_type in subjects.keys():
    fig = plt.figure(figsize=(4.4, 5), dpi=500)
    gs1 = gridspec.GridSpec(5, 4, figure=fig)
    gs1.update(wspace=0.05, hspace=0.05)  # set the spacing between axes.

    for sub, i in zip(subjects[sub_type], range(20)):
        # load connectivity matrix
        sub_file = f"CorrMatrix_Subject{str(sub).zfill(3)}.mat"
        fc_file = DATA_DIR / f"symetrical_corr_mat/{sub_file}"
        fc = io.loadmat(fc_file)

        # Xfc = abs(fc['CorrMatrix'])
        Xfc = fc["CorrMatrix"]
        np.fill_diagonal(Xfc, 0)

        # Plot
        ax = plt.subplot(gs1[i])
        im = ax.imshow(Xfc, cmap="RdYlBu", vmax=np.max(Xfc), vmin=-np.max(Xfc))
        ax.set_title("Sub{0}".format(str(sub).zfill(3)), fontsize=5, y=0.8)

        # thick line between the large cells
        # use clip_on=False and hide the spines to avoid that the border cells look different
        ax.axvline(123 - 0.5, color="black", lw=0.4, clip_on=False, linestyle="--")
        ax.axhline(123 - 0.5, color="black", lw=0.4, clip_on=False, linestyle="--")

        # add the column names as labels
        ax.set_xticks([])
        ax.set_yticks([])
        if i == 13 or i == 14 or i == 15 or i == 16:
            plt.xticks(x, labels, fontsize=5)
        if (i / 4).is_integer():
            plt.yticks(x, labels, fontsize=5)
        [
            ax.spines[axis].set_linewidth(0.6)
            for axis in ["top", "bottom", "left", "right"]
        ]

    # borders
    fig.subplots_adjust(bottom=0.01, top=0.99, right=0.99, left=0.07)

    # colorbar
    cbar_ax = fig.add_axes([0.32, 0.09, 0.66, 0.05])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
    cbar.ax.tick_params(labelsize=6)
    cbar.outline.set_linewidth(0.6)
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.text(0.5, 1.2, sub_type, horizontalalignment="center")
    # plt.savefig(DATA_DIR / f"/plots/connectivity_matrices_{sub_type}", transparent=True)
    plt.show()
