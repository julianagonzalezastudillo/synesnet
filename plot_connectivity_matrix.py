import os.path
from scipy import io
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import DivergingNorm
import matplotlib.colors
import matplotlib.gridspec as gridspec


path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/'
net_path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/net_metrics/'
if not os.path.exists(net_path):
    os.makedirs(net_path)

n_name = []
with open(path + "BN_Atlas_246_LUT_reoriented.txt", "r") as filestream:
    for line in filestream:
        n_name.append(line.split(",")[0])

subjects = {'synes': np.arange(1, 18),
            'ctr': np.arange(18, 35)}
labels = ['Left', 'Right']
x = [round(len(n_name)/4), round(len(n_name)/4*3)]

# cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["red", "white", "blue"] )
# cmap = matplotlib.colors.ListedColormap(["white", [1., .8, 0.], [1., .4, 0.], (1., 0., 0.)])
# bounds = [0., .25, 0.5, .75, 1.]
# norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
cmap = plt.cm.get_cmap('plasma')
cmap = cmap.reversed()
for sub_type in subjects.keys():
    fig = plt.figure(figsize=(4.4, 5), dpi=500)
    gs1 = gridspec.GridSpec(5, 4, figure=fig)
    gs1.update(wspace=0.05, hspace=0.05)  # set the spacing between axes.

    for sub, i in zip(subjects[sub_type], range(20)):
        sub_file = 'CorrMatrix_Subject{0}.mat'.format(str(sub).zfill(3))
        print(sub_file)

        # load correlation
        fc_file = path + 'symetrical_corr_mat/' + sub_file
        fc = io.loadmat(fc_file)

        # connectivity matrix
        Xfc = abs(fc['CorrMatrix'])
        np.fill_diagonal(Xfc, 0)

        ax = plt.subplot(gs1[i])
        im = ax.imshow(Xfc, cmap=cmap)
        ax.set_title('Sub{0}'.format(str(sub).zfill(3)), fontsize=5, y=0.8)

        # thick line between the large cells
        # use clip_on=False and hide the spines to avoid that the border cells look different
        ax.axvline(123 - 0.5, color='black', lw=0.4, clip_on=False, linestyle='--')
        ax.axhline(123 - 0.5, color='black', lw=0.4, clip_on=False, linestyle='--')

        # add the column names as labels
        ax.set_xticks([])
        ax.set_yticks([])
        if i == 13 or i == 14 or i == 15 or i == 16:
            plt.xticks(x, labels, fontsize=5)
        if (i/4).is_integer():
            plt.yticks(x, labels, fontsize=5)
        [ax.spines[axis].set_linewidth(0.6) for axis in ['top', 'bottom', 'left', 'right']]

    # borders
    fig.subplots_adjust(bottom=0.01, top=0.99, right=0.99, left=0.07)

    # colorbar
    cbar_ax = fig.add_axes([0.32, 0.09, 0.66, 0.05])
    cbar = fig.colorbar(im, cax = cbar_ax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=6)
    cbar.outline.set_linewidth(0.6)
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.text(0.5, 1.2, sub_type, horizontalalignment='center')
    plt.savefig(os.getcwd() + '/plots/connectivity_matrices_{0}'.format(sub_type), transparent=True)
    plt.show()
