import os.path
from scipy import io
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
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

# num_sub = len(os.listdir(path + 'symetrical_corr_mat'))
subjects = {'synes': np.arange(1, 18),
            'ctr': np.arange(18, 35)}
for sub_type in subjects.keys():
    fig = plt.figure(figsize=(4, 5), dpi=500)
    gs1 = gridspec.GridSpec(5, 4)
    gs1.update(wspace=0.05, hspace=0.05)  # set the spacing between axes.
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["red", "white", "blue"])
    divnorm = DivergingNorm(vmin=-1, vcenter=0, vmax=1)
    for sub, i in zip(subjects[sub_type], range(20)):
        sub_file = 'CorrMatrix_Subject{0}.mat'.format(str(sub).zfill(3))
        print(sub_file)

        # load correlation
        fc_file = path + 'symetrical_corr_mat/' + sub_file
        fc = io.loadmat(fc_file)

        # connectivity matrix
        Xfc = fc['CorrMatrix']
        np.fill_diagonal(Xfc, 0)

        print(np.max(abs(Xfc)))
        ax = plt.subplot(gs1[i])
        # ax = sns.heatmap(Xfc, square = True, cmap = cmap, norm=divnorm)
        im = ax.imshow(Xfc, cmap=cmap, norm=divnorm)
        ax.set_title('Sub{0}'.format(str(sub).zfill(3)), fontsize=5, y=0.8)

        # thick line between the large cells
        # use clip_on=False and hide the spines to avoid that the border cells look different
        ax.axvline(123 - 0.5, color='black', lw=0.4, clip_on=False, linestyle='--')
        ax.axhline(123 - 0.5, color='black', lw=0.4, clip_on=False, linestyle='--')

        #add the column names as labels
        ax.set_xticks([])
        ax.set_yticks([])
        [ax.spines[axis].set_linewidth(0.6) for axis in ['top', 'bottom', 'left', 'right']]

    # borders
    fig.subplots_adjust(bottom=0.01)
    fig.subplots_adjust(top=0.99)
    fig.subplots_adjust(right=0.99)
    fig.subplots_adjust(left=0.01)
    
    #colorbar
    cbar_ax = fig.add_axes([0.28, 0.09, 0.69, 0.05])
    cbar = fig.colorbar(im, cax = cbar_ax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=6)
    cbar.outline.set_linewidth(0.6)
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.text(0, 1.2, sub_type, horizontalalignment='center')
    plt.show()