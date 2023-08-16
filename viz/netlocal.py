import numpy as np
import matplotlib.pyplot as plt


def plot_3d_local_metric(X, xyz, n_name, cmap='RdYlBu', return_scatter=False, **kwargs):
    """
    Generate a 3D scatter plot with variable point sizes and colors.

    Arguments:
        X: local metric vector
        xyz: numpy array of xyz coordinates
        n_names: list of strings representing the names of each point
        cmap: colormap to use (default is 'RdYlBu')
    """

    fig = plt.figure(figsize=(8, 8), dpi=400)
    ax = fig.add_subplot(111, projection='3d')

    x = xyz[:, 0]
    y = xyz[:, 1]
    z = xyz[:, 2]
    scatter = ax.scatter(x, y, z, s=abs(X/max(abs(X)))*80, c=X, cmap=cmap,
                         alpha=0.8, linewidth=0.5, **kwargs)
    for i in range(len(n_name)):
        ax.text(x[i], y[i], z[i], n_name[i], fontdict={'size': 3, 'color': 'grey'})

    # Show the plot
    ax.view_init(elev=90, azim=-90)
    ax.grid(False)
    ax.axis('off')

    # Set the background color to transparent
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    # Adjust the spacing between the plot and the colorbar
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.83, 0.25, 0.03, 0.5])

    # Add a colorbar to the plot
    cbar = fig.colorbar(scatter, cax=cbar_ax, )
    cbar.ax.tick_params(labelsize=10)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(labelsize=8)
    fig.tight_layout()
    # Set the title of the plot
    # title = ax.set_title( '{0}{1}'.format( X_name, corr_type), fontsize=10)
    # title.set_position( [0.5, -0.9] )

    if return_scatter:
        return fig, ax, scatter, cbar

    return fig, ax
