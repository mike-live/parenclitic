import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
from .transform_data import calc_edge_weight

def plot_kde(kde, p, I, threshold_p, mask, data, ages, x_label, y_label, x_min, x_max):
    x_min = np.min(data[0, :])
    x_max = np.max(data[0, :])
    y_min = np.min(data[1, :])
    y_max = np.max(data[1, :])
    X, Y = np.mgrid[x_min:x_max:500j, y_min:y_max:500j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    Z = np.reshape(kde(positions).T, X.shape) # score_samples
    ZI = np.zeros(X.shape[0] * X.shape[1])

    for i, pt in enumerate(positions.T):
        print(pt.shape)
        ZI[i] = calc_edge_weight(pt[0], pt[1], kde, p, I, x_min, x_max)
        if ZI[i] > threshold_p:
            ZI[i] = 1

    ZI = np.reshape(ZI, X.shape)
    params = {'legend.fontsize': 'x-large',
          'figure.figsize': (16, 9),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
    pylab.rcParams.update(params)

    fig, (ax1, ax2) = plt.subplots(1, 2)
    im = ax1.imshow(np.rot90(Z), cmap="viridis",
              extent=[x_min, x_max, y_min, y_max])
    fig.colorbar(im, ax=ax1)
#    ax.plot(data[0, mask], data[1, mask], 'b.', markersize=2, alpha=.5)
#    ax.plot(data[0, np.invert(mask)], data[1, np.invert(mask)], 'r.', markersize=2, alpha=.5)
#    ax.plot(data[0, :], data[1, :], cmap=plt.cm.gist_earth_r, markersize=2, alpha=.5)

    im = ax2.imshow(np.rot90(ZI), cmap="hot",
              extent=[x_min, x_max, y_min, y_max])
    sc = ax2.scatter(data[0, :], data[1, :], c = ages, alpha = .5, cmap='cool')
    ax1.axis('equal')
    ax2.axis('equal')
    ax1.set_xlim([x_min, x_max])
    ax1.set_ylim([y_min, y_max])
    ax2.set_xlim([x_min, x_max])
    ax2.set_ylim([y_min, y_max])
    fig.colorbar(sc, ax=ax2)

    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label)
    ax2.set_xlabel(x_label)
    ax2.set_ylabel(y_label)
    ax1.set_title('KDE for age < 65')
    
    plt.show()
    return fig


def plot_parenclitic_kdes(kde_path, X, ages, mask, kdes, p, I, genes_names, threshold_p = 0.5, min_x = 0, max_x = 1):
    k = X.shape[1]
    for i in range(k):
        for j in range(k):
            if (i == j): continue
            data = np.array([X[:, i], X[:, j]])
            print(data.shape)
            kde = kdes[i, j]
            gene_i = genes_names[i]
            gene_j = genes_names[j]
            fig = plot_kde(kde, p[i, j], I[i, j], threshold_p, mask, data, ages, gene_i, gene_j, min_x, max_x)
            fig.savefig(kde_path + 'kde_' + gene_i + '_' + gene_j + '.png')
            plt.close(fig)
