import numpy as np

def make_synthetic_move_data(config, X, start_group, finish_group, alpha, relative_noise):
    Xs = []
    for id_start in start_group:
        Xs.append(X[id_start, :])
    for id_finish in finish_group:
        Xs.append(X[id_finish, :])

    for id_start in start_group:
        for id_finish in finish_group:
            sample_sigma = X[id_start, :].std()
            xcur = X[id_start, :] * (1 - alpha) + X[id_finish, :] * alpha
            noise = np.random.randn(*xcur.shape) * relative_noise * sample_sigma
            xcur += noise
            Xs.append(xcur)

    Xs = np.array(Xs)
    mask = np.ones(Xs.shape[0], dtype = np.int8) * 2
    mask[:len(start_group)] = -1
    mask[len(start_group):len(start_group) + len(finish_group)] = +1

    y = mask.copy()
    y[mask == -1] = 0
    y[mask == +1] = 1

    return Xs, y, mask
