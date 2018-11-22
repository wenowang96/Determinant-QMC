from glob import glob
import sys
import numpy as np
import util


def get_mu_n(path):
    n_sample, sign, density = \
        util.load(path, "meas_eqlt/n_sample", "meas_eqlt/sign",
                        "meas_eqlt/density")
    mask = (n_sample == n_sample.max())
    print(f"complete: {mask.sum()}/{len(n_sample)}")
    sign, density = sign[mask], density[mask]
    nj = util.jackknife(sign, density.sum(1))
    return util.load_firstfile(path, "metadata/mu")[0], nj[0], nj[1]


def main(argv):
    n = 3
    target = float(argv[1])
    data = np.array([get_mu_n(path) for path in argv[2:]])
    data = data[np.argsort(data[:, 0])]
    diffs = data[:, 1].copy()
    diffs[0] = 0
    diffs[1:] -= data[:-1, 1]
    print(np.hstack((data, diffs[:, None])))

    data[:, 1] -= target
    c = np.abs(data[:, 1]).argmin()
    if c == 0:
        c = 1
    elif c == data.shape[0] - 1:
        c = data.shape[0] - 2
    p = np.polyfit(data[c-1:c+2, 0], data[c-1:c+2, 1], 2)
    r = np.roots(p)
    print(r[np.abs(r - data[c, 0]).argmin()])

if __name__ == "__main__":
    main(sys.argv)