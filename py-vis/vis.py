import numpy as np
import matplotlib.pyplot as plt
import struct


def read_floats_from_bin(path: str, num_values: int) -> np.array:
    with open(path, "rb") as fin:
        return np.array(struct.unpack("f" * num_values, fin.read(4 * num_values)))


def corr_plot(
    b1: str,
    b2: str,
    num_markers: int,
    # outpath: str
):
    num_values = num_markers * (num_markers - 1) / 2
    v1 = read_floats_from_bin(b1, num_values)
    v2 = read_floats_from_bin(b2, num_values)
    corr = np.corrcoef(v1, v2)[0, 1]
    diag = np.linspace(-1, 1, 10)
    props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)

    plt.figure()
    plt.plot(diag, diag, "k--")
    plt.plot(v1, v2, "x")
    plt.xlabel(r"$\rho$")
    plt.ylabel(r"$\sin(\pi / 2 \tau_B)$")
    ax = plt.gca()
    ax.text(
        0.05,
        0.95,
        r"$\rho={}$".format(corr),
        transform=ax.transAxes,
        fontsize=14,
        verticalalignment="top",
        bbox=props,
    )
    plt.tight_layout()
    # plt.savefig(outpath)
