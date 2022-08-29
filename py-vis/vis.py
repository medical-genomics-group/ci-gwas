import numpy as np
import matplotlib.pyplot as plt
import struct

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc("font", size=SMALL_SIZE)  # controls default text sizes
plt.rc("axes", titlesize=SMALL_SIZE)  # fontsize of the axes title
plt.rc("axes", labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("ytick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("legend", fontsize=SMALL_SIZE)  # legend fontsize
plt.rc("figure", titlesize=BIGGER_SIZE)  # fontsize of the figure title


def read_floats_from_bin(path: str, num_values: int) -> np.array:
    with open(path, "rb") as fin:
        return np.array(struct.unpack("f" * num_values, fin.read(4 * num_values)))


def corr_plot(
    b1: str,
    b2: str,
    num_markers: int,
    title="",
):
    num_values = int(num_markers * (num_markers - 1) / 2)
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
    if title:
        plt.title(title)
    plt.tight_layout()
