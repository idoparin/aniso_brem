"""Compare electron-ion and electron-electron bremsstrahlung cross sections."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from aniso_brem import cspe, eebls, eebapx, eeblser

try:
    from utils import bremsstrahlung_cross_section

    HAS_EI_APPROX = True
except ImportError:
    HAS_EI_APPROX = False


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--electron-energy",
        type=float,
        default=1.0e4, # 1.5e4,
        help="Electron kinetic energy in keV (default: 15000)",
    )
    parser.add_argument(
        "--nsteps",
        type=int,
        default=100,
        help="Number of photon energy samples (default: 100)",
    )
    parser.add_argument(
        "--xsteps",
        type=int,
        default=500,
        help="Integration steps for eebls (default: 500)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("cross_section_comparison.pdf"),
        help="Output plot path",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the plot interactively",
    )
    args = parser.parse_args()

    electron_energy = args.electron_energy
    photon_energies = np.logspace(1, 4, args.nsteps)

    eics = []
    eiappx = []
    eecs = []
    eeappx = []
    eeultrarel = []

    for photon_energy in photon_energies:
        eics.append(cspe(electron_energy, photon_energy))
        eecs.append(eebls(electron_energy, photon_energy, args.xsteps))
        eeappx.append(eebapx(electron_energy, photon_energy))
        eeultrarel.append(eeblser(electron_energy, photon_energy))

        if HAS_EI_APPROX:
            eiappx.append(
                bremsstrahlung_cross_section(
                    np.array([electron_energy]),
                    np.array([photon_energy]),
                    z=1.0,
                )[0]
            )

    eics = np.asarray(eics)
    eecs = np.asarray(eecs)
    eeappx = np.asarray(eeappx)
    eeultrarel = np.asarray(eeultrarel)

    fig, ax = plt.subplots()
    lw = 1.5

    ax.plot(photon_energies, eics, label="e-i CS", linewidth=lw)
    if HAS_EI_APPROX:
        ax.plot(
            photon_energies,
            np.asarray(eiappx),
            label="e-i APPROX CS",
            linewidth=lw,
            linestyle="--",
            color="red",
        )
    ax.plot(photon_energies, eecs, label="e-e CS", linewidth=lw, color="green")
    ax.plot(
        photon_energies,
        eeappx,
        label="e-e APPROX CS",
        linewidth=lw,
        color="green",
        linestyle="--",
    )
    ax.plot(
        photon_energies,
        eeultrarel,
        label="e-e ULTRAREL CS",
        linewidth=lw,
        color="green",
        linestyle=":",
    )

    ax.set_xlabel("Photon energy, keV")
    ax.set_ylabel("Cross-section per photon energy, cm$^2$/keV")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(frameon=False)
    ax.set_title(
        f"Bremsstrahlung cross-sections (electron energy {electron_energy / 1e3:.2f} MeV)"
    )

    fig.savefig(args.output)
    print(f"Wrote {args.output}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
