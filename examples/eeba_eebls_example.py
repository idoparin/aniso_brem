"""Compare isotropic e-e bremsstrahlung: analytic eebls vs eeba integrated over angles.

``eebls`` returns d(sigma)/dk (Haug 1998 eq. 1). ``eeba`` is the double-differential
lab-frame cross section; ``eeba_integrated`` applies 2 pi int eeba d mu times m_e c^2,
matching the ``cs1`` convention and ``eebapx`` / ``eebls``.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from aniso_brem import eebls, eeba_integrated_array, eeb_emin


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--electron-energy",
        type=float,
        default=1500.0,
        help="Electron kinetic energy in keV (default: 1500 = 1.5 MeV)",
    )
    parser.add_argument(
        "--xsteps",
        type=int,
        default=500,
        help="Integration steps for eebls (default: 500)",
    )
    parser.add_argument(
        "--nsteps",
        type=int,
        default=200,
        help="Number of photon energy samples (default: 200)",
    )
    parser.add_argument(
        "--n-mu",
        type=int,
        default=200,
        help="Cos(pitch angle) quadrature points for eeba integration (default: 200)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("eeba_eebls_comparison.pdf"),
        help="Output plot path",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the plot interactively",
    )
    args = parser.parse_args()

    ee = args.electron_energy
    ep_max = min(1.0e4, 0.99 * ee)
    photon_energies = np.logspace(1, np.log10(ep_max), args.nsteps)

    cs_analytic = np.array(
        [eebls(ee, ep, args.xsteps) for ep in photon_energies]
    )
    cs_numeric = eeba_integrated_array(
        np.ascontiguousarray(photon_energies, dtype=np.float64),
        ee,
        args.n_mu,
    )

    fig, axes = plt.subplots(2, 1, figsize=(7, 8), sharex=True,
                              gridspec_kw={"height_ratios": [3, 1]})
    ax, ax_ratio = axes
    lw = 1.8

    valid_a = cs_analytic > 0
    valid_n = cs_numeric > 0
    ax.plot(
        photon_energies[valid_a],
        cs_analytic[valid_a],
        label=r"Analytic eebls (Haug 1998 eq. 1, d$\sigma$/d$k$)",
        lw=lw,
        color="C0",
    )
    ax.plot(
        photon_energies[valid_n],
        cs_numeric[valid_n],
        label=rf"Numeric $\int$ eeba d$\Omega$ ($N_\mu$ = {args.n_mu})",
        lw=lw,
        linestyle="--",
        color="C1",
    )

    ax.axvline(ee, color="gray", lw=0.8, linestyle=":", label=f"Kinematic limit ({ee:.0f} keV)")
    ax.set_ylabel(r"d$\sigma$/d$\epsilon_\gamma$  [cm$^2$ keV$^{-1}$]")
    ax.set_yscale("log")
    ax.legend(frameon=False, fontsize=9)
    ax.set_title(
        f"Isotropic e-e bremsstrahlung  ($E_e$ = {ee/1e3:.2f} MeV)"
    )

    both_valid = valid_a & valid_n & (ee >= 1.05 * np.array([eeb_emin(ep) for ep in photon_energies]))
    ratio = np.where(
        both_valid,
        cs_numeric / np.where(cs_analytic > 0, cs_analytic, np.nan),
        np.nan,
    )
    ax_ratio.plot(photon_energies, ratio, color="C2", lw=lw)
    ax_ratio.axhline(1.0, color="k", lw=0.8, linestyle="--")
    ax_ratio.fill_between(photon_energies, 0.99, 1.01, color="k", alpha=0.08, label="±1%")
    ax_ratio.set_ylabel("Numeric / eebls")
    ax_ratio.set_ylim(
        max(0.85, np.nanmin(ratio[both_valid]) * 0.98),
        min(1.15, np.nanmax(ratio[both_valid]) * 1.02),
    )
    ax_ratio.legend(frameon=False, fontsize=8)
    ax_ratio.set_xlabel(r"Photon energy $\epsilon_\gamma$ [keV]")
    ax_ratio.set_xscale("log")

    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"Wrote {args.output}")

    mask = both_valid & (photon_energies < 0.95 * ee)
    r = ratio[mask]
    print(
        f"Agreement (ee >= 1.05*E_min, ep < 0.95*E_e): "
        f"max dev = {np.abs(r - 1).max()*100:.2f}%, "
        f"mean dev = {np.abs(r - 1).mean()*100:.4f}%"
    )

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
