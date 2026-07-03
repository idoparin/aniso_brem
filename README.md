# aniso-brem

Routines to calculate Hard X-ray and continuum gamma-ray flux using electron-ion and electron-electron bremsstrahlung cross sections (Fortran + Python interface).

## Install

Requires Python >= 3.9, NumPy, and gfortran:

```bash
pip install -e .
# or
bash compex
```

On Python >= 3.12, the legacy `f2py -c` path (`bash compex --inplace`) needs
`meson` and `ninja` (`pip install meson ninja` or `pip install -e ".[dev]"`).
The default `compex` / `pip install -e .` build uses CMake and does not need meson.

Optional test dependencies:

```bash
pip install -e ".[test,examples]"
```

## Third-party code

`fortran/dilog490.f` is ACM Algorithm 490 (Ginsberg & Zaborowski, 1975). See `THIRD_PARTY_NOTICES.md`.
