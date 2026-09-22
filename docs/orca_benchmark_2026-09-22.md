# Choosing the level of theory: ORCA benchmark, 2026-09-22

Why this pipeline runs `wB97X-D3 / def2-TZVP / DefGrid3 / RIJCOSX`, measured rather
than assumed. The keyword lines these numbers justify live in
[`scripts/orca_settings.py`](../scripts/orca_settings.py); nothing else defines them.

## Setup

| | |
|---|---|
| molecule | `Cc1ccc(C)c(Oc2ccccc2)c1` — a diaryl ether, 29 atoms, 549 basis functions |
| geometry | one shared geometry for every frequency run (B3LYP D3BJ/def2-TZVP TightOpt), so the comparison isolates the integral treatment |
| hardware | 16 cores, one node, **node-local scratch** |
| ORCA | 6.1.0 |

The shared geometry is correct for comparing timings and integral error. It is *not*
valid for comparing imaginary-frequency counts between functionals, since a geometry
optimised with B3LYP is not a stationary point for ωB97X-D3 — the ωB97X rows below
come from a separate paired optimisation + frequency run.

## Results

| configuration | wall | rel. | imag. | Gibbs error vs exact |
|---|---|---|---|---|
| B3LYP RIJCOSX DefGrid2 | 6.2 min | 1.0× | 0 | −0.127 kcal/mol |
| B3LYP RIJCOSX DefGrid3 | 14.6 min | 2.4× | 0 | −0.162 kcal/mol |
| **B3LYP exact 4-centre** | **109.0 min** | **17.7×** | 0 | reference |
| B3LYP RIJK | — | — | — | **refused by ORCA** |
| ωB97X-D3 RIJCOSX DefGrid2 | 7.5 min | 1.2× | **1** | — |
| ωB97X-D3 RIJCOSX DefGrid3 | 17.2 min | 2.8× | **0** | — |

Frequency error against the exact-integral reference:

| | all modes (median / max) | modes > 100 cm⁻¹ (median / max) |
|---|---|---|
| COSX DefGrid2 | 0.32 / 10.1 cm⁻¹ | 0.35 / 5.5 cm⁻¹ |
| COSX DefGrid3 | 0.28 / 9.2 cm⁻¹ | 0.34 / 1.9 cm⁻¹ |

Thermochemical components (kcal/mol, vs exact):

| | ZPE | entropy term | Gibbs |
|---|---|---|---|
| COSX DefGrid2 | +0.0420 | −0.2121 | −0.1270 |
| COSX DefGrid3 | +0.0270 | −0.1767 | −0.1616 |

## What the numbers decide

**RIJCOSX is not a preference, it is the only option.** ORCA refuses RIJK outright:

```
WARNING: Analytical Hessian not available with RIJK approximation
If you want approximation for Coulomb AND Exchange please choose RIJCOSX!
```

and exact 4-centre integrals cost 17.7× at 29 atoms, which is hopeless for the
55–140 atom molecules in this campaign. The price of COSX is 0.13–0.16 kcal/mol in
Gibbs and ~0.3 cm⁻¹ median in the frequencies — an order of magnitude below both the
DFT error itself and the ~3 kcal/mol that counts as a meaningful difference in the
residual target.

**DefGrid3 buys Hessian stability, not energy accuracy.** It barely moves Gibbs
(−0.127 → −0.162 kcal/mol, i.e. slightly *further* from exact, well inside the noise)
but it cuts the worst frequency error on real vibrational modes from 5.5 to 1.9 cm⁻¹,
and it is what removed the spurious imaginary mode from the ωB97X-D3 run. COSX's grid
noise shows up in the low-frequency modes, and the low-frequency modes are what the
entropy term is most sensitive to.

**ωB97X-D3 is affordable and — the thing that had to be verified before committing —
ORCA does support an analytic Hessian for it.** It costs 1.2× B3LYP at the same grid.
Its absolute energies sit ~90 kcal/mol below B3LYP's for this molecule; that is a
functional offset, not an error, and it is exactly why `atom_ref.csv` has to be
recomputed with ωB97X-D3 too.

**def2-TZVP, not def2-TZVPPD.** The diffuse shells in TZVPPD matter for anions,
Rydberg states and polarisabilities. These are neutral closed-shell monomers, so the
diffuse functions buy essentially nothing while adding basis functions and making the
SCF harder. (The atom references previously ran at TZVPPD while the molecules ran at
TZVP — an inconsistency that would have silently corrupted every ΔG. It is now
impossible: both sides import the same string, and `parse_atom_ref.py` re-reads the
level of theory echoed in each output and refuses to write `atom_ref.csv` on a
mismatch.)

## Separate result: node-local scratch

Same input, same 8 cores, same node, differing only in the working directory:

| working directory | progress |
|---|---|
| node-local `/tmp` | 25 optimisation cycles, converged, 32 min (1.30 min/cycle) |
| `/groups` (NFS) | 17 cycles, still running at 175 min (10.33 min/cycle) |

**7.9×.** ORCA rewrites its integral, density and Hessian scratch continuously, and on
NFS every one of those reads and writes crosses the network. The energies agree, so
this costs nothing in accuracy. It also keeps ~145 GB of scratch per 3915 molecules
off `/groups`, since `.gbw`, `.densities` and `.property.txt` are ~95% of what ORCA
writes and nothing downstream reads any of them. Implemented in
[`submit_free_energy.csh`](../submit_free_energy.csh).

## Net effect on the campaign

| | factor |
|---|---|
| ωB97X-D3 + DefGrid3 + tightened convergence | 2.8× slower |
| node-local scratch | 7.9× faster |
| **combined** | **≈2.8× faster** |

For 3915 molecules: ~14,300 CPU-hours → ~5,000 CPU-hours, and ~3 TB of `/groups`
traffic → ~17 GB of retained output.

## Raw data

Archived at `../../orca_bench_2026-09-22.tar.gz` (inputs, job scripts and ORCA
outputs; the multi-GB scratch was node-local and was never retained).
