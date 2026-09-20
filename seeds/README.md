# Iteration-0 seed set

`seed_v9.csv` — 3,915 monomers selected from a 2,057,755-molecule pool
(PI1M + OMG + omics + PolyInfo, merged and filtered). This is the input list for the
first DFT campaign: ~14,300 CPU-hours, about 6 days at 100 concurrent jobs.

## Columns

| column | meaning |
|---|---|
| `rank` | row id, no ordering significance (the file is shuffled) |
| `db` | source database |
| `smiles` | polymer SMILES with `*` connection points |
| `canon` | canonical SMILES after `*` → C end-capping |
| `component` | which selection component chose this molecule (see below) |
| `n_atoms` | atom count including hydrogens |
| `est_cost_h` | predicted DFT wall time, from `t ≈ 5.82e-5 · N^2.69` fitted on 220 finished jobs |
| `is_validation` | **True for the 1,000 `random` rows — hold these out of GNN training** |

## Why four components

Iteration 0 has to satisfy four requirements at once, and measurement showed no single
strategy satisfies more than one of them, so the seed is built as four blocks and every
row records which block it came from.

| component | n | job | evidence |
|---|---|---|---|
| `coverage` | 1,999 | cover chemical space so later AL rounds extend the frontier instead of rediscovering basics | greedy max-coverage at τ=0.4 reaches 90.9% of the pool; 2,000 random molecules reach ~79% |
| `random` | 1,000 | the only unbiased yardstick, and it keeps the model calibrated to the pool distribution the AL acquisition depends on | every other block is selected by a criterion correlated with the model or the coverage |
| `dopt` | 500 | make the composition baseline's coefficients well-determined | sequential D-optimality; at N=500 the baseline's held-out MAE was 11.73 vs 13.16 (random) and 17.03 (coverage-greedy) |
| `feature` | 416 | guarantee every FS5 feature has enough data to pin its coefficient | a feature with 3 supporting molecules leaks 3.26 kcal/mol into its neighbours when one DFT run is wrong by 10; at 12 supporters it leaks <0.5 |

`coverage` and `random` are poor at determining the baseline (leverage efficiency 0.22x
and 0.35x), while `dopt` and `feature` are poor at covering chemical space (the feature
block alone covers 2.2%). That is why both exist.

## Pool filtering

Element whitelist `C H O N F Cl Br S P Si I`, plus eleven substructure rules removing
generative-database artefacts: Si/Ge/Sn multiple bonds, C=I, hypervalent iodine, I-P,
I-N, aromatic C-P, aromatic peroxide, and 2H/3H isotopes. 17,398 molecules (0.84%) were
dropped. See `scripts/build_clean_pool.py`, which self-tests every rule against a known
positive example — the isotope rule was originally written `[2H,3H]`, which parses
cleanly and never matches, and had been silently passing deuterated structures through.

Ge and Sn were removed deliberately: the downstream generative model targets organic
polymers and will never emit them, while the two highest-leverage components were
spending 50-121x their pool share chasing them. Si (polysiloxanes) and P
(polyphosphazenes) were kept — both are backbone elements of real polymer families.
The tell for an artefact is the source split: Ge appears 2,020 times in PI1M
(RNN-generated) against 14 in PolyInfo (real, synthesised polymers).

## Known limits

Three features still sit below 12 supporters — `BP_N-N_TRIPLE` (3), `BP_I-O_SINGLE` (7),
`BP_Br-S_SINGLE` (11) — because the pool itself has too few carriers. Molecules
containing them should be treated as out-of-vocabulary downstream; they account for
~0.01% of the pool.
