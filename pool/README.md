# Candidate pool

`pool_molecules.csv.gz` — 2,057,755 monomers, the space active learning selects from.
Merged from PI1M, OMG, omics and PolyInfo, deduplicated on SMILES, then filtered.

Two columns, `db` and `smiles`. Everything else the pipeline uses — canonical SMILES,
atom counts, Morgan fingerprints — is a pure function of `smiles` and is regenerated
rather than stored. That is what keeps this file at 14 MB instead of 4.2 GB: the
fingerprint array alone is 3.92 GB and carries no information the SMILES do not.
Verified on a 200,000-molecule sample: regenerating `canon` and `n_atoms` from `smiles`
reproduces the stored values for 100.0000% of rows.

## Rebuilding the caches

```bash
# fingerprints + metadata (~40 min, 2M molecules)
ALLOW_UNFILTERED_REBUILD=1 python ml/active_learning_v3.py   # writes al_v2_cache/
python ml/build_clean_pool.py                                # applies the filter
```

Reproducibility depends on the RDKit version, because canonical SMILES and aromaticity
perception can change between releases and the difference would be silent. This pool
was built with:

    RDKit 2025.09.2   numpy 2.2.6   pandas 2.3.3

## What was filtered out

Starting from 2,144,924 raw molecules across the four databases:

| removed | why |
|---|---|
| non-whitelisted elements | 27 elements appear in the raw data. The whitelist is `C H O N F Cl Br S P Si I` |
| net-charged species | the DFT pipeline runs every molecule as a neutral singlet, so their ΔG is invalid |
| 17,398 further molecules (0.84%) | eleven substructure rules, below |

The substructure rules target generative-database artefacts: Si/Ge/Sn multiple bonds
(disilenes, germenes, stannenes), C=I, hypervalent iodine, I-P, I-N, aromatic C-P,
aromatic peroxide, and 2H/3H isotopes. `ml/build_clean_pool.py` self-tests every rule
against a known positive example and aborts if one fails to match — the isotope rule
was originally written `[2H,3H]`, which parses cleanly, never matches, and had been
silently passing 12,617 deuterated structures while reporting zero hits.

Deuterium deserves its own note, since deuterated polymers are real. It was removed for
three independent reasons: the FS5 featuriser cannot see it (`GetSymbol()` returns "H",
so a deuterated molecule and its protiated twin have identical feature vectors but
different ΔG); `generate_xyz.py` writes "H" into the .xyz, so ORCA would compute the
protiated molecule anyway; and the source split is 1,215 PI1M against 2 PolyInfo.

## Elements kept and dropped

| element | pool share | decision |
|---|---|---|
| C, O, N, S, F, Cl, Br | 11-100% | keep |
| Si | 4.3% | **keep** — polysiloxanes are a real polymer family |
| P | 2.1% | **keep** — polyphosphazenes have a P=N backbone |
| I | 0.2% | keep |
| Ge, Sn | 0.08%, 0.11% | **drop** |
| Na, Se, B, Fe, As, Zn, Ca, Pb, K, Ni, Te, Cd, Co, Li | <0.26% each | drop |

Ge and Sn were the interesting case. They were inside the original element whitelist,
so the two highest-leverage seed components went hunting for molecules carrying them:
Ge is 0.08% of the pool but 4.0% of the D-optimal block and 9.7% of the feature block,
a 50-121x enrichment, spent on chemistry an organic-polymer generative model will never
emit. The tell that they are artefacts rather than rare-but-real chemistry is the source
split — Ge appears 2,020 times in RNN-generated PI1M against 14 times in PolyInfo, the
database of polymers that have actually been synthesised. Si and P show no such split
and were kept, even though their enrichment is also high: a coefficient needs a roughly
fixed number of samples regardless of how rare the feature is, so over-representation is
the correct price for chemistry you intend to keep.
