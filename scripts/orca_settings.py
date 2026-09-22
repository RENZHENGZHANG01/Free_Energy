"""
Single source of truth for the ORCA level of theory used by this pipeline.

Every .inp generator imports its keyword line from here. They MUST agree, because

    Delta_G(molecule) = G(molecule) - SUM_i n_i * G(atom_i)

is a difference of two absolute energies. If the molecules and the isolated-atom
references are computed with different functionals, basis sets, grids or integral
approximations, the subtraction does not cancel anything and Delta_G is garbage --
not slightly wrong, but wrong by tens of kcal/mol per atom. Keeping the strings in
one file is what stops that from happening silently when one script is edited and
the others are not.

-------------------------------------------------------------------------------
WHY THESE SETTINGS  (benchmark, 2026-09-22)
-------------------------------------------------------------------------------
Measured on Cc1ccc(C)c(Oc2ccccc2)c1 (29 atoms, 549 basis functions), 16 cores, one
node, node-local scratch, identical geometry for every frequency run:

    frequency job                     wall     rel    imag   Gibbs err vs exact
    RIJCOSX  DefGrid2 (B3LYP)         6.2 min   1.0x    0      -0.127 kcal/mol
    RIJCOSX  DefGrid3 (B3LYP)        14.6 min   2.4x    0      -0.162 kcal/mol
    exact 4-centre    (B3LYP)       109.0 min  17.7x    0      (reference)
    RIJK              (B3LYP)        REFUSED -- see below
    RIJCOSX  DefGrid2 (wB97X-D3)      7.5 min   1.2x    1
    RIJCOSX  DefGrid3 (wB97X-D3)     17.2 min   2.8x    0

RIJCOSX: not a preference, the only option. ORCA refuses RIJK outright --
    "Analytical Hessian not available with RIJK approximation. If you want
     approximation for Coulomb AND Exchange please choose RIJCOSX!"
-- and exact 4-centre integrals cost 17.7x at 29 atoms, which is hopeless for the
55-140 atom molecules in this campaign. The price of COSX is 0.13-0.16 kcal/mol in
Gibbs and ~0.3 cm^-1 median in the frequencies, an order of magnitude below both
the DFT error itself and the ~3 kcal/mol that counts as a meaningful difference in
the residual target.

DefGrid3: bought for Hessian stability, not for energy. It barely moves Gibbs
(-0.127 -> -0.162 kcal/mol, i.e. slightly further from exact, well inside the
noise) but it cuts the worst frequency error on real vibrational modes from 5.5 to
1.9 cm^-1, and it is what removed the spurious imaginary mode from the wB97X-D3
run. Low-frequency modes are where COSX's grid noise shows up, and low-frequency
modes are exactly what the entropy term is most sensitive to.

wB97X-D3: a range-separated hybrid, better for the charge-transfer and long-chain
systems in this pool than B3LYP, and -- the thing that had to be verified before
committing -- ORCA does support an ANALYTIC Hessian for it. It costs only 1.2x
B3LYP at the same grid. Its absolute energies sit ~90 kcal/mol below B3LYP's for
this molecule, which is a functional offset, not an error: it is why atom_ref MUST
be recomputed with wB97X-D3 as well.

def2-TZVP, not def2-TZVPPD: the diffuse functions in TZVPPD matter for anions,
Rydberg states and polarisabilities. These are neutral closed-shell monomers, so
the diffuse shells buy essentially nothing while adding basis functions and making
the SCF harder to converge. (The atom references previously used TZVPPD while the
molecules used TZVP -- an inconsistency that would have corrupted every Delta_G.)

TightOpt on the geometry: a loose minimum leaves residual gradients that turn into
small imaginary frequencies in the Hessian that follows. VeryTightSCF on the
frequency job: the Hessian differentiates the SCF solution, so SCF noise is
amplified into the very low-frequency modes that dominate the entropy.

Step 1 stays PBE/def2-SVP. It is a cheap pre-optimiser whose only job is to hand
step 2 a sane geometry; step 2 re-optimises from scratch, so the functional used
here cannot affect the final numbers.

-------------------------------------------------------------------------------
ONE CHANGE PROPAGATES
-------------------------------------------------------------------------------
Changing anything below invalidates every previously computed Delta_G. The atom
references have to be recomputed with the same string (scripts/generate_atom_ref_inp.py
-> submit_atom_ref.csh -> parse_atom_ref.py --write-production --use G) and the whole
molecular dataset re-run. Do not edit one line "just to try something".
"""
import os

# ---------------------------------------------------------------------------
# Level of theory
# ---------------------------------------------------------------------------
# Cheap pre-optimisation. Output geometry only; never used for energies.
STEP1_OPT = "! PBE def2-SVP TightSCF Opt RIJCOSX def2/J"

# Production geometry optimisation.
STEP2_OPT = "! wB97X-D3 def2-TZVP TightSCF TightOpt DefGrid3 RIJCOSX def2/J"

# Production frequencies / thermochemistry. Same functional, basis and grid as
# STEP2_OPT -- a Hessian evaluated at a different level than the geometry it sits
# on is not a Hessian at a stationary point, and produces imaginary frequencies
# that are pure artefact.
FREQ = "! wB97X-D3 def2-TZVP VeryTightSCF Freq DefGrid3 RIJCOSX def2/J"

# Isolated-atom references. MUST be byte-identical in method to FREQ.
ATOM_REF = FREQ

# Thermochemistry conditions (ORCA defaults, stated explicitly so they are visible
# in the .inp and cannot drift between the molecule and atom runs).
TEMP_K = 298.15

# ---------------------------------------------------------------------------
# Parallelism
# ---------------------------------------------------------------------------
# ORCA only uses the cores named in %pal; the scheduler request alone does nothing.
# Reading NSLOTS (set by SGE from "-pe smp N") makes the two impossible to
# disagree, which is the failure this pipeline has already hit once.
def nprocs():
    for var in ("ORCA_NPROCS", "NSLOTS"):
        v = os.environ.get(var)
        if v and v.strip().isdigit() and int(v) > 0:
            return int(v)
    return 8


def pal_block(n=None):
    """'%pal nprocs N end' block, or '' for a serial run (n <= 1)."""
    n = nprocs() if n is None else n
    return "" if n <= 1 else f"%pal nprocs {n} end\n\n"


# ---------------------------------------------------------------------------
# File naming
# ---------------------------------------------------------------------------
# Method-neutral on purpose. These files used to be called "_step2_b3lyp_opt" and
# kept that name after the functional changed, which is how a directory ends up
# full of wB97X results labelled B3LYP.
STEP1_TAG = "step1_pbe_opt"   # pre-optimisation; PBE is in the name because it is
                              # genuinely fixed and never the production method
STEP2_TAG = "step2_opt"
FREQ_TAG = "freq"
