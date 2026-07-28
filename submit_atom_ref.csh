#!/bin/csh
#$ -M rzhang4@nd.edu
#$ -m ae
#$ -pe smp 2
#$ -q long
#$ -N atom_ref
#$ -o atom_ref_job.log
#$ -j y
#
# Isolated-atom reference energies for BOTH method sets, run serially in one job
# (single atoms are cheap -- SECONDS each serial; 26 jobs is not worth queueing).
# The ORCA inputs deliberately have NO %pal: an atom has ~10-30 basis functions, so
# MPI overhead dominates -- an earlier "%pal nprocs 8" attempt spun 20 min on MPI
# "Read -1" errors and completed only 2/13. Serial is both correct and far faster.
#
#   b3lyp : B3LYP D3BJ def2-TZVP        <- CONSISTENT with the existing 1277 molecules
#   wb97x : wB97X-D3 def2-TZVPPD        <- future use only; do NOT mix with B3LYP data
#
# Inputs come from scripts/generate_atom_ref_inp.py (ground-state multiplicities baked in).

cd $SGE_O_WORKDIR

module load orca/6.1.0
set ORCA = /software/o/orca/6.1.0/orca_6_1_0_linux_x86-64_shared_openmpi418_avx2/orca

echo "======================================"
echo "atom_ref: isolated-atom references"
echo "Host: `hostname`   Date: `date`"
echo "======================================"

foreach tag ( b3lyp wb97x )
  set d = data/atom_ref/$tag
  echo ""
  echo "########## $tag ##########"
  if ( ! -d $d ) then
    echo "  missing $d -- run scripts/generate_atom_ref_inp.py first"
    continue
  endif
  foreach inp ( $d/*.inp )
    # use csh's :r (root) modifier -- a sed '$' inside backticks would be parsed
    # by csh as a variable name and abort with "Illegal variable name".
    set base = $inp:r
    set out = ${base}.out
    echo "--- `basename $inp` : `head -1 $inp` ---"
    $ORCA $inp > $out
    grep -E "FINAL SINGLE POINT ENERGY|Final Gibbs free energy" $out | tail -2
  end
end

echo ""
echo "======================================"
echo "Parsing results"
echo "======================================"
module load conda
source /afs/crc.nd.edu/x86_64_linux/c/conda/24.7.1/etc/profile.d/conda.csh
conda activate rdkit_env
python scripts/parse_atom_ref.py

echo ""
echo "Review the table above, then write the production file with:"
echo "   python scripts/parse_atom_ref.py --write-production --use G"
echo "atom_ref job finished: `date`"
