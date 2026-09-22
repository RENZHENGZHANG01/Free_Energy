#!/bin/csh
#$ -M rzhang4@nd.edu
#$ -m ae
#$ -pe smp 2
#$ -q long
#$ -N atom_ref
#$ -o atom_ref_job.log
#$ -j y
#
# Isolated-atom reference energies: G(atom) for every element in the pool vocabulary.
#
#   Delta_G = G(molecule) - SUM_i n_i * G(atom_i)
#
# so these MUST be at the same level of theory as the molecular frequency runs. The
# keyword line is not written here or in the .inp generator -- both come from
# scripts/orca_settings.py, and scripts/parse_atom_ref.py re-reads the level of
# theory echoed in each .out and refuses to write atom_ref.csv if it disagrees.
#
# The ORCA inputs deliberately have NO %pal: an atom has ~10-30 basis functions, so
# MPI overhead dominates -- an earlier "%pal nprocs 8" attempt spun 20 min on MPI
# "Read -1" errors and completed only 2 of 13. Serial is both correct and far faster
# (seconds per atom), which is why all 13 run serially inside this one job.
#
# Inputs come from scripts/generate_atom_ref_inp.py (ground-state multiplicities
# baked in -- see its docstring; a wrong multiplicity is silent and poisons every
# molecule containing that element).

cd $SGE_O_WORKDIR

module load orca/6.1.0
set ORCA = /software/o/orca/6.1.0/orca_6_1_0_linux_x86-64_shared_openmpi418_avx2/orca

echo "======================================"
echo "atom_ref: isolated-atom references"
echo "Host: `hostname`   Date: `date`"
echo "======================================"

# Every subdirectory holding .inp files. Normally just "production"; "legacy_b3lyp"
# appears only if generate_atom_ref_inp.py was run with --legacy-b3lyp.
foreach d ( data/atom_ref/*/ )
  set tag = `basename $d`
  set n = `ls $d/*.inp |& grep -c '\.inp'`
  if ( $n == 0 ) continue
  echo ""
  echo "########## $tag ($n atoms) ##########"
  foreach inp ( $d/*.inp )
    # csh's :r (root) modifier -- a sed '$' inside backticks would be read by csh
    # as a variable name and abort with "Illegal variable name".
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

# --use G: the molecular side of Delta_G is a Gibbs free energy, so the atom side
# must be one too. --write-production also runs the keyword-consistency check and
# backs up the previous atom_ref.csv before overwriting.
python scripts/parse_atom_ref.py --write-production --use G

echo ""
echo "atom_ref job finished: `date`"
