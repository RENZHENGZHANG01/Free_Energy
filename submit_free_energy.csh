#!/bin/csh
#$ -M rzhang4@nd.edu
#$ -m abe
#$ -pe smp 8
#$ -q  long
#$ -N FE_$MOL

# ORCA runs on the compute node's LOCAL disk, not on /groups.
#
# Measured on one 29-atom molecule, identical input, identical 8 cores, same node,
# the only difference being the working directory:
#     node-local /tmp   25 optimisation cycles, converged, 32 min  (1.30 min/cycle)
#     /groups (NFS)     17 cycles, still running at 175 min       (10.33 min/cycle)
# i.e. 7.9x. ORCA rewrites its integral, density and Hessian scratch continuously, and
# on NFS every one of those reads and writes crosses the network. The energies agree,
# so this costs nothing in accuracy.
#
# It also stops ~145 GB of scratch per 3915 molecules from ever reaching /groups:
# .gbw, .densities and .property.txt are 95% of what ORCA writes and nothing
# downstream reads any of them.

cd $SGE_O_WORKDIR
set WORK = `pwd`
set SCRATCH = /tmp/orca_${MOL}_$$
mkdir -p $SCRATCH

# Always clean up, including on failure -- otherwise a crashed job leaves hundreds of
# MB on a compute node that nothing will ever collect.
onintr cleanup

#################################
# 1. LOAD MODULES
#################################
module load orca/6.1.0
module load conda
source /afs/crc.nd.edu/x86_64_linux/c/conda/24.7.1/etc/profile.d/conda.csh
conda activate rdkit_env

set ORCA = /software/o/orca/6.1.0/orca_6_1_0_linux_x86-64_shared_openmpi418_avx2/orca

# ORCA only parallelises over the cores named in "%pal nprocs" inside the .inp;
# the "-pe smp" request above allocates them but does not tell ORCA about them.
# Passing NSLOTS through to the generators makes the two impossible to disagree.
if ( ! $?NSLOTS ) setenv NSLOTS 8
setenv ORCA_NPROCS $NSLOTS

echo "======================================"
echo "Running FREE ENERGY for MOL = $MOL"
echo "Scratch : $SCRATCH"
echo "Project : $WORK"
echo "======================================"

#################################
# 2. STEP 1: PBE OPT
#################################
cp $WORK/data/opt_inp/${MOL}_step1_pbe_opt.inp $SCRATCH/
cd $SCRATCH
echo "STEP1 PBE OPT: $MOL"
$ORCA ${MOL}_step1_pbe_opt.inp > ${MOL}_step1_pbe_opt.out

# Keep the log and the optimised geometry; everything else is scratch.
cp ${MOL}_step1_pbe_opt.out $WORK/data/opt_out/
if ( -s ${MOL}_step1_pbe_opt.xyz ) cp ${MOL}_step1_pbe_opt.xyz $WORK/data/opt_inp/

grep -q "ORCA TERMINATED NORMALLY" ${MOL}_step1_pbe_opt.out
if ( $status != 0 ) then
    echo "ERROR: PBE OPT failed for $MOL"
    goto cleanup
endif

#################################
# 3. GENERATE STEP2 INPUT (ONLY THIS MOL)
#################################
cd $WORK
python scripts/generate_step2_opt_inp_from_xyz.py --mol $MOL

#################################
# 4. STEP 2: PRODUCTION OPT
#################################
cp $WORK/data/opt_inp/${MOL}_step2_opt.inp $SCRATCH/
cd $SCRATCH
echo "STEP2 wB97X-D3 OPT: $MOL"
$ORCA ${MOL}_step2_opt.inp > ${MOL}_step2_opt.out

cp ${MOL}_step2_opt.out $WORK/data/opt_out/
if ( -s ${MOL}_step2_opt.xyz ) cp ${MOL}_step2_opt.xyz $WORK/data/opt_inp/

grep -q "ORCA TERMINATED NORMALLY" ${MOL}_step2_opt.out
if ( $status != 0 ) then
    echo "ERROR: STEP2 OPT failed for $MOL"
    goto cleanup
endif

#################################
# 5. GENERATE FREQ INPUT (ONLY THIS MOL)
#################################
cd $WORK
python scripts/generate_freq_inp.py --mol $MOL

#################################
# 6. FREQ
#################################
cp $WORK/data/freq_inp/${MOL}_freq.inp $SCRATCH/
cd $SCRATCH
echo "FREQ: $MOL"
$ORCA ${MOL}_freq.inp > ${MOL}_freq.out

cp ${MOL}_freq.out $WORK/data/freq_out/
# The .hess is the one large file worth keeping: it allows the thermochemistry to be
# recomputed (a different temperature, or a quasi-harmonic entropy correction) without
# repeating the Hessian, which is the most expensive step in the whole pipeline.
if ( -s ${MOL}_freq.hess ) cp ${MOL}_freq.hess $WORK/data/freq_out/

grep -q "ORCA TERMINATED NORMALLY" ${MOL}_freq.out
if ( $status != 0 ) then
    echo "ERROR: FREQ failed for $MOL"
    goto cleanup
endif

echo "FINISHED FREE ENERGY FOR $MOL"

cleanup:
cd $WORK
rm -rf $SCRATCH
echo "scratch cleaned: $SCRATCH"
