#!/bin/csh
#$ -M rzhang4@nd.edu
#$ -m abe
#$ -pe mpi-48 48 
#$ -q  hpc
#$ -N FE_$MOL


cd $SGE_O_WORKDIR

#################################
# 1. LOAD MODULES
#################################
module load orca/6.1.0
module load conda
source /afs/crc.nd.edu/x86_64_linux/c/conda/24.7.1/etc/profile.d/conda.csh
conda activate rdkit_env

echo "======================================"
echo "Running FREE ENERGY for MOL = $MOL"
echo "Conda env:"
conda info --envs | grep '*'
echo "======================================"

#################################
# 2. STEP 1: PBE OPT
#################################
set inp1 = data/opt_inp/${MOL}_step1_pbe_opt.inp
set out1 = data/opt_out/${MOL}_step1_pbe_opt.out

echo "STEP1 PBE OPT: $MOL"
orca $inp1 > $out1

# HARD SYNC: wait for xyz
set xyz1 = data/opt_inp/${MOL}_step1_pbe_opt.xyz
while ( ! -s $xyz1 )
    echo "Waiting for $xyz1 ..."
    sleep 5
end

#################################
# 3. GENERATE STEP2 INPUT (ONLY THIS MOL)
#################################
python scripts/generate_step2_opt_inp_from_xyz.py --mol $MOL

#################################
# 4. STEP 2: B3LYP OPT
#################################
set inp2 = data/opt_inp/${MOL}_step2_b3lyp_opt.inp
set out2 = data/opt_out/${MOL}_step2_b3lyp_opt.out

echo "STEP2 B3LYP OPT: $MOL"
orca $inp2 > $out2

grep -q "ORCA TERMINATED NORMALLY" $out2
if ( $status != 0 ) then
    echo "ERROR: B3LYP OPT failed for $MOL"
    exit 1
endif

set xyz2 = data/opt_inp/${MOL}_step2_b3lyp_opt.xyz
while ( ! -s $xyz2 )
    echo "Waiting for $xyz2 ..."
    sleep 5
end

#################################
# 5. GENERATE FREQ INPUT (ONLY THIS MOL)
#################################
python scripts/generate_freq_inp.py --mol $MOL

#################################
# 6. FREQ
#################################
set finp = data/freq_inp/${MOL}_freq.inp
set fout = data/freq_out/${MOL}_freq.out

echo "FREQ: $MOL"
orca $finp > $fout

#################################
# 7. EXTRACT THERMO
#################################
python scripts/extract_thermo.py --mol $MOL

#################################
# 8. COMPUTE DELTA G
#################################
python scripts/compute_deltaG.py --mol $MOL

echo "🎉 FINISHED FREE ENERGY FOR $MOL"
