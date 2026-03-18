#!/bin/bash
#SBATCH --job-name=lovo_burden_test_MAF001
#SBATCH --mail-type=END
#SBATCH --partition=cpuq
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G

INPUT_FILE="/path/to/lovo_jobs.tsv"
OUTPUT_DIR="/path/to/lovo_results"
PGEN_PREFIX="/path/to/pgen_prefix"
ANNO_FILE="/path/to/annot"
SET_LIST="/path/to/setlist"
MASK_DEF="/path/to/masks.txt"
CONTAINER="/path/to/regenie_v4.0.sif"
PGEN_DATASET="/path/to/pgen_dataset"
RAREVAR_FILES="/path/to/rarevar_files"

mkdir -p log_files
mkdir -p "${OUTPUT_DIR}"

tail -n +2 "$INPUT_FILE" | while IFS=$'\t' read -r Gene Mask Threshold trait phenoFile covarFile covarColList phenoColList pred path_validated_input path_regenie_step1_preds; do

    mask_lovo="${Gene},${Mask},${Threshold}"
    output_file="${OUTPUT_DIR}/${Gene}-${Mask}-${Threshold}-${trait}-MAC1"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=lovo_${Gene}_${Mask}_${Threshold}
#SBATCH --mail-type=END
#SBATCH --partition=cpuq
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --output=log_files/regenie_lovo_${Gene}_${Mask}_${Threshold}.log

mkdir -p log_files

module load singularity

singularity exec --cleanenv \
  -B \$path_validated_input \
  -B ${PGEN_DATASET} \
  -B \$path_regenie_step1_preds \
  -B ${RAREVAR_FILES} \
  -B ${OUTPUT_DIR} \
  ${CONTAINER} \
  regenie \
  --step 2 \
  --pgen ${PGEN_PREFIX} \
  --phenoFile \$phenoFile \
  --covarFile \$covarFile \
  --covarColList \$covarColList \
  --phenoColList \$phenoColList \
  --minMAC 1 \
  --bsize 400 \
  --threads 8 \
  --pred \$pred \
  --anno-file ${ANNO_FILE} \
  --set-list ${SET_LIST} \
  --mask-def ${MASK_DEF} \
  --vc-tests skat,skato,acatv,acato \
  --out ${output_file} \
  --mask-lovo ${mask_lovo}

EOF

done
