##!/usr/bin/env bash
set -euo pipefail

BASE="/path/to/project_dir" 
TMPDIR="$BASE/tmp"
REF="/path/to/project_dir/reference" 
REF_IMG="/path/to/project_dir/reference/img/file"

# Select full or filtered bed file and change the output names _full or _filtered 
BED="/path/to/bed"

# Resources known sites hg38
KNOWN_SITES_indels=$BASE/known_sites_hg38/hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz
KNOWN_SITES_mills=$BASE/known_sites_hg38/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
KNOWN_SITES_dbsnp=$BASE/known_sites_hg38/hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz
KNOWN_SITES_hapmap=/$BASE/known_sites_hg38/hapmap_3.3.hg38.vcf.gz
KNOWN_SITES_omni=$BASE/known_sites_hg38/1000G_omni2.5.hg38.vcf.gz
KNOWN_SITES_1000G=$BASE/known_sites_hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz

######################################### Preprocessing ##################################################################

# -------------------- preprocessing steps for WES pipeline --------------------------------------------
# ----------- SAMPLE INFORMATION --------------------------------------------
# Define Patient and Sample IDs
patient="patientID"
SAMPLE="sample_name"
FASTQ_R1=$BASE/$patient/${SAMPLE}_1.fastq.gz
FASTQ_R2=$BASE/$patient/${SAMPLE}_2.fastq.gz
TRIM_R1=$BASE/$patient/${SAMPLE}_1_trimmed.fastq.gz
TRIM_R2=$BASE/$patient/${SAMPLE}_2_trimmed.fastq.gz

# ----------- Conda environments -----------------------------------------
# conda activate gatk

# ----------- Enable logging --------------------------------------------
exec > >(tee -a "$BASE/$patient/${SAMPLE}_pipeline.log") 2>&1
echo "=== Pipeline started for sample $SAMPLE at $(date) ==="
THREADS=48


# ----------- 1. FastQC  ------------------------------------------
echo "Step 1: FastQC on raw reads"
mkdir -p $BASE/$patient/01_gatk_qc
mkdir -p $BASE/$patient/01_gatk_qc/1_fastQC_raw_reads
mkdir -p $BASE/$patient/01_gatk_qc/2_trimmed
mkdir -p $BASE/$patient/01_gatk_qc/3_fastQC_on_trimmed
fastqc -t $THREADS -o $BASE/$patient/01_gatk_qc/1_fastQC_raw_reads $FASTQ_R1 $FASTQ_R2

# ----------- 2. Trimming ------------------------------------------
mkdir -p $BASE/$patient/01_gatk_qc/2_trimmed

fastp \
  -i $FASTQ_R1 \
  -I $FASTQ_R2 \
  -o $TRIM_R1 \
  -O $TRIM_R2 \
  --detect_adapter_for_pe \
  --trim_front1 15 \
  --trim_front2 15 \
  --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
  --length_required 50 \
  --html $BASE/WES_germline/01_gatk_WES_qc/2_trimmed/fastp_${SAMPLE}.html \
  --json $BASE/WES_germline/01_gatk_WES_qc/2_trimmed/fastp_${SAMPLE}.json

# ----------- 3. FastQC on Trimmed ------------------------------------------
echo "Step 3: FastQC on trimmed reads for sample ${SAMPLE}"
fastqc -t $THREADS -o $BASE/$patient/01_gatk_qc/3_fastQC_on_trimmed $TRIM_R1 $TRIM_R2


# ----------- 4. FASTQ to uBAM  ------------------------------------------
mkdir -p $BASE/$patient/02_gatk_unaligned
gatk --java-options "-Xmx32G -Djava.io.tmpdir=$TMPDIR" FastqToSam \
  -F1 $TRIM_R1 \
  -F2 $TRIM_R2 \
  -O $BASE/$patient/02_gatk_unaligned/${SAMPLE}_unaligned_read_pairs.bam \
  -SM $SAMPLE -RG rg11 -PL ILLUMINA -LB lib11 -PU run0011


# ----------- 5. Alignment BWA + sort --------------------------------
mkdir -p $BASE/$patient/03_gatk_aligned
gatk --java-options "-Xmx386G -Djava.io.tmpdir=$TMPDIR" BwaSpark \
  -I $BASE/$patient/02_gatk_unaligned/${SAMPLE}_unaligned_read_pairs.bam \
  -O $BASE/$patient/03_gatk_aligned/${SAMPLE}_aligned_read_pairs.bam \
  -R $REF \
  -image $REF_IMG
 
# ----------- 6. Mark Duplicates ----------------------------------------
mkdir -p $BASE/$patient/04_gatk_duplMarked
gatk --java-options "-Xmx386G -Djava.io.tmpdir=$TMPDIR" MarkDuplicates \
  -I $BASE/$patient/03_gatk_aligned/${SAMPLE}_aligned_read_pairs.bam \
  -O $BASE/$patient/04_gatk_duplMarked/${SAMPLE}_aligned_duplMarked.bam \
  -M $BASE/$patient/04_gatk_duplMarked/${SAMPLE}_duplMarked_metrics.txt

# ----------- 7. Estimate library complexity -------------------------------
gatk --java-options "-Xmx386G -Djava.io.tmpdir=$TMPDIR" EstimateLibraryComplexity \
  -I $BASE/$patient/03_gatk_aligned/${SAMPLE}_aligned_read_pairs.bam \
  -O $BASE/$patient/04_gatk_duplMarked/${SAMPLE}_libComplex_metrics.txt


# ----------- 8. Sorting ------------------------------------
mkdir -p $BASE/$patient/05_gatk_coordSorted
gatk --java-options "-Xmx386G -Djava.io.tmpdir=$TMPDIR" SortSamSpark \
  -I $BASE/$patient/04_gatk_duplMarked/${SAMPLE}_aligned_duplMarked.bam \
  -O $BASE/W$patient/05_gatk_coordSorted/${SAMPLE}_aligned_duplMarked_coordSorted.bam \
  --sort-order coordinate


# ----------- 9. Base Quality Score Recalibration -------------------
mkdir -p $BASE/$patient/06_gatk_baseRecal
gatk --java-options "-Xmx426G -Djava.io.tmpdir=$TMPDIR" BQSRPipelineSpark \
  -I $BASE/$patient/05_gatk_coordSorted/${SAMPLE}_aligned_duplMarked_coordSorted.bam \
  -R $REF \
  --known-sites $KNOWN_SITES_dbsnp \
  --known-sites $KNOWN_SITES_indels \
  --known-sites $KNOWN_SITES_mills \
  -O $BASE/$patient/06_gatk_baseRecal/${SAMPLE}_aligned_duplMarked_coordSorted_BQSR.bam


####################################################################################################

# -------------------- preprocessing steps for Targeted Panel Pipeline ----------------------------
# fastq files preprocesing with GeneGlobe

# ----------- 1. SAMPLE INFORMATION --------------------------------------------
# Define Patient and Sample IDs
patient="patientID"
SAMPLE="sample_name"

# Input BAM file (generated from GeneGlobe)
RAW_BAM_DIR="/path/to/raw_bam_files"
RAW_SAMPLE="$RAW_BAM_DIR/${SAMPLE}.bam"
BAM_RG="$RAW_BAM_DIR/${SAMPLE}_RG.bam"

# ----------- 2. Add Reading Groups --------------------------------------------
echo "Running GATK AddOrReplaceReadGroups for: $SAMPLE"
gatk AddOrReplaceReadGroups \
    -I "$RAW_SAMPLE" \
    -O "$BAM_RG" \
    -ID "$SAMPLE" \
    -LB QIAseq_Targeted_DNA_Pro \
    -PL QIAseq \
    -SM "$SAMPLE" \
    -PU unit1 \
    --SORT_ORDER coordinate
###########################################################################################


########################## Variant Calling ################################################

# ----------- 10. HaplotypeCaller (WES) ----------------------
mkdir -p $BASE/$patient/07_gatk_germVariants
gatk --java-options "-Xmx38G -Djava.io.tmpdir=$TMPDIR" HaplotypeCaller \
  -R $REF \
  -I $BASE/$patient/06_gatk_baseRecal/${SAMPLE}_aligned_duplMarked_coordSorted_BQSR.bam \
  -L $BED \
  -O $BASE/$patient/07_gatk_germVariants/${SAMPLE}_germVariants.vcf.gz \
  --annotation-group StandardAnnotation \
  --annotation-group StandardHCAnnotation \
  -A DepthPerAlleleBySample -A StrandBiasBySample \
  --min-base-quality-score 20 \
  --standard-min-confidence-threshold-for-calling 30.0 \
  --native-pair-hmm-threads 8 \
  --verbosity INFO


# ----------- 10. HaplotypeCaller (Targeted Panel) ----------------------
mkdir -p $BASE/07_gatk_germVariants
VCF_OUT=$BASE/07_gatk_germVariants/${SAMPLE}_germVariants_after_variantCalling.vcf

gatk --java-options "-Xmx38G -Djava.io.tmpdir=$TMPDIR" HaplotypeCaller \
  -R $REF \
  -I $BAM_RG \
  -L $BED \
  -O $BASE/07_gatk_germVariants/${SAMPLE}_germVariants.vcf.gz \
  --annotation-group StandardAnnotation \
  --annotation-group StandardHCAnnotation \
  --pcr-indel-model AGGRESSIVE \
  --min-base-quality-score 20 \
  --standard-min-confidence-threshold-for-calling 30.0 \
  --native-pair-hmm-threads 48 \
  --verbosity INFO
 

########################## Filtering (WES) ##################################################

# ------------------------ VQSR (WES) ------------------------------------------------ 
#  ----------------------- 11A. SNP recalibration --------------------
mkdir -p $BASE/$patient/08_gatk_germVariantsFiltered
echo "Step 11: VariantRecalibration SNP"

gatk --java-options "-Xmx38G -Djava.io.tmpdir=$TMPDIR" VariantRecalibrator \
  -R $REF \
  -V $BASE/$patient/07_gatk_germVariants/${SAMPLE}_germVariants.vcf.gz \
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $KNOWN_SITES_hapmap \
  --resource:omni,known=false,training=true,truth=false,prior=12.0 $KNOWN_SITES_omni \
  --resource:1000G,known=false,training=true,truth=false,prior=10.0 $KNOWN_SITES_1000G \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $KNOWN_SITES_dbsnp \
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
  -mode SNP \
  -O $BASE/$patient/08_gatk_germVariantsFiltered/${SAMPLE}_snps.recal \
  --tranches-file $BASE/$patient/08_gatk_germVariantsFiltered/${SAMPLE}_snps.tranches 
 
# ---------------- 11Β.  ApplyVQSR for SNPs ---------------------------
echo "Step 12a: Apply SNP Recalibration"

gatk --java-options "-Xmx38G -Djava.io.tmpdir=$TMPDIR" ApplyVQSR \
  -R $REF \
  -V $BASE/$patient/07_gatk_germVariants/${SAMPLE}_germVariants.vcf.gz \
  -O $BASE/$patient/08_gatk_germVariantsFiltered/${SAMPLE}_recalibrated_snps.vcf.gz \
  --truth-sensitivity-filter-level 99.9 \
  --tranches-file $BASE/$patient/08_gatk_germVariantsFiltered/${SAMPLE}_snps.tranches \
  --recal-file $BASE/$patient/08_gatk_germVariantsFiltered/${SAMPLE}_snps.recal \
  -mode SNP

# ----------------------- 12A. INDEL recalibration --------------------
mkdir -p $BASE/$patient/08_gatk_germVariantsFiltered
echo "Step 12: VariantRecalibration INDEL"

gatk --java-options "-Xmx38G -Djava.io.tmpdir=$TMPDIR" VariantRecalibrator \
  -R $REF \
  -V $BASE/$patient/07_gatk_germVariants/${SAMPLE}_germVariants.vcf.gz \
  --resource:mills,known=false,training=true,truth=true,prior=12.0 $KNOWN_SITES_mills \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $KNOWN_SITES_dbsnp \
  -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
  -mode INDEL \
  -O $BASE/$patient/08_gatk_germVariantsFiltered/${SAMPLE}_indels.recal \
  --tranches-file $BASE/$patient/08_gatk_germVariantsFiltered/${SAMPLE}_indels.tranches
 

# ---------------- 12B.  ApplyVQSR for INDELs ---------------------------
echo "Step 12b: Apply INDEL Recalibration"
gatk --java-options "-Xmx38G -Djava.io.tmpdir=$TMPDIR" ApplyVQSR \
  -R $REF \
  -V $BASE/$patient/08_gatk_germVariantsFiltered/${SAMPLE}_recalibrated_snps.vcf.gz \
  -O $BASE/$patient/08_gatk_germVariantsFiltered/${SAMPLE}_germVariants_filtered.vcf.gz \
  --truth-sensitivity-filter-level 99.9 \
  --tranches-file $BASE/$patient/08_gatk_germVariantsFiltered/${SAMPLE}_indels.tranches \
  --recal-file $BASE/$patient/08_gatk_germVariantsFiltered/${SAMPLE}_indels.recal \
  -mode INDEL

############################# Filtering (Targeted Panel)#############################

# ------------------------ 11. Hard Filtering (for targeted panels) ------------------------ 
echo "Step 11: Hard filtering (SNPs + INDELs)..."
mkdir -p $BASE/$patient/08_gatk_germVariantsFiltered
gatk VariantFiltration \
  -V $BASE/$patient/07_gatk_germVariants/${SAMPLE}_germVariants.vcf.gz \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -O $BASE/$patient/08_gatk_germVariantsFiltered/${SAMPLE}_germVariants_filtered.vcf.gz

################################### Select PASS in ROI #####################################

# ----------- 13. Extract PASS variants ---------------------------
echo "Step 13: Final filtering for PASS variants only"

mkdir -p $BASE/$patient/09_gatk_PASS
gatk --java-options "-Xmx38G -Djava.io.tmpdir=$TMPDIR" SelectVariants \
  -R $REF \
  -V $BASE/$patient/08_gatk_germVariantsFiltered/${SAMPLE}_germVariants_filtered.vcf.gz \
  -L $BED \
  --exclude-filtered \
  -O $BASE/$patient/09_gatk_PASS/${SAMPLE}_germVariants_PASS_inROI.vcf.gz

echo "Pipeline completed for sample $SAMPLE"
