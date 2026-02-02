# HiC_icenorm_for_AggregatePeakAnalysis
Aggregate_Peak_Analysis_for_icenorm_HiC

# Aggregate peak analysis (APA) is from https://github.com/ay-lab/Utilities/tree/main/APA


# I have a list of target genes and Hi-C data from different samples 
# I need to test if the promoter region of these target genes is very active in specific samples

# Step 1. Generate bins of BED files covering the promoter regions of the target genes -> get Total Unique Bins at TSS = 289

Run Generate_Promoter_BED_Bins_for_HiC_APA_Analysis_hg19.R

Get Promoter_Bins_5Kb_hg19.bed 

# Step 2. Convert bed file into bedpe file

awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $1, $2, $3}’
 /share/lab_teng/trainee/JingjingWu/EBV/HiC_icenorm_AggregatePeakAnalysis/Promoter_Bins_5Kb_hg19.bed > /share/lab_teng/trainee/JingjingWu/EBV/HiC_icenorm_AggregatePeakAnalysis/Promoter_Anchors.bedpe

# Step 3. Run Analysis from Juicer

juicer_tools.1.6.2_linux_jcuda.0.8.jar

Individual_Promoter_Results_LCL_HiC_icenorm_AggregatePeakAnalysis_1.28.2026.sh

Core steps:

JAR="/share/lab_teng/trainee/JingjingWu/EBV/HiC_icenorm_AggregatePeakAnalysis/juicer_tools.1.6.2_linux_jcuda.0.8.jar"
HIC="/share/lab_teng/trainee/JingjingWu/EBV/hicpro2juice_out/LCL_5000_iced.hic"
BEDPE="/share/lab_teng/trainee/JingjingWu/EBV/HiC_icenorm_AggregatePeakAnalysis/hiccups_style.bedpe"
OUT_DIR="/share/lab_teng/trainee/JingjingWu/EBV/APA_Results/Individual_Promoter_Results_1.28.2026"

    java -Xmx8g -jar  apa \
        -n 0 \
        -u \
        -r 5000 \
        -k NONE \
        $HIC temp.bedpe "$OUT_DIR/$GENE_ID"

# Step 4. Collect data in tables

Generate_APA_Results_in_Tables.py
