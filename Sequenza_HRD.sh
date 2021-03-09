#! /usr/bin/bash

# ----------------------------------- Description ------------------------------------- #
# Perform:                                                                              #
#   Trimming (fastp) +                                                                  #
#   Alignment (Sentieon-bwa) +                                                          #
#   statistic analysis +                                                                #
#   Deduplication (optional) +                                                          #
#   Allelic specific CNV calling (Sequencza) +                                          #
#   Calculate and plot HRD genomic scars (HRDecipher)                                   #
# ------------------------------------------------------------------------------------- #
# Usage:                                                                                #
#   [admin@kai]$bash Sequenza_HRD.sh [input_folder] [output_folder] [BED]               #
#                                                                                       #
# input_folder should contain pairs of matched tumor/normal pair-end sequenced Fastqs   #
#                                                                                       #
# Naming convention:                                                                    #
# tumor:   ${sampleID}_tumor_R[1|2].fastq.gz                                            #
# normal:  ${sampleID}_normal_R[1|2].fastq.gz                                           #
#                                                                                       #
# 1. Depends on DNA capture methods (i.e. targeted amplicon based or hybrid capture     #
#    based), user can choose to either perform or not perform deduplication step.       #
# 2. As this assay designed to sequence on backbone region and gene region at different #
#    depths, two BED files must be prepared to indicate those two regions separately.   #
# ------------------------------------------------------------------------------------- #


# --------------------------- set parameters --------------------------- #
# ---------------------------------------------------------------------- #
sentieon_license="192.168.1.186:8990"
thread=8

# whether perform deduplicate step (true || false)
dedup=true

# path to software 
fastp="/data/ngs/softs/fastp/fastp"
sentieon="/data/ngs/softs/sentieon/sentieon-genomics-201808.08/bin/sentieon"
bcftools="/public/software/bcftools-1.9/bcftools"
samtools="/public/software/samtools-1.9/samtools"
bamdst="/public/software/bamdst/bamdst"
HRDecipher="/public/home/kai/softwares/HRDecipher/HRDecipher.py"

ref="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"
gcwig="/public/database/GATK_Resource_Bundle/hg19/hg19.gc50.wig.gz"

# switch (on||off)
do_trim="on"
do_align="on"
do_dedup="on"
do_qc="on"
do_seqz="on"
do_cnv="on"
do_hrd="on"

# ------------------------------ argparser ----------------------------- #
# ---------------------------------------------------------------------- #
if [[  $1 == '-h'  ]]; then
    echo "Usage: ./Sequenza_HRD.sh [input_folder] [output_folder] [BED]"
    echo "-------------------------------------------------------------------------"
    echo "[input_folder] should contain fastq files with following naming system:"
    echo "\${sampleID}_[tumor|normal]_R[1|2].fastq.gz"
    exit 0
fi

input_folder=`realpath $1`
output_folder=`realpath $2`
bed=`realpath $3`

if [[ ! -d $input_folder  ]];
then
    echo "Error: input_folder does not Found!"
    exit 1
fi

if [[ ! -f $bed  ]];
then
    echo "Error: BED file does not Found!"
    exit 1
fi

# ----------------------  orgnise output dir  -------------------------- #
# ---------------------------------------------------------------------- #
if [[ ! -d $output_folder ]]; 
then
    mkdir $output_folder
fi

trim_dir=$output_folder/trim/;
if [[ ! -d $trim_dir ]]; then
    mkdir $trim_dir
fi

align_dir=$output_folder/align/;
if [[ ! -d $align_dir ]]; then
    mkdir $align_dir
fi

cnv_dir=$output_folder/sequenza_cnv/;
if [[ ! -d $cnv_dir ]]; then
    mkdir $cnv_dir
fi

hrd_dir=${output_folder}/HRD/;
if [[  ! -d $hrd_dir  ]];
then 
    mkdir $hrd_dir
fi

qc_dir=$output_folder/qc/;
if [[ ! -d $qc_dir ]]; then
    mkdir $qc_dir
fi


# ---------------------------  LOGGING  -------------------------------- #
# ---------------------------------------------------------------------- #
echo "LOGGING: `date --rfc-3339=seconds` -- Analysis started"
echo "LOGGING: This is the Sequenza_HRD.sh pipeline"
echo "========================================================"
echo "LOGGING: -- settings -- input folder -- ${input_folder}"
echo "LOGGING: -- settings -- output folder -- ${output_folder}"
echo "LOGGING: -- settings -- BED file -- ${bed}"
echo "========================================================"

echo "sampleID,fastq_size,raw_reads,raw_bases,clean_reads,clean_bases,\
qc30_rate,mapping_rate(%),on-target_percent(%),\
mean_depth,mean_dedup_depth,dup_rate(%),\
average_insert_size,std_insert_size,\
Uniformity_0.1X(%),Uniformity_0.2X(%),\
Uniformity_0.5X(%),Uniformity_1X(%),\
50x_depth_percent(%),100x_depth_percent(%),\
150x_depth_percent(%),200x_depth_percent(%),\
300x_depth_percent(%),400x_depth_percent(%),\
500x_depth_percent(%)" > $qc_dir/QC_summary.csv



# ---------------------------------------------------------------------- #
# ---------------------------  Pipeline  ------------------------------- #
# ---------------------------------------------------------------------- #
export SENTIEON_LICENSE=${sentieon_license};

for ifile in $input_folder/*_tumor_R1.fastq.gz;
do 
    sampleID=`basename ${ifile%%"_tumor"*}`

    # step1 - trim reads
    if [[  $do_trim == "on"  ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- trimming reads";

        $fastp --in1 $input_folder/${sampleID}_tumor_R1.fastq.gz \
        --in2 $input_folder/${sampleID}_tumor_R2.fastq.gz \
        --out1 $trim_dir/${sampleID}_trim_R1.tumor.fastq.gz \
        --out2 $trim_dir/${sampleID}_trim_R2.tumor.fastq.gz \
        -c --length_required 3 --detect_adapter_for_pe -p \
        --thread ${thread} \
        --html $trim_dir/${sampleID}.tumor.trim.html \
        --json $trim_dir/${sampleID}.tumor.trim.json;

        $fastp --in1 $input_folder/${sampleID}_normal_R1.fastq.gz \
        --in2 $input_folder/${sampleID}_normal_R2.fastq.gz \
        --out1 $trim_dir/${sampleID}_trim_R1.normal.fastq.gz \
        --out2 $trim_dir/${sampleID}_trim_R2.normal.fastq.gz \
        -c --length_required 3 --detect_adapter_for_pe -p \
        --thread ${thread} \
        --html $trim_dir/${sampleID}.normal.trim.html \
        --json $trim_dir/${sampleID}.normal.trim.json;
    fi

    # step2 - align & sort
    if [[  $do_align == "on"  ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- alignment & sorting";

        ($sentieon bwa mem -M -R "@RG\tID:${sampleID}.tumor\tSM:${sampleID}.tumor\tPL:illumina" \
        -t ${thread} -K 10000000 ${ref} \
        $trim_dir/${sampleID}_trim_R1.tumor.fastq.gz \
        $trim_dir/${sampleID}_trim_R2.tumor.fastq.gz \
        || echo -n 'error' ) \
        | ${sentieon} util sort -r ${ref} -o ${align_dir}/${sampleID}.tumor.sorted.bam \
        -t ${thread} --sam2bam -i -;

        ($sentieon bwa mem -M -R "@RG\tID:${sampleID}.normal\tSM:${sampleID}.normal\tPL:illumina" \
        -t ${thread} -K 10000000 ${ref} \
        $trim_dir/${sampleID}_trim_R1.normal.fastq.gz \
        $trim_dir/${sampleID}_trim_R2.normal.fastq.gz \
        || echo -n 'error' ) \
        | ${sentieon} util sort -r ${ref} -o ${align_dir}/${sampleID}.normal.sorted.bam \
        -t ${thread} --sam2bam -i -;
    fi

    # step3 (optional) - remove duplicates
    if [[  $do_dedup == "on"  ]];
    then
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- deduplication";
        
        if [[ $dedup == true ]]; then
            ${sentieon} driver -t ${thread} \
            -i ${align_dir}/${sampleID}.tumor.sorted.bam \
            --algo LocusCollector \
            --fun score_info ${align_dir}/${sampleID}.tumor.score.txt;

            ${sentieon} driver -t ${thread} \
            -i ${align_dir}/${sampleID}.tumor.sorted.bam \
            --algo Dedup \
            --score_info ${align_dir}/${sampleID}.tumor.score.txt \
            --metrics ${align_dir}/${sampleID}.tumor.dedup_metrics.txt \
            ${align_dir}/${sampleID}.tumor.sorted.dedup.bam;

            ${sentieon} driver -t ${thread} \
            -i ${align_dir}/${sampleID}.normal.sorted.bam \
            --algo LocusCollector \
            --fun score_info ${align_dir}/${sampleID}.normal.score.txt;

            ${sentieon} driver -t ${thread} \
            -i ${align_dir}/${sampleID}.normal.sorted.bam \
            --algo Dedup \
            --score_info ${align_dir}/${sampleID}.normal.score.txt \
            --metrics ${align_dir}/${sampleID}.normal.dedup_metrics.txt \
            ${align_dir}/${sampleID}.normal.sorted.dedup.bam;
        fi
    fi


    # step4 - quality control
    if [[  $do_qc == "on"  ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- qc";

        if [[ ! -d $qc_dir/${sampleID}.normal ]]; then
            mkdir $qc_dir/${sampleID}.normal
        fi

        if [[ ! -d $qc_dir/${sampleID}.tumor ]]; then
            mkdir $qc_dir/${sampleID}.tumor
        fi

        $bamdst -p $bed -o $qc_dir/${sampleID}.normal \
        ${align_dir}/${sampleID}.normal.sorted.dedup.bam;

        $bamdst -p $bed -o $qc_dir/${sampleID}.tumor \
        ${align_dir}/${sampleID}.tumor.sorted.dedup.bam;

        $samtools stats -@ ${thread} ${align_dir}/${sampleID}.normal.sorted.dedup.bam > ${qc_dir}/${sampleID}.normal.stats.txt;
        $samtools stats -@ ${thread} ${align_dir}/${sampleID}.tumor.sorted.dedup.bam > ${qc_dir}/${sampleID}.tumor.stats.txt;

        tumor_r1=$(du $input_folder/${sampleID}_tumor_R1.fastq.gz -shL |awk '{print $1}');
        tumor_r2=$(du $input_folder/${sampleID}_tumor_R2.fastq.gz -shL |awk '{print $1}');
        normal_r1=$(du $input_folder/${sampleID}_normal_R1.fastq.gz -shL |awk '{print $1}');
        normal_r2=$(du $input_folder/${sampleID}_normal_R2.fastq.gz -shL |awk '{print $1}');

        normal_raw_reads=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.normal.trim.json', 'r')); \
        print(fh['summary']['before_filtering']['total_reads'])"`

        normal_clean_reads=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.normal.trim.json', 'r')); \
        print(fh['summary']['after_filtering']['total_reads'])"`

        normal_raw_bases=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.normal.trim.json', 'r')); \
        print(fh['summary']['before_filtering']['total_bases'])"`

        normal_clean_bases=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.normal.trim.json', 'r')); \
        print(fh['summary']['after_filtering']['total_bases'])"`

        tumor_raw_reads=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.tumor.trim.json', 'r')); \
        print(fh['summary']['before_filtering']['total_reads'])"`

        tumor_clean_reads=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.tumor.trim.json', 'r')); \
        print(fh['summary']['after_filtering']['total_reads'])"`

        tumor_raw_bases=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.tumor.trim.json', 'r')); \
        print(fh['summary']['before_filtering']['total_bases'])"`

        tumor_clean_bases=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.tumor.trim.json', 'r')); \
        print(fh['summary']['after_filtering']['total_bases'])"`

        normal_qc_rate=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.normal.trim.json', 'r')); \
        print(fh['summary']['before_filtering']['q30_rate'])"`

        tumor_qc_rate=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.tumor.trim.json', 'r')); \
        print(fh['summary']['before_filtering']['q30_rate'])"`

        normal_mapping_rate=$(grep "Fraction of Mapped Reads" $qc_dir/${sampleID}.normal/coverage.report | awk -F"\t" '{print $2}');
        tumor_mapping_rate=$(grep "Fraction of Mapped Reads" $qc_dir/${sampleID}.tumor/coverage.report | awk -F"\t" '{print $2}');

        normal_mean_depth=$(grep "Average depth" $qc_dir/${sampleID}.normal/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
        normal_mean_dedup_depth=$(grep "Average depth(rmdup)" $qc_dir/${sampleID}.normal/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
        normal_dup_rate=$(grep "Fraction of PCR duplicate reads" $qc_dir/${sampleID}.normal/coverage.report |awk -F"\t" '{print $2}');
        
        tumor_mean_depth=$(grep "Average depth" $qc_dir/${sampleID}.tumor/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
        tumor_mean_dedup_depth=$(grep "Average depth(rmdup)" $qc_dir/${sampleID}.tumor/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
        tumor_dup_rate=$(grep "Fraction of PCR duplicate reads" $qc_dir/${sampleID}.tumor/coverage.report |awk -F"\t" '{print $2}');

        normal_on_target=$(grep "Fraction of Target Reads in all reads" $qc_dir/${sampleID}.normal/coverage.report |awk -F"\t" '{print $2}');
        tumor_on_target=$(grep "Fraction of Target Reads in all reads" $qc_dir/${sampleID}.tumor/coverage.report |awk -F"\t" '{print $2}');

        normal_insert_size=$(awk -F"\t" '$2 == "insert size average:" {print $3}' ${qc_dir}/${sampleID}.normal.stats.txt);
        normal_insert_std=$(awk -F"\t" '$2 == "insert size standard deviation:" {print $3}' ${qc_dir}/${sampleID}.normal.stats.txt);

        tumor_insert_size=$(awk -F"\t" '$2 == "insert size average:" {print $3}' ${qc_dir}/${sampleID}.tumor.stats.txt);
        tumor_insert_std=$(awk -F"\t" '$2 == "insert size standard deviation:" {print $3}' ${qc_dir}/${sampleID}.tumor.stats.txt);
        
        normal_50x=$(less -S $qc_dir/${sampleID}.normal/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 50) count+=1} END {print count/NR*100}');
        normal_100x=$(less -S $qc_dir/${sampleID}.normal/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 100) count+=1} END {print count/NR*100}');
        normal_150x=$(less -S $qc_dir/${sampleID}.normal/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 150) count+=1} END {print count/NR*100}');
        normal_200x=$(less -S $qc_dir/${sampleID}.normal/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 200) count+=1} END {print count/NR*100}');
        normal_300x=$(less -S $qc_dir/${sampleID}.normal/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 300) count+=1} END {print count/NR*100}');
        normal_400x=$(less -S $qc_dir/${sampleID}.normal/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 400) count+=1} END {print count/NR*100}');
        normal_500x=$(less -S $qc_dir/${sampleID}.normal/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 500) count+=1} END {print count/NR*100}');

        tumor_50x=$(less -S $qc_dir/${sampleID}.tumor/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 50) count+=1} END {print count/NR*100}');
        tumor_100x=$(less -S $qc_dir/${sampleID}.tumor/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 100) count+=1} END {print count/NR*100}');
        tumor_150x=$(less -S $qc_dir/${sampleID}.tumor/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 150) count+=1} END {print count/NR*100}');
        tumor_200x=$(less -S $qc_dir/${sampleID}.tumor/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 200) count+=1} END {print count/NR*100}');
        tumor_300x=$(less -S $qc_dir/${sampleID}.tumor/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 300) count+=1} END {print count/NR*100}');
        tumor_400x=$(less -S $qc_dir/${sampleID}.tumor/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 400) count+=1} END {print count/NR*100}');
        tumor_500x=$(less -S $qc_dir/${sampleID}.tumor/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 500) count+=1} END {print count/NR*100}');

        normal_01x=$(less -S $qc_dir/${sampleID}.normal/depth.tsv.gz | awk -v depth=${normal_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.1) count+=1} END {print count/NR*100}');
        normal_02x=$(less -S $qc_dir/${sampleID}.normal/depth.tsv.gz | awk -v depth=${normal_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.2) count+=1} END {print count/NR*100}');
        normal_05x=$(less -S $qc_dir/${sampleID}.normal/depth.tsv.gz | awk -v depth=${normal_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.5) count+=1} END {print count/NR*100}');
        normal_1x=$(less -S $qc_dir/${sampleID}.normal/depth.tsv.gz | awk -v depth=${normal_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth) count+=1} END {print count/NR*100}');

        tumor_01x=$(less -S $qc_dir/${sampleID}.tumor/depth.tsv.gz | awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.1) count+=1} END {print count/NR*100}');
        tumor_02x=$(less -S $qc_dir/${sampleID}.tumor/depth.tsv.gz | awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.2) count+=1} END {print count/NR*100}');
        tumor_05x=$(less -S $qc_dir/${sampleID}.tumor/depth.tsv.gz | awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.5) count+=1} END {print count/NR*100}');
        tumor_1x=$(less -S $qc_dir/${sampleID}.tumor/depth.tsv.gz | awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth) count+=1} END {print count/NR*100}');

        echo "${sampleID}.normal,${normal_r1}/${normal_r2},${normal_raw_reads},${normal_raw_bases},${normal_clean_reads},${normal_clean_bases},\
        ${normal_qc_rate},${normal_mapping_rate},${normal_on_target},${normal_mean_depth},${normal_mean_dedup_depth},${normal_dup_rate},\
        ${normal_insert_size},${normal_insert_std},${normal_01x},${normal_02x},${normal_05x},${normal_1x},\
        ${normal_50x},${normal_100x},${normal_150x},${normal_200x},${normal_300x},${normal_400x},${normal_500x}" \
        >> ${qc_dir}/QC_summary.csv

        echo "${sampleID}.tumor,${tumor_r1}/${tumor_r2},${tumor_raw_reads},${tumor_raw_bases},${tumor_clean_reads},${tumor_clean_bases},\
        ${tumor_qc_rate},${tumor_mapping_rate},${tumor_on_target},${tumor_mean_depth},${tumor_mean_dedup_depth},${tumor_dup_rate},\
        ${tumor_insert_size},${tumor_insert_std},${tumor_01x},${tumor_02x},${tumor_05x},${tumor_1x},\
        ${tumor_50x},${tumor_100x},${tumor_150x},${tumor_200x},${tumor_300x},${tumor_400x},${tumor_500x}" \
        >> ${qc_dir}/QC_summary.csv
    fi

    # step5 - generate seqz file
    if [[  $do_seqz == "on"  ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- generating seqz file";

        if [[ $dedup == true ]]; then
            normal_bam=${align_dir}/${sampleID}.normal.sorted.dedup.bam
            tumor_bam=${align_dir}/${sampleID}.tumor.sorted.dedup.bam
        else
            normal_bam=${align_dir}/${sampleID}.normal.sorted.bam
            tumor_bam=${align_dir}/${sampleID}.tumor.sorted.bam
        fi

        sequenza-utils bam2seqz \
        --chromosome chr1 chr2 chr3 chr4 chr5 chr6 chr7 \
        chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 \
        chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
        --parallel ${thread} \
        -n $normal_bam \
        -t $tumor_bam \
        --fasta ${ref} -gc ${gcwig} -f illumina \
        -o $cnv_dir/${sampleID}.seqz.gz;

        ls $cnv_dir/${sampleID}*.seqz.gz | \
        xargs -i -P ${thread} sh -c \
        "sequenza-utils seqz_binning \
        -s {} -w 50 -o {}.bin.gz"

        zcat $cnv_dir/${sampleID}_chr{1..22}.seqz.gz.bin.gz \
        $cnv_dir/${sampleID}_chrX.seqz.gz.bin.gz \
        $cnv_dir/${sampleID}_chrY.seqz.gz.bin.gz \
        > $cnv_dir/${sampleID}.bin.seqz;

        grep -v "chromosome" $cnv_dir/${sampleID}.bin.seqz > $cnv_dir/tmp.contents;
        head -n 1 $cnv_dir/${sampleID}.bin.seqz > $cnv_dir/tmp.header;
        cat $cnv_dir/tmp.header $cnv_dir/tmp.contents > $cnv_dir/${sampleID}.bin.seqz;
        rm $cnv_dir/tmp.header $cnv_dir/tmp.contents;

        bgzip $cnv_dir/${sampleID}.bin.seqz;
        tabix -f -s 1 -b 2 -e 2 -S 1 $cnv_dir/${sampleID}.bin.seqz.gz;

        rm $cnv_dir/${sampleID}_chr*.gz;
    fi

    # step6 - process seqz data, normalization, segmentation and estimate cellularity and ploidy
    if [[  $do_cnv == "on"  ]];
    then 
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- CNV calling ";

        Rscript -e "library(sequenza); \
        extract <- sequenza.extract('$cnv_dir/${sampleID}.bin.seqz.gz'); \
        cp <- sequenza.fit(extract); \
        sequenza.results(sequenza.extract = extract, \
                        cp.table = cp, \
                        sample.id = '${sampleID}', \
                        out.dir = '$cnv_dir/${sampleID}')"
    fi

    # step7 - HRDecipher
    if [[  $do_hrd == "on" ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- run HRDecipher ";

        ploidy=$(awk 'NR ==2 {print $2}' $cnv_dir/${sampleID}/${sampleID}_alternative_solutions.txt);

        awk -F"\t" -v ploidy=${ploidy} -v sample=${sampleID} \
        'BEGIN {OFS="\t"; print "SampleID","Chromosome","Start_position","End_position","total_cn","A_cn","B_cn","ploidy"};
        NR!=1 {OFS="\t"; print sample,$1,$2,$3,$10,$11,$12,ploidy}' \
        $cnv_dir/${sampleID}/${sampleID}_segments.txt > $hrd_dir/${sampleID}.pre_hrd.tsv;

        python3 $HRDecipher $hrd_dir/${sampleID}.pre_hrd.tsv $cnv_dir/${sampleID}.bin.seqz.gz;
    fi
done


if [[  $do_hrd == "on" ]];
then
    cat $hrd_dir/*.hrd.tsv | awk 'NR==1 || !/(HRD-sum)/' > $hrd_dir/HRD_results.tsv
fi


echo "LOGGING: `date --rfc-3339=seconds` -- Analysis finished";