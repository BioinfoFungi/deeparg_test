# Trimming and QC using Trimmomatic
trimmomatic PE -phred33  \
	F.fq.gz R.fq.gz  \
	./output/trimmomatic/F.fq.gz.paired ./output/trimmomatic/F.fq.gz.unpaired \
	./output/trimmomatic/R.fq.gz.paired ./output/trimmomatic/R.fq.gz.unpaired \
	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Merging paired end reads using Vsearch
vsearch \
    --fastq_mergepairs  ./output/trimmomatic/F.fq.gz.paired \
    --reverse  ./output/trimmomatic/R.fq.gz.paired \
    --fastaout  ./output/vsearch/F.fq.gz.paired.merged \
    --fastaout_notmerged_fwd ./output/vsearch/F.fq.gz.paired.unmerged \
    --fastaout_notmerged_rev ./output/vsearch/R.fq.gz.paired.unmerged



cat  ./output/vsearch/F.fq.gz.paired.merged \
    ./output/vsearch/F.fq.gz.paired.unmerged \
    ./output/vsearch/R.fq.gz.paired.unmerged > ./output/clean/test.clean

# Run DeepARG-SS to identify ARG-like reads
deeparg predict \
    --type nucl \
    --model SS -d /home/wangyang/workspace/gusphdproj-deeparg-ss-fbe063e24cf7/database \
    -i ./output/clean/test.clean \
    -o ./output/deeparg/test.clean.deeparg \
    --arg-alignment-identity 80 \
    --min-prob 0.8 \
    --arg-alignment-evalue 1e-10

# Quantification of ARG-like counts
sort -k1,1 -k2,2n \
    ./output/deeparg/test.clean.deeparg.mapping.ARG  \
    | bedtools merge -c 12,5 -o sum,distinct >./output/deeparg/test.clean.deeparg.mapping.ARG.merged


python merge.py  ./output/deeparg/test.clean.deeparg.mapping.ARG /home/wangyang/workspace/gusphdproj-deeparg-ss-fbe063e24cf7/database


# Normalize to 16S rRNAs - this may take a while
bowtie2 -f \
    --fast-local \
    --no-unal \
    -x /home/wangyang/workspace/gusphdproj-deeparg-ss-fbe063e24cf7/database/data/gg13/dataset \
    -U ./output/clean/test.clean \
    -S ./output/normalize/test.clean.sam

samtools view -bS ./output/normalize/test.clean.sam > ./output/normalize/test.clean.bam
samtools sort ./output/normalize/test.clean.bam -o ./output/normalize/test.clean.sorted.bam

bedtools merge -i ./output/normalize/test.clean.sorted.bam -c 1 -o count > ./output/normalize/test.clean.sorted.bam.merged

python mapping.py output/normalize/test.clean /home/wangyang/workspace/gusphdproj-deeparg-ss-fbe063e24cf7/database/data/gg13/dataset


python normalize.py ./output/normalize/test.clean.sorted.bam.merged.quant ./output/deeparg/test.clean.deeparg.mapping.ARG.merged.quant


