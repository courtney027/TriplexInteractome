#!/bin/bash

rootpath=/share/home/Triplex/CUTTAG_DNMT1
Datapath=/share/home/Triplex/CUTTAG_DNMT1/Raw

cd ${rootpath}
mkdir -p {qc, qc_clean, align, stat, part, rmdup, flt, bw, peak}

for name in D7 D8 D9 D10
do
        fastqc -t 8 -o ${rootpath}/qc/ ${Datapath}/${name}/${name}_1.fq.gz
        fastqc -t 8 -o ${rootpath}/qc/ ${Datapath}/${name}/${name}_2.fq.gz &
done
wait

for name in D7 D8 D9 D10
do
        fastp -i ${Datapath}/${name}/${name}_1.fq.gz -o ${rootpath}/qc_clean/${name}_R1.clean.fq.gz \
        -I ${Datapath}/${name}/${name}_2.fq.gz -O ${rootpath}/qc_clean/${name}_R2.clean.fq.gz &
done
wait

cd ${rootpath}/qc_clean
ls *R1.clean.fq.gz | while read id; do bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant -p 10 -x /share/home/Triplex/bowtie2_index/hg38_spike/hg38_spike -1 ${id} -2 ${id%_R1*}_R2.clean.fq.gz | samtools sort -@ 10 -o - > ${rootpath}/align/${id%_R1*}.sorted.bam; done

cd ${rootpath}/align
ls *.bam | xargs -i samtools index {}
ls *.bam | while read id; do samtools flagstat ${id} > ${rootpath}/stat/$(basename $id ".bam").stat ; done

ls *sorted.bam | while read id ; do samtools view -b -L /share/home/Triplex/Ref_genome/hg38/hg38.bed ${id} > ${rootpath}/part/${id%%.*}.hg38.bam; done
ls *sorted.bam | while read id ; do samtools view -b -L /share/home/Triplex/Ref_genome/hg38/CnT_SpikeIn.bed ${id} > ${rootpath}/part/${id%%.*}.spike.bam; done

cd ${rootpath}/part
ls *.bam | while read id; do echo ${id}; /share/apps/anaconda3/envs/python3/bin/picard -Xmx4g  MarkDuplicates REMOVE_DUPLICATES=true Input=${id} OUTPUT=${id%.*}.rmdup.bam METRICS_FILE=${id%.*}.metrics; done
mv *.rmdup.bam ${rootpath}/rmdup

cd ${rootpath}/rmdup
ls *.rmdup.bam | while read id ; do samtools view -F 4 -f 3 -bh -q 20 ${id} > ${rootpath}/flt/${id%.rmdup*}.flt.bam ; done

cd ${rootpath}/flt
# Samples were also merged by samtools.
ls *.bam | while read id
do
        samtools sort -n ${id} > ${id%.*}.sortn.bam
done

ls *spike.flt.sortn.bam | while read id
do

        num=`samtools view ${id} | wc -l`
        num=`expr ${num} / 2`
        a=1000
        sc=`perl -e "print sprintf('%.4f',$a/$num)"`
        echo "${id}-${sc}" >> scalefactor.txt

        samtools view -H ${id%spike.flt.sortn*}hg38.flt.sortn.bam | grep ^@SQ | cut -d : -f 2- | sed 's/LN://g' | sort -k 1,1 -k 2,2n > ref
        bedtools bamtobed -i ${id%spike.flt.sortn*}hg38.flt.sortn.bam -bedpe | \
        perl -ne 'use List::Util qw(min max); @t=split; print join("\t",$t[0],min($t[1],$t[4]),max($t[2],$t[5]));print "\n"' | sort -k 1,1 -k 2,2n -k 3,3n | \
        bedtools genomecov -i - -g ref -bga -split -scale ${sc} | \
        sort -k 1,1 -k 2,2n > ${id%spike.flt.sortn*}hg38.bg
        ls *.bg | while read x; do LC_COLLATE=C sort -k1,1 -k2,2n ${x} > ${x}.sorted ; done
        bedGraphToBigWig ${id%spike.flt.sortn*}hg38.bg.sorted ref ${id%spike.flt.sortn*}hg38.bw && rm ref
done
mv *.bw ${rootpath}/bw

macs2 callpeak -t D_DMSO_BQQ.hg38.flt.bam -f BAMPE -g hs -n DNMT1_DMSO_BQQ_q001 -q 0.01 --keep-dup=all --outdir ${rootpath}/peak
cat ${rootpath}/peak/DNMT1_DMSO_BQQ_q001_peaks.narrowPeak | bedtools intersect -v -a - -b /share/home/Triplex/Ref_genome/hg38-blacklist.v2.bed > ${rootpath}/peak/DNMT1_DMSO_BQQ_q001_noBlacklist.bed
