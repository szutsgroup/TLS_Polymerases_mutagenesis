#!/bin/bash

# generating a GATK Panel-of-Normals from all the samples to obtained common heterozygous loci

# run MuTect2
for i in *bam;
do
java -jar GenomeAnalysisTK.jar -T MuTect2 \
-I:tumor ${i} \
--artifact_detection_mode \
-R reference/GRCh38_EBV.fasta
-o $(basename $i .bam).vcf.gz
done

# combine variants
java -jar GenomeAnalysisTK.jar -T CombineVariants \
-minN 10
$(ls *vcf.gz | awk '{print "-V "$1'}) | tr "\n" " ") \
--setKey "null" \
--filteredAreUncalled \
--filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
-o pon_combinevariants.vcf.gz \
-R reference/GRCh38_EBV.fasta \
--genotypemergeoption UNIQUIFY

# get loci

java -jar picard.jar MakeSitesOnlyVcf \
    I=pon_combinevariants.vcf.gz \
    O=pon_siteonly.vcf.gz

# restrict the position list to have at least 10kb intermutational distance
gzip -cd pon_siteonly.vcf.gz | awk 'BEGIN{a=0; b=0}{a = $2-b; if (a > 10000) {print $1"\t"$2"\t"$4"\t"$5"\tOK"; b = $2; a = 0}} > loci_list.txt

