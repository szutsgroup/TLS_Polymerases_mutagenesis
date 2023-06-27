#!/bin/bash

for i in *bam;
do
python3 getaf.py loci_list.txt $i $(basename $i .bam)_getaf.txt
done