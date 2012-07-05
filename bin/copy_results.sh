#!/bin/bash

SAMPLE=$1
DESTINATION=$2

if [ -z "$SAMPLE" ]
then 
    echo "No sample given!"
    exit 1;
fi

if [ -z "$DESTINATION" ]
then 
    echo "No destination given!"
    exit 1;
fi

if [ ! -x "$DESTINATION" ]
then
    mkdir $DESTINATION
fi

cp *${SAMPLE}_1*vcf $DESTINATION/
cp *${SAMPLE}_1*rare*exonic_variant_function $DESTINATION/
cp *${SAMPLE}_1*rare*splicing $DESTINATION/
cp *${SAMPLE}_1*rare*splicing.indisp $DESTINATION/
cp *${SAMPLE}_1*rare*exonic_variant_function.not_syn.indisp $DESTINATION/

