#!/bin/bash

# Inputs
INPUT_VCF="28samples_ref1_snps_pass_bial.vcf.gz"
OUTPUT_ROH="28samples_ref1_snps_pass_bial.roh.txt"
EDITED_OUTPUT="28samples_ref1_snps_pass_bial.roh.edited.txt"

# bcftools roh
bcftools roh --AF-dflt 0.4 -I -G30 --rec-rate 1.4e-9 "$INPUT_VCF" > "$OUTPUT_ROH"

# Extract columns
grep "RG" "$OUTPUT_ROH" | cut -f 2,3,6 > "$EDITED_OUTPUT"

echo "Done"