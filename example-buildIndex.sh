#/bin/bash
INPUT_FILE="../data/HS22.fasta"
APPEND_RC="0"
OUTPUT_FILE="../data/HS22.fasta.indexed"

./buildIndex ${INPUT_FILE} ${APPEND_RC} ${OUTPUT_FILE}