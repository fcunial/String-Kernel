#!/bin/bash

# IO
INPUT_FILE="./data/fna/problem0.fna"
INPUT_FORMAT="1"
APPEND_RC="0"
OUTPUT_FILE="./data/bwt/problem0.bwt"

echo $PWD
./buildBWT ${INPUT_FILE} ${INPUT_FORMAT} ${APPEND_RC} ${OUTPUT_FILE}
