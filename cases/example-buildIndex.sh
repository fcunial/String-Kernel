#!/bin/bash


# IO
INPUT_FILE="./data/bwt/problem0.bwt"
SHARP_POSITION="61933857"
OUTPUT_FILE="./data/problem0.indexed"

echo $PWD
./buildIndex ${INPUT_FILE} ${SHARP_POSITION} ${OUTPUT_FILE}
