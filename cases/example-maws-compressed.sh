#!/bin/bash

# IO
INPUT_FILE="./data/problem0.indexed"
WRITE_MAWS="1"
OUTPUT_FILE="./data/problem0.maws.compressed"

# Number of threads
N_THREADS="1"

# Statistics
MIN_HISTOGRAM_LENGTH="0"
MAX_HISTOGRAM_LENGTH="0"

# Filtering by length
MIN_LENGTH="2"
MAX_LENGTH="100000"

# Filtering by score
SELECT_SCORE="0"
SCORE_THRESHOLD="0"
COMPUTE_SCORES="0"

# Compression
COMPRESS_OUTPUT="1"

# Number of workpakages
N_WORKPACKAGES_RATE="2"


./run_MAWs_single ${INPUT_FILE} ${N_THREADS} \
	${MIN_LENGTH} ${MAX_LENGTH} \
	${MIN_HISTOGRAM_LENGTH} ${MAX_HISTOGRAM_LENGTH} \
	${COMPUTE_SCORES} ${SELECT_SCORE} ${SCORE_THRESHOLD} \
	${WRITE_MAWS} ${OUTPUT_FILE} \
	${COMPRESS_OUTPUT} ${N_WORKPACKAGES_RATE}
