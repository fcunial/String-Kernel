//#include "simple_parse.c"
#include <jansson.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    json_t *json;
    json_error_t error;
    if (argc != 2) {
        fprintf(stderr, "Usage: %s config.json\n", argv[0]);
        exit(-1);
    }
    json_t *root;
    root = json_load_file(argv[1], 0, &error);

    /* root is the JSON object {"foo": "bar", "quux": true} */
    //const char *str;
    //int boolean;

    char *INPUT_FILE_PATH;
    const uint8_t WRITE_MAWS;
    char *OUTPUT_FILE_PATH;
    const uint8_t N_THREADS;
    const uint64_t MIN_MAW_LENGTH;
    const uint64_t MAX_MAW_LENGTH;
    const uint64_t MIN_HISTOGRAM_LENGTH;
    const uint64_t MAX_HISTOGRAM_LENGTH;
    const uint8_t COMPUTE_SCORES;
    unsigned char SELECTED_SCORE;
    double SELECTED_SCORE_THRESHOLD;
    uint8_t COMPRESS_OUTPUT;
    
    //json_unpack(root, "{s:s, s:b}", "foo", &str, "quux", &boolean);


    json_unpack(root, "{s:s, s:i, s:s, s:i}",//, s:b, s:b, s:b, s:b, s:b, s:b, s:b, s:b
	 "INPUT_FILE", &INPUT_FILE_PATH, "WRITE_MAWS", &WRITE_MAWS, "OUTPUT_FILE", &OUTPUT_FILE_PATH, "N_THREADS", &N_THREADS);
    printf("%s %i %s %i",INPUT_FILE_PATH, WRITE_MAWS, OUTPUT_FILE_PATH, N_THREADS);

    /*
    	"INPUT_FILE": "./data/HS22.fasta.indexed",
    	"WRITE_MAWS": 1,
    	"OUTPUT_FILE": "./data/HS22.fasta.maws.basic",
	"N_THREADS": 1,
	"MIN_LENGTH": 2,
	"MAX_LENGTH": 100000,
	"MIN_HISTOGRAM_LENGTH": 0,
	"MAX_HISTOGRAM_LENGTH": 0,
	"COMPUTE_SCORES": 0,
	"SELECT_SCORE": 0,
	"SCORE_THRESHOLD": 0,
	"COMPRESS_OUTPUT": 0
    */

    return 0;
}
