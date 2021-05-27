/**
 * @author Fabio Cunial
 */
#include "./malloc_count/malloc_count.h"  // For measuring memory usage
#include "./iterator/DNA5_Basic_BWT.h"
#include "./io/io.h"
#include "./io/bufferedFileWriter.h"
#include <jansson.h>

int main(int argc, char **argv) {

	if (argc != 2) {
        	fprintf(stderr, "Usage: %s config.json\n", argv[0]);
        	exit(-1);
    	}

    	json_error_t error;
	json_t *root;

	char *INPUT_FILE_PATH;
	const uint8_t IS_FASTA;
	char *OUTPUT_FILE_PATH;
	const uint8_t APPEND_RC;	

	root = json_load_file(argv[1], 0, &error);

	if(!root){
		fprintf(stderr, "error: on line %d: %s\n", error.line, error.text);
		return 1;
	}

	json_unpack(root, "{s:s, s:I, s:I, s:s}", "OUTPUT_FILE", &OUTPUT_FILE_PATH, "INPUT_FORMAT", &IS_FASTA, "APPEND_RC", &APPEND_RC, "INPUT_FILE", &INPUT_FILE_PATH);
	
	printf("INPUT_FILE_PATH %s  \n", INPUT_FILE_PATH);
	printf("IS_FASTA %i  \n", IS_FASTA);
	printf("OUTPUT_FILE_PATH %s  \n", OUTPUT_FILE_PATH);
	printf("APPEND_RC %i  \n", APPEND_RC);


	uint64_t nBytes;
	double t, loadingTime, indexingTime, serializationTime;
	Concatenation sequence;
	BwtIndex_t *index;
	
	t=getTime();
	sequence=IS_FASTA?loadFASTA(INPUT_FILE_PATH,APPEND_RC):loadPlainText(INPUT_FILE_PATH,APPEND_RC);
	loadingTime=getTime()-t;	
	t=getTime();
	index=buildBwtIndex(sequence.buffer,sequence.length,Basic_bwt_free_text);
	indexingTime=getTime()-t;
	t=getTime();
	nBytes=serializeBwtIndex(index,OUTPUT_FILE_PATH);
	if (nBytes==0) {
		printf("ERROR in serializeBwtIndex().");
		return 1;
	}
	serializationTime=getTime()-t;
	printf( "%llu,%llu,%d|%lf,%lf,%lf|%llu \n",
    		(long long unsigned int)(sequence.inputLength),
    		(long long unsigned int)(sequence.length),
			sequence.hasRC,
			
	        loadingTime,
			indexingTime,
			serializationTime,
			
			(long long unsigned int)malloc_count_peak()
	      );
	freeBwtIndex(index);
	return 0;
}
