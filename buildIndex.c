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
	const uint64_t SHARP_POSITION;
	char *OUTPUT_FILE_PATH;	

	root = json_load_file(argv[1], 0, &error);

	if(!root){
		fprintf(stderr, "error: on line %d: %s\n", error.line, error.text);
		exit(-1);
	}

	json_unpack(root, "{s:s, s:I, s:s}", "OUTPUT_FILE", &OUTPUT_FILE_PATH, "SHARP_POSITION",&SHARP_POSITION, "INPUT_FILE", &INPUT_FILE_PATH);
	
	printf("INPUT_FILE_PATH %s  \n", INPUT_FILE_PATH);
	printf("OUTPUT_FILE_PATH %s  \n", OUTPUT_FILE_PATH);
	printf("SHARP_POSITION %ld  \n", SHARP_POSITION);


	uint64_t nBytes;
	double t, loadingTime, indexingTime, serializationTime;
	Concatenation sequence;
	BwtIndex_t *index;
	
	t=getTime();
	sequence=loadBWT(INPUT_FILE_PATH);
	loadingTime=getTime()-t;	
	

	t=getTime();
	index=buildBwtIndex(sequence.buffer,sequence.length-1,SHARP_POSITION,Basic_bwt_free_text);//-1 because we already add sharp
	indexingTime=getTime()-t;
	
	t=getTime();
	nBytes=serializeBwtIndex(index,OUTPUT_FILE_PATH);
	if (nBytes==0) {
		printf("ERROR in serializeBwtIndex().");
		return 1;
	}
	serializationTime=getTime()-t;

	printf( "%llu|%lf,%lf,%lf|%llu \n",
    		(long long unsigned int)(sequence.length),
			
	        loadingTime,
			indexingTime,
			serializationTime,
			
			(long long unsigned int)malloc_count_peak()
	      );
	freeBwtIndex(index);
	return 0;
}
