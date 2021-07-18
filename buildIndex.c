/**
 * @author Fabio Cunial
 */
#include "./malloc_count/malloc_count.h"  // For measuring memory usage
#include "./iterator/DNA5_Basic_BWT.h"
#include "./io/io.h"
#include "./io/bufferedFileWriter.h"


/**
 * 1: input file path;
 * 2: position of $;
 * 3: output file path. If the file already exists, its content is overwritten.
 */
int main(int argc, char **argv) {
	char *INPUT_FILE_PATH = argv[1];
	const uint64_t SHARP_POSITION = atoi(argv[2]);
	char *OUTPUT_FILE_PATH = argv[3];

	printf("INPUT_FILE_PATH %s  \n", INPUT_FILE_PATH);
	printf("SHARP_POSITION %li  \n", SHARP_POSITION);
	printf("OUTPUT_FILE_PATH %s  \n", OUTPUT_FILE_PATH);

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
