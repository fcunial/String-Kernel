/**
 * @author Saad Alowayyed
 */
#include "./malloc_count/malloc_count.h"  // For measuring memory usage
#include "./BWT/BwtBuilders.h"
#include "./io/io.h"
#include "./io/bufferedFileWriter.h"


/**
 * 1: input file path;
 * 2: input format;
 * 3: append reverse-complement
 * 4: output file path. If the file already exists, its content is overwritten.
 */

int main(int argc, char **argv) {

	char *INPUT_FILE_PATH = argv[1];
	const uint8_t IS_FASTA = atoi(argv[2]);
	const uint8_t APPEND_RC = atoi(argv[3]);
	char *OUTPUT_FILE_PATH = argv[4];
	
	printf("INPUT_FILE_PATH %s  \n", INPUT_FILE_PATH);
	printf("IS_FASTA %i  \n", IS_FASTA);
	printf("OUTPUT_FILE_PATH %s  \n", OUTPUT_FILE_PATH);
	printf("APPEND_RC %i  \n\n", APPEND_RC);

	
	double t, loadingTime, indexingTime;
	Concatenation sequence;
	uint64_t sharpPosition;
	uint8_t *bwt;

	t=getTime();
	printf("Loading  %s...\n",INPUT_FILE_PATH);
	sequence=IS_FASTA?loadFASTA(INPUT_FILE_PATH,APPEND_RC):loadPlainText(INPUT_FILE_PATH,APPEND_RC);
	loadingTime=getTime()-t;

	t=getTime();
	printf("Building BWT ...\n");
	bwt=useDivsufsort(sequence.buffer,sequence.length, &sharpPosition);
	indexingTime=getTime()-t;

	printf("Finish building BWT, the parameters are \n");
	printf("sharpPosition %ld\n",sharpPosition);

	printf("Writing BWT to %s ...\n",OUTPUT_FILE_PATH); 
	///*
	FILE *fp;
	fp = fopen(OUTPUT_FILE_PATH, "w" );
	fwrite(bwt , sequence.length+1 , sizeof(uint8_t) , fp );// length+1 is the RIGHT one
	fclose(fp);
	//*/
	printf("Done ... \n");
	
	
	printf( "%llu,%llu,%d|%lf,%lf|%llu \n",
    		(long long unsigned int)(sequence.inputLength),
    		(long long unsigned int)(sequence.length),
			sequence.hasRC,
			
	        	loadingTime,
			indexingTime,
			
			(long long unsigned int)malloc_count_peak()
	      );
	free(bwt);
	return 0;
}
