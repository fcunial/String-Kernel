#include <stdio.h>
#include <stdlib.h>
#include "./malloc_count/malloc_count.h"  // For measuring memory usage
#include "./io/io.h"
#include "./iterator/DNA5_Basic_BWT.h"
#include "./callbacks/MRWs_single.h"


/** 
 * 1: input file path;
 * 2: append reverse-complement (1/0);
 * 3: minimum MRW length;
 * 4: minFreq;
 * 5: maxFreq;
 * 6: write MRWs to a file (1/0);
 * 7: output file path (read only if the previous argument is 1).
 */
int main(int argc, char **argv) {
	char *INPUT_FILE_PATH = argv[1];
	const unsigned char APPEND_RC = atoi(argv[2]);
	const unsigned long MIN_MRW_LENGTH = atoi(argv[3]);
	const unsigned long MIN_FREQ = atoi(argv[4]);
	const unsigned long MAX_FREQ = atoi(argv[5]);
	const unsigned char WRITE_MRWS = atoi(argv[6]);
	char *OUTPUT_FILE_PATH = NULL;
	if (WRITE_MRWS==1) OUTPUT_FILE_PATH=argv[7];
	double t, tPrime, loadingTime, indexingTime, processingTime;
	unsigned long nMRWs;
	Concatenation sequence;
	Basic_BWT_t *bbwt;
	FILE *file;

	t=getTime();
	sequence=loadFASTA(INPUT_FILE_PATH,APPEND_RC);
	tPrime=getTime();
	loadingTime=tPrime-t;
	
	t=tPrime;
	bbwt=Build_BWT_index_from_text(sequence.buffer,sequence.length,Basic_bwt_free_text);
	indexingTime=getTime()-t;
	
	// Erasing output file
	if (WRITE_MRWS==1) {
		file=fopen(OUTPUT_FILE_PATH,"w");
		fclose(file);
	}
	
	t=getTime();
	nMRWs=find_MRWs_single(bbwt,MIN_MRW_LENGTH,MIN_FREQ,MAX_FREQ,WRITE_MRWS,OUTPUT_FILE_PATH);
	processingTime=getTime()-t;
	free_Basic_BWT(bbwt);
	
	if (WRITE_MRWS==1) fclose(file);
	
	printf( "%ld,%ld,%d,%ld,%ld,%ld,%f,%f,%f,%ld,%ld \n", 
	        sequence.inputLength,
	        sequence.length,
			sequence.hasRC,
			MIN_MRW_LENGTH,
			MIN_FREQ,
			MAX_FREQ,
			loadingTime,
			indexingTime,
			processingTime,
			malloc_count_peak(),
			nMRWs
	      );
	return 0;
}