#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "./malloc_count/malloc_count.h"  // For measuring memory usage
#include "./io/io.h"
#include "./iterator/DNA5_Basic_BWT.h"
#include "./callbacks/MAWs_single.h"


/** 
 * 0: input file path;
 * 1: append reverse-complement (1/0);
 * 2: minimum MAW length;
 * 3: write MAWs to a file (1/0).
 */
int main(int argc, char **argv) {
	char *INPUT_FILE_PATH = argv[0];
	const unsigned char APPEND_RC = atoi(argv[1]);
	const unsigned char MIN_MAW_LENGTH = atoi(argv[2]);
	const unsigned char WRITE_MAWS = atoi(argv[3]);
	char *OUTPUT_FILE_PATH = argv[4];
	double t, tPrime, loadingTime, indexingTime, processingTime;
	unsigned long long nMAWs;
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
	if (WRITE_MAWS==1) {
		file=fopen(OUTPUT_FILE_PATH,"w");
		fclose(file);
	}
	
	t=getTime();
	nMAWs=find_MAWs_single(bbwt,MIN_MAW_LENGTH,WRITE_MAWS,OUTPUT_FILE_PATH);
	processingTime=getTime()-t;
	free_Basic_BWT(bbwt);
	
	if (WRITE_MAWS==1) fclose(file);
	
	printf( "%lld,%lld,%d,%f,%f,%f,%ld,%lld", 
	        sequence.inputLength,
	        sequence.length,
			sequence.hasRC,
			loadingTime,
			indexingTime,
			processingTime,
			malloc_count_peak(),
			nMAWs
	      );
	return 0;
}