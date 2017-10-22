#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "./malloc_count-master/malloc_count.h"  // For measuring memory usage
#include "io.h"
#include "DNA5_Basic_BWT.h"
#include "MAWs_single.h"


/** 
 * 
 */
int main(int argc, char **argv) {
	const unsigned char *INPUT_FILE_PATH = argv[0];
	const unsigned char APPEND_RC = atoi(argv[1]);
	const unsigned char MIN_MAW_LENGTH = atoi(argv[2]);
	const unsigned char N_CORES = atoi(argv[3]);
	const unsigned char STORE_MAWS = atoi(argv[4]);
	const unsigned char *OUTPUT_FILE_PATH = argv[5];
	double t, loadingTime, indexingTime, processingTime;
	unsigned long long nMAWs;
	Concatenation *sequence;
	Basic_BWT_t *bbwt;
	
	t=gettime();
	sequence=loadFASTA(INPUT_FILE_PATH,APPEND_RC);
	loadingTime=gettime()-t;
	
	t=gettime();
	bbwt=Build_BWT_index_from_text(sequence->buffer,sequence->length,Basic_bwt_free_text);
	indexingTime=gettime()-t;
	
	// Erasing output file
	if (STORE_MAWS==1) {
		file=fopen(OUTPUT_FILE_PATH,"w");
		fclose(file);
	}
	
	t=gettime();
	nMAWs=find_MAWs_single(bbwt,MIN_MAW_LENGTH,STORE_MAWS,OUTPUT_FILE_PATH);
	processingTime=gettime()-t;
	free_Basic_BWT(bbwt);
	
	if (STORE_MAWS==1) fclose(OUTPUT_FILE_PATH);
	
	printf( "%d,%d,%d,%d,%d,%d,%d,%d,%d", 
	        sequence->inputLength,
	        sequence->length,
			sequence->hasRC,
			N_CORES,
			loadingTime,
			indexingTime,
			processingTime,
			malloc_count_peak(),
			nMAWs
	      );
	return 0;
}


/**
 * In microseconds
 */
static double gettime() {
	struct timeval ttime;
	gettimeofday(&ttime,0);
	return ttime.tv_sec+ttime.tv_usec*0.000001;
}
