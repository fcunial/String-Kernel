#include <stdio.h>
#include <stdlib.h>
#include "./malloc_count/malloc_count.h"  // For measuring memory usage
#include "./io/io.h"
#include "./iterator/DNA5_Basic_BWT.h"
#include "./callbacks/MAWs_single.h"


/** 
 * 1: input file path;
 * 2: append reverse-complement (1/0);
 * 3: minimum MAW length;
 * 4: write MAWs to a file (1/0);
 * 5: assigns scores to each MAW (1/0); used only if MAWs are written to a file;
 * 6: output file path (read only if the previous argument is 1).
 */
int main(int argc, char **argv) {
	char *INPUT_FILE_PATH = argv[1];
	const unsigned char APPEND_RC = atoi(argv[2]);
	const unsigned int MIN_MAW_LENGTH = atoi(argv[3]);
	const unsigned char WRITE_MAWS = atoi(argv[4]);
	const unsigned char COMPUTE_SCORES = atoi(argv[5]);
	char *OUTPUT_FILE_PATH = NULL;
	if (WRITE_MAWS==1) OUTPUT_FILE_PATH=argv[6];
	double t, tPrime, loadingTime, indexingTime, processingTime;
	FILE *file;
	Concatenation sequence;
	Basic_BWT_t *bbwt;
	MAWs_callback_state_t MAWs_state;
	SLT_iterator_t_single_string SLT_iterator;

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
	
	MAWs_initialize(&MAWs_state,sequence.length,MIN_MAW_LENGTH,WRITE_MAWS,COMPUTE_SCORES,OUTPUT_FILE_PATH);
	SLT_iterator=new_SLT_iterator(MAWs_callback,&MAWs_state,bbwt,SLT_stack_trick);
	t=getTime();
	SLT_execute_iterator(&SLT_iterator);
	processingTime=getTime()-t;
	MAWs_finalize(&MAWs_state);
	free_Basic_BWT(bbwt);
	
	printf( "%lu,%lu,%u,%u,%lf,%lf,%lf,%llu,%u \n", 
	        sequence.inputLength,
	        sequence.length,
			sequence.hasRC,
			MIN_MAW_LENGTH,
			loadingTime,
			indexingTime,
			processingTime,
			(unsigned long long)malloc_count_peak(),
			MAWs_state.nMAWs
	      );
	return 0;
}