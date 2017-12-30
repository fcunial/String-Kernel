#include <stdio.h>
#include <stdlib.h>
#include "./malloc_count/malloc_count.h"  // For measuring memory usage
#include "./io/io.h"
#include "./iterator/DNA5_Basic_BWT.h"
#include "./callbacks/MAWs_single.h"
#include "./callbacks/lengthScores.h"


/** 
 * 1: input file path;
 * 2: append reverse-complement (1/0);
 * 3: minimum MRW length;
 * 4: minFreq;
 * 5: maxFreq;
 * 6: min histogram length;
 * 7: max histogram length;
 * 8: write MRWs to a file (1/0);
 * 9: assigns a score to each MRW (1/0); used only if MRWs are written to a file;
 * 10: output file path (read only if the previous argument is 1).
 */
int main(int argc, char **argv) {
	char *INPUT_FILE_PATH = argv[1];
	const unsigned char APPEND_RC = atoi(argv[2]);
	const unsigned int MIN_MRW_LENGTH = atoi(argv[3]);
	const unsigned int MIN_FREQ = atoi(argv[4]);
	const unsigned int MAX_FREQ = atoi(argv[5]);
	const unsigned int MIN_HISTOGRAM_LENGTH = atoi(argv[6]);
	const unsigned int MAX_HISTOGRAM_LENGTH = atoi(argv[7]);
	const unsigned char WRITE_MRWS = atoi(argv[8]);
	const unsigned char COMPUTE_SCORES = atoi(argv[9]);
	char *OUTPUT_FILE_PATH = NULL;
	if (WRITE_MRWS==1) OUTPUT_FILE_PATH=argv[10];
	double t, tPrime, loadingTime, indexingTime, processingTime;
	FILE *file;
	Concatenation sequence;
	Basic_BWT_t *bbwt;
	MAWs_callback_state_t MRWs_state;
	SLT_iterator_t_single_string SLT_iterator;

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
	
	MRWs_initialize(&MRWs_state,sequence.length,MIN_MRW_LENGTH,MIN_FREQ,MAX_FREQ,MIN_HISTOGRAM_LENGTH,MAX_HISTOGRAM_LENGTH,WRITE_MRWS,COMPUTE_SCORES,OUTPUT_FILE_PATH);
	MRWs_state.lengthScoreCallback=lengthScore2;
	SLT_iterator=new_SLT_iterator(MRWs_callback,&MRWs_state,bbwt,SLT_stack_trick);
	t=getTime();
	SLT_execute_iterator(&SLT_iterator);
	processingTime=getTime()-t;
	MRWs_finalize(&MRWs_state);
	free_Basic_BWT(bbwt);
	
	printf( "%lu,%lu,%u,%u,%u,%u,%lf,%lf,%lf,%llu,%u \n", 
	        sequence.inputLength,
	        sequence.length,
			sequence.hasRC,
			MIN_MRW_LENGTH,
			MIN_FREQ,
			MAX_FREQ,
			loadingTime,
			indexingTime,
			processingTime,
			(unsigned long long)malloc_count_peak(),
			MRWs_state.nMAWs
	      );
	if (MIN_HISTOGRAM_LENGTH>0) printLengthHistogram(&MRWs_state);
	return 0;
}