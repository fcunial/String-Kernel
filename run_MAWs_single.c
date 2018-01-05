#include <stdio.h>
#include "./malloc_count/malloc_count.h"  // For measuring memory usage
#include "./iterator/DNA5_Basic_BWT.h"
#include "./callbacks/MAWs_single.h"
#include "./io/io.h"
#include "./io/bufferedFileWriter.h"
#include "scores.h"


extern unsigned char SELECTED_SCORE;
extern double SELECTED_SCORE_THRESHOLD;


/** 
 * 1: input file path;
 * 2: append reverse-complement (1/0);
 * 3: minimum MAW length;
 * 4: min histogram length;
 * 5: max histogram length;
 * 6: write MAWs to a file (1/0);
 * 7: assigns scores to each MAW (1/0); used only if MAWs are written to a file;
 * 8: score ID for selecting MAWs;
 * 9: min absolute value of a score for a MAW to be selected;
 * 10: compresses output (1/0); used only if MAWs are written to a file and scores are not
 *     computed;
 * 11: output file path (read only if the previous argument is 1).
 */
int main(int argc, char **argv) {
	char *INPUT_FILE_PATH = argv[1];
	const unsigned char APPEND_RC = atoi(argv[2]);
	const unsigned int MIN_MAW_LENGTH = atoi(argv[3]);
	const unsigned int MIN_HISTOGRAM_LENGTH = atoi(argv[4]);
	const unsigned int MAX_HISTOGRAM_LENGTH = atoi(argv[5]);
	const unsigned char WRITE_MAWS = atoi(argv[6]);
	const unsigned char COMPUTE_SCORES = atoi(argv[7]);
	SELECTED_SCORE = atoi(argv[8]);
	SELECTED_SCORE_THRESHOLD = atof(argv[9]);
	const unsigned char COMPRESS_OUTPUT = atoi(argv[10]);
	char *OUTPUT_FILE_PATH = NULL;
	if (WRITE_MAWS==1) OUTPUT_FILE_PATH=argv[11];
	double t, tPrime, loadingTime, indexingTime, processingTime;
	Concatenation sequence;
	Basic_BWT_t *bbwt;
	SLT_iterator_t_single_string SLT_iterator;
	MAWs_callback_state_t MAWs_state;
	FILE *file;
	buffered_file_writer_t bufferedFileWriter;
	score_state_t scoreState;

	// Building the BWT
	t=getTime();
	sequence=loadFASTA(INPUT_FILE_PATH,APPEND_RC);
	tPrime=getTime();
	loadingTime=tPrime-t;
	t=tPrime;
	bbwt=Build_BWT_index_from_text(sequence.buffer,sequence.length,Basic_bwt_free_text);
	indexingTime=getTime()-t;
	
	// Initializing application state
	if (WRITE_MAWS!=0) {
		file=fopen(OUTPUT_FILE_PATH,"w");
		fclose(file);
		initializeBufferedFileWriter(&bufferedFileWriter,OUTPUT_FILE_PATH);
	}
	MAWs_initialize(&MAWs_state,sequence.length,MIN_MAW_LENGTH,MIN_HISTOGRAM_LENGTH,MAX_HISTOGRAM_LENGTH,WRITE_MAWS==0?NULL:&bufferedFileWriter,COMPRESS_OUTPUT);
	if (COMPUTE_SCORES!=0) {
		scoreInitialize(&scoreState);
		MAWs_state.scoreState=&scoreState;
	}
	
	// Running the iterator
	SLT_iterator=new_SLT_iterator(MAWs_callback,&MAWs_state,bbwt,SLT_stack_trick);
	t=getTime();
	SLT_execute_iterator(&SLT_iterator);
	processingTime=getTime()-t;
	printf( "%lu,%lu,%u,%u,%lf,%lf,%lf,%llu,%u,%u,%lf \n", 
	        sequence.inputLength,
	        sequence.length,
			sequence.hasRC,
			MIN_MAW_LENGTH,
			loadingTime,
			indexingTime,
			processingTime,
			(unsigned long long)malloc_count_peak(),
			MAWs_state.nMAWs,
			MAWs_state.maxLength,
			((double)MAWs_state.nMAWMaxreps)/MAWs_state.nMaxreps
	      );
	if (MIN_HISTOGRAM_LENGTH>0) printLengthHistogram(&MAWs_state);

	// Finalizing application state
	MAWs_finalize(&MAWs_state);
	if (WRITE_MAWS!=0) finalizeBufferedFileWriter(&bufferedFileWriter);
	if (COMPUTE_SCORES!=0) scoreFinalize(&scoreState);
	free_Basic_BWT(bbwt);
	
	return 0;
}