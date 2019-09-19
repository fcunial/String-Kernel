/**
 * @author Fabio Cunial, Filippo Gambarotto
 */
#include "./malloc_count/malloc_count.h"  // For measuring memory usage
#include "./iterator/DNA5_Basic_BWT.h"
#include "./callbacks/MAWs_single.h"
#include "./io/io.h"
#include "./io/bufferedFileWriter.h"
#include "scores.h"

/**
 * For communicating with $scores.c$.
 */
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
 * 11: output file path (used only if argument 6 equals 1). If the file already exists, 
 *     its content is overwritten.
 * 12: max length of a MAW.
 */
int main(int argc, char **argv) {
	char *INPUT_FILE_PATH = argv[1];
	const uint8_t APPEND_RC = atoi(argv[2]);
	const uint64_t MIN_MAW_LENGTH = atoi(argv[3]);
	const uint64_t MIN_HISTOGRAM_LENGTH = atoi(argv[4]);
	const uint64_t MAX_HISTOGRAM_LENGTH = atoi(argv[5]);
	const uint8_t WRITE_MAWS = atoi(argv[6]);
	const uint8_t COMPUTE_SCORES = atoi(argv[7]);
	SELECTED_SCORE = atoi(argv[8]);
	SELECTED_SCORE_THRESHOLD = atof(argv[9]);
	const uint8_t COMPRESS_OUTPUT = atoi(argv[10]);
	char *OUTPUT_FILE_PATH = NULL;
	if (WRITE_MAWS==1) OUTPUT_FILE_PATH=argv[11];
	uint64_t MAX_LENGTH = atoi(argv[12]);
	double t, tPrime, loadingTime, indexingTime, processingTime;
	Concatenation sequence;
	BwtIndex_t *bbwt;
	UnaryIterator_t iterator;
	MAWs_callback_state_t MAWs_state;
	FILE *file;
	BufferedFileWriter_t bufferedFileWriter;
	ScoreState_t scoreState;

	// Building the BWT
	t=getTime();
	sequence=loadFASTA(INPUT_FILE_PATH,APPEND_RC);
	tPrime=getTime();
	loadingTime=tPrime-t;
	t=tPrime;
	bbwt=buildBwtIndex(sequence.buffer,sequence.length,Basic_bwt_free_text);
	indexingTime=getTime()-t;
	
	// Initializing application state
	if (WRITE_MAWS!=0) {
		file=fopen(OUTPUT_FILE_PATH,"w");  // Cleaning the old content of the file
		fclose(file);
		initializeBufferedFileWriter(&bufferedFileWriter,OUTPUT_FILE_PATH);
	}
	MAWs_initialize(&MAWs_state,sequence.length,MIN_MAW_LENGTH,MIN_HISTOGRAM_LENGTH,MAX_HISTOGRAM_LENGTH,WRITE_MAWS==0?NULL:&bufferedFileWriter,COMPRESS_OUTPUT);
	if (COMPUTE_SCORES!=0) {
		scoreInitialize(&scoreState);
		MAWs_state.scoreState=&scoreState;
	}
	
	// Running the iterator
	iterator=newIterator(MAWs_callback,&MAWs_state,bbwt,MAX_LENGTH-2);
	t=getTime();
	iterate(&iterator);
	processingTime=getTime()-t;
	printf( "%llu,%llu,%u,%llu|%lf,%lf,%lf|%llu|%llu,%llu,%llu,%lf \n", 
	        (long long unsigned int)(sequence.inputLength),
	        (long long unsigned int)(sequence.length),
			sequence.hasRC,
			(long long unsigned int)(MIN_MAW_LENGTH),
			
			loadingTime,
			indexingTime,
			processingTime,
			
			(long long unsigned int)malloc_count_peak(),
			
			(long long unsigned int)(MAWs_state.nMAWs),
			(long long unsigned int)(MAWs_state.minObservedLength),
			(long long unsigned int)(MAWs_state.maxObservedLength),
			((double)MAWs_state.nMAWMaxreps)/MAWs_state.nMaxreps
	      );
	if (MIN_HISTOGRAM_LENGTH>0) printLengthHistogram(&MAWs_state);

	// Finalizing application state
	MAWs_finalize(&MAWs_state);
	if (WRITE_MAWS!=0) finalizeBufferedFileWriter(&bufferedFileWriter);
	if (COMPUTE_SCORES!=0) scoreFinalize(&scoreState);
	freeBwtIndex(bbwt);
	
	return 0;
}