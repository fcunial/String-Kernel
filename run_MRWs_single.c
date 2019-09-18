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
 * 3: minimum MRW length;
 * 4: minFreq;
 * 5: maxFreq;
 * 6: min histogram length;
 * 7: max histogram length;
 * 8: write MRWs to a file (1/0);
 * 9: assigns a score to each MRW (1/0); used only if MRWs are written to a file;
 * 10: score ID for selecting MRWs;
 * 11: min absolute value of a score for a MRW to be selected;
 * 12: compresses output (1/0); used only if MRWs are written to a file and scores are not
 *     computed;
 * 13: output file path (read only if the previous argument is 1);
 * 14: max length of a MRW.
 */
int main(int argc, char **argv) {
	char *INPUT_FILE_PATH = argv[1];
	const uint8_t APPEND_RC = atoi(argv[2]);
	const uint64_t MIN_MRW_LENGTH = atoi(argv[3]);
	const uint64_t MIN_FREQ = atoi(argv[4]);
	const uint64_t MAX_FREQ = atoi(argv[5]);
	const uint64_t MIN_HISTOGRAM_LENGTH = atoi(argv[6]);
	const uint64_t MAX_HISTOGRAM_LENGTH = atoi(argv[7]);
	const uint8_t WRITE_MRWS = atoi(argv[8]);
	const uint8_t COMPUTE_SCORES = atoi(argv[9]);
	SELECTED_SCORE = atoi(argv[10]);
	SELECTED_SCORE_THRESHOLD = atof(argv[11]);
	const uint8_t COMPRESS_OUTPUT = atoi(argv[12]);
	char *OUTPUT_FILE_PATH = NULL;
	if (WRITE_MRWS==1) OUTPUT_FILE_PATH=argv[13];
	uint64_t MAX_LENGTH = atoi(argv[14]);
	double t, tPrime, loadingTime, indexingTime, processingTime;
	Concatenation sequence;
	BwtIndex_t *bbwt;
	UnaryIterator_t SLT_iterator;
	MAWs_callback_state_t MRWs_state;
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
	if (WRITE_MRWS!=0) {
		file=fopen(OUTPUT_FILE_PATH,"w");
		fclose(file);
		initializeBufferedFileWriter(&bufferedFileWriter,OUTPUT_FILE_PATH);
	}
	MRWs_initialize(&MRWs_state,sequence.length,MIN_MRW_LENGTH,MIN_FREQ,MAX_FREQ,MIN_HISTOGRAM_LENGTH,MAX_HISTOGRAM_LENGTH,WRITE_MRWS==0?NULL:&bufferedFileWriter,COMPRESS_OUTPUT);
	if (COMPUTE_SCORES!=0) {
		scoreInitialize(&scoreState);
		MRWs_state.scoreState=&scoreState;
	}
	
	// Running the iterator
	SLT_iterator=newIterator(MRWs_callback,&MRWs_state,bbwt,MAX_LENGTH-2);
	t=getTime();
	iterate(&SLT_iterator);
	processingTime=getTime()-t;
	printf( "%llu,%llu,%u,%llu,%llu,%llu|%lf,%lf,%lf|%llu|%llu,%llu,%lf \n", 
	        sequence.inputLength,
	        sequence.length,
			sequence.hasRC,
			MIN_MRW_LENGTH,
			MIN_FREQ,
			MAX_FREQ,
			loadingTime,
			indexingTime,
			processingTime,
			(uint64_t)malloc_count_peak(),
			MRWs_state.nMAWs,
			MRWs_state.maxLength,
			((double)MRWs_state.nMAWMaxreps)/MRWs_state.nMaxreps
	      );
	if (MIN_HISTOGRAM_LENGTH>0) printLengthHistogram(&MRWs_state);
	
	// Finalizing application state
	MRWs_finalize(&MRWs_state);
	if (WRITE_MRWS!=0) finalizeBufferedFileWriter(&bufferedFileWriter);
	if (COMPUTE_SCORES!=0) scoreFinalize(&scoreState);
	freeBwtIndex(bbwt);
	
	return 0;
}