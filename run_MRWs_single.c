/**
 * @author Fabio Cunial, Filippo Gambarotto
 */
#include <limits.h>
#include "./malloc_count/malloc_count.h"  // For measuring memory usage
#include "./iterator/DNA5_Basic_BWT.h"
#include "./callbacks/MAWs_single.h"
#include "./io/io.h"
#include "./io/bufferedFileWriter.h"
#include "scores.h"
#include <jansson.h>

/**
 * For communicating with $scores.c$.
 */
extern unsigned char SELECTED_SCORE;
extern double SELECTED_SCORE_THRESHOLD;


int main(int argc, char **argv) {
	if (argc != 2) {
        	fprintf(stderr, "Usage: %s config.json\n", argv[0]);
        	exit(-1);
    	}

    	json_error_t error;
	json_t *root;

	char *INPUT_FILE_PATH;
	const uint8_t WRITE_MRWS;
	char *OUTPUT_FILE_PATH;
	const uint8_t N_THREADS;
	const uint8_t N_WORKPACKAGES_RATE;
	const uint64_t MIN_MRW_LENGTH;
	const uint64_t MAX_MRW_LENGTH;
	const uint64_t LOW_FREQ;
	const uint64_t HIGH_FREQ;
	const uint64_t MIN_HISTOGRAM_LENGTH;
	const uint64_t MAX_HISTOGRAM_LENGTH;
	const uint8_t COMPUTE_SCORES;
	uint8_t COMPRESS_OUTPUT;
	

	root = json_load_file(argv[1], 0, &error);

	if(!root){
		fprintf(stderr, "error: on line %d: %s\n", error.line, error.text);
		return 1;
	}

	json_unpack(root, "{s:s, s:I, s:I, s:I, s:I , s:I, s:I, s:I s:I, s:I, s:I, s:I, s:F, s:I, s:s}", "OUTPUT_FILE", &OUTPUT_FILE_PATH, "WRITE_MRWS", &WRITE_MRWS, "N_THREADS", &N_THREADS, "N_WORKPACKAGES_RATE", &N_WORKPACKAGES_RATE, "MIN_LENGTH", &MIN_MRW_LENGTH, "MAX_LENGTH", &MAX_MRW_LENGTH, "LOW_FREQUENCY", &LOW_FREQ, "HIGH_FREQUENCY", &HIGH_FREQ, "MIN_HISTOGRAM_LENGTH", &MIN_HISTOGRAM_LENGTH, "MAX_HISTOGRAM_LENGTH", &MAX_HISTOGRAM_LENGTH, "COMPUTE_SCORES", &COMPUTE_SCORES, "SELECT_SCORE", &SELECTED_SCORE, "SCORE_THRESHOLD", &SELECTED_SCORE_THRESHOLD, "COMPRESS_OUTPUT", &COMPRESS_OUTPUT, "INPUT_FILE", &INPUT_FILE_PATH);
	
	//atoi()

	printf("INPUT_FILE_PATH %s  \n", INPUT_FILE_PATH);
	printf("WRITE_MRWS %i  \n", WRITE_MRWS);
	printf("OUTPUT_FILE_PATH %s  \n", OUTPUT_FILE_PATH);
	printf("N_THREADS %i  \n", N_THREADS);
	printf("MIN_MRW_LENGTH %li \n", MIN_MRW_LENGTH);
	printf("MAX_MRW_LENGTH %li \n", MAX_MRW_LENGTH);
	printf("LOW_FREQUENCY %li \n", LOW_FREQ);
	printf("HIGH_FREQUENCY %li \n", HIGH_FREQ);
	printf("MIN_HISTOGRAM_LENGTH %li \n", MIN_HISTOGRAM_LENGTH);
	printf("MAX_HISTOGRAM_LENGTH %li \n", MAX_HISTOGRAM_LENGTH);
	printf("COMPUTE_SCORES %i  \n", COMPUTE_SCORES);
	printf("SELECTED_SCORE %i  \n", SELECTED_SCORE);
	printf("SELECTED_SCORE_THRESHOLD %f  \n", SELECTED_SCORE_THRESHOLD);
	printf("COMPRESS_OUTPUT %i  \n", COMPRESS_OUTPUT);

	uint64_t nBytes;
	double t, loadingTime, processingTime;
	BwtIndex_t *bbwt;
	MAWs_callback_state_t MRWs_state;
	ScoreState_t scoreState;
	
	// Loading the index
	t=getTime();
	bbwt=newBwtIndex();
	nBytes=deserializeBwtIndex(bbwt,INPUT_FILE_PATH);
	if (nBytes==0) {
		printf("ERROR while reading the index \n");
		return 1;
	}
	loadingTime=getTime()-t;
	
	// Initializing application state
	MRWs_initialize(&MRWs_state,bbwt->textLength,MIN_MRW_LENGTH,LOW_FREQ,HIGH_FREQ,MIN_HISTOGRAM_LENGTH,MAX_HISTOGRAM_LENGTH,WRITE_MRWS==0?NULL:OUTPUT_FILE_PATH,COMPRESS_OUTPUT);
	if (COMPUTE_SCORES!=0) {
		scoreInitialize(&scoreState,bbwt->dnaProbabilities,bbwt->logDnaProbabilities);
		MRWs_state.scoreState=&scoreState;
	}
	
	// Running the iterator
	t=getTime();
	if (N_THREADS==1) iterate_sequential( bbwt, 
	                                      MIN_MRW_LENGTH>=2?MIN_MRW_LENGTH-2:MIN_MRW_LENGTH,MAX_MRW_LENGTH-2,HIGH_FREQ,ULLONG_MAX,1,0,
                             			  MRWs_callback,cloneMAWState,mergeMAWState,MRWs_finalize,&MRWs_state,sizeof(MAWs_callback_state_t)
				                        );
	else iterate_parallel( bbwt,
				           MIN_MRW_LENGTH>=2?MIN_MRW_LENGTH-2:MIN_MRW_LENGTH,MAX_MRW_LENGTH-2,HIGH_FREQ,ULLONG_MAX,1,0,
						   N_THREADS, N_WORKPACKAGES_RATE,
					       MRWs_callback,cloneMAWState,mergeMAWState,MRWs_finalize,&MRWs_state,sizeof(MAWs_callback_state_t)
					     );
	processingTime=getTime()-t;
	printf( "%llu,%llu,%llu,%llu,%llu|%lf,%lf|%llu|%llu,%llu,%llu,%lf \n", 
		    (long long unsigned int)(bbwt->textLength),
			(long long unsigned int)(MIN_MRW_LENGTH),
			(long long unsigned int)(MAX_MRW_LENGTH),
			(long long unsigned int)(LOW_FREQ),
			(long long unsigned int)(HIGH_FREQ),
			
			loadingTime,
			processingTime,
			
			(long long unsigned int)malloc_count_peak(),
			
			(long long unsigned int)(MRWs_state.nMAWs),
			(long long unsigned int)(MRWs_state.minObservedLength),
			(long long unsigned int)(MRWs_state.maxObservedLength),
			((double)MRWs_state.nMAWMaxreps)/MRWs_state.nMaxreps
	      );
	if (MIN_HISTOGRAM_LENGTH>0) printLengthHistogram(&MRWs_state);
	
	// Finalizing application state
	if (COMPUTE_SCORES!=0) scoreFinalize(&scoreState);
	freeBwtIndex(bbwt);
	
	return 0;
}
