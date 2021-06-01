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
	const uint8_t WRITE_MAWS;
	char *OUTPUT_FILE_PATH;
	const uint8_t N_THREADS;
	const uint8_t N_WORKPACKAGES_RATE;
	const uint64_t MIN_MAW_LENGTH;
	const uint64_t MAX_MAW_LENGTH;
	const uint64_t MIN_HISTOGRAM_LENGTH;
	const uint64_t MAX_HISTOGRAM_LENGTH;
	const uint8_t COMPUTE_SCORES;
	uint8_t COMPRESS_OUTPUT;
	

	root = json_load_file(argv[1], 0, &error);

	if(!root){
		fprintf(stderr, "error: on line %d: %s\n", error.line, error.text);
		return 1;
	}

	json_unpack(root, "{s:s, s:I, s:I, s:I, s:I ,s:I, s:I, s:I, s:I, s:I, s:F, s:I, s:s}", "OUTPUT_FILE", &OUTPUT_FILE_PATH, "WRITE_MAWS", &WRITE_MAWS, "N_THREADS", &N_THREADS, "N_WORKPACKAGES_RATE", &N_WORKPACKAGES_RATE ,"MIN_LENGTH", &MIN_MAW_LENGTH, "MAX_LENGTH", &MAX_MAW_LENGTH, "MIN_HISTOGRAM_LENGTH", &MIN_HISTOGRAM_LENGTH, "MAX_HISTOGRAM_LENGTH", &MAX_HISTOGRAM_LENGTH, "COMPUTE_SCORES", &COMPUTE_SCORES, "SELECT_SCORE", &SELECTED_SCORE, "SCORE_THRESHOLD", &SELECTED_SCORE_THRESHOLD, "COMPRESS_OUTPUT", &COMPRESS_OUTPUT, "INPUT_FILE", &INPUT_FILE_PATH);
	
	//atoi()

	printf("INPUT_FILE_PATH %s  \n", INPUT_FILE_PATH);
	printf("WRITE_MAWS %i  \n", WRITE_MAWS);
	printf("OUTPUT_FILE_PATH %s  \n", OUTPUT_FILE_PATH);
	printf("N_THREADS %i  \n", N_THREADS);
	printf("N_WORKPACKAGES_RATE %i  \n", N_WORKPACKAGES_RATE);
	printf("MIN_MAW_LENGTH %li \n", MIN_MAW_LENGTH);
	printf("MAX_MAW_LENGTH %li \n", MAX_MAW_LENGTH);
	printf("MIN_HISTOGRAM_LENGTH %li \n", MIN_HISTOGRAM_LENGTH);
	printf("MAX_HISTOGRAM_LENGTH %li \n", MAX_HISTOGRAM_LENGTH);
	printf("COMPUTE_SCORES %i  \n", COMPUTE_SCORES);
	printf("SELECTED_SCORE %i  \n", SELECTED_SCORE);
	printf("SELECTED_SCORE_THRESHOLD %f  \n", SELECTED_SCORE_THRESHOLD);
	printf("COMPRESS_OUTPUT %i  \n", COMPRESS_OUTPUT);

	uint64_t nBytes;
	double t, loadingTime, processingTime;
	BwtIndex_t *bbwt;
	MAWs_callback_state_t MAWs_state;
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
	MAWs_initialize(&MAWs_state,bbwt->textLength,MIN_MAW_LENGTH,MIN_HISTOGRAM_LENGTH,MAX_HISTOGRAM_LENGTH,WRITE_MAWS==0?NULL:OUTPUT_FILE_PATH,COMPRESS_OUTPUT);
	if (COMPUTE_SCORES!=0) {
		scoreInitialize(&scoreState,bbwt->dnaProbabilities,bbwt->logDnaProbabilities);
		MAWs_state.scoreState=&scoreState;
	}
	
	// Running the iterator
	t=getTime();
	if (N_THREADS==1) iterate_sequential( bbwt, 
	                                      MIN_MAW_LENGTH>=2?MIN_MAW_LENGTH-2:MIN_MAW_LENGTH,MAX_MAW_LENGTH-2,0,ULLONG_MAX,1,0,
                             			  MAWs_callback,cloneMAWState,mergeMAWState,MAWs_finalize,&MAWs_state,sizeof(MAWs_callback_state_t)
				                        );
	else iterate_parallel( bbwt,
				           MIN_MAW_LENGTH>=2?MIN_MAW_LENGTH-2:MIN_MAW_LENGTH,MAX_MAW_LENGTH-2,0,ULLONG_MAX,1,0,
						   N_THREADS,N_WORKPACKAGES_RATE,
					       MAWs_callback,cloneMAWState,mergeMAWState,MAWs_finalize,&MAWs_state,sizeof(MAWs_callback_state_t)
					     );
	processingTime=getTime()-t;
	printf( "%llu,%llu,%llu|%lf,%lf|%llu|%llu,%llu,%llu,%lf \n", 
	        (long long unsigned int)(bbwt->textLength),
			(long long unsigned int)(MIN_MAW_LENGTH),
			(long long unsigned int)(MAX_MAW_LENGTH),
			
			loadingTime,
			processingTime,
			
			(long long unsigned int)malloc_count_peak(),
			
			(long long unsigned int)(MAWs_state.nMAWs),
			(long long unsigned int)(MAWs_state.minObservedLength),
			(long long unsigned int)(MAWs_state.maxObservedLength),
			((double)MAWs_state.nMAWMaxreps)/MAWs_state.nMaxreps
	      );
	if (MIN_HISTOGRAM_LENGTH>0) printLengthHistogram(&MAWs_state);

	// Finalizing application state
	if (COMPUTE_SCORES!=0) scoreFinalize(&scoreState);
	freeBwtIndex(bbwt);
	
	return 0;
}
