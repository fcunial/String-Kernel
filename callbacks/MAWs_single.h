#ifndef MAWs_single_h
#define MAWs_single_h


#include <stdio.h>
#include "../iterator/SLT_single_string.h"


/** 
 * Function of string length used to score each MAW. Must be defined by the user.
 */
typedef double (*lengthScore_callback_t)(unsigned int length);


typedef struct {
	unsigned int textLength;
	unsigned int minLength;  // Minimum length of a MAW to be reported
	unsigned int nMAWs;  // Total number of reported MAWs
	
	// Character stack
	unsigned char *char_stack;  // Indexed from zero
	unsigned int char_stack_capacity;  // Number of characters in the stack
	
	// Output buffer
	unsigned char writeMAWs;  // 0 iff MAWs should not be written to the output
	char *MAWs_buffer;
	unsigned int MAWs_buffer_capacity;  // Maximum number of chars in the buffer
	unsigned int MAWs_buffer_size;  // Number of chars currently in the buffer
	FILE *file;
	
	// Scores
	unsigned char computeScores;
	unsigned int *leftFreqs, *rightFreqs;  // Frequency of each left/right extension. Only for ACGT, indexed from zero.
	char *scoreBuffer;  // Temporary space for the string representation of a score
	lengthScore_callback_t lengthScoreCallback;  // To compute a length-based score
	double *score_stack;  // Indexed from zero
	
	// Histograms
	unsigned int lengthHistogramMin, lengthHistogramMax, lengthHistogramSize;
	unsigned int *lengthHistogram;
	
	// Minimal rare words
	unsigned int minFreq, maxFreq;
} MAWs_callback_state_t;


void MAWs_callback(const SLT_params_t SLT_params, void *intern_state);


/**
 * @param minLength (>=2) considers only MAWs of length at least $minLength$;
 *
 * @param lengthHistogramMin,lengthHistogramMax computes the number of MAWs with length
 * $i$ for all $i \in [lengthHistogramMin..lengthHistogramMax]$; the first (respectively,
 * last) cell of the histogram contains the number of MAWs with length at most (at least)
 * equal to the corresponding length; no histogram is computed if $lengthHistogramMin==0$;
 *
 * @param writeMAWs 0 iff MAWs should not be written to the output; otherwise, MAWs are  
 * appended to file $filePath$;
 *
 * @param computeScores 0 iff scores should not be added to $filePath$; used only if
 * $writeMAWs$ is nonzero.
 */
void MAWs_initialize( MAWs_callback_state_t *state,
			    	  unsigned int textLength, 
					  unsigned int minLength, 
					  unsigned int lengthHistogramMin,
					  unsigned int lengthHistogramMax,
					  unsigned char writeMAWs, 
					  unsigned char computeScores, 
					  char *filePath );


/**
 * Flushes the output buffers one more time, if any, and frees space.
 */
void MAWs_finalize(MAWs_callback_state_t *state);


void printLengthHistogram(MAWs_callback_state_t *state);


void MRWs_callback(const SLT_params_t SLT_params, void *intern_state);


/**
 * Detects minimal rare words $W$ such that $minFreq \leq f(W) < maxFreq$ and 
 * $f(V) \geq maxFreq$ for every substring $V$ of $W$. 
 * See $MAWs_initialize$ for details on the arguments.
 */
void MRWs_initialize( MAWs_callback_state_t *state,
			    	  unsigned int textLength, 
					  unsigned int minLength, 
					  unsigned int minFreq, 
					  unsigned int maxFreq, 
					  unsigned int lengthHistogramMin,
					  unsigned int lengthHistogramMax,
					  unsigned char writeMRWs, 
					  unsigned char computeScores, 
					  char *filePath );


/**
 * Flushes the output buffers one more time, if any, and frees space.
 */
void MRWs_finalize(MAWs_callback_state_t *state);


#endif