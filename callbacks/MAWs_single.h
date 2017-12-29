#ifndef MAWs_single_h
#define MAWs_single_h


#include <stdio.h>
#include "../iterator/SLT_single_string.h"


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
	
	// Minimal rare words
	unsigned int minFreq, maxFreq;
} MAWs_callback_state_t;


void MAWs_callback(const SLT_params_t SLT_params, void *intern_state);


/**
 * @param minLength (>=2) considers only MAWs of length at least $minLength$;
 * @param writeMAWs 0 iff MAWs should not be written to the output; otherwise, MAWs are  
 * appended to file $filePath$;
 * @param computeScores 0 iff scores should not be added to $filePath$; used only if
 * $writeMAWs$ is nonzero.
 */
void MAWs_initialize( MAWs_callback_state_t *state,
			    	  unsigned int textLength, 
					  unsigned int minLength, 
					  unsigned char writeMAWs, 
					  unsigned char computeScores, 
					  char *filePath );


/**
 * Flushes the output buffers one more time, if any, and frees space.
 */
void MAWs_finalize(MAWs_callback_state_t *state);


void MRWs_callback(const SLT_params_t SLT_params, void *intern_state);


/**
 * Detects minimal rare words $W$ such that $minFreq \leq f(W) < maxFreq$ and 
 * $f(V) \geq maxFreq$ for every substring $V$ of $W$.
 *
 * @param minLength (>=2) considers only MRWs of length at least $minLength$;
 * @param writeMAWs 0 iff MAWs should not be written to the output; otherwise, MAWs are  
 * appended to file $filePath$;
 * @param computeScores 0 iff scores should not be added to $filePath$; used only if
 * $writeMAWs$ is nonzero.
 */
void MRWs_initialize( MAWs_callback_state_t *state,
			    	  unsigned int textLength, 
					  unsigned int minLength, 
					  unsigned int minFreq, 
					  unsigned int maxFreq, 
					  unsigned char writeMRWs, 
					  unsigned char computeScores, 
					  char *filePath );


/**
 * Flushes the output buffers one more time, if any, and frees space.
 */
void MRWs_finalize(MAWs_callback_state_t *state);


#endif