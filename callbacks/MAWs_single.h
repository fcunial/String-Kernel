#ifndef MAWs_single_h
#define MAWs_single_h


#include "../iterator/SLT_single_string.h"
#include "../io/bufferedFileWriter.h"
#include "../scores.h"


typedef struct {
	unsigned int textLength;
	unsigned int minLength;  // Minimum length of a MAW to be reported
	unsigned int nMAWs;  // Total number of reported MAWs
	unsigned int maxLength;  // Maximum observed length of a MAW
	unsigned int nMaxreps;  // Number of visited maximal repeats
	unsigned int nMAWMaxreps;  // N. of visited maxreps that are the infix of a MAW
	
	// Character stack
	unsigned long *char_stack;  // We push numbers in [0..3] of two bits each.
	unsigned int char_stack_capacity;  // Number of characters that can fit in the stack
	
	// Output buffer
	buffered_file_writer_t *outputFile;
	
	// Scores
	unsigned int *leftFreqs, *rightFreqs;  // Frequency of each left/right extension. Only for ${A,C,G,T}$, indexed from zero.
	score_state_t *scoreState;
	
	// Histograms
	unsigned int lengthHistogramMin, lengthHistogramMax, lengthHistogramSize;
	unsigned int *lengthHistogram;
	
	// Compressed output
	unsigned char compressOutput;  // 0 iff MAWs should not be represented in compressed form in the output.
	unsigned long *compressionBuffers[4][4][4];  // Bits inside each long in the buffer are assumed to be stored from LSB to MSB.
	unsigned int compressionBuffersLength[4][4][4];
	unsigned int compressionBuffersCapacity[4][4][4];
	unsigned long *runs_stack;  // Tells whether a node of the ST is a run of a single character or not (1/0)
	
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
 * @param file NULL iff MAWs should not be written to the output; otherwise, MAWs are  
 * appended to $file$.
 */
void MAWs_initialize( MAWs_callback_state_t *state,
			    	  unsigned int textLength, 
					  unsigned int minLength, 
					  unsigned int lengthHistogramMin,
					  unsigned int lengthHistogramMax,
					  buffered_file_writer_t *file,
					  unsigned char compressOutput );


/**
 * Flushes the output buffers one more time, if any, and frees up space.
 */
void MAWs_finalize(MAWs_callback_state_t *state);


void printLengthHistogram(MAWs_callback_state_t *state);


void MRWs_callback(const SLT_params_t SLT_params, void *intern_state);


/**
 * Detects minimal rare words $W$ such that $minFreq \leq f(W) < maxFreq$ and 
 * $f(V) \geq maxFreq$ for every substring $V$ of $W$.
 * See $MAWs_initialize$ for details on the input arguments.
 */
void MRWs_initialize( MAWs_callback_state_t *state,
			    	  unsigned int textLength, 
					  unsigned int minLength, 
					  unsigned int minFreq, 
					  unsigned int maxFreq, 
					  unsigned int lengthHistogramMin,
					  unsigned int lengthHistogramMax,
					  buffered_file_writer_t *file,
					  unsigned char compressOutput );


/**
 * Flushes the output buffers one more time, if any, and frees up space.
 */
void MRWs_finalize(MAWs_callback_state_t *state);


/**
 * Creates a clone of $from$, with all statistics reset to zero.
 * If $from$ has an output file, the path of the output file of the clone has prefix  
 * $pathPrefix$ and is followed by $id$.
 */
void cloneMAWState(MAWs_callback_state_t *from, MAWs_callback_state_t *to, char *pathPrefix, char id, char *cloneBuffer);


/**
 * Merges the statistics of $from$ into those of $to$.
 */
void mergeMAWState(MAWs_callback_state_t *from, MAWs_callback_state_t *to);


#endif