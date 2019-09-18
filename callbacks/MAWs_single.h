/**
 * Unary-iterator callback for computing the minimal absent and the minimal rare words of
 * a single string.
 *
 * @author Fabio Cunial
 */
#ifndef MAWs_single_h
#define MAWs_single_h


#include "../iterator/SLT_single_string.h"
#include "../io/bufferedFileWriter.h"
#include "../scores.h"


typedef struct {
	uint64_t textLength;
	uint64_t minLength;  // Minimum length of a MAW to be reported
	uint64_t nMAWs;  // Total number of reported MAWs
	uint64_t maxLength;  // Maximum observed length of a MAW
	uint64_t nMaxreps;  // Number of visited maximal repeats
	uint64_t nMAWMaxreps;  // N. of visited maxreps that are the infix of a MAW
	
	// Character stack
	uint64_t *char_stack;  // We push numbers in [0..3] of two bits each.
	uint64_t char_stack_capacity;  // Number of characters that can fit in the stack
	
	// Output buffer
	BufferedFileWriter_t *outputFile;
	
	// Scores
	uint64_t *leftFreqs, *rightFreqs;  // Frequency of each left/right extension. Only for ${A,C,G,T}$, indexed from zero.
	ScoreState_t *scoreState;
	
	// Histograms
	uint64_t lengthHistogramMin, lengthHistogramMax, lengthHistogramSize;
	uint64_t *lengthHistogram;
	
	// Compressed output
	uint8_t compressOutput;  // 0 iff MAWs should not be represented in compressed form in the output.
	uint64_t *compressionBuffers[4][4][4];  // Bits inside each long in the buffer are assumed to be stored from LSB to MSB.
	uint64_t compressionBuffersLength[4][4][4];
	uint64_t compressionBuffersCapacity[4][4][4];
	uint64_t *runs_stack;  // Tells whether a node of the ST is a run of a single character or not (1/0)
	
	// Minimal rare words
	uint64_t minFreq, maxFreq;
} MAWs_callback_state_t;


void MAWs_callback(const RightMaximalString_t RightMaximalString, void *applicationData);


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
			    	  uint64_t textLength, 
					  uint64_t minLength, 
					  uint64_t lengthHistogramMin,
					  uint64_t lengthHistogramMax,
					  BufferedFileWriter_t *file,
					  uint8_t compressOutput );


/**
 * Flushes the output buffers one more time, if any, and frees up space.
 */
void MAWs_finalize(MAWs_callback_state_t *state);


void printLengthHistogram(MAWs_callback_state_t *state);


void MRWs_callback(const RightMaximalString_t RightMaximalString, void *applicationData);


/**
 * Detects minimal rare words $W$ such that $minFreq \leq f(W) < maxFreq$ and 
 * $f(V) \geq maxFreq$ for every substring $V$ of $W$.
 * See $MAWs_initialize$ for details on the input arguments.
 */
void MRWs_initialize( MAWs_callback_state_t *state,
			    	  uint64_t textLength, 
					  uint64_t minLength, 
					  uint64_t minFreq, 
					  uint64_t maxFreq, 
					  uint64_t lengthHistogramMin,
					  uint64_t lengthHistogramMax,
					  BufferedFileWriter_t *file,
					  uint8_t compressOutput );


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