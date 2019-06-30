/**
 * 
 *
 * @author Fabio Cunial
 */
#ifndef scores_h
#define scores_h

#include "./iterator/SLT_single_string.h"
#include "./io/bufferedFileWriter.h"


typedef struct {
	double *scores;  // List of scores
	double *score_stack;
	unsigned int score_stack_capacity;  // In elements
	char *scoreBuffer;  // String representation of a score
} score_state_t;


void scoreInitialize(score_state_t *scoreState);


void scoreFinalize(score_state_t *scoreState);


void scoreClone(score_state_t *from, score_state_t *to);


/**
 * Called for each MAW $W=aVb$ where $V$ is described by $SLT_params$.
 *
 * @param leftCharID,rightCharID (in [0..3]) position of characters $a$ and $b$ in the 
 * alphabet;
 * @param leftFreq,rightFreq frequency of $aV$ and $Vb$ in the text.
 */
void scoreCallback(unsigned int leftCharID, unsigned int rightCharID, unsigned int leftFreq, unsigned int rightFreq, unsigned int textLength, SLT_params_t *SLT_params, score_state_t *scoreState);


/**
 * Called whenever a character is pushed on the character stack.
 *
 * @param charID position of the character in the alphabet;
 * @param stringDepth depth of the stack after the character has been pushed.
 */
void scorePush(unsigned char charID, unsigned int stringDepth, score_state_t *scoreState);


/**
 * Prints all scores to $file$.
 */
void scorePrint(score_state_t *scoreState, buffered_file_writer_t *file);


/**
 * Returns a number different from 0 iff the scores in $scoreState$ have an
 * implementation-defined property.
 */
char scoreSelect(score_state_t *scoreState);


#endif