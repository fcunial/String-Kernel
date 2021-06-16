/**
 * @author Fabio Cunial
 */
#include "scores.h"
#include "./io/io.h"
#include "./io/bits.h"
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <string.h>

#ifndef SCORE_BUFFER_CAPACITY
#define SCORE_BUFFER_CAPACITY 50  // In characters. The buffer does not grow.
#endif
#ifndef N_SCORES
#define N_SCORES 8
#endif
#ifndef LENGTH_SCORE_EPSILON
#define LENGTH_SCORE_EPSILON 0.1
#endif
#ifndef INITIAL_SCORE_STACK_CAPACITY
#define INITIAL_SCORE_STACK_CAPACITY 128  // In doubles.
#endif


/**
 * For $scoreSelect$. Values are assumed to be set elsewhere.
 */
unsigned char SELECTED_SCORE;
double SELECTED_SCORE_THRESHOLD;


inline void scoreInitialize(ScoreState_t *scoreState, double *dnaProbabilities, double *logDnaProbabilities) {
	scoreState->scores=(double *)calloc(N_SCORES,sizeof(double));
	scoreState->scoreStackCapacity=INITIAL_SCORE_STACK_CAPACITY;
	scoreState->scoreStack=(double *)malloc(scoreState->scoreStackCapacity*sizeof(double));
	scoreState->scoreBuffer=(char *)malloc(SCORE_BUFFER_CAPACITY*sizeof(char));
	scoreState->dnaProbabilities=dnaProbabilities;
	scoreState->logDnaProbabilities=logDnaProbabilities;
}


inline void scoreFinalize(ScoreState_t *scoreState) {
	free(scoreState->scores);
	free(scoreState->scoreStack);
	free(scoreState->scoreBuffer);
}


inline void scoreClone(ScoreState_t *from, ScoreState_t *to) {
	to->scores=(double *)malloc(N_SCORES*sizeof(double));
	to->scoreStackCapacity=from->scoreStackCapacity;
	to->scoreStack=(double *)malloc(to->scoreStackCapacity*sizeof(double));
	memcpy(to->scoreStack,from->scoreStack,to->scoreStackCapacity*sizeof(double));
	to->scoreBuffer=(char *)malloc(SCORE_BUFFER_CAPACITY*sizeof(char));
	to->dnaProbabilities=from->dnaProbabilities;
	to->logDnaProbabilities=from->logDnaProbabilities;
}


/**
 * Length score used in \cite{crochemore2016linear}.     
 */
static inline double lengthScore1(uint64_t length) {
	return 1.0/(length*length);
}


/**
 * Length score used in \cite{smola2003fast}.     
 */
static inline double lengthScore2(uint64_t length) {
	return pow(LENGTH_SCORE_EPSILON,length);
}


/**
 * Stores in  $scoreState->scores$ the following scores for a MAW $W=aVb$:
 * 
 * 0. the expected frequency of $W$ if the text is generated by an IID source;
 * 1. the probability of observing $W$ according to the model in (1);
 * 2. a z-score based on (1) (see \cite{apostolico2000efficient,apostolico2003monotony});
 *
 * 3. the expected frequency of $W$, if the text is generated by a Markov chain of order 
 * at most $|W|-2$ (see \cite{almirantis2016optimal,brendel1986linguistics});
 * 4. an estimate of the probability of observing $W$ according to the model in (3)
 * (see \cite{qi2004whole,apostolico2008fast});
 * 5. a z-score based on (4);
 *
 * 6. the length-based score defined by $lengthScore1()$;
 * 7. the length-based score defined by $lengthScore2()$.
 */
inline void scoreCallback(uint8_t leftCharID, uint8_t rightCharID, uint64_t leftFreq, uint64_t rightFreq, uint64_t textLength, RightMaximalString_t *RightMaximalString, ScoreState_t *scoreState) {
	const uint64_t STRING_LENGTH = RightMaximalString->length+2;
	register double tmp;
	register double expectedFrequencyIID, probabilityIID, zScoreIID;
	register double expectedFrequencyMarkov, probabilityMarkov, zScoreMarkov;
	register double ls1, ls2;

	// IID
	probabilityIID=pow(M_E,scoreState->logDnaProbabilities[leftCharID]+scoreState->scoreStack[RightMaximalString->length-1]+scoreState->logDnaProbabilities[rightCharID]);
	expectedFrequencyIID=probabilityIID*(textLength-STRING_LENGTH+1);
	zScoreIID=-expectedFrequencyIID/sqrt(expectedFrequencyIID*(1-probabilityIID));

	// Markov
	expectedFrequencyMarkov=((double)(leftFreq*rightFreq))/RightMaximalString->frequency;
	tmp=textLength-STRING_LENGTH+2;
	tmp=(tmp+1)/(tmp*tmp);
	probabilityMarkov=expectedFrequencyMarkov*tmp;
	zScoreMarkov=-expectedFrequencyMarkov/fmax(sqrt(expectedFrequencyMarkov),1.0);

	// Length
	ls1=lengthScore1(STRING_LENGTH);
	ls2=lengthScore2(STRING_LENGTH);
	
	// Output
	scoreState->scores[0]=expectedFrequencyIID;
	scoreState->scores[1]=probabilityIID;
	scoreState->scores[2]=zScoreIID;
	scoreState->scores[3]=expectedFrequencyMarkov;
	scoreState->scores[4]=probabilityMarkov;
	scoreState->scores[5]=zScoreMarkov;
	scoreState->scores[6]=ls1;
	scoreState->scores[7]=ls2;
}


/**
 * Updates just the log of the product of character probabilities of the string in the 
 * stack.
 */
inline void scorePush(uint8_t charID, uint64_t stringDepth, ScoreState_t *scoreState) {
	const uint64_t CAPACITY = scoreState->scoreStackCapacity;
	register double probabilityIID;
	
	if (stringDepth>CAPACITY) {
		scoreState->scoreStackCapacity+=MY_CEIL(scoreState->scoreStackCapacity*ALLOC_GROWTH_NUM,ALLOC_GROWTH_DENOM);
		scoreState->scoreStack=(double *)realloc(scoreState->scoreStack,scoreState->scoreStackCapacity*sizeof(double));
	}
	probabilityIID=scoreState->logDnaProbabilities[charID];
	if (stringDepth>1) probabilityIID+=scoreState->scoreStack[stringDepth-2];
	scoreState->scoreStack[stringDepth-1]=probabilityIID;
}


/**
 * Scores are separated by $OUTPUT_SEPARATOR_1$.
 * The procedure assumes $scoreState->scoreBuffer$ to be large enough to contain a score.
 */
inline void scorePrint(ScoreState_t *scoreState, BufferedFileWriter_t *file) {
	register uint8_t i, nCharacters;

	nCharacters=sprintf(scoreState->scoreBuffer,"%c",OUTPUT_SEPARATOR_1);
	writeChars(scoreState->scoreBuffer,nCharacters-1,file);
	for (i=0; i<N_SCORES; i++) {
		nCharacters=sprintf(scoreState->scoreBuffer,"%g%c",scoreState->scores[i],OUTPUT_SEPARATOR_1);
		writeChars(scoreState->scoreBuffer,nCharacters-1,file);
	}
}


/**
 * Returns 1 iff the absolute value of $scores[SELECTED_SCORE]$ is at least 
 * $SELECTED_SCORE_THRESHOLD$.
 */
inline uint8_t scoreSelect(ScoreState_t *scoreState) {
	return fabs(scoreState->scores[SELECTED_SCORE])>=SELECTED_SCORE_THRESHOLD?1:0;
}