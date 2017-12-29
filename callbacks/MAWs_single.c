#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../io/io.h"
#include "MAWs_single.h"


void MAWs_initialize( MAWs_callback_state_t *state, 
	                  unsigned int textLength, 
				      unsigned int minLength, 
					  unsigned int lengthHistogramMin,
					  unsigned int lengthHistogramMax,
				      unsigned char writeMAWs, 
				      unsigned char computeScores, 
				      char *filePath ) {
	unsigned int i;
	
	state->textLength=textLength;
	state->minLength=minLength;
	state->nMAWs=0;
	
	// Character stack
	state->char_stack_capacity=BUFFER_CHUNK;  // Arbitrary choice
	state->char_stack=(unsigned char *)malloc(state->char_stack_capacity*sizeof(unsigned char));
	
	// Output buffer
	state->writeMAWs=writeMAWs;
	if (writeMAWs==0) return;
	state->MAWs_buffer_capacity=BUFFER_CHUNK;  // Arbitrary choice
	state->MAWs_buffer=(char *)malloc(state->MAWs_buffer_capacity*sizeof(char));
	state->MAWs_buffer_size=0;
	state->file=fopen(filePath,"a");
	
	// Scores
	state->computeScores=computeScores;
	if (computeScores==0) return;
	state->leftFreqs=(unsigned int *)malloc(strlen(DNA_ALPHABET)*sizeof(unsigned int));
	state->rightFreqs=(unsigned int *)malloc(strlen(DNA_ALPHABET)*sizeof(unsigned int));
	state->scoreBuffer=(char *)malloc(50*sizeof(char));  // Arbitrary choice
	
	// Histograms
	state->lengthHistogramMin=lengthHistogramMin;
	state->lengthHistogramMax=lengthHistogramMax;
	if (lengthHistogramMin==0) return;
	state->lengthHistogramSize=lengthHistogramMax-lengthHistogramMin+1;
	state->lengthHistogram=(unsigned int *)malloc(state->lengthHistogramSize*sizeof(unsigned int));
	for (i=0; i<state->lengthHistogramSize; i++) state->lengthHistogram[i]=0;
}


void MAWs_finalize(MAWs_callback_state_t *state) {
	// Character stack
	free(state->char_stack);
	
	// Output buffer
	if (state->writeMAWs==0) return;
	if (state->MAWs_buffer_size>0) {
		// Flushing the buffer one more time
		fwrite(state->MAWs_buffer,state->MAWs_buffer_size,sizeof(char),state->file);
	}
	fclose(state->file);
	free(state->MAWs_buffer);
	
	// Scores
	if (state->computeScores==0) return;
	free(state->leftFreqs);
	free(state->rightFreqs);
	free(state->scoreBuffer);
	
	// Histograms
	if (state->lengthHistogramMin==-1) return;
	free(state->lengthHistogram);
}


/**
 * Pushes to $state->char_stack$ the label of the last Weiner link, i.e. the first 
 * character of the right-maximal string described by $SLT_params$.
 */
static void pushChar(SLT_params_t SLT_params, MAWs_callback_state_t *state) {
	const unsigned int CAPACITY = state->char_stack_capacity;
	
	if (SLT_params.string_depth>CAPACITY) {
		state->char_stack_capacity+=(CAPACITY*ALLOC_GROWTH_NUM)/ALLOC_GROWTH_DENOM;
		state->char_stack=(unsigned char *)realloc(state->char_stack,state->char_stack_capacity*sizeof(unsigned char));
	}
	state->char_stack[SLT_params.string_depth-1]=SLT_params.WL_char;
}


/**
 * Sets just the cells of $state->{left,right}Freqs$ that correspond to ones in 
 * $SLT_params.{left,right}_extension_bitmap$.
 */
static void initLeftRightFreqs(SLT_params_t SLT_params, MAWs_callback_state_t *state) {
	unsigned int i, j;
	unsigned char char_mask;
	unsigned int frequency;
	
	char_mask=1;
	for (i=1; i<=4; i++) {
		char_mask<<=1;
		if (!(SLT_params.left_extension_bitmap & char_mask)) continue;
		frequency=0;
		for (j=0; j<=5; j++) frequency+=SLT_params.left_right_extension_freqs[i][j];
		state->leftFreqs[i-1]=frequency;
	}
	char_mask=1;
	for (j=1; j<=4; j++) {
		char_mask<<=1;
		if (!(SLT_params.right_extension_bitmap & char_mask)) continue;
		frequency=0;
		for (i=0; i<=5; i++) frequency+=SLT_params.left_right_extension_freqs[i][j];
		state->rightFreqs[j-1]=frequency;
	}
}


/**
 * Prints to $state->file$ a string $aWb$, where $W$ is the maximal repeat described by
 * $SLT_params$, and $a,b$ are characters that correspond to its left- and right-
 * extensions in the text. The string is terminated by $OUTPUT_SEPARATOR$.
 */
static void printMAW(SLT_params_t SLT_params, unsigned char a, unsigned char b, MAWs_callback_state_t *state) {
	const unsigned int STRING_LENGTH = SLT_params.string_depth+2;
	unsigned int i;
	
	if (state->MAWs_buffer_size+STRING_LENGTH+1 > state->MAWs_buffer_capacity) {
		fwrite(state->MAWs_buffer,state->MAWs_buffer_size,sizeof(char),state->file);
		state->MAWs_buffer_size=0;
	}
	state->MAWs_buffer[state->MAWs_buffer_size++]=a;
	for (i=0; i<SLT_params.string_depth; i++) state->MAWs_buffer[state->MAWs_buffer_size++]=DNA_ALPHABET[state->char_stack[SLT_params.string_depth-1-i]-1];
	state->MAWs_buffer[state->MAWs_buffer_size++]=b;
	state->MAWs_buffer[state->MAWs_buffer_size++]=OUTPUT_SEPARATOR;
}


/**
 * Prints the first $nCharacters$ characters of $state->scoreBuffer$ to $state.file$,
 * followed by $OUTPUT_SEPARATOR$.
 */
static void printScore(double score, MAWs_callback_state_t *state) {
	unsigned int i;
	unsigned int nCharacters;

	nCharacters=sprintf(state->scoreBuffer,"%g %c",score,OUTPUT_SEPARATOR);
	if (state->MAWs_buffer_size+nCharacters > state->MAWs_buffer_capacity) {
		fwrite(state->MAWs_buffer,state->MAWs_buffer_size,sizeof(char),state->file);
		state->MAWs_buffer_size=0;
	}
	for (i=0; i<nCharacters; i++) state->MAWs_buffer[state->MAWs_buffer_size++]=state->scoreBuffer[i];
}


/**
 * Prints to $state->file$ the following scores for each MAW $W=aVb$ of length $k$:
 *
 * 1. the expected frequency of $W$, if the text is generated by a Markov chain of order 
 * at most $k-2$ (see \cite{almirantis2016optimal,brendel1986linguistics});
 * 2. an estimate of the probability of observing $W$ according to the model in (1)
 * (see \cite{qi2004whole,apostolico2008fast});
 * 3. a z-score based on (1);
 * 4. the length-based score used in \cite{crochemore2016linear}.
 *
 * @param leftCharID,rightCharID (in [0..3]) position of characters $a$ and $b$ in the 
 * alphabet.
 */
static void printScores(unsigned int leftCharID, unsigned int rightCharID, SLT_params_t SLT_params, MAWs_callback_state_t *state) {
	const unsigned int STRING_LENGTH = SLT_params.string_depth+2;
	double tmp, expectedFrequency, probability, zScore, lengthScore;
	
	// Expected frequency Markov
	expectedFrequency=((double)(state->leftFreqs[leftCharID]*state->rightFreqs[rightCharID]))/SLT_params.interval_size;
	printScore(expectedFrequency,state);
	
	// Probability Markov
	tmp=state->textLength-STRING_LENGTH+2;
	tmp=(tmp+1)/(tmp*tmp);
	probability=expectedFrequency*tmp;
	printScore(probability,state);
	
	// Z-score Markov
	zScore=-expectedFrequency/fmax(sqrt(expectedFrequency),1.0);
	printScore(zScore,state);
	
	// Length-based
	lengthScore=1.0/(STRING_LENGTH*STRING_LENGTH);
	printScore(lengthScore,state);
}


static void incrementHistogram(SLT_params_t SLT_params, MAWs_callback_state_t *state) {
	unsigned int length, position;
	
	length=SLT_params.string_depth+2;
	if (length>=state->lengthHistogramMax) position=state->lengthHistogramSize-1;
	else if (length<=state->lengthHistogramMin) position=0;
	else position=length-state->lengthHistogramMin;
	state->lengthHistogram[position]++;
}


inline void printLengthHistogram(MAWs_callback_state_t *state) {
	printf("Histogram of lengths [%d..%d]:",state->lengthHistogramMin,state->lengthHistogramMax);
	for (unsigned int i=0; i<state->lengthHistogramSize; i++) printf("%d,%d \n",state->lengthHistogramMin+i,state->lengthHistogram[i]);
}


void MAWs_callback(const SLT_params_t SLT_params, void *intern_state) {
	unsigned char i, j;
	unsigned char char_mask1, char_mask2;
	MAWs_callback_state_t *state = (MAWs_callback_state_t *)(intern_state);

	if (state->writeMAWs!=0 && SLT_params.string_depth!=0) pushChar(SLT_params,state);
	if (SLT_params.nleft_extensions<2 || SLT_params.string_depth+2<state->minLength) return;
	if (state->writeMAWs!=0 && state->computeScores!=0) initLeftRightFreqs(SLT_params,state);
	char_mask1=1;
	for (i=1; i<=4; i++) {
		char_mask1<<=1;
		if (!(SLT_params.left_extension_bitmap & char_mask1)) continue;
		char_mask2=1;
		for (j=1; j<=4; j++) {
			char_mask2<<=1;
			if ( !(SLT_params.right_extension_bitmap & char_mask2) ||
				 (SLT_params.left_right_extension_freqs[i][j]>0)
			   ) continue;
			state->nMAWs++;
			if (state->lengthHistogramMin>0) incrementHistogram(SLT_params,state);
			if (state->writeMAWs==0) continue;
			printMAW(SLT_params,DNA_ALPHABET[i-1],DNA_ALPHABET[j-1],state);
			if (state->computeScores!=0) printScores(i-1,j-1,SLT_params,state);
		}
	}
}


void MRWs_initialize( MAWs_callback_state_t *state,
			    	  unsigned int textLength, 
					  unsigned int minLength, 
					  unsigned int minFreq, 
					  unsigned int maxFreq, 
					  unsigned int lengthHistogramMin,
					  unsigned int lengthHistogramMax,
					  unsigned char writeMRWs, 
					  unsigned char computeScores, 
					  char *filePath ) {
	MAWs_initialize(state,textLength,minLength,lengthHistogramMin,lengthHistogramMax,writeMRWs,computeScores,filePath);
	state->minFreq=minFreq;
	state->maxFreq=maxFreq;
}


void MRWs_finalize(MAWs_callback_state_t *state) {
	MAWs_finalize(state);
}


void MRWs_callback(const SLT_params_t SLT_params, void *intern_state) {
	unsigned char i, j;
	unsigned char char_mask1, char_mask2;
	MAWs_callback_state_t *state = (MAWs_callback_state_t *)(intern_state);

	if (state->writeMAWs!=0 && SLT_params.string_depth!=0) pushChar(SLT_params,state);
	if (SLT_params.nleft_extensions<2 || SLT_params.string_depth+2<state->minLength || SLT_params.interval_size<state->maxFreq) return;
	if (state->writeMAWs!=0 && state->computeScores!=0) initLeftRightFreqs(SLT_params,state);
	char_mask1=1;
	for (i=1; i<=4; i++) {
		char_mask1<<=1;
		if ( !(SLT_params.left_extension_bitmap & char_mask1) ||
			 state->leftFreqs[i-1]<state->maxFreq
		   ) continue;
		char_mask2=1;
		for (j=1; j<=4; j++) {
			char_mask2<<=1;
			if ( !(SLT_params.right_extension_bitmap & char_mask2) ||
				 state->rightFreqs[j-1]<state->maxFreq ||
				 (SLT_params.left_right_extension_freqs[i][j]>=state->maxFreq) ||
				 (SLT_params.left_right_extension_freqs[i][j]<state->minFreq)
			   ) continue;
			state->nMAWs++;
			if (state->lengthHistogramMin>0) incrementHistogram(SLT_params,state);
			if (state->writeMAWs==0) continue;
			printMAW(SLT_params,DNA_ALPHABET[i-1],DNA_ALPHABET[j-1],state);
			if (state->computeScores!=0) printScores(i-1,j-1,SLT_params,state);
		}
	}
}