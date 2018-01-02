#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../io/io.h"
#include "MAWs_single.h"


#ifndef INITIAL_CHAR_STACK_CAPACITY
#define INITIAL_CHAR_STACK_CAPACITY 8  // In characters. The stack can grow.
#endif
#ifndef INITIAL_MAWS_BUFFER_CAPACITY
#define INITIAL_MAWS_BUFFER_CAPACITY 1000  // In characters. The buffer can grow.
#endif
#ifndef SCORE_BUFFER_CAPACITY
#define SCORE_BUFFER_CAPACITY 50  // In characters. The buffer does not grow.
#endif


void MAWs_initialize( MAWs_callback_state_t *state, 
	                  unsigned int textLength, 
				      unsigned int minLength, 
					  unsigned int lengthHistogramMin,
					  unsigned int lengthHistogramMax,
				      unsigned char writeMAWs, 
				      unsigned char computeScores, 
					  unsigned char compressOutput,
				      char *filePath ) {
	unsigned int i, j, k;
	
	state->textLength=textLength;
	state->minLength=minLength;
	state->lengthHistogramMin=lengthHistogramMin;
	state->lengthHistogramMax=lengthHistogramMax;
	state->writeMAWs=writeMAWs;
	state->computeScores=computeScores;
	state->compressOutput=compressOutput;
	state->nMAWs=0;
	state->maxLength=0;
	state->nMaxreps=0;
	state->nMAWMaxreps=0;
	
	// Character stack
	if (state->writeMAWs!=0) {
		state->char_stack_capacity=INITIAL_CHAR_STACK_CAPACITY;  // In characters
		state->char_stack=(unsigned long *)malloc(state->char_stack_capacity>>2);  // In bytes
	}
	else {
		state->char_stack_capacity=0;
		state->char_stack=NULL;
	}
	
	// Output buffer
	state->MAWs_buffer_size=0;
	if (state->writeMAWs!=0) {
		state->MAWs_buffer_capacity=INITIAL_MAWS_BUFFER_CAPACITY;  // In characters
		state->MAWs_buffer=(char *)malloc(state->MAWs_buffer_capacity*sizeof(char));
		state->file=fopen(filePath,"a");
	}
	else {
		state->MAWs_buffer_capacity=0;
		state->MAWs_buffer=NULL;
		state->file=NULL;
	}
	
	// Scores
	state->lengthScoreCallback=NULL;
	if (state->writeMAWs!=0 && state->computeScores!=0) {
		state->leftFreqs=(unsigned int *)malloc(strlen(DNA_ALPHABET)*sizeof(unsigned int));
		state->rightFreqs=(unsigned int *)malloc(strlen(DNA_ALPHABET)*sizeof(unsigned int));
		state->scoreBuffer=(char *)malloc(SCORE_BUFFER_CAPACITY*sizeof(char));  // Arbitrary choice
		state->score_stack=(double *)malloc(state->char_stack_capacity*sizeof(double));
	}
	else {
		state->leftFreqs=NULL;
		state->rightFreqs=NULL;
		state->scoreBuffer=NULL;
		state->score_stack=NULL;
	}
	
	// Histograms
	if (state->lengthHistogramMin!=0) {
		state->lengthHistogramSize=lengthHistogramMax-lengthHistogramMin+1;
		state->lengthHistogram=(unsigned int *)malloc(state->lengthHistogramSize*sizeof(unsigned int));
		for (i=0; i<state->lengthHistogramSize; i++) state->lengthHistogram[i]=0;
	}
	else {
		state->lengthHistogramSize=0;
		state->lengthHistogram=NULL;
	}
	
	// Compressed output
	for (i=0; i<4; i++) {
		for (j=0; j<4; j++) {
			for (k=0; k<4; k++) state->compressionBuffersLength[i][j][k]=0;
		}
	}
	if (state->writeMAWs!=0 && state->computeScores==0 && state->compressOutput!=0) {
		for (i=0; i<4; i++) {
			for (j=0; j<i; j++) {
				for (k=0; k<j; k++) state->compressionBuffers[i][j][k]=(unsigned long *)malloc(BUFFER_CHUNK);
				for (k=j+1; k<4; k++) state->compressionBuffers[i][j][k]=(unsigned long *)malloc(BUFFER_CHUNK);
			}
			for (j=i+1; j<4; j++) {
				for (k=0; k<j; k++) state->compressionBuffers[i][j][k]=(unsigned long *)malloc(BUFFER_CHUNK);
				for (k=j+1; k<4; k++) state->compressionBuffers[i][j][k]=(unsigned long *)malloc(BUFFER_CHUNK);
			}
		}
		for (i=0; i<4; i++) {
			for (j=0; j<4; j++) {
				for (k=0; k<4; k++) state->compressionBuffersCapacity[i][j][k]=BUFFER_CHUNK<<3;  // In bits
			}
		}
		state->runs_stack=(unsigned char *)malloc(state->char_stack_capacity*sizeof(unsigned char));
	}
	else {
		for (i=0; i<4; i++) {
			for (j=0; j<i; j++) {
				for (k=0; k<j; k++) state->compressionBuffers[i][j][k]=NULL;
			}
		}
		for (i=0; i<4; i++) {
			for (j=0; j<4; j++) {
				for (k=0; k<4; k++) state->compressionBuffersCapacity[i][j][k]=0;
			}
		}
		state->runs_stack=NULL;
	}
}


/**
 * Prints to $state->file$ all MAWs stored in $state->compressionBuffers$.
 */
static void printCompressedMAWs(MAWs_callback_state_t *state) {
	unsigned char i, j, k, p, q;
	unsigned char cell, rem;
	unsigned int infixLength, stringLength;
	unsigned long mask;
	
	for (i=0; i<4; i++) {
		for (j=0; j<4; j++) {
			for (k=0; k<4; k++) {
				infixLength=state->compressionBuffersLength[i][j][k];
				if (infixLength==0) continue;
				stringLength=infixLength+2;
				if (state->MAWs_buffer_size+(stringLength+1+infixLength+1) > state->MAWs_buffer_capacity) {
					fwrite(state->MAWs_buffer,state->MAWs_buffer_size,sizeof(char),state->file);
					state->MAWs_buffer_size=0;
				}
				// Printing MAW
				state->MAWs_buffer[state->MAWs_buffer_size++]=DNA_ALPHABET[i];
				for (p=0; p<infixLength; p++) state->MAWs_buffer[state->MAWs_buffer_size++]=DNA_ALPHABET[j];
				state->MAWs_buffer[state->MAWs_buffer_size++]=DNA_ALPHABET[k];
				state->MAWs_buffer[state->MAWs_buffer_size++]=OUTPUT_SEPARATOR_1;
				// Printing bitmap
				p=state->compressionBuffersLength[i][j][k]-1;
				cell=p/BITS_PER_LONG; rem=p%BITS_PER_LONG;
				for (p=0; p<cell; p++) {
					mask=1L;
					for (q=0; q<BITS_PER_LONG; q++) {
						state->MAWs_buffer[state->MAWs_buffer_size++]=(state->compressionBuffers[i][j][k][p]&mask)==0?'0':'1';
						mask<<=1;
					}
				}
				mask=1L;
				for (q=0; q<=rem; q++) {
					state->MAWs_buffer[state->MAWs_buffer_size++]=(state->compressionBuffers[i][j][k][cell]&mask)==0?'0':'1';
					mask<<=1;
				}
				state->MAWs_buffer[state->MAWs_buffer_size++]=OUTPUT_SEPARATOR_2;
			}
		}
	}	
}


void MAWs_finalize(MAWs_callback_state_t *state) {
	unsigned int i, j, k;

	// Character stack
	if (state->writeMAWs!=0) free(state->char_stack);
	
	// Output buffer
	if (state->writeMAWs!=0) {
		if (state->computeScores==0 && state->compressOutput!=0) printCompressedMAWs(state);
		if (state->MAWs_buffer_size>0) {
			// Flushing the buffer one more time
			fwrite(state->MAWs_buffer,state->MAWs_buffer_size,sizeof(char),state->file);
		}
		fclose(state->file);
		free(state->MAWs_buffer);
	}
	
	// Scores
	if (state->computeScores!=0) {
		free(state->leftFreqs);
		free(state->rightFreqs);
		free(state->scoreBuffer);
		free(state->score_stack);
	}
	
	// Histograms
	if (state->lengthHistogramMin!=0) free(state->lengthHistogram);
	
	// Compressed output
	if (state->writeMAWs!=0 && state->computeScores==0 && state->compressOutput!=0) {
		for (i=0; i<4; i++) {
			for (j=0; j<i; j++) {
				for (k=0; k<j; k++) free(state->compressionBuffers[i][j][k]);
				for (k=j+1; k<4; k++) free(state->compressionBuffers[i][j][k]);
			}
			for (j=i+1; j<4; j++) {
				for (k=0; k<j; k++) free(state->compressionBuffers[i][j][k]);
				for (k=j+1; k<4; k++) free(state->compressionBuffers[i][j][k]);
			}
		}
		free(state->runs_stack);
	}
}


/**
 * Pushes to $state->char_stack$ the ID of the character of the last Weiner link, i.e. of
 * the first character of the nonempty right-maximal string described by $SLT_params$.
 * $state->char_stack$ contains numbers in $[0..3]$ represented with two bits.
 *
 * If $state->computeScores$ is nonzero, pushes to $state->score_stack$ the log of the 
 * product of character probabilities of the string described by $SLT_params$.
 * If $state->compressOutput$ is nonzero, pushes to $state->runs_stack$ a one if the
 * right-maximal string is $a^n$ for some character $a$, and pushes a zero otherwise.
 */
static void pushChar(SLT_params_t SLT_params, MAWs_callback_state_t *state) {
	const unsigned int CAPACITY = state->char_stack_capacity;
	unsigned char b, c, rem, flag;
	unsigned int bit, cell;
	const unsigned long MASK = 3L;
	double value;
	
	if (SLT_params.string_depth>CAPACITY) {
		state->char_stack_capacity+=MY_CEIL(state->char_stack_capacity*ALLOC_GROWTH_NUM,ALLOC_GROWTH_DENOM);
		state->char_stack=(unsigned long *)realloc(state->char_stack,MY_CEIL(state->char_stack_capacity<<1,8));
		if (state->computeScores!=0) state->score_stack=(double *)realloc(state->score_stack,state->char_stack_capacity*sizeof(double));
		if (state->computeScores==0 && state->compressOutput) state->runs_stack=(unsigned char *)realloc(state->runs_stack,state->char_stack_capacity*sizeof(unsigned char));
	}
	c=SLT_params.WL_char-1;
	bit=(SLT_params.string_depth-1)<<1; cell=bit/BITS_PER_LONG; rem=bit%BITS_PER_LONG;
	state->char_stack[cell]&=~(MASK<<rem);
	state->char_stack[cell]|=c<<rem;
	if (state->computeScores!=0) {
		value=LOG_DNA_ALPHABET_PROBABILITIES[SLT_params.WL_char-1];
		if (SLT_params.string_depth>1) value+=state->score_stack[SLT_params.string_depth-2];
		state->score_stack[SLT_params.string_depth-1]=value;
	}
	else if (state->compressOutput) {
		if (SLT_params.string_depth<=1) flag=1;
		else {
			if (state->runs_stack[SLT_params.string_depth-2]==0) flag=0;
			else {
				bit=(SLT_params.string_depth-2)<<1; cell=bit/BITS_PER_LONG; rem=bit%BITS_PER_LONG;
				b=(state->char_stack[cell]&(MASK<<rem))>>rem;
				flag=c==b?1:0;
			}
		}
		state->runs_stack[SLT_params.string_depth-1]=flag;
	}
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
 * extensions in the text. The string is terminated by $OUTPUT_SEPARATOR_1$.
 *
 * Remark: the procedure enlarges $state->MAWs_buffer$ if necessary.
 */
static void printMAW(SLT_params_t SLT_params, char a, char b, MAWs_callback_state_t *state) {
	const unsigned int STRING_LENGTH = SLT_params.string_depth+2;
	const unsigned int CAPACITY = state->MAWs_buffer_capacity;
	unsigned long mask;
	unsigned char rem;
	unsigned int bit;
	int cell;
	
	if (STRING_LENGTH+1>CAPACITY) {
		state->MAWs_buffer_capacity=(STRING_LENGTH+1)<<1;
		state->MAWs_buffer=(char *)realloc(state->MAWs_buffer,state->MAWs_buffer_capacity*sizeof(char));
	}
	if (state->MAWs_buffer_size+STRING_LENGTH+1 > CAPACITY) {
		fwrite(state->MAWs_buffer,state->MAWs_buffer_size,sizeof(char),state->file);
		state->MAWs_buffer_size=0;
	}
	state->MAWs_buffer[state->MAWs_buffer_size++]=a;
	bit=SLT_params.string_depth<<1; cell=bit/BITS_PER_LONG; 
	rem=bit%BITS_PER_LONG; mask=3L<<rem;
	while (cell>=0) {
		state->MAWs_buffer[state->MAWs_buffer_size++]=DNA_ALPHABET[(state->char_stack[cell]&mask)>>rem];
		if (rem==0) { cell--; rem=INITIAL_REST; mask=INITIAL_MASK; }
		else { rem-=2; mask>>=2; }
	}
	state->MAWs_buffer[state->MAWs_buffer_size++]=b;
	state->MAWs_buffer[state->MAWs_buffer_size++]=OUTPUT_SEPARATOR_1;
}


/**
 * Prints $score$ to $state->file$, followed by $OUTPUT_SEPARATOR_1$.
 * The procedure assumes both $state->scoreBuffer$ and $state->MAWs_buffer$ to be large
 * enough to contain a score.
 */
static void printScore(double score, MAWs_callback_state_t *state) {
	unsigned int i;
	unsigned int nCharacters;

	nCharacters=sprintf(state->scoreBuffer,"%g%c",score,OUTPUT_SEPARATOR_1);
	if (state->MAWs_buffer_size+nCharacters > state->MAWs_buffer_capacity) {
		fwrite(state->MAWs_buffer,state->MAWs_buffer_size,sizeof(char),state->file);
		state->MAWs_buffer_size=0;
	}
	for (i=0; i<nCharacters; i++) state->MAWs_buffer[state->MAWs_buffer_size++]=state->scoreBuffer[i];
}


/**
 * Prints $separator$ to $state->file$.
 */
static void printSeparator(MAWs_callback_state_t *state, char separator) {
	if (state->MAWs_buffer_size+1 > state->MAWs_buffer_capacity) {
		fwrite(state->MAWs_buffer,state->MAWs_buffer_size,sizeof(char),state->file);
		state->MAWs_buffer_size=0;
	}
	state->MAWs_buffer[state->MAWs_buffer_size++]=separator;
}


/**
 * Prints to $state->file$ the following scores for each MAW $W=aVb$ of length $k$:
 *
 * 1. the expected frequency of $W$, if the text is generated by a Markov chain of order 
 * at most $k-2$ (see \cite{almirantis2016optimal,brendel1986linguistics});
 * 2. an estimate of the probability of observing $W$ according to the model in (1)
 * (see \cite{qi2004whole,apostolico2008fast});
 * 3. a z-score based on (1);
 * 4. a length-based score defined by the user, if any;
 * 5. the probability of observing $W$ if the text is generated by an IID source;
 * 6. the expected frequency of $W$ according to the model in (1);
 * 7. a z-score based on (1) (see \cite{apostolico2000efficient,apostolico2003monotony}).
 *
 * @param leftCharID,rightCharID (in [0..3]) position of characters $a$ and $b$ in the 
 * alphabet.
 */
static void printScores(unsigned int leftCharID, unsigned int rightCharID, SLT_params_t SLT_params, MAWs_callback_state_t *state) {
	const unsigned int STRING_LENGTH = SLT_params.string_depth+2;
	double tmp, zScore, lengthScore;
	double expectedFrequencyMarkov, expectedFrequencyIID, probabilityMarkov, probabilityIID;

	// Expected frequency Markov
	expectedFrequencyMarkov=((double)(state->leftFreqs[leftCharID]*state->rightFreqs[rightCharID]))/SLT_params.interval_size;
	printScore(expectedFrequencyMarkov,state);
		
	// Probability Markov
	tmp=state->textLength-STRING_LENGTH+2;
	tmp=(tmp+1)/(tmp*tmp);
	probabilityMarkov=expectedFrequencyMarkov*tmp;
	printScore(probabilityMarkov,state);
	
	// Z-score Markov
	zScore=-expectedFrequencyMarkov/fmax(sqrt(expectedFrequencyMarkov),1.0);
	printScore(zScore,state);
	
	// Length-based score, if any.
	if (state->lengthScoreCallback!=NULL) {
		lengthScore=state->lengthScoreCallback(STRING_LENGTH);
		printScore(lengthScore,state);
	}
		
	// IID probability
	probabilityIID=pow(M_E,LOG_DNA_ALPHABET_PROBABILITIES[leftCharID]+state->score_stack[SLT_params.string_depth-1]+LOG_DNA_ALPHABET_PROBABILITIES[rightCharID]);
	printScore(probabilityIID,state);
	
	// Expected frequency IID
	expectedFrequencyIID=probabilityIID*(state->textLength-STRING_LENGTH+1);
	printScore(expectedFrequencyIID,state);
	
	// Z-score IID
	zScore=-expectedFrequencyIID/sqrt(expectedFrequencyIID*(1-probabilityIID));
	printScore(zScore,state);
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
	printf("Histogram of lengths [%d..%d]:\n",state->lengthHistogramMin,state->lengthHistogramMax);
	for (unsigned int i=0; i<state->lengthHistogramSize; i++) printf("%d,%d \n",state->lengthHistogramMin+i,state->lengthHistogram[i]);
}


/**
 * Stores in compressed form a MAW $a b^n c$, where $a=DNA_ALPHABET[i]$, 
 * $b=DNA_ALPHABET[j]$, $c=DNA_ALPHABET[k]$, $a \neq b$, $b \neq c$, $n \geq 1$.
 */
static void compressMAW(unsigned char i, unsigned char j, unsigned char k, unsigned int n, MAWs_callback_state_t *state) {
	unsigned int bit, cell;
	unsigned long mask;
	
	if (n>state->compressionBuffersLength[i][j][k]) {
		state->compressionBuffersLength[i][j][k]=n;
		if (n>state->compressionBuffersCapacity[i][j][k]) {
			state->compressionBuffersCapacity[i][j][k]=n<<1;  // In bits
			state->compressionBuffers[i][j][k]=(unsigned long *)realloc(state->compressionBuffers[i][j][k],MY_CEIL(n<<1,8));  // In bytes
		}
	}	
	bit=n-1; cell=bit/BITS_PER_LONG; mask=1L<<(bit%BITS_PER_LONG);
	state->compressionBuffers[i][j][k][cell]&=~mask;
	state->compressionBuffers[i][j][k][cell]|=mask;
	// The bit with ID $bit$ is set at most once during the whole traversal.
}


void MAWs_callback(const SLT_params_t SLT_params, void *intern_state) {
	unsigned char i, j;
	unsigned char found, char_mask1, char_mask2;
	MAWs_callback_state_t *state = (MAWs_callback_state_t *)(intern_state);

	if (state->writeMAWs!=0 && SLT_params.string_depth!=0) pushChar(SLT_params,state);
	if (SLT_params.nleft_extensions<2 || SLT_params.string_depth+2<state->minLength) return;
	state->nMaxreps++;
	if (state->writeMAWs!=0 && state->computeScores!=0) initLeftRightFreqs(SLT_params,state);
	char_mask1=1; found=0;
	for (i=1; i<=4; i++) {
		char_mask1<<=1;
		if ((SLT_params.left_extension_bitmap&char_mask1)==0) continue;
		char_mask2=1;
		for (j=1; j<=4; j++) {
			char_mask2<<=1;
			if ( (SLT_params.right_extension_bitmap&char_mask2)==0 ||
				 (SLT_params.left_right_extension_freqs[i][j]>0)
			   ) continue;
			state->nMAWs++;
			if (found==0) found=1;
			if (SLT_params.string_depth+2>state->maxLength) state->maxLength=SLT_params.string_depth+2;
			if (state->lengthHistogramMin>0) incrementHistogram(SLT_params,state);
			if (state->writeMAWs==0) continue;
			if ( state->compressOutput!=0 && 
			     i!=SLT_params.WL_char && j!=SLT_params.WL_char && 
			     state->runs_stack[SLT_params.string_depth-1]!=0 
			   ) compressMAW(i-1,SLT_params.WL_char-1,j-1,SLT_params.string_depth,state);
			else {
				printMAW(SLT_params,DNA_ALPHABET[i-1],DNA_ALPHABET[j-1],state);
				if (state->computeScores!=0) printScores(i-1,j-1,SLT_params,state);
				printSeparator(state,OUTPUT_SEPARATOR_2);
			}
		}
	}
	if (found!=0) state->nMAWMaxreps++;
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
					  unsigned char compressOutput,
					  char *filePath ) {
	MAWs_initialize(state,textLength,minLength,lengthHistogramMin,lengthHistogramMax,writeMRWs,computeScores,compressOutput,filePath);
	state->minFreq=minFreq;
	state->maxFreq=maxFreq;
	if (state->leftFreqs==NULL) state->leftFreqs=(unsigned int *)malloc(strlen(DNA_ALPHABET)*sizeof(unsigned int));
	if (state->rightFreqs==NULL) state->rightFreqs=(unsigned int *)malloc(strlen(DNA_ALPHABET)*sizeof(unsigned int));
}


void MRWs_finalize(MAWs_callback_state_t *state) {
	MAWs_finalize(state);
	if (state->leftFreqs!=NULL) free(state->leftFreqs);
	if (state->rightFreqs!=NULL) free(state->rightFreqs);
}


void MRWs_callback(const SLT_params_t SLT_params, void *intern_state) {
	unsigned char i, j;
	unsigned char found, char_mask1, char_mask2;
	MAWs_callback_state_t *state = (MAWs_callback_state_t *)(intern_state);

	if (state->writeMAWs!=0 && SLT_params.string_depth!=0) pushChar(SLT_params,state);
	if (SLT_params.nleft_extensions<2 || SLT_params.string_depth+2<state->minLength || SLT_params.interval_size<state->maxFreq) return;
	state->nMaxreps++;
	initLeftRightFreqs(SLT_params,state);
	char_mask1=1; found=0;
	for (i=1; i<=4; i++) {
		char_mask1<<=1;	
		if ( (SLT_params.left_extension_bitmap&char_mask1)==0 ||
			 state->leftFreqs[i-1]<state->maxFreq
		   ) continue;
		char_mask2=1;
		for (j=1; j<=4; j++) {
			char_mask2<<=1;
			if ( (SLT_params.right_extension_bitmap&char_mask2)==0 ||
				 state->rightFreqs[j-1]<state->maxFreq ||
				 (SLT_params.left_right_extension_freqs[i][j]>=state->maxFreq) ||
				 (SLT_params.left_right_extension_freqs[i][j]<state->minFreq)
			   ) continue;
			state->nMAWs++;
			if (found==0) found=1;
			if (SLT_params.string_depth+2>state->maxLength) state->maxLength=SLT_params.string_depth+2;
			if (state->lengthHistogramMin>0) incrementHistogram(SLT_params,state);
			if (state->writeMAWs==0) continue;
			if ( state->compressOutput!=0 && 
			     i!=SLT_params.WL_char && j!=SLT_params.WL_char && 
			     state->runs_stack[SLT_params.string_depth-1]!=0 
			   ) compressMAW(i-1,SLT_params.WL_char-1,j-1,SLT_params.string_depth,state);
			else {
				printMAW(SLT_params,DNA_ALPHABET[i-1],DNA_ALPHABET[j-1],state);
				if (state->computeScores!=0) printScores(i-1,j-1,SLT_params,state);
				printSeparator(state,OUTPUT_SEPARATOR_2);
			}
		}
	}
	if (found!=0) state->nMAWMaxreps++;
}