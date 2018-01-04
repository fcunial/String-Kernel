#include "MAWs_single.h"
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "../io/io.h"


#ifndef INITIAL_CHAR_STACK_CAPACITY
#define INITIAL_CHAR_STACK_CAPACITY 100  // In characters. The stack can grow.
#endif


void MAWs_initialize( MAWs_callback_state_t *state, 
	                  unsigned int textLength, 
				      unsigned int minLength, 
					  unsigned int lengthHistogramMin,
					  unsigned int lengthHistogramMax,
					  buffered_file_writer_t *file,
					  unsigned char compressOutput ) {
	unsigned int i, j, k;
	
	state->textLength=textLength;
	state->minLength=minLength;
	state->lengthHistogramMin=lengthHistogramMin;
	state->lengthHistogramMax=lengthHistogramMax;
	state->compressOutput=compressOutput;
	state->nMAWs=0;
	state->maxLength=0;
	state->nMaxreps=0;
	state->nMAWMaxreps=0;
	
	// Output buffer
	state->outputFile=file;
	
	// Character stack
	if (file!=NULL) {
		state->char_stack_capacity=INITIAL_CHAR_STACK_CAPACITY;  // In characters
		state->char_stack=(unsigned long *)malloc(state->char_stack_capacity>>2);  // In bytes
	}
	else {
		state->char_stack_capacity=0;
		state->char_stack=NULL;
	}
	
	// Scores
	state->leftFreqs=(unsigned int *)malloc(strlen(DNA_ALPHABET)*sizeof(unsigned int));
	state->rightFreqs=(unsigned int *)malloc(strlen(DNA_ALPHABET)*sizeof(unsigned int));
	state->scoreState=NULL;
	
	// Histograms
	if (state->lengthHistogramMin!=0) {
		state->lengthHistogramSize=lengthHistogramMax-lengthHistogramMin+1;
		state->lengthHistogram=(unsigned int *)calloc(state->lengthHistogramSize,sizeof(unsigned int));
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
	if (state->outputFile!=NULL && state->compressOutput!=0) {
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
		state->runs_stack=(unsigned long *)malloc(MY_CEIL(state->char_stack_capacity,8));
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
 * Prints to $state->outputFile$ all MAWs stored in $state->compressionBuffers$.
 */
static inline void printCompressedMAWs(MAWs_callback_state_t *state) {
	unsigned char i, j, k, p;
	unsigned int infixLength;
	
	for (i=0; i<4; i++) {
		for (j=0; j<4; j++) {
			for (k=0; k<4; k++) {
				infixLength=state->compressionBuffersLength[i][j][k];
				if (infixLength==0) continue;
				writeChar(DNA_ALPHABET[i],state->outputFile);
				for (p=0; p<infixLength; p++) writeChar(DNA_ALPHABET[j],state->outputFile);
				writeChar(DNA_ALPHABET[k],state->outputFile);
				writeChar(OUTPUT_SEPARATOR_1,state->outputFile);
				writeBits(state->compressionBuffers[i][j][k],state->compressionBuffersLength[i][j][k]-1,state->outputFile);
				writeChar(OUTPUT_SEPARATOR_2,state->outputFile);
			}
		}
	}	
}


void MAWs_finalize(MAWs_callback_state_t *state) {
	unsigned int i, j, k;

	// Character stack
	if (state->outputFile!=NULL) free(state->char_stack);
	
	// Output buffer
	if (state->outputFile!=NULL && state->compressOutput!=0) printCompressedMAWs(state);
	
	// Histograms
	if (state->lengthHistogramMin!=0) free(state->lengthHistogram);
	
	// Scores
	free(state->leftFreqs);
	free(state->rightFreqs);
	
	// Compressed output
	if (state->outputFile!=NULL && state->compressOutput!=0) {
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
 * If $state->compressOutput$ is nonzero, pushes to $state->runs_stack$ a one if the
 * right-maximal string is $a^n$ for some character $a$, and pushes a zero otherwise.
 */
static void pushChar(SLT_params_t SLT_params, MAWs_callback_state_t *state) {
	const unsigned int CAPACITY = state->char_stack_capacity;
	unsigned char c, flag;
	
	if (SLT_params.string_depth>CAPACITY) {
		state->char_stack_capacity+=MY_CEIL(state->char_stack_capacity*ALLOC_GROWTH_NUM,ALLOC_GROWTH_DENOM);
		state->char_stack=(unsigned long *)realloc(state->char_stack,MY_CEIL(state->char_stack_capacity<<1,8));
		if (state->compressOutput) state->runs_stack=(unsigned long *)realloc(state->runs_stack,MY_CEIL(state->char_stack_capacity,8));
	}
	c=SLT_params.WL_char-1;
	writeTwoBits(state->char_stack,SLT_params.string_depth-1,c);
	if (state->scoreState!=NULL) scorePush(c,SLT_params.string_depth,state->scoreState);
	else if (state->compressOutput) {
		if (SLT_params.string_depth<=1) flag=1;
		else {
			if (readBit(state->runs_stack,SLT_params.string_depth-2)==0) flag=0;
			else flag=c==readTwoBits(state->char_stack,SLT_params.string_depth-2)?1:0;
		}
		writeBit(state->runs_stack,SLT_params.string_depth-1,flag);
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
 * Prints to $state->outputFile$ a string $aWb$, where $W$ is the maximal repeat described
 * by $SLT_params$, and $a,b$ are characters that correspond to its left- and right-
 * extensions in the text. The string is terminated by $OUTPUT_SEPARATOR_1$.
 */
static inline void printMAW(SLT_params_t SLT_params, char a, char b, MAWs_callback_state_t *state) {
	writeChar(a,state->outputFile);
	writeTwoBitsReversed(state->char_stack,SLT_params.string_depth-1,state->outputFile,DNA_ALPHABET);
	writeChar(b,state->outputFile);
	writeChar(OUTPUT_SEPARATOR_1,state->outputFile);
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
	if (n>state->compressionBuffersLength[i][j][k]) {
		state->compressionBuffersLength[i][j][k]=n;
		if (n>state->compressionBuffersCapacity[i][j][k]) {
			state->compressionBuffersCapacity[i][j][k]=n<<1;  // In bits
			state->compressionBuffers[i][j][k]=(unsigned long *)realloc(state->compressionBuffers[i][j][k],MY_CEIL(n<<1,8));  // In bytes
		}
	}
	writeBit(state->compressionBuffers[i][j][k],n-1,1);
	// The bit with ID $bit$ is set at most once during the whole traversal.
}


void MAWs_callback(SLT_params_t SLT_params, void *intern_state) {
	unsigned char i, j;
	unsigned char found, char_mask1, char_mask2;
	MAWs_callback_state_t *state = (MAWs_callback_state_t *)(intern_state);

	if (state->outputFile!=NULL && SLT_params.string_depth!=0) pushChar(SLT_params,state);
	if (SLT_params.nleft_extensions<2 || SLT_params.string_depth+2<state->minLength) return;
	state->nMaxreps++;
	if (state->outputFile!=NULL && state->scoreState!=NULL) initLeftRightFreqs(SLT_params,state);
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
			if (state->outputFile==NULL) continue;
			if ( state->compressOutput!=0 && 
			     i!=SLT_params.WL_char && j!=SLT_params.WL_char && 
				 readBit(state->runs_stack,SLT_params.string_depth-1)!=0
			   ) compressMAW(i-1,SLT_params.WL_char-1,j-1,SLT_params.string_depth,state);
			else {
				printMAW(SLT_params,DNA_ALPHABET[i-1],DNA_ALPHABET[j-1],state);
				if (state->scoreState!=NULL) {
					scoreCallback(i-1,j-1,state->leftFreqs[i-1],state->rightFreqs[j-1],state->textLength,&SLT_params,state->scoreState);
					scorePrint(state->scoreState,state->outputFile);
				}
				writeChar(OUTPUT_SEPARATOR_2,state->outputFile);
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
					  buffered_file_writer_t *file,
					  unsigned char compressOutput ) {
	MAWs_initialize(state,textLength,minLength,lengthHistogramMin,lengthHistogramMax,file,compressOutput);
	state->minFreq=minFreq;
	state->maxFreq=maxFreq;
}


void MRWs_finalize(MAWs_callback_state_t *state) {
	MAWs_finalize(state);
}


void MRWs_callback(SLT_params_t SLT_params, void *intern_state) {
	unsigned char i, j;
	unsigned char found, char_mask1, char_mask2;
	MAWs_callback_state_t *state = (MAWs_callback_state_t *)(intern_state);

	if (state->outputFile!=NULL && SLT_params.string_depth!=0) pushChar(SLT_params,state);
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
			if (state->outputFile==NULL) continue;
			if ( state->compressOutput!=0 && 
			     i!=SLT_params.WL_char && j!=SLT_params.WL_char && 
				 readBit(state->runs_stack,SLT_params.string_depth-1)!=0
			   ) compressMAW(i-1,SLT_params.WL_char-1,j-1,SLT_params.string_depth,state);
			else {
				printMAW(SLT_params,DNA_ALPHABET[i-1],DNA_ALPHABET[j-1],state);
				if (state->scoreState!=NULL) {
					scoreCallback(i-1,j-1,state->leftFreqs[i-1],state->rightFreqs[j-1],state->textLength,&SLT_params,state->scoreState);
					scorePrint(state->scoreState,state->outputFile);
				}
				writeChar(OUTPUT_SEPARATOR_2,state->outputFile);
			}
		}
	}
	if (found!=0) state->nMAWMaxreps++;
}