/**
 * @author Fabio Cunial
 */
#include "MAWs_single.h"
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "../io/io.h"
#include "../io/bits.h"


#ifndef INITIAL_CHAR_STACK_CAPACITY
#define INITIAL_CHAR_STACK_CAPACITY 128  // In characters. The stack can grow.
#endif
#ifndef MY_CEIL
#define MY_CEIL(N,D) (1+((N)-1)/(D))  // ceil(N/D) where N and D are integers.
#endif


static void initCompressedOutput(MAWs_callback_state_t *state) {
	uint8_t i, j, k;
	
	for (i=0; i<4; i++) {
		for (j=0; j<4; j++) {
			for (k=0; k<4; k++) state->compressionBuffersLength[i][j][k]=0;
		}
	}
	if (state->outputFile!=NULL && state->compressOutput!=0) {
		for (i=0; i<4; i++) {
			for (j=0; j<4; j++) {
				for (k=0; k<4; k++) state->compressionBuffersCapacity[i][j][k]=BUFFER_CHUNK<<3;  // In bits
			}
		}
		for (i=0; i<4; i++) {
			for (j=0; j<i; j++) {
				for (k=0; k<j; k++) state->compressionBuffers[i][j][k]=(uint64_t *)malloc(MY_CEIL(BUFFER_CHUNK,BYTES_PER_LONG));
				for (k=j+1; k<4; k++) state->compressionBuffers[i][j][k]=(uint64_t *)malloc(MY_CEIL(BUFFER_CHUNK,BYTES_PER_LONG));
			}
			for (j=i+1; j<4; j++) {
				for (k=0; k<j; k++) state->compressionBuffers[i][j][k]=(uint64_t *)malloc(MY_CEIL(BUFFER_CHUNK,BYTES_PER_LONG));
				for (k=j+1; k<4; k++) state->compressionBuffers[i][j][k]=(uint64_t *)malloc(MY_CEIL(BUFFER_CHUNK,BYTES_PER_LONG));
			}
		}
	}
}


void MAWs_initialize( MAWs_callback_state_t *state, 
	                  uint64_t textLength, 
				      uint64_t minLength, 
					  uint64_t lengthHistogramMin,
					  uint64_t lengthHistogramMax,
					  BufferedFileWriter_t *file,
					  uint8_t compressOutput ) {	
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
		state->char_stack=(uint64_t *)malloc(MY_CEIL(state->char_stack_capacity<<1,BITS_PER_LONG)*BYTES_PER_LONG);  // In bytes
	}
	else {
		state->char_stack_capacity=0;
		state->char_stack=NULL;
	}
	
	// Scores
	state->leftFreqs=(uint64_t *)malloc(strlen(DNA_ALPHABET)*sizeof(uint64_t));
	state->rightFreqs=(uint64_t *)malloc(strlen(DNA_ALPHABET)*sizeof(uint64_t));
	state->scoreState=NULL;
	
	// Histograms
	if (state->lengthHistogramMin!=0) {
		state->lengthHistogramSize=lengthHistogramMax-lengthHistogramMin+1;
		state->lengthHistogram=(uint64_t *)calloc(state->lengthHistogramSize,sizeof(uint64_t));
	}
	else {
		state->lengthHistogramSize=0;
		state->lengthHistogram=NULL;
	}
	
	// Compressed output
	initCompressedOutput(state);
	if (state->outputFile!=NULL && state->compressOutput!=0) state->runs_stack=(uint64_t *)malloc(MY_CEIL(state->char_stack_capacity,8));
	else state->runs_stack=NULL;
}


static void mergeCompressedOutput(MAWs_callback_state_t *from, MAWs_callback_state_t *to) {
	uint8_t i, j, k;
	uint64_t p, nBits, nLongs;
	uint64_t *tmp;
	
	for (i=0; i<4; i++) {
		for (j=0; j<i; j++) {
			for (k=0; k<j; k++) {
				nBits=to->compressionBuffersLength[i][j][k];
				if (from->compressionBuffersLength[i][j][k]>nBits) nBits=from->compressionBuffersLength[i][j][k];
				if (nBits==0) continue;
				tmp=(uint64_t *)calloc(MY_CEIL(nBits,8),1);
				if (from->compressionBuffersLength[i][j][k]!=0) {
					nLongs=MY_CEIL(from->compressionBuffersLength[i][j][k],BITS_PER_LONG);
					for (p=0; p<nLongs; p++) tmp[p]|=from->compressionBuffers[i][j][k][p];
				}
				if (to->compressionBuffersLength[i][j][k]!=0) {
					nLongs=MY_CEIL(to->compressionBuffersLength[i][j][k],BITS_PER_LONG);
					for (p=0; p<nLongs; p++) tmp[p]|=to->compressionBuffers[i][j][k][p];
				}
				to->compressionBuffersLength[i][j][k]=nBits;
				to->compressionBuffers[i][j][k]=tmp;
			}
			for (k=j+1; k<4; k++) {
				nBits=to->compressionBuffersLength[i][j][k];
				if (from->compressionBuffersLength[i][j][k]>nBits) nBits=from->compressionBuffersLength[i][j][k];
				if (nBits==0) continue;
				tmp=(uint64_t *)calloc(MY_CEIL(nBits,8),1);
				if (from->compressionBuffersLength[i][j][k]!=0) {
					nLongs=MY_CEIL(from->compressionBuffersLength[i][j][k],BITS_PER_LONG);
					for (p=0; p<nLongs; p++) tmp[p]|=from->compressionBuffers[i][j][k][p];
				}
				if (to->compressionBuffersLength[i][j][k]!=0) {
					nLongs=MY_CEIL(to->compressionBuffersLength[i][j][k],BITS_PER_LONG);
					for (p=0; p<nLongs; p++) tmp[p]|=to->compressionBuffers[i][j][k][p];
				}
				to->compressionBuffersLength[i][j][k]=nBits;
				to->compressionBuffers[i][j][k]=tmp;
			}
		}
		for (j=i+1; j<4; j++) {
			for (k=0; k<j; k++) {
				nBits=to->compressionBuffersLength[i][j][k];
				if (from->compressionBuffersLength[i][j][k]>nBits) nBits=from->compressionBuffersLength[i][j][k];
				if (nBits==0) continue;
				tmp=(uint64_t *)calloc(MY_CEIL(nBits,8),1);
				if (from->compressionBuffersLength[i][j][k]!=0) {
					nLongs=MY_CEIL(from->compressionBuffersLength[i][j][k],BITS_PER_LONG);
					for (p=0; p<nLongs; p++) tmp[p]|=from->compressionBuffers[i][j][k][p];
				}
				if (to->compressionBuffersLength[i][j][k]!=0) {
					nLongs=MY_CEIL(to->compressionBuffersLength[i][j][k],BITS_PER_LONG);
					for (p=0; p<nLongs; p++) tmp[p]|=to->compressionBuffers[i][j][k][p];
				}
				to->compressionBuffersLength[i][j][k]=nBits;
				to->compressionBuffers[i][j][k]=tmp;
			}
			for (k=j+1; k<4; k++) {
				nBits=to->compressionBuffersLength[i][j][k];
				if (from->compressionBuffersLength[i][j][k]>nBits) nBits=from->compressionBuffersLength[i][j][k];
				if (nBits==0) continue;
				tmp=(uint64_t *)calloc(MY_CEIL(nBits,8),1);
				if (from->compressionBuffersLength[i][j][k]!=0) {
					nLongs=MY_CEIL(from->compressionBuffersLength[i][j][k],BITS_PER_LONG);
					for (p=0; p<nLongs; p++) tmp[p]|=from->compressionBuffers[i][j][k][p];
				}
				if (to->compressionBuffersLength[i][j][k]!=0) {
					nLongs=MY_CEIL(to->compressionBuffersLength[i][j][k],BITS_PER_LONG);
					for (p=0; p<nLongs; p++) tmp[p]|=to->compressionBuffers[i][j][k][p];
				}
				to->compressionBuffersLength[i][j][k]=nBits;
				to->compressionBuffers[i][j][k]=tmp;
			}
		}
	}
}


void cloneMAWState(MAWs_callback_state_t *from, MAWs_callback_state_t *to, char *pathPrefix, char id, char *cloneBuffer) {
	uint64_t nBytes;
	
	to->textLength=from->textLength;
	to->minLength=from->minLength;
	to->lengthHistogramMin=from->lengthHistogramMin;
	to->lengthHistogramMax=from->lengthHistogramMax;
	to->compressOutput=from->compressOutput;
	to->nMAWs=0;
	to->maxLength=0;
	to->nMaxreps=0;
	to->nMAWMaxreps=0;
	to->minFreq=from->minFreq;
	to->maxFreq=from->maxFreq;
	
	// Output buffer
	if (from->outputFile!=NULL) {
		to->outputFile=(BufferedFileWriter_t *)malloc(sizeof(BufferedFileWriter_t));
		sprintf(cloneBuffer,"%s-%d.out",pathPrefix,id);
		initializeBufferedFileWriter(to->outputFile,cloneBuffer);
	}
	
	// Character stack
	if (from->char_stack!=NULL) {
		to->char_stack_capacity=from->char_stack_capacity;
		nBytes=MY_CEIL(to->char_stack_capacity<<1,8);
		to->char_stack=(uint64_t *)malloc(nBytes);
		memcpy(to->char_stack,from->char_stack,nBytes);
	}
	
	// Scores
	to->leftFreqs=(uint64_t *)malloc(strlen(DNA_ALPHABET)*sizeof(uint64_t));
	to->rightFreqs=(uint64_t *)malloc(strlen(DNA_ALPHABET)*sizeof(uint64_t));
	if (from->scoreState!=NULL) {
		to->scoreState=(ScoreState_t *)malloc(sizeof(ScoreState_t));
		scoreClone(from->scoreState,to->scoreState);
	}
	
	// Histograms
	if (from->lengthHistogramMin!=0) {
		to->lengthHistogramMin=from->lengthHistogramMin;
		to->lengthHistogramMax=from->lengthHistogramMax;
		to->lengthHistogramSize=from->lengthHistogramSize;
		to->lengthHistogram=(uint64_t *)calloc(to->lengthHistogramSize,sizeof(uint64_t));
	}
	else {
		to->lengthHistogramMin=0;
		to->lengthHistogramMax=0;
		to->lengthHistogramSize=0;
		to->lengthHistogram=NULL;
	}
	
	// Compressed output
	initCompressedOutput(to);
	if (to->outputFile!=NULL && to->compressOutput!=0) {
		nBytes=MY_CEIL(to->char_stack_capacity,8);
		to->runs_stack=(uint64_t *)malloc(nBytes);
		memcpy(to->runs_stack,from->runs_stack,nBytes);
	}
	else to->runs_stack=NULL;
}


void mergeMAWState(MAWs_callback_state_t *from, MAWs_callback_state_t *to) {
	uint64_t i;
	
	to->nMAWs+=from->nMAWs;
	to->maxLength=from->maxLength>to->maxLength?from->maxLength:to->maxLength;
	to->nMaxreps+=from->nMaxreps;
	to->nMAWMaxreps+=from->nMAWMaxreps;
	
	// Histograms, assumed to be of the same length.
	if (from->lengthHistogramMin!=0) {
		for (i=0; i<from->lengthHistogramSize; i++) to->lengthHistogram[i]+=from->lengthHistogram[i];
	}
	
	// Compressed output
	if (from->outputFile!=NULL && from->compressOutput!=0) mergeCompressedOutput(from,to);
}


/**
 * Prints to $state->outputFile$ all MAW encodings stored in $state->compressionBuffers$.
 *
 * Remark: the last bit of a compressed buffer is not printed, since it is always one.
 * If a bitvector has just its last bit to one, it is not printed.
 */
static void printCompressedMAWs(MAWs_callback_state_t *state) {
	uint8_t i, j, k;
	uint64_t p, infixLength;
	
	for (i=0; i<4; i++) {
		for (j=0; j<4; j++) {
			for (k=0; k<4; k++) {
				infixLength=state->compressionBuffersLength[i][j][k];
				if (infixLength==0) continue;
				writeChar(DNA_ALPHABET[i],state->outputFile);
				for (p=0; p<infixLength; p++) writeChar(DNA_ALPHABET[j],state->outputFile);
				writeChar(DNA_ALPHABET[k],state->outputFile);
				writeChar(OUTPUT_SEPARATOR_1,state->outputFile);
				if (infixLength==1 || hasOneBit(state->compressionBuffers[i][j][k],infixLength-2)==1) writeBits(state->compressionBuffers[i][j][k],infixLength-2,state->outputFile);
				writeChar(OUTPUT_SEPARATOR_2,state->outputFile);
			}
		}
	}	
}


void MAWs_finalize(MAWs_callback_state_t *state) {
	uint8_t i, j, k;

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
 * the first character of the nonempty right-maximal string described by $rightMaximalString$.
 * $state->char_stack$ contains numbers in $[0..3]$ represented with two bits.
 *
 * If $state->compressOutput$ is nonzero, the procedure pushes to $state->runs_stack$ a 
 * one if the right-maximal string is $a^n$ for some character $a$, and it pushes a zero
 * otherwise.
 */
static void pushChar(RightMaximalString_t rightMaximalString, MAWs_callback_state_t *state) {
	const uint64_t CAPACITY = state->char_stack_capacity;
	uint8_t c, flag;
	
	if (rightMaximalString.length>CAPACITY) {
		state->char_stack_capacity+=MY_CEIL(state->char_stack_capacity*ALLOC_GROWTH_NUM,ALLOC_GROWTH_DENOM);
		state->char_stack=(uint64_t *)realloc(state->char_stack,MY_CEIL(state->char_stack_capacity<<1,8));
		if (state->compressOutput) state->runs_stack=(uint64_t *)realloc(state->runs_stack,MY_CEIL(state->char_stack_capacity,8));
	}
	c=rightMaximalString.firstCharacter-1;
	writeTwoBits(state->char_stack,rightMaximalString.length-1,c);
	if (state->scoreState!=NULL) scorePush(c,rightMaximalString.length,state->scoreState);
	else if (state->compressOutput) {
		if (rightMaximalString.length<=1) flag=1;
		else {
			if (readBit(state->runs_stack,rightMaximalString.length-2)==0) flag=0;
			else flag=c==readTwoBits(state->char_stack,rightMaximalString.length-2)?1:0;
		}
		writeBit(state->runs_stack,rightMaximalString.length-1,flag);
	}
}


/**
 * Sets just the cells of $state->{left,right}Freqs$ that correspond to ones in 
 * $rightMaximalString.{left,right}_extension_bitmap$.
 */
static void initLeftRightFreqs(RightMaximalString_t rightMaximalString, MAWs_callback_state_t *state) {
	uint8_t i, j, char_mask;
	uint64_t frequency;
	
	char_mask=1;
	for (i=1; i<=4; i++) {
		char_mask<<=1;
		if (!(rightMaximalString.leftExtensionBitmap & char_mask)) continue;
		frequency=0;
		for (j=0; j<=5; j++) frequency+=rightMaximalString.frequency_leftRight[i][j];
		state->leftFreqs[i-1]=frequency;
	}
	char_mask=1;
	for (j=1; j<=4; j++) {
		char_mask<<=1;
		if (!(rightMaximalString.rightExtensionBitmap & char_mask)) continue;
		frequency=0;
		for (i=0; i<=5; i++) frequency+=rightMaximalString.frequency_leftRight[i][j];
		state->rightFreqs[j-1]=frequency;
	}
}


/**
 * Prints to $state->outputFile$ a string $aWb$, where $W$ is the maximal repeat described
 * by $rightMaximalString$, and $a,b$ are characters that correspond to its left- and right-
 * extensions in the text. The string is terminated by $OUTPUT_SEPARATOR_1$.
 */
static inline void printMAW(RightMaximalString_t rightMaximalString, uint8_t a, uint8_t b, MAWs_callback_state_t *state) {
	writeChar(a,state->outputFile);
	if (rightMaximalString.length!=0) writeTwoBitsReversed(state->char_stack,rightMaximalString.length-1,state->outputFile,DNA_ALPHABET);
	writeChar(b,state->outputFile);
	writeChar(OUTPUT_SEPARATOR_1,state->outputFile);
}


static void incrementLengthHistogram(RightMaximalString_t rightMaximalString, MAWs_callback_state_t *state) {
	uint64_t length, position;
	
	length=rightMaximalString.length+2;
	if (length>=state->lengthHistogramMax) position=state->lengthHistogramSize-1;
	else if (length<=state->lengthHistogramMin) position=0;
	else position=length-state->lengthHistogramMin;
	state->lengthHistogram[position]++;
}


inline void printLengthHistogram(MAWs_callback_state_t *state) {
	printf("Histogram of lengths [%llu..%llu]:\n",state->lengthHistogramMin,state->lengthHistogramMax);
	for (uint64_t i=0; i<state->lengthHistogramSize; i++) printf("%llu,%llu \n",state->lengthHistogramMin+i,state->lengthHistogram[i]);
}


/**
 * Stores, in compressed form, a MAW $a b^n c$, where $a=DNA_ALPHABET[i]$, 
 * $b=DNA_ALPHABET[j]$, $c=DNA_ALPHABET[k]$, $a \neq b$, $b \neq c$, $n \geq 1$.
 *
 * Remark: a bit of the buffer is set to one at most once during the whole traversal.
 */
static void compressMAW(uint8_t i, uint8_t j, uint8_t k, uint64_t n, MAWs_callback_state_t *state) {
	if (n>state->compressionBuffersLength[i][j][k]) {
		state->compressionBuffersLength[i][j][k]=n;
		if (n>state->compressionBuffersCapacity[i][j][k]) {
			state->compressionBuffersCapacity[i][j][k]=n<<1;  // In bits
			state->compressionBuffers[i][j][k]=(uint64_t *)realloc(state->compressionBuffers[i][j][k],MY_CEIL(n<<1,8));  // In bytes
		}
	}
	writeBit(state->compressionBuffers[i][j][k],n-1,1);
}


void MAWs_callback(RightMaximalString_t rightMaximalString, void *applicationData) {
	uint8_t i, j;
	uint8_t found, char_mask1, char_mask2;
	MAWs_callback_state_t *state = (MAWs_callback_state_t *)(applicationData);

	if (state->outputFile!=NULL && rightMaximalString.length!=0) pushChar(rightMaximalString,state);
	if (rightMaximalString.nLeftExtensions<2 || rightMaximalString.length+2<state->minLength) return;
	state->nMaxreps++;
	if (state->outputFile!=NULL && state->scoreState!=NULL) initLeftRightFreqs(rightMaximalString,state);
	char_mask1=1; found=0;
	for (i=1; i<=4; i++) {
		char_mask1<<=1;
		if ((rightMaximalString.leftExtensionBitmap&char_mask1)==0) continue;
		char_mask2=1;
		for (j=1; j<=4; j++) {
			char_mask2<<=1;
			if ( (rightMaximalString.rightExtensionBitmap&char_mask2)==0 ||
				 (rightMaximalString.frequency_leftRight[i][j]>0)
			   ) continue;
			if (state->scoreState!=NULL) {
				scoreCallback(i-1,j-1,state->leftFreqs[i-1],state->rightFreqs[j-1],state->textLength,&rightMaximalString,state->scoreState);
				if (scoreSelect(state->scoreState)==0) continue;
			}
			state->nMAWs++;
			if (found==0) found=1;
			if (rightMaximalString.length+2>state->maxLength) state->maxLength=rightMaximalString.length+2;
			if (state->lengthHistogramMin>0) incrementLengthHistogram(rightMaximalString,state);
			if (state->outputFile==NULL) continue;
			if ( state->compressOutput!=0 && 
			     i!=rightMaximalString.firstCharacter && j!=rightMaximalString.firstCharacter && 
				 readBit(state->runs_stack,rightMaximalString.length-1)!=0
			   ) compressMAW(i-1,rightMaximalString.firstCharacter-1,j-1,rightMaximalString.length,state);
			else printMAW(rightMaximalString,DNA_ALPHABET[i-1],DNA_ALPHABET[j-1],state);
			if (state->scoreState!=NULL) scorePrint(state->scoreState,state->outputFile);
			writeChar(OUTPUT_SEPARATOR_2,state->outputFile);
		}
	}
	if (found!=0) state->nMAWMaxreps++;
}


void MRWs_initialize( MAWs_callback_state_t *state,
			    	  uint64_t textLength, 
					  uint64_t minLength, 
					  uint64_t minFreq, 
					  uint64_t maxFreq, 
					  uint64_t lengthHistogramMin,
					  uint64_t lengthHistogramMax,
					  BufferedFileWriter_t *file,
					  uint8_t compressOutput ) {
	MAWs_initialize(state,textLength,minLength,lengthHistogramMin,lengthHistogramMax,file,compressOutput);
	state->minFreq=minFreq;
	state->maxFreq=maxFreq;
}


void MRWs_finalize(MAWs_callback_state_t *state) {
	MAWs_finalize(state);
}


void MRWs_callback(RightMaximalString_t rightMaximalString, void *applicationData) {
	uint8_t i, j;
	uint8_t found, char_mask1, char_mask2;
	MAWs_callback_state_t *state = (MAWs_callback_state_t *)(applicationData);

	if (state->outputFile!=NULL && rightMaximalString.length!=0) pushChar(rightMaximalString,state);
	if (rightMaximalString.nLeftExtensions<2 || rightMaximalString.length+2<state->minLength || rightMaximalString.frequency<state->maxFreq) return;
	state->nMaxreps++;
	initLeftRightFreqs(rightMaximalString,state);
	char_mask1=1; found=0;
	for (i=1; i<=4; i++) {
		char_mask1<<=1;	
		if ( (rightMaximalString.leftExtensionBitmap&char_mask1)==0 ||
			 state->leftFreqs[i-1]<state->maxFreq
		   ) continue;
		char_mask2=1;
		for (j=1; j<=4; j++) {
			char_mask2<<=1;
			if ( (rightMaximalString.rightExtensionBitmap&char_mask2)==0 ||
				 state->rightFreqs[j-1]<state->maxFreq ||
				 (rightMaximalString.frequency_leftRight[i][j]>=state->maxFreq) ||
				 (rightMaximalString.frequency_leftRight[i][j]<state->minFreq)
			   ) continue;
			if (state->scoreState!=NULL) {
				scoreCallback(i-1,j-1,state->leftFreqs[i-1],state->rightFreqs[j-1],state->textLength,&rightMaximalString,state->scoreState);
				if (scoreSelect(state->scoreState)==0) continue;
			}
			state->nMAWs++;
			if (found==0) found=1;
			if (rightMaximalString.length+2>state->maxLength) state->maxLength=rightMaximalString.length+2;
			if (state->lengthHistogramMin>0) incrementLengthHistogram(rightMaximalString,state);
			if (state->outputFile==NULL) continue;
			if ( state->compressOutput!=0 && 
			     i!=rightMaximalString.firstCharacter && j!=rightMaximalString.firstCharacter && 
				 readBit(state->runs_stack,rightMaximalString.length-1)!=0
			   ) compressMAW(i-1,rightMaximalString.firstCharacter-1,j-1,rightMaximalString.length,state);
			else printMAW(rightMaximalString,DNA_ALPHABET[i-1],DNA_ALPHABET[j-1],state);
			if (state->scoreState!=NULL) scorePrint(state->scoreState,state->outputFile);
			writeChar(OUTPUT_SEPARATOR_2,state->outputFile);
		}
	}
	if (found!=0) state->nMAWMaxreps++;
}