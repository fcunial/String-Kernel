#include "MAWs_single.h"
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "../io/io.h"
#include "../io/bits.h"


#ifndef INITIAL_CHAR_STACK_CAPACITY
#define INITIAL_CHAR_STACK_CAPACITY 128  // In characters. The stack can grow.
#endif


static const unsigned char BYTES_PER_LONG = sizeof(unsigned long);
static const unsigned char BITS_PER_LONG = BYTES_PER_LONG<<3;


static void initCompressedOutput(MAWs_callback_state_t *state) {
	unsigned int i, j, k;
	
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
				for (k=0; k<j; k++) state->compressionBuffers[i][j][k]=(unsigned long *)malloc(MY_CEIL(BUFFER_CHUNK,BYTES_PER_LONG));
				for (k=j+1; k<4; k++) state->compressionBuffers[i][j][k]=(unsigned long *)malloc(MY_CEIL(BUFFER_CHUNK,BYTES_PER_LONG));
			}
			for (j=i+1; j<4; j++) {
				for (k=0; k<j; k++) state->compressionBuffers[i][j][k]=(unsigned long *)malloc(MY_CEIL(BUFFER_CHUNK,BYTES_PER_LONG));
				for (k=j+1; k<4; k++) state->compressionBuffers[i][j][k]=(unsigned long *)malloc(MY_CEIL(BUFFER_CHUNK,BYTES_PER_LONG));
			}
		}
	}
}


void MAWs_initialize( MAWs_callback_state_t *state, 
	                  unsigned int textLength, 
				      unsigned int minLength, 
					  unsigned int lengthHistogramMin,
					  unsigned int lengthHistogramMax,
					  buffered_file_writer_t *file,
					  unsigned char compressOutput ) {	
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
		state->char_stack=(unsigned long *)malloc(MY_CEIL(state->char_stack_capacity<<1,BITS_PER_LONG)*BYTES_PER_LONG);  // In bytes
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
	initCompressedOutput(state);
	if (state->outputFile!=NULL && state->compressOutput!=0) state->runs_stack=(unsigned long *)malloc(MY_CEIL(state->char_stack_capacity,8));
	else state->runs_stack=NULL;
}


static void mergeCompressedOutput(MAWs_callback_state_t *from, MAWs_callback_state_t *to) {
	unsigned int i, j, k, p, nBits, nLongs;
	unsigned long *tmp;
	
	for (i=0; i<4; i++) {
		for (j=0; j<i; j++) {
			for (k=0; k<j; k++) {
				nBits=to->compressionBuffersLength[i][j][k];
				if (from->compressionBuffersLength[i][j][k]>nBits) nBits=from->compressionBuffersLength[i][j][k];
				if (nBits==0) continue;
				tmp=(unsigned long *)calloc(MY_CEIL(nBits,8),1);
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
				tmp=(unsigned long *)calloc(MY_CEIL(nBits,8),1);
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
				tmp=(unsigned long *)calloc(MY_CEIL(nBits,8),1);
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
				tmp=(unsigned long *)calloc(MY_CEIL(nBits,8),1);
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
	unsigned int nBytes;
	
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
		to->outputFile=(buffered_file_writer_t *)malloc(sizeof(buffered_file_writer_t));
		sprintf(cloneBuffer,"%s-%d.out",pathPrefix,id);
		initializeBufferedFileWriter(to->outputFile,cloneBuffer);
	}
	
	// Character stack
	if (from->char_stack!=NULL) {
		to->char_stack_capacity=from->char_stack_capacity;
		nBytes=MY_CEIL(to->char_stack_capacity<<1,8);
		to->char_stack=(unsigned long *)malloc(nBytes);
		memcpy(to->char_stack,from->char_stack,nBytes);
	}
	
	// Scores
	to->leftFreqs=(unsigned int *)malloc(strlen(DNA_ALPHABET)*sizeof(unsigned int));
	to->rightFreqs=(unsigned int *)malloc(strlen(DNA_ALPHABET)*sizeof(unsigned int));
	if (from->scoreState!=NULL) {
		to->scoreState=(score_state_t *)malloc(sizeof(score_state_t));
		scoreClone(from->scoreState,to->scoreState);
	}
	
	// Histograms
	if (from->lengthHistogramMin!=0) {
		to->lengthHistogramMin=from->lengthHistogramMin;
		to->lengthHistogramMax=from->lengthHistogramMax;
		to->lengthHistogramSize=from->lengthHistogramSize;
		to->lengthHistogram=(unsigned int *)calloc(to->lengthHistogramSize,sizeof(unsigned int));
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
		to->runs_stack=(unsigned long *)malloc(nBytes);
		memcpy(to->runs_stack,from->runs_stack,nBytes);
	}
	else to->runs_stack=NULL;
}


void mergeMAWState(MAWs_callback_state_t *from, MAWs_callback_state_t *to) {
	unsigned int i;
	
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
				if (infixLength==1 || hasOneBit(state->compressionBuffers[i][j][k],infixLength-2)==1) writeBits(state->compressionBuffers[i][j][k],infixLength-2,state->outputFile);
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
 * If $state->compressOutput$ is nonzero, the procedure pushes to $state->runs_stack$ a 
 * one if the right-maximal string is $a^n$ for some character $a$, and it pushes a zero
 * otherwise.
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


static void incrementLengthHistogram(SLT_params_t SLT_params, MAWs_callback_state_t *state) {
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
 * Stores, in compressed form, a MAW $a b^n c$, where $a=DNA_ALPHABET[i]$, 
 * $b=DNA_ALPHABET[j]$, $c=DNA_ALPHABET[k]$, $a \neq b$, $b \neq c$, $n \geq 1$.
 *
 * Remark: a bit of the buffer is set to one at most once during the whole traversal.
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
			if (state->scoreState!=NULL) {
				scoreCallback(i-1,j-1,state->leftFreqs[i-1],state->rightFreqs[j-1],state->textLength,&SLT_params,state->scoreState);
				if (scoreSelect(state->scoreState)==0) continue;
			}
			state->nMAWs++;
			if (found==0) found=1;
			if (SLT_params.string_depth+2>state->maxLength) state->maxLength=SLT_params.string_depth+2;
			if (state->lengthHistogramMin>0) incrementLengthHistogram(SLT_params,state);
			if (state->outputFile==NULL) continue;
			if ( state->compressOutput!=0 && 
			     i!=SLT_params.WL_char && j!=SLT_params.WL_char && 
				 readBit(state->runs_stack,SLT_params.string_depth-1)!=0
			   ) compressMAW(i-1,SLT_params.WL_char-1,j-1,SLT_params.string_depth,state);
			else printMAW(SLT_params,DNA_ALPHABET[i-1],DNA_ALPHABET[j-1],state);
			if (state->scoreState!=NULL) scorePrint(state->scoreState,state->outputFile);
			writeChar(OUTPUT_SEPARATOR_2,state->outputFile);
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
			if (state->scoreState!=NULL) {
				scoreCallback(i-1,j-1,state->leftFreqs[i-1],state->rightFreqs[j-1],state->textLength,&SLT_params,state->scoreState);
				if (scoreSelect(state->scoreState)==0) continue;
			}
			state->nMAWs++;
			if (found==0) found=1;
			if (SLT_params.string_depth+2>state->maxLength) state->maxLength=SLT_params.string_depth+2;
			if (state->lengthHistogramMin>0) incrementLengthHistogram(SLT_params,state);
			if (state->outputFile==NULL) continue;
			if ( state->compressOutput!=0 && 
			     i!=SLT_params.WL_char && j!=SLT_params.WL_char && 
				 readBit(state->runs_stack,SLT_params.string_depth-1)!=0
			   ) compressMAW(i-1,SLT_params.WL_char-1,j-1,SLT_params.string_depth,state);
			else printMAW(SLT_params,DNA_ALPHABET[i-1],DNA_ALPHABET[j-1],state);
			if (state->scoreState!=NULL) scorePrint(state->scoreState,state->outputFile);
			writeChar(OUTPUT_SEPARATOR_2,state->outputFile);
		}
	}
	if (found!=0) state->nMAWMaxreps++;
}