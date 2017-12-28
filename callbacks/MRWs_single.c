#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "MRWs_single.h"


/**
 * Space used by the application
 */
typedef struct {
	unsigned int textLength;
	unsigned int minLength;  // Minimum desired length of a MRW
	unsigned int minFreq, maxFreq;
	unsigned int nMRWs;  // Total number of MRWs
	
	// Character stack
	unsigned char *char_stack;  // Indexed from zero
	unsigned int char_stack_capacity;  // Number of characters in the stack
	
	// Output buffer
	unsigned char writeMRWs;  // 0 iff MRWs should not be written to the output
	unsigned char *MRW_buffer;
	unsigned int MRW_buffer_capacity;  // Maximum number of chars in the buffer
	unsigned int MRW_buffer_size;  // Number of chars currently in the buffer
	FILE *file;
	
	// Score
	unsigned char computeScore;
	unsigned int *leftFreqs, *rightFreqs;  // Frequency of each left/right extension. Only for ACGT, indexed from zero.
	char *probabilityBuffer;  // Temporary space for the string representation of the probability
	char *scoreBuffer;  // Temporary space for the string representation of the score
} MRWs_callback_state_t;


static unsigned char alpha4_to_ACGT[4] = {'A','C','G','T'};


static void SLT_MRWs_callback(const SLT_params_t SLT_params, void *intern_state) {
	unsigned char i, j;
	unsigned int k;
	unsigned char char_mask1, char_mask2;
	unsigned int capacity, frequency;
	unsigned int nCharacters, nCharactersProbability, nCharactersScore;
	double multiplier, probability, score;
	MRWs_callback_state_t *state = (MRWs_callback_state_t *)(intern_state);

	// Pushing to $char_stack$ the label of the last Weiner link
	if (state->writeMRWs!=0 && SLT_params.string_depth!=0) {
		capacity=state->char_stack_capacity;
		if (SLT_params.string_depth>capacity) {
			state->char_stack_capacity+=(capacity*ALLOC_GROWTH_NUM)/ALLOC_GROWTH_DENOM;
			state->char_stack=(unsigned char *)realloc(state->char_stack,state->char_stack_capacity);
		}
		state->char_stack[SLT_params.string_depth-1]=SLT_params.WL_char;
	}
	if ( SLT_params.interval_size<state->maxFreq || 
		 SLT_params.string_depth+2<state->minLength ||
		 SLT_params.nleft_extensions<2
	   ) return;

	// Initializing $leftFreqs$ and $rightFreqs$
	char_mask1=1;
	for (i=1; i<=4; i++) {
		char_mask1<<=1;
		if (!(SLT_params.left_extension_bitmap & char_mask1)) continue;
		frequency=0;
		for (j=0; j<=5; j++) frequency+=SLT_params.left_right_extension_freqs[i][j];
		state->leftFreqs[i-1]=frequency;
	}
	char_mask1=1;
	for (j=1; j<=4; j++) {
		char_mask1<<=1;
		if (!(SLT_params.right_extension_bitmap & char_mask1)) continue;
		frequency=0;
		for (i=0; i<=5; i++) frequency+=SLT_params.left_right_extension_freqs[i][j];
		state->rightFreqs[j-1]=frequency;
	}
	
	// Detecting MRWs
	char_mask1=1;
	for (i=1; i<5; i++) {
		char_mask1<<=1;
		if ( !(SLT_params.left_extension_bitmap & char_mask1) ||
			 state->leftFreqs[i-1]<state->maxFreq
		   ) continue;
		char_mask2=1;
		for (j=1; j<5; j++) {
			char_mask2<<=1;
			if ( !(SLT_params.right_extension_bitmap & char_mask2) ||
				 state->rightFreqs[j-1]<state->maxFreq ||
				 (SLT_params.left_right_extension_freqs[i][j]>=state->maxFreq) ||
				 (SLT_params.left_right_extension_freqs[i][j]<state->minFreq)
			   ) continue;
			state->nMRWs++;
			if (state->writeMRWs!=0) {
				nCharacters=SLT_params.string_depth+2+1;
				if (state->computeScore!=0) {
					multiplier=state->textLength-(SLT_params.string_depth+2)+2;
					multiplier=(multiplier+1)/(multiplier*multiplier);
					probability=(multiplier*(state->leftFreqs[i-1]*state->rightFreqs[j-1]))/SLT_params.interval_size;
					nCharactersProbability=sprintf(state->probabilityBuffer,"%g",probability);
					nCharacters+=nCharactersProbability+1;
					score=((double)(SLT_params.left_right_extension_freqs[i][j]))/(state->textLength-(SLT_params.string_depth+2)+1);
					score=(score-probability)/probability;
					nCharactersScore=sprintf(state->scoreBuffer,"%g",score);
					nCharacters+=nCharactersScore+1;
				}
				// Flushing the buffer if it gets full
				if (state->MRW_buffer_size+nCharacters > state->MRW_buffer_capacity) {
					fwrite(state->MRW_buffer,state->MRW_buffer_size,sizeof(char),state->file);
					state->MRW_buffer_size=0;
				}
				state->MRW_buffer[state->MRW_buffer_size++]=alpha4_to_ACGT[i-1];
				for (k=0; k<SLT_params.string_depth; k++) state->MRW_buffer[state->MRW_buffer_size++]=alpha4_to_ACGT[state->char_stack[SLT_params.string_depth-1-k]-1];
				state->MRW_buffer[state->MRW_buffer_size++]=alpha4_to_ACGT[j-1];
				state->MRW_buffer[state->MRW_buffer_size++]=OUTPUT_SEPARATOR;
				if (state->computeScore!=0) {
					for (k=0; k<nCharactersProbability; k++) state->MRW_buffer[state->MRW_buffer_size++]=state->probabilityBuffer[k];
					state->MRW_buffer[state->MRW_buffer_size++]=OUTPUT_SEPARATOR;
					for (k=0; k<nCharactersScore; k++) state->MRW_buffer[state->MRW_buffer_size++]=state->scoreBuffer[k];
					state->MRW_buffer[state->MRW_buffer_size++]=OUTPUT_SEPARATOR;
				}
			}
		}
	}
}


unsigned int find_MRWs_single(Basic_BWT_t *BBWT, unsigned int textLength, unsigned int minLength, unsigned int minFreq, unsigned int maxFreq, unsigned char writeMRWs, unsigned char computeScore, char *filePath) {
	SLT_iterator_t_single_string SLT_iterator;
	MRWs_callback_state_t state;
	
	// Running the iterator
	state.textLength=textLength;
	state.minLength=minLength;
	state.minFreq=minFreq;
	state.maxFreq=maxFreq;
	state.nMRWs=0;
	state.char_stack_capacity=BUFFER_CHUNK;  // Arbitrary choice
	state.char_stack=(unsigned char *)malloc(state.char_stack_capacity);
	state.writeMRWs=writeMRWs;
	state.computeScore=computeScore;
	state.MRW_buffer_capacity=BUFFER_CHUNK;  // Arbitrary choice
	state.MRW_buffer=(unsigned char *)malloc(state.MRW_buffer_capacity);
	state.MRW_buffer_size=0;
	if (writeMRWs!=0) state.file=fopen(filePath,"a");
	if (computeScore!=0) {
		state.leftFreqs=(unsigned int *)malloc(4);
		state.rightFreqs=(unsigned int *)malloc(4);
		state.probabilityBuffer=(char *)malloc(50);
		state.scoreBuffer=(char *)malloc(50);
	}
	
	SLT_iterator=new_SLT_iterator(SLT_MRWs_callback,&state,BBWT,SLT_stack_trick);
	SLT_execute_iterator(&SLT_iterator);
	
	// Flushing the buffer one more time
	if (writeMRWs!=0 && state.MRW_buffer_size>0) fwrite(state.MRW_buffer,state.MRW_buffer_size,sizeof(char),state.file);
	
	// Cleaning
	if (writeMRWs!=0) fclose(state.file);
	free(state.char_stack);
	free(state.MRW_buffer);
	if (computeScore!=0) {
		free(state.leftFreqs);
		free(state.rightFreqs);
		free(state.probabilityBuffer);
		free(state.scoreBuffer);
	}
	return state.nMRWs;
}
