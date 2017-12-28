#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "MRWs_single.h"


// Character stack reallocation rate
#ifndef ALLOC_GROWTH_NUM
#define ALLOC_GROWTH_NUM 4
#endif
#ifndef ALLOC_GROWTH_DENOM
#define ALLOC_GROWTH_DENOM 3
#endif


/**
 * Space used by the application
 */
typedef struct {
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
} MRWs_callback_state_t;


static unsigned char alpha4_to_ACGT[4] = {'A','C','G','T'};


static void SLT_MRWs_callback(const SLT_params_t SLT_params, void *intern_state) {
	unsigned char i, j;
	unsigned int k;
	unsigned char char_mask1, char_mask2, markedLeft, markedRight;
	unsigned int capacity, frequency;
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

	// Initializing $markedLeft$ and $markedRight$
	markedLeft=0;
	char_mask1=1;
	for (i=1; i<=4; i++) {
		char_mask1<<=1;
		if (!(SLT_params.left_extension_bitmap & char_mask1)) continue;
		frequency=0;
		for (j=0; j<=5; j++) frequency+=SLT_params.left_right_extension_freqs[i][j];
		if (frequency<state->maxFreq) continue;
		markedLeft|=char_mask1;
	}
	markedRight=0;
	char_mask1=1;
	for (j=1; j<=4; j++) {
		char_mask1<<=1;
		if (!(SLT_params.right_extension_bitmap & char_mask1)) continue;
		frequency=0;
		for (i=0; i<=5; i++) frequency+=SLT_params.left_right_extension_freqs[i][j];
		if (frequency<state->maxFreq) continue;
		markedRight|=char_mask1;
	}
	
	// Detecting MRWs
	char_mask1=1;
	for (i=1; i<5; i++) {
		char_mask1<<=1;
		if (!(markedLeft & char_mask1)) continue;
		char_mask2=1;
		for (j=1; j<5; j++) {
			char_mask2<<=1;
			if ( !(markedRight & char_mask2) ||
				 (SLT_params.left_right_extension_freqs[i][j]>=state->maxFreq)
			   ) continue;
			state->nMRWs++;
			if (state->writeMRWs!=0) {
				// Flushing the buffer if it gets full
				if (state->MRW_buffer_size+SLT_params.string_depth+3 > state->MRW_buffer_capacity) {
					fwrite(state->MRW_buffer,state->MRW_buffer_size,sizeof(char),state->file);
					state->MRW_buffer_size=0;
				}
				state->MRW_buffer[state->MRW_buffer_size++]=alpha4_to_ACGT[i-1];
				for (k=0; k<SLT_params.string_depth; k++) state->MRW_buffer[state->MRW_buffer_size++]=alpha4_to_ACGT[state->char_stack[SLT_params.string_depth-1-k]-1];
				state->MRW_buffer[state->MRW_buffer_size++]=alpha4_to_ACGT[j-1];
				state->MRW_buffer[state->MRW_buffer_size++]=OUTPUT_SEPARATOR;
			}
		}
	}
}


unsigned int find_MRWs_single(Basic_BWT_t *BBWT, unsigned int minLength, unsigned int minFreq, unsigned int maxFreq, unsigned char writeMRWs, char *filePath) {
	SLT_iterator_t_single_string SLT_iterator;
	MRWs_callback_state_t state;
	
	state.minLength=minLength;
	state.minFreq=minFreq;
	state.maxFreq=maxFreq;
	state.nMRWs=0;
	state.char_stack_capacity=minLength;
	state.char_stack=(unsigned char *)malloc(state.char_stack_capacity);
	state.writeMRWs=writeMRWs;
	state.MRW_buffer_capacity=minLength*10;  // Arbitrary choice
	state.MRW_buffer=(unsigned char *)malloc(state.MRW_buffer_capacity);
	state.MRW_buffer_size=0;
	if (writeMRWs!=0) state.file=fopen(filePath,"a");
	
	SLT_iterator=new_SLT_iterator(SLT_MRWs_callback,&state,BBWT,SLT_stack_trick);
	SLT_execute_iterator(&SLT_iterator);
	
	if (writeMRWs!=0) fclose(state.file);
	free(state.char_stack);
	free(state.MRW_buffer);
	return state.nMRWs;
}
