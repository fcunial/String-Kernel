#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "MAWs_single.h"


// Character stack reallocation rate
#ifndef ALLOC_GROWTH_NUM
#define ALLOC_GROWTH_NUM 4
#endif
#ifndef ALLOC_GROWTH_DENOM
#define ALLOC_GROWTH_DENOM 3
#endif
#ifndef MAW_SEPARATOR
#define MAW_SEPARATOR '\n'
#endif


/**
 * Space used by the application
 */
typedef struct {
	unsigned int minLength;  // Minimum desired length of a MAW
	unsigned int nMAWs;  // Total number of MAWs
	
	// Character stack
	unsigned char *char_stack;  // Indexed from zero
	unsigned int char_stack_capacity;  // Number of characters in the stack
	
	// Output buffer
	unsigned char writeMAWs;  // 0 iff MAWs should not be written to the output
	unsigned char *MAW_buffer;
	unsigned int MAW_buffer_capacity;  // Maximum number of chars in the buffer
	unsigned int MAW_buffer_size;  // Number of chars currently in the buffer
	FILE *file;
} MAWs_callback_state_t;


static unsigned char alpha4_to_ACGT[4] = {'A','C','G','T'};


static void SLT_MAWs_callback(const SLT_params_t SLT_params, void *intern_state) {
	unsigned int i, j, k;
	unsigned int capacity, char_mask1, char_mask2;
	MAWs_callback_state_t *state = (MAWs_callback_state_t *)(intern_state);

	// Pushing to $char_stack$ the label of the last Weiner link
	if (state->writeMAWs!=0 && SLT_params.string_depth!=0) {
		capacity=state->char_stack_capacity;
		if (SLT_params.string_depth>capacity) {
			state->char_stack_capacity+=(capacity*ALLOC_GROWTH_NUM)/ALLOC_GROWTH_DENOM;
			state->char_stack=(unsigned char *)realloc(state->char_stack,state->char_stack_capacity);
		}
		state->char_stack[SLT_params.string_depth-1]=SLT_params.WL_char;
	}

	// Detecting MAWs
	if (SLT_params.nleft_extensions<2 || SLT_params.string_depth+2<state->minLength) return;
	char_mask1=1;
	for (i=1; i<5; i++) {
		char_mask1<<=1;
		if (!(SLT_params.left_extension_bitmap & char_mask1)) continue;
		char_mask2=1;
		for (j=1; j<5; j++) {
			char_mask2<<=1;
			if ( !(SLT_params.right_extension_bitmap & char_mask2) ||
				 (SLT_params.left_right_extension_freqs[i][j]>0)
			   ) continue;
			state->nMAWs++;
			if (state->writeMAWs!=0) {
				// Flushing the buffer if it gets full
				if (state->MAW_buffer_size+SLT_params.string_depth+3 > state->MAW_buffer_capacity) {
					fwrite(state->MAW_buffer,state->MAW_buffer_size,sizeof(char),state->file);
					state->MAW_buffer_size=0;
				}
				state->MAW_buffer[state->MAW_buffer_size++]=alpha4_to_ACGT[i-1];
				for (k=0; k<SLT_params.string_depth; k++) state->MAW_buffer[state->MAW_buffer_size++]=alpha4_to_ACGT[state->char_stack[SLT_params.string_depth-1-k]-1];
				state->MAW_buffer[state->MAW_buffer_size++]=alpha4_to_ACGT[j-1];
				state->MAW_buffer[state->MAW_buffer_size++]=MAW_SEPARATOR;
			}
		}
	}
}


unsigned int find_MAWs_single(Basic_BWT_t *BBWT, unsigned int minLength, unsigned char writeMAWs, char *filePath) {
	SLT_iterator_t_single_string SLT_iterator;
	MAWs_callback_state_t state;
	
	state.minLength=minLength;
	state.nMAWs=0;
	state.char_stack_capacity=minLength;
	state.char_stack=(unsigned char *)malloc(state.char_stack_capacity);
	state.writeMAWs=writeMAWs;
	state.MAW_buffer_capacity=minLength*10;  // Arbitrary choice
	state.MAW_buffer=(unsigned char *)malloc(state.MAW_buffer_capacity);
	state.MAW_buffer_size=0;
	if (writeMAWs!=0) state.file=fopen(filePath,"a");
	
	SLT_iterator=new_SLT_iterator(SLT_MAWs_callback,&state,BBWT,SLT_stack_trick);
	SLT_execute_iterator(&SLT_iterator);
	
	if (writeMAWs!=0) fclose(state.file);
	free(state.char_stack);
	free(state.MAW_buffer);
	return state.nMAWs;
}
