#include"stdlib.h"
#include"stdio.h"
#include"SLT_MAWs_single_string.h"
#include <string.h>




#define alloc_growth_num 4
#define alloc_growth_denom 3


typedef struct 
{
	unsigned int minlen;
	unsigned int nMAW_capacity;
	unsigned char * MAW_buffer;
	unsigned int MAW_buffer_idx;
	unsigned int char_stack_capacity;
	unsigned char * char_stack;
	unsigned int nMAWs;
	double * KL;
	unsigned int KL_capacity;

	FILE *file;
} MAWs_callback_state_t;


static inline unsigned int is_powof2(unsigned int x)
{
	return ((x&(x-1))==0);

};
static unsigned char alpha4_to_ACGT[4]={'A','C','G','T'};

void SLT_MAWs_callback(const SLT_params_t * SLT_params,void * intern_state, unsigned int mem)
{
	MAWs_callback_state_t * state= (MAWs_callback_state_t*)(intern_state);
	unsigned int char_mask1;
	unsigned int char_mask2;
	unsigned int i,j,k,h;

	if(state->nMAW_capacity==0 && mem) {
		state->nMAW_capacity=1<<16;
		state->MAW_buffer=(unsigned char *) malloc(state->nMAW_capacity);
	}
	if(SLT_params->string_depth!=0)
	{
		if(SLT_params->string_depth>state->char_stack_capacity)
		{
			state->char_stack_capacity=(state->char_stack_capacity+1)*alloc_growth_num/alloc_growth_denom;
			state->char_stack=(unsigned char *)realloc(state->char_stack,state->char_stack_capacity);
		};
		state->char_stack[SLT_params->string_depth-1]=SLT_params->WL_char;
	}

	// Check that we are at a maximal repeat of length at least minlen-2
	if(SLT_params->string_depth + 2 < state->minlen || 
		SLT_params->nleft_extensions < 2 || SLT_params->nright_extensions < 2)
		return;

	char_mask1=1;
	for(i=1;i<5;i++) {
		char_mask1<<=1;
		char_mask2=1;
		for(j=1;j<5;j++) {
			char_mask2<<=1;
			if(((SLT_params->right_extension_bitmap&char_mask2)
					&& (SLT_params->left_extension_bitmap&char_mask1)
							&& SLT_params->left_right_extension_freqs[i][j]==0)) {
				// We have a MAW. Write it to the output
				state->nMAWs++;
				if(mem) {
					if(state->MAW_buffer_idx+SLT_params->string_depth+3>state->nMAW_capacity) {
						fwrite(state->MAW_buffer, state->MAW_buffer_idx, sizeof(char), state->file);
						state->MAW_buffer_idx=0;
					}
					state->MAW_buffer[state->MAW_buffer_idx++]=alpha4_to_ACGT[i-1];
					for(k=0;k<SLT_params->string_depth;k++){
						if (state->char_stack[SLT_params->string_depth-k-1]-1 > 3)
							printf("%d ", state->char_stack[SLT_params->string_depth-k-1]-1);
						state->MAW_buffer[state->MAW_buffer_idx++]=
								alpha4_to_ACGT[state->char_stack[SLT_params->string_depth-k-1]-1];
					}
					state->MAW_buffer[state->MAW_buffer_idx++]=alpha4_to_ACGT[j-1];
					state->MAW_buffer[state->MAW_buffer_idx++]='\n';
				}
			}
		}
	}

};
unsigned int SLT_find_MAWs_single_string(Basic_BWT_t * BBWT1, unsigned int minlen, unsigned int * _nMAWs1,
		double * _output_result, unsigned int mem)
{
	SLT_iterator_t_single_string * SLT_iterator;
	MAWs_callback_state_t state;

	unsigned int i;
	state.nMAWs=0;
	state.MAW_buffer=0;
	state.MAW_buffer_idx=0;
	state.nMAW_capacity=0;
	state.minlen=minlen;
	state.char_stack_capacity=4;
	state.char_stack=(unsigned char *) malloc(state.char_stack_capacity);
	FILE *f;
	if(mem) {
		f=fopen("output.txt", "a");
		state.file=f;
	}
	SLT_iterator=new_SLT_iterator(SLT_MAWs_callback,&state,BBWT1,SLT_stack_trick, mem);
	SLT_execute_iterator(SLT_iterator);

};
