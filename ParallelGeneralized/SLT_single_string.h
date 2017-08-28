#ifndef SLT_h
#define SLT_h
#include"indexed_DNA5_seq.h"
#include"DNA5_Basic_BWT.h"

typedef struct 
{
	unsigned int string_depth;
	unsigned char nright_extensions;
	unsigned char nleft_extensions;
	unsigned char right_extension_bitmap;
	unsigned char left_extension_bitmap;
	unsigned int left_right_extension_freqs[6][6];
	unsigned int left_ext_bwt_start[5];
	unsigned int bwt_start;
	unsigned int revbwt_start;
	unsigned int interval_size;
	unsigned int WL_char;
} SLT_params_t;
typedef void (*SLT_MAWs_callback_t)(const SLT_params_t * SLT_params,void * intern_state, unsigned int mem);

#ifndef SLT_arbitrary_order
#define SLT_arbitrary_order 0
#endif
#ifndef SLT_lex_order
#define SLT_lex_order 1
#endif
#ifndef SLT_stack_trick 
#define SLT_stack_trick 2
#endif


typedef struct 
{
	SLT_MAWs_callback_t SLT_callback;
	unsigned int options;
	void * intern_state;
	unsigned int mem;
	Basic_BWT_t * BBWT;
}SLT_iterator_t_single_string;


static inline SLT_iterator_t_single_string * new_SLT_iterator(SLT_MAWs_callback_t SLT_callback,
		void * intern_state,Basic_BWT_t * BBWT,unsigned int options, unsigned int mem)
{
	SLT_iterator_t_single_string * SLT_iterator=(SLT_iterator_t_single_string *)malloc(sizeof(SLT_iterator_t_single_string));
	SLT_iterator->SLT_callback=SLT_callback;
	SLT_iterator->intern_state=intern_state;
	SLT_iterator->BBWT=BBWT;
	SLT_iterator->options=options;
	SLT_iterator->mem=mem;
	return SLT_iterator;
};

static inline void free_SLT_iterator(SLT_iterator_t_single_string * SLT_iterator)
{
	free(SLT_iterator);
};

void SLT_execute_iterator(SLT_iterator_t_single_string * SLT_iterator); 


#endif
