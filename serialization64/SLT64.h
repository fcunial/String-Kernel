#ifndef SLT64_h
#define SLT64_h
#include"indexed_DNA5_seq64.h"
#include"DNA5_Basic_BWT64.h"

typedef struct 
{
	unsigned int string_depth;
	unsigned char nright_extensions;
	unsigned char nleft_extensions;
	unsigned char right_extension_bitmap;
	unsigned char left_extension_bitmap;
	unsigned long long left_right_extension_freqs[6][6];
	unsigned long long left_ext_bwt_start[5];
	unsigned long long bwt_start;
	unsigned long  revbwt_start;
	unsigned long interval_size;
	unsigned int WL_char;
} SLT64_params_t;

typedef void (*SLT64_callback_t)(const SLT64_params_t * SLT64_params,void * intern_state);

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
	SLT64_callback_t SLT_callback;
	unsigned int options;
	void * intern_state;
	Basic_BWT64_t * BBWT;
}SLT64_iterator_t;


static inline SLT64_iterator_t * new_SLT64_iterator(SLT64_callback_t SLT_callback,
		void * intern_state,Basic_BWT64_t * BBWT,unsigned int options)
{
	SLT64_iterator_t * SLT_iterator=(SLT64_iterator_t *)malloc(sizeof(SLT64_iterator_t));
	SLT_iterator->SLT_callback=SLT_callback;
	SLT_iterator->intern_state=intern_state;
	SLT_iterator->BBWT=BBWT;
	SLT_iterator->options=options;
	return SLT_iterator;
};

static inline void free_SLT64_iterator(SLT64_iterator_t * SLT_iterator)
{
	free(SLT_iterator);
};

void SLT64_execute_iterator(SLT64_iterator_t * SLT_iterator); 


#endif
