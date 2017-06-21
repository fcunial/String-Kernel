#ifndef SLT_h
#define SLT_h
#include"indexed_DNA5_seq.h"
#include"DNA5_Basic_BWT.h"

typedef struct 
{
	unsigned int string_depth;
	unsigned char nright_extensions1;
	unsigned char nleft_extensions1;
	unsigned char nright_extensions2;
	unsigned char nleft_extensions2;
	unsigned char right_extension_bitmap1;
	unsigned char right_extension_bitmap2;
	unsigned char left_extension_bitmap1;
	unsigned char left_extension_bitmap2;
	unsigned int left_right_extension_freqs1[6][6];
	unsigned int left_right_extension_freqs2[6][6];
	unsigned int interval_size1;
	unsigned int interval_size2;
	unsigned int WL_char;
} SLT_joint_params_t;

typedef void (*SLT_joint_callback_t)(const SLT_joint_params_t * SLT_params,void * intern_state, unsigned int memory);
typedef void* (*SLT_cloner_t)(void* p, unsigned int t);
typedef void (*SLT_combiner_t)(void** intern_state, void* state, unsigned int t, unsigned int mem);
typedef void (*SLT_free_t)(void* intern_state, unsigned int mem);

#ifndef SLT_arbitrary_order
#define SLT_arbitrary_order 0
#endif
#ifndef SLT_lex_order
#define SLT_lex_order 1
#endif
#ifndef SLT_stack_trick
#define SLT_stack_trick 4
#endif
#ifndef SLT_joint_and_enum
#define SLT_joint_and_enum 0
#endif
#ifndef SLT_joint_or_enum
#define SLT_joint_or_enum 4
#endif


typedef struct
{
	SLT_joint_callback_t SLT_callback;
	SLT_cloner_t SLT_cloner;
	SLT_combiner_t SLT_combiner;
	SLT_free_t SLT_free;
	unsigned int options;
	void * intern_state;
	Basic_BWT_t * BBWT1;
	Basic_BWT_t * BBWT2;
	unsigned int mem;
	unsigned int cores;
} SLT_joint_iterator_t;

typedef struct
{
	unsigned int string_depth;
	unsigned int interval_start1;
	unsigned int interval_start2;
	unsigned int child_freqs1[6];
	unsigned int child_freqs2[6];
	unsigned int interval_size1;
	unsigned int interval_size2;
	unsigned char WL_char;
} SLT_stack_item_t;


static inline SLT_joint_iterator_t * new_SLT_joint_iterator(SLT_joint_callback_t SLT_callback,SLT_cloner_t SLT_cloner,
		SLT_combiner_t SLT_combiner, SLT_free_t SLT_free, void * intern_state,
		Basic_BWT_t * BBWT1,Basic_BWT_t * BBWT2,unsigned int options, unsigned int mem,  unsigned int cores)
{
	SLT_joint_iterator_t * SLT_iterator=(SLT_joint_iterator_t *)malloc(sizeof(SLT_joint_iterator_t));
	SLT_iterator->SLT_callback=SLT_callback;
	SLT_iterator->SLT_cloner=SLT_cloner;
	SLT_iterator->SLT_combiner=SLT_combiner;
	SLT_iterator->SLT_free=SLT_free;
	SLT_iterator->intern_state=intern_state;
	SLT_iterator->BBWT1=BBWT1;
	SLT_iterator->BBWT2=BBWT2;
	SLT_iterator->options=options;
	SLT_iterator->mem=mem;
	SLT_iterator->cores= cores;
	return SLT_iterator;
};

static inline void free_SLT_joint_iterator(SLT_joint_iterator_t * SLT_iterator)
{
	free(SLT_iterator);
};

void SLT_joint_execute_iterator(SLT_joint_iterator_t * SLT_iterator);
void SLT_slave(SLT_joint_iterator_t * SLT_iterator, SLT_stack_item_t stack_item, void* intern_state);

#endif
