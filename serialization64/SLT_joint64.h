#ifndef SLT_joint_h
#define SLT_joint_h
#include"indexed_DNA5_seq64.h"
#include"DNA5_Basic_BWT64.h"

typedef struct 
{
	unsigned long long string_depth;
	unsigned char nright_extensions1;
	unsigned char nleft_extensions1;
	unsigned char nright_extensions2;
	unsigned char nleft_extensions2;
	unsigned char right_extension_bitmap1;
	unsigned char right_extension_bitmap2;
	unsigned char left_extension_bitmap1;
	unsigned char left_extension_bitmap2;
	unsigned long long left_right_extension_freqs1[6][6];
	unsigned long long left_right_extension_freqs2[6][6];
//	unsigned int left_ext_bwt_start1[5];
//	unsigned int left_ext_bwt_start2[5];
//	unsigned int bwt_start1;
//	unsigned int bwt_start2;
//	unsigned int revbwt_start1;
//	unsigned int revbwt_start2;
	unsigned long long interval_size1;
	unsigned long long interval_size2;
	unsigned int WL_char;
} SLT64_joint_params_t;
typedef void (*SLT64_joint_callback_t)(const SLT64_joint_params_t * SLT_params,void * intern_state);

#ifndef SLT_arbitrary_order
#define SLT_arbitrary_order 0
#endif
#ifndef SLT_lex_order
#define SLT_lex_order 1
#endif
#ifndef SLT_stack_trick 
#define SLT_stack_trick 2
#endif
#ifndef SLT_joint_and_enum 
#define SLT_joint_and_enum 0
#endif
#ifndef SLT_joint_or_enum 
#define SLT_joint_or_enum 4
#endif


typedef struct 
{
	SLT64_joint_callback_t SLT_callback;
	unsigned int options;
	void * intern_state;
	Basic_BWT64_t * BBWT1;
	Basic_BWT64_t * BBWT2;
} SLT64_joint_iterator_t;


static inline SLT64_joint_iterator_t * new_SLT64_joint_iterator(SLT64_joint_callback_t SLT_callback,
		void * intern_state, Basic_BWT64_t * BBWT1, Basic_BWT64_t * BBWT2,
		unsigned int options)
{
	SLT64_joint_iterator_t * SLT_iterator=(SLT64_joint_iterator_t *)malloc(sizeof(SLT64_joint_iterator_t));
	SLT_iterator->SLT_callback=SLT_callback;
	SLT_iterator->intern_state=intern_state;
	SLT_iterator->BBWT1=BBWT1;
	SLT_iterator->BBWT2=BBWT2;
	SLT_iterator->options=options;
	return SLT_iterator;
};

static inline void free_SLT64_joint_iterator(SLT64_joint_iterator_t * SLT_iterator)
{
	free(SLT_iterator);
};

void SLT64_joint_execute_iterator(SLT64_joint_iterator_t * SLT_iterator); 


#endif
