/**
 * Interface of the one-string iterator
 */
#ifndef SLT_single_string_h
#define SLT_single_string_h


#ifndef SLT_arbitrary_order
#define SLT_arbitrary_order 0
#endif
#ifndef SLT_lex_order
#define SLT_lex_order 1
#endif
#ifndef SLT_stack_trick 
#define SLT_stack_trick 2
#endif


#include "DNA5_Basic_BWT.h"


/** 
 * Function called by the one-string iterator for every string it enumerates.
 *
 * @param intern_state pointer to space used by the application that implements the 
 * callback.
 */
typedef void (*SLT_callback_t)(const SLT_params_t *SLT_params, void *intern_state);


/**
 * The one-string iterator
 */
typedef struct {
	SLT_callback_t SLT_callback;
	void *intern_state;
	Basic_BWT_t *BBWT;
	unsigned int options;
} SLT_iterator_t_single_string;


/** 
 * Data returned by the one-string iterator
 */
typedef struct {
	// Properties of the current right-maximal string W
	unsigned int string_depth;
	unsigned int bwt_start;
	unsigned int revbwt_start;
	unsigned int interval_size;
	unsigned int WL_char;  // ID of the label of the last Weiner link (i.e. of the first character of W). 0=#. 1..4: ACGT. Remark: this variable can never be 0.
	
	// Right extensions
	unsigned char nright_extensions;
	unsigned char right_extension_bitmap;  // The i-th LSB is one iff the i-th character is a right extension. '#' is the 0-th LSB.
	
	// Left extensions
	unsigned char nleft_extensions;
	unsigned char left_extension_bitmap;  // The i-th LSB is one iff the i-th character is a left extension. '#' is the 0-th LSB.
	unsigned int left_ext_bwt_start[5];
	
	// Frequency of every pair of left (rows) and right (columns) extension.
	unsigned int left_right_extension_freqs[6][6];
} SLT_params_t;


static inline SLT_iterator_t_single_string *new_SLT_iterator(SLT_callback_t SLT_callback, void *intern_state, Basic_BWT_t *BBWT, unsigned int options) {
	SLT_iterator_t_single_string *SLT_iterator;
	SLT_iterator->SLT_callback=SLT_callback;
	SLT_iterator->intern_state=intern_state;
	SLT_iterator->BBWT=BBWT;
	SLT_iterator->options=options;
	return SLT_iterator;
};


static inline void free_SLT_iterator(SLT_iterator_t_single_string *SLT_iterator) {
	free(SLT_iterator);
};


void SLT_execute_iterator(SLT_iterator_t_single_string *SLT_iterator);


#endif