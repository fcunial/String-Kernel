/**
 * Interface of the one-string iterator of suffix-link tree nodes.
 */
#ifndef SLT_single_string_h
#define SLT_single_string_h

#include "DNA5_Basic_BWT.h"

#ifndef SLT_arbitrary_order
#define SLT_arbitrary_order 0
#endif
#ifndef SLT_lex_order
#define SLT_lex_order 1
#endif
#ifndef SLT_stack_trick 
#define SLT_stack_trick 2
#endif


/** 
 * Data returned by the one-string iterator to the callback function
 */
typedef struct {
	// Properties of the current right-maximal string W
	unsigned int string_depth;
	unsigned int bwt_start;
	unsigned int revbwt_start;
	unsigned int interval_size;
	unsigned char WL_char;  // ID of the label of the last Weiner link (i.e. of the first character of W). 0=#; 1..4=ACGT. Remark: this variable can never be 0.
	
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


/** 
 * Function called by the one-string iterator for every string it enumerates.
 */
typedef void (*SLT_callback_t)(const SLT_params_t SLT_params, void *intern_state);


/**
 * The one-string iterator
 */
typedef struct {
	SLT_callback_t SLT_callback;  // Callback function
	void *intern_state;  // Space handled by the function that implements the callback
	Basic_BWT_t *BBWT;
	unsigned short options;
} SLT_iterator_t_single_string;


/**
 * 
 */
SLT_iterator_t_single_string new_SLT_iterator(SLT_callback_t SLT_callback, void *intern_state, Basic_BWT_t *BBWT, unsigned short options);


/**
 *
 */
void SLT_execute_iterator(SLT_iterator_t_single_string *SLT_iterator);


#endif