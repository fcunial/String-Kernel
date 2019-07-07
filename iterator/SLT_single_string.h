/**
 * @author Djamal Belazzougui, Fabio Cunial
 */
#ifndef SLT_single_string_h
#define SLT_single_string_h

#include "DNA5_Basic_BWT.h"

#ifndef SLT_lex_order
#define SLT_lex_order 1
#endif
#ifndef SLT_stack_trick 
#define SLT_stack_trick 2
#endif


/** 
 * The representation of a right-maximal string W returned to the callback function.
 */
typedef struct {
	unsigned int length;  // Length of W
	unsigned int bwtStart;
	unsigned int frequency;  // Number of occurrences of W in the text.
	unsigned char firstCharacter;  // The first character of W. Can only be one of the following: 1=A, 2=C, 3=G, 4=T.
	
	unsigned char nRightExtensions;  // Number of distinct characters to the right of W, including # and N.
	unsigned char rightExtensionBitmap;  // LSBs: 0=#, 1=A, 2=C, 3=G, 4=T, 5=N.
	unsigned char nLeftExtensions;  // Number of distinct characters to the left of W, including # and N.
	unsigned char leftExtensionBitmap;  // LSBs: 0=#, 1=A, 2=C, 3=G, 4=T, 5=N.
	unsigned int bwtStart_left[5];  // 0=A, 1=C, 2=G, 3=T, 4=N.
	
	// Frequency of every pair of left- (rows) and right- (columns) extension.
	unsigned int frequency_leftRight[6][6];  // 0=#, 1=A, 2=C, 3=G, 4=T, 5=N.
} RightMaximalString_t;


/** 
 * Function called by the iterator for every string it enumerates.
 * 
 * @param applicationData pointer to a memory area maintained by the program that 
 * implements the callback function. The iterator does not touch this area.
 */
typedef void (*SLT_callback_t)(const RightMaximalString_t RightMaximalString, void *applicationData);


/**
 * The iterator of right-maximal substrings that works on one input string. 
 *
 * Remark: it uses just the BWT of the forward string.
 */
typedef struct {
	SLT_callback_t SLT_callback;  // Callback function
	void *applicationData;  // Memory area managed by the callback function
	Basic_BWT_t *BBWT;
	unsigned short flags;
} UnaryIterator_t;


/**
 * Allocates a new iterator.
 */
UnaryIterator_t newIterator(SLT_callback_t SLT_callback, void *applicationData, Basic_BWT_t *BBWT, unsigned short flags);


/**
 * Runs the iterator.
 */
void run(UnaryIterator_t *SLT_iterator);


#endif