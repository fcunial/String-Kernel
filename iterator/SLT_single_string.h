/**
 * Iterator of all right-maximal substrings of one input string, based just on the BWT of 
 * the forward string (i.e. it does not use a bidirectional index).
 *
 * @author Djamal Belazzougui, Fabio Cunial
 */
#ifndef SLT_single_string_h
#define SLT_single_string_h

#include "DNA5_Basic_BWT.h"

/** 
 * Order in which nodes are pushed on the iterator stack.
 * 0: no specification; 
 * 1: no specification, but with the stack trick; 
 * 2: lexicographic, without the stack trick.
 */
#ifndef TRAVERSAL_ORDER
#define TRAVERSAL_ORDER 1
#endif

/** 
 * A substring is considered right- (respectively, left-) maximal iff it is followed 
 * (respectively, preceded) by:
 * 0: at least two distinct characters in {#,A,C,G,T,N};
 * 1: at least two distinct characters in {#,A,C,G,T,N}, or at least two Ns (i.e. any two 
 * occurrences of N are considered as distinct characters);
 * 2: at least two distinct characters in {A,C,G,T}.
 */
#ifndef TRAVERSAL_MAXIMALITY
#define TRAVERSAL_MAXIMALITY 0
#endif


/** 
 * The representation of a right-maximal string W sent to the callback function.
 */
typedef struct {
	uint64_t length;  // Length of W
	uint64_t bwtStart;
	uint64_t frequency;  // Number of occurrences of W in the text.
	uint8_t firstCharacter;  // The first character of W. Can only be one of the following: 1=A, 2=C, 3=G, 4=T.
	
	uint8_t nRightExtensions;  // Number of distinct characters to the right of W, including # and N.
	uint8_t rightExtensionBitmap;  // LSBs: 0=#, 1=A, 2=C, 3=G, 4=T, 5=N.
	uint8_t nLeftExtensions;  // Number of distinct characters to the left of W, including # and N.
	uint8_t leftExtensionBitmap;  // LSBs: 0=#, 1=A, 2=C, 3=G, 4=T, 5=N.
	uint64_t bwtStart_left[5];  // 0=A, 1=C, 2=G, 3=T, 4=N.
	
	// Frequency of every pair of left- (rows) and right- (columns) extension.
	uint64_t frequency_leftRight[6][6];  // 0=#, 1=A, 2=C, 3=G, 4=T, 5=N.
} RightMaximalString_t;


/** 
 * Callback function issued by the iterator on every string it enumerates.
 * 
 * @param applicationData pointer to a memory area maintained by the program that 
 * implements the callback function. The iterator does not touch this area.
 */
typedef void (*SLT_callback_t)(const RightMaximalString_t RightMaximalString, void *applicationData);


typedef struct {
	SLT_callback_t SLT_callback;  // Callback function
	void *applicationData;  // Memory area managed by the callback function
	BwtIndex_t *BBWT;  // The string is assumed to contain at least one character, followed by the sharp.
	uint64_t maxLength;  // Maximum length of a substring to be enumerated
} UnaryIterator_t;


UnaryIterator_t newIterator(SLT_callback_t SLT_callback, void *applicationData, BwtIndex_t *BBWT, uint64_t maxLength);


void iterate(UnaryIterator_t *SLT_iterator);


#endif