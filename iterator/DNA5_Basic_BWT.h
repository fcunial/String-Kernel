/**
 * @author Djamal Belazzougui, Fabio Cunial
 */
#ifndef DNA5_Basic_BWT_h
#define DNA5_Basic_BWT_h


#define Basic_bwt_no_free_text 0
#define Basic_bwt_free_text 0


#include <stdlib.h>
#include "indexed_DNA5_seq.h"


/**
 * A simple BWT index that supports just rank and access.
 */
typedef struct {
	unsigned int *indexedBWT;
	unsigned int size;  // Size of $indexedBWT$
	
	unsigned int cArray[5];  // C array. 0=A, 1=C, 2=G, 3=T, 4=N.
	unsigned int sharpPosition;  // Position of the sharp in the BWT
	unsigned int textLength;  // Length of the text, excluding the sharp.
} BwtIndex_t;


/**
 * Allocates the memory for the index, without creating it.
 */
BwtIndex_t *newBwtIndex();


void freeBwtIndex(BwtIndex_t *bwtIndex);


/**
 * Creates the index.
 * 
 * @param options for dbwt.
 */
BwtIndex_t *buildBwtIndex(unsigned char *text, unsigned int textLength, unsigned int options);


#endif