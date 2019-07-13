/**
 * @author Djamal Belazzougui, Fabio Cunial
 */
#include <stdio.h>
#include <stdlib.h>
#include "DNA5_Basic_BWT.h"

/**
 * External algorithm used for building the BWT.
 * 0=dbwt; 1=divsufsort.
 */
#ifndef CONSTRUCTION_ALGORITHM
	#define CONSTRUCTION_ALGORITHM 1
#endif
#if CONSTRUCTION_ALGORITHM == 0
	#include "../dbwt/dbwt.h"
#elif CONSTRUCTION_ALGORITHM == 1
	#include "divsufsort.h"
#endif


BwtIndex_t *newBwtIndex() {
	return (BwtIndex_t *)calloc(1,sizeof(BwtIndex_t));
}


void freeBwtIndex(BwtIndex_t *Basic_BWT) {
	if (Basic_BWT!=NULL) {
		if (Basic_BWT->indexedBWT!=NULL) {
			free_basic_DNA5_seq(Basic_BWT->indexedBWT);
			Basic_BWT->indexedBWT=NULL;
		}
		free(Basic_BWT);
	}
}


/**
 * Builds the BWT of T# from the BWT of T# built by dbwt.
 * 
 * @param text the string T (without the final sharp), of length $length$;
 * @param Basic_BWT to set the position of the sharp;
 * @return a pointer to the BWT, or NULL if construction failed.
 */
#if CONSTRUCTION_ALGORITHM == 0
static inline unsigned char *useDbwt(unsigned char *text, unsigned int length, BwtIndex_t *Basic_BWT, unsigned int options) {
	unsigned char *bwt;
	
	bwt=dbwt_bwt(text,length,&Basic_BWT->sharpPosition,options);
	bwt[Basic_BWT->sharpPosition]='A';
	return bwt;
}
#endif


/**
 * Builds the BWT of T# from the suffix array of T built by divsufsort.
 * 
 * @param text the string T (without the final sharp), of length $length$;
 * @param Basic_BWT to set the position of the sharp;
 * @return a pointer to the BWT, or NULL if construction failed.
 */
#if CONSTRUCTION_ALGORITHM == 1
static inline unsigned char *useDivsufsort(unsigned char *text, unsigned int length, BwtIndex_t *Basic_BWT) {
	unsigned int i;
	unsigned int error, textPosition;
	int *suffixArray;
	unsigned char *bwt = NULL;
	
	suffixArray=(int *)malloc(length*sizeof(int));
	error=divsufsort(text,suffixArray,length);
	if (error) {
		free(suffixArray);
		return NULL;
	}
	bwt=(unsigned char *)malloc(length+1);
	bwt[0]=text[length-1];
	for (i=0; i<length; i++) {
		textPosition=suffixArray[i];
		if (textPosition==0) {
			Basic_BWT->sharpPosition=i+1;		
			bwt[i+1]='A';
		}
		else bwt[i+1]=text[textPosition-1];
	}
	free(suffixArray);
	return bwt;
}
#endif


BwtIndex_t *buildBwtIndex(unsigned char *text, unsigned int length, unsigned int options) {
	unsigned int i;
	unsigned char *bwt;
	BwtIndex_t *bwtIndex = newBwtIndex();
	unsigned int tmpArray[4];
	
	// Building the BWT
	#if CONSTRUCTION_ALGORITHM == 0
		bwt=useDbwt(text,length,bwtIndex,options);
	#elif CONSTRUCTION_ALGORITHM == 1
		bwt=useDivsufsort(text,length,bwtIndex);
	#endif
	if (bwt==NULL) {
		freeBwtIndex(bwtIndex);
		return NULL;
	}
	
	// Building the index
	bwtIndex->indexedBWT=build_basic_DNA5_seq(bwt,length+1,&bwtIndex->size,tmpArray);
	if (bwtIndex->indexedBWT==NULL) {
		freeBwtIndex(bwtIndex);
		return NULL;
	}
	bwtIndex->cArray[0]=0;
	bwtIndex->cArray[1]=tmpArray[0]-1;  // Since # is replaced by an A in the BWT.
	for (i=2; i<=4; i++) bwtIndex->cArray[i]=bwtIndex->cArray[i-1]+tmpArray[i-1];
	bwtIndex->textLength=length;
	return bwtIndex;
}
