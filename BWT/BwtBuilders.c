#include <math.h>
#include "BwtBuilders.h"
#include "../io/bits.h"
#include "divsufsort64.h"

/**
 * Builds the BWT of T# from the suffix array of T built by divsufsort.
 * 
 * @param text the string T (without the final sharp), of length $length$;
 * @return a pointer to the BWT, or NULL if construction failed.
 */
inline uint8_t *useDivsufsort(char *text, uint64_t length, uint64_t *sharpPosition) {
  	register uint32_t error;
	register uint64_t i, textPosition;
	register uint8_t *bwt = NULL;
	register int64_t *suffixArray = NULL;
	
	suffixArray=(int64_t *)malloc(length*sizeof(int64_t));
	error=divsufsort64((uint8_t *)text,suffixArray,length);
	if (error) {
		free(suffixArray);
		return NULL;
		}
	bwt=(uint8_t *)malloc(length+1);
	bwt[0]=text[length-1];
	for (i=0; i<length; i++) {
		textPosition=suffixArray[i];
		if (textPosition==0) {
			*sharpPosition=i+1;		
			bwt[i+1]=DNA_ALPHABET[0];
		}
		else bwt[i+1]=text[textPosition-1];
	}
	free(suffixArray);
	return bwt;
}
