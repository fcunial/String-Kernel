/**
 * @author Djamal Belazzougui, Fabio Cunial
 */
#include <math.h>
#include "DNA5_Basic_BWT.h"
#include "../io/bits.h"
#include <string.h>

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

//** static / extern inline
static void computeProbabilities(BwtIndex_t *index) {
	register uint8_t i;//
	index->textLengthDNA=index->cArray[4];
	for (i=0; i<4; i++) {
		index->dnaProbabilities[i]=((double)(index->cArray[i+1]-index->cArray[i]))/index->textLength;
		index->logDnaProbabilities[i]=log(index->dnaProbabilities[i]);
	}
}

BwtIndex_t *buildBwtIndex(char *text, uint64_t length, uint64_t sharpPosition, uint32_t options) {
	register uint8_t i;//
	register uint8_t *bwt = (uint8_t *)text;//
	BwtIndex_t *bwtIndex = newBwtIndex();
	uint64_t tmpArray[4];
	
 
	// Indexing the BWT
	bwtIndex->sharpPosition = sharpPosition;

	bwtIndex->indexedBWT=build_basic_DNA5_seq(bwt,length+1,&bwtIndex->size,tmpArray);
	if (bwtIndex->indexedBWT==NULL) {
		freeBwtIndex(bwtIndex);
		return NULL;
	}

	bwtIndex->cArray[0]=0;
	bwtIndex->cArray[1]=tmpArray[0]-1;  // Since # is replaced by an A in the BWT.
	for (i=2; i<=4; i++) bwtIndex->cArray[i]=bwtIndex->cArray[i-1]+tmpArray[i-1];
	bwtIndex->textLength=length;
	computeProbabilities(bwtIndex);
	
	return bwtIndex;
}



/**
 * Remark: the procedure stores just $size$, $sharpPosition$, $textLength$ and $cArray$,
 * since the other values of $BwtIndex_t$ can be derived from them.
 */
uint64_t serializeBwtIndex(BwtIndex_t *index, const char *path) {
	register uint8_t i;//
	register uint64_t tmp;//
	register FILE *file;//
	uint64_t tmpArray[8];
	
	file=fopen(path,"w");
	if (file==NULL) return 0;
	tmpArray[0]=index->size;
	tmpArray[1]=index->sharpPosition;
	tmpArray[2]=index->textLength;
	for (i=0; i<5; i++) tmpArray[3+i]=index->cArray[i];
	tmp=8-fwrite(&tmpArray,BYTES_PER_LONG,8,file);
	if (tmp) {
		fclose(file);
		return 0;
	}

	tmp=serialize(index->indexedBWT,index->textLength,file);
	fclose(file);
	return tmp==0?0:8*BYTES_PER_LONG+tmp;
}


uint64_t deserializeBwtIndex(BwtIndex_t *index, const char *path) {
	register uint8_t i;//
	register uint64_t tmp, nAllocatedBytes;//
	register uint32_t *pointer;//
	register FILE *file;//
	uint64_t tmpArray[8];
	
	file=fopen(path,"r");
	if (file==NULL) return 0;
	tmp=8-fread(&tmpArray,BYTES_PER_LONG,8,file);
	if (tmp) {
		fclose(file);
		return 0;
	}
	index->size=tmpArray[0];
	index->sharpPosition=tmpArray[1];
	index->textLength=tmpArray[2];
	for (i=0; i<5; i++) index->cArray[i]=tmpArray[3+i];
	computeProbabilities(index);
	
	index->indexedBWT=NULL;
	nAllocatedBytes=getIndexSize(index->textLength);
	pointer=(uint32_t *)calloc(1,nAllocatedBytes);

	if (pointer==NULL) {
		fclose(file);
		return 0;
	}
	index->indexedBWT=alignIndex(pointer);
	tmp=deserialize(index->indexedBWT,index->textLength,file);
	
	fclose(file);
	return tmp==0?0:8*BYTES_PER_LONG+tmp;
}
