#include <stdlib.h>
#include "bufferedFileWriter.h"


#ifndef TWO_BIT_MASK
#define TWO_BIT_MASK 3L  // 2-bit selector
#endif
#ifndef INITIAL_BUFFER_CAPACITY
#define INITIAL_BUFFER_CAPACITY 1000  // In characters
#endif


static const unsigned char BITS_PER_LONG = sizeof(unsigned long)<<3;
static const unsigned int BYTES_PER_CHAR = sizeof(char);
static const unsigned char INITIAL_REM = BITS_PER_LONG-2;
static const unsigned long INITIAL_MASK = 3L<<INITIAL_REM;


inline void initializeBufferedFileWriter(buffered_file_writer_t *file, char *path) {
	file->size=0;
	file->capacity=INITIAL_BUFFER_CAPACITY;  // In characters
	file->buffer=(char *)malloc(file->capacity*sizeof(char));
	file->file=fopen(path,"a");
}


inline void finalizeBufferedFileWriter(buffered_file_writer_t *file) {
	if (file->size>0) {
		// Flushing the buffer one more time
		fwrite(file->buffer,BYTES_PER_CHAR,file->size,file->file);
	}
	fclose(file->file);
	free(file->buffer);
}


/**
 * Prepares $buffer$ to host $nCharacters$ additional characters.
 */
inline static void resize(unsigned int nCharacters, buffered_file_writer_t *file) {	
	if (nCharacters>file->capacity) {
		file->capacity=nCharacters<<1;
		file->buffer=(char *)realloc(file->buffer,file->capacity*BYTES_PER_CHAR);
	}
	if (file->size+nCharacters > file->capacity) {
		fwrite(file->buffer,BYTES_PER_CHAR,file->size,file->file);
		file->size=0;
	}
}


inline void writeChar(char c, buffered_file_writer_t *to) {
	resize(1,to);
	to->buffer[to->size++]=c;
}


inline void writeChars(char *from, unsigned int last, buffered_file_writer_t *to) {
	unsigned int i;
	
	resize(last+1,to);
	for (i=0; i<=last; i++) to->buffer[to->size++]=from[i];
}


inline void writeBits(unsigned long *from, unsigned int lastBit, buffered_file_writer_t *to) {
	unsigned char i, j;
	unsigned char cell, rem;
	unsigned long mask;
	
	resize(lastBit+1,to);
	cell=lastBit/BITS_PER_LONG; rem=lastBit%BITS_PER_LONG;
	for (i=0; i<cell; i++) {
		mask=1L;
		for (j=0; j<BITS_PER_LONG; j++) {
			to->buffer[to->size++]=(from[i]&mask)==0?'0':'1';
			mask<<=1;
		}
	}
	mask=1L;
	for (i=0; i<=rem; i++) {
		to->buffer[to->size++]=(from[cell]&mask)==0?'0':'1';
		mask<<=1;
	}
}


inline void writeTwoBitsReversed(unsigned long *from, unsigned int last, buffered_file_writer_t *to, char *alphabet) {
	unsigned char rem;
	int cell;
	unsigned int bit;
	unsigned long mask;
		
	resize(last+1,to);
	bit=last<<1; cell=bit/BITS_PER_LONG;
	rem=bit%BITS_PER_LONG; mask=TWO_BIT_MASK<<rem;
	while (cell>=0) {
		to->buffer[to->size++]=alphabet[(from[cell]&mask)>>rem];
		if (rem==0) { cell--; rem=INITIAL_REM; mask=INITIAL_MASK; }
		else { rem-=2; mask>>=2; }
	}
}