#include <stddef.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/time.h>
#include "bits.h"


#ifndef BIT_MASK
#define BIT_MASK 1L  // 1-bit selector
#endif
#ifndef TWO_BIT_MASK
#define TWO_BIT_MASK 3L  // 2-bit selector
#endif


static const unsigned char BITS_PER_LONG = sizeof(unsigned long)<<3;


void printLong(unsigned long number) {
	unsigned char i;
	unsigned long mask = 1L<<(BITS_PER_LONG-1);
	
	for (i=0; i<BITS_PER_LONG; i++) {
		printf("%c",(number&mask)==0?'0':'1');
		mask>>=1;
	}
	printf("\n");
}


char readTwoBits(unsigned long *buffer, unsigned int i) {
	unsigned int bit, cell;
	unsigned char rem;
	
	bit=i<<1; cell=bit/BITS_PER_LONG; rem=bit%BITS_PER_LONG;
	return (buffer[cell]&(TWO_BIT_MASK<<rem))>>rem;
}


void writeTwoBits(unsigned long *buffer, unsigned int i, unsigned char value) {
	unsigned int bit, cell;
	unsigned char rem;
	
	bit=i<<1; cell=bit/BITS_PER_LONG; rem=bit%BITS_PER_LONG;
	buffer[cell]&=~(TWO_BIT_MASK<<rem);
	buffer[cell]|=value<<rem;
}


char readBit(unsigned long *buffer, unsigned int i) {
	unsigned int bit, cell;
	unsigned char rem;
	
	bit=i; cell=bit/BITS_PER_LONG; rem=bit%BITS_PER_LONG;			
	return (buffer[cell]&(BIT_MASK<<rem))==0?0:1;
}


void writeBit(unsigned long *buffer, unsigned int i, unsigned char value) {
	unsigned int bit, cell;
	unsigned char rem;
	
	bit=i; cell=bit/BITS_PER_LONG; rem=bit%BITS_PER_LONG;
	buffer[cell]&=~(BIT_MASK<<rem);
	if (value==1) buffer[cell]|=BIT_MASK<<rem;
}