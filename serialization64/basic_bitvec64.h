#ifndef basic_bitvec64_h
#define basic_bitvec64_h
#include<stdlib.h>

#ifndef bits_per_word
#define bits_per_word ((sizeof(unsigned int)*bits_per_byte))
#endif
#ifndef bits_per_byte
#define bits_per_byte 8
#endif

static inline unsigned int ismarkedbit64(unsigned long long bitpos,unsigned int * bitvec)
{
	return (bitvec[bitpos/bits_per_word]>>(bitpos%bits_per_word))&1;
};




static inline void mark_bit64(unsigned long long bitpos,unsigned int * bitvec)
{
	bitvec[bitpos/bits_per_word]|=(1<<(bitpos%bits_per_word));
};

static inline unsigned int test_and_mark_bit64(unsigned long long bitpos,unsigned int * bitvec)
{
	unsigned int bitpos_in_word=(unsigned int) (bitpos%bits_per_word);
	unsigned long long word_pos=bitpos/bits_per_word;
	unsigned int oldbit=(bitvec[word_pos]>>bitpos_in_word)&1;
	bitvec[word_pos]|=1<<bitpos_in_word;
	return oldbit;
};

static inline unsigned int * new_bitvec64(unsigned long long size)
{
	return calloc((size+bits_per_word-1)/(bits_per_word),
			bits_per_word/bits_per_byte);
};





#endif
