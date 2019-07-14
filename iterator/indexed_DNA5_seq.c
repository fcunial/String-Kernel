/**
 * @author Djamal Belazzougui, Fabio Cunial
 */
#include <stdlib.h>
#include <stdio.h>
#include "indexed_DNA5_seq.h"

#define DNA5_alphabet_size 5


#define DNA5_bits_per_byte 8
#define DNA5_bytes_per_word ((sizeof(unsigned int)))
#define DNA5_bits_per_word (((DNA5_bytes_per_word)*(DNA5_bits_per_byte)))

/**
 * Since the alphabet has size 5, rather than packing each character into a word using 3
 * bits, it is more space-efficient to encode a sequence of X consecutive characters as a 
 * block of Y bits, where $Y=ceil(log2(5^X))$, which represents the sequence as a number 
 * in base 5. We call \emph{miniblock} such a block of Y bits that encodes X characters. 
 * Function Y(X) is represented in <img src="miniblock-5.pdf">, and its first values are 
 * 3,5,7,10,12,14,17,19,21,24. We use X=3, Y=7 in the code, since this already achieves 
 * 2.3333 bits per character.
 *
 * Figure <img src="miniblock-3-14.pdf"> shows Y(X) for other alphabet sizes.
 */
#define DNA5_chars_per_miniblock 3
#define BITS_PER_MINIBLOCK 7
#define MINIBLOCK_MASK 127  // The seven LSBs set to all ones
#define MINIBLOCKS_PER_WORD ((DNA5_bits_per_word)/(BITS_PER_MINIBLOCK))

/**
 * We assume that memory is allocated in blocks called \emph{pages}, and that a page 
 * is big enough to contain a pointer to memory.
 */
#define BYTES_PER_PAGE 8
#define BITS_PER_PAGE (BYTES_PER_PAGE*8)


#define DNA5_words_per_block 32
#define DNA5_bytes_per_block (((DNA5_words_per_block)*(DNA5_bytes_per_word)))
#define DNA5_bits_per_block (((DNA5_words_per_block)*(DNA5_bits_per_word)))
#define DNA5_header_size_in_words 4
#define DNA5_header_size_in_bytes (((DNA5_header_size_in_words)*(DNA5_bytes_per_word)))
#define DNA5_header_size_in_bits (((DNA5_header_size_in_words)*(DNA5_bits_per_word)))
#define DNA5_useful_words_per_block (((DNA5_words_per_block)-(DNA5_header_size_in_words)))
#define DNA5_useful_bits_per_block (DNA5_useful_words_per_block*DNA5_bits_per_word)
#define DNA5_miniblocks_per_block (((DNA5_useful_bits_per_block)/BITS_PER_MINIBLOCK))
#define DNA5_chars_per_block (((DNA5_miniblocks_per_block)*(DNA5_chars_per_miniblock)))
#define DNA5_ceildiv(x,y) ((((x)+(y)-1)/(y)))
#define DNA5_ceilround(x,y) (((DNA5_ceildiv(x,y))*y))
#define DNA5_floordiv(x,y) ((x)/(y))






// -------------------------- UNDERSTOOD

/**
 * Writes the seven LSBs of $value$ into the $miniblockID$-th miniblock.
 * The procedure assumes that we have a sequence of 7-bit miniblocks packed consecutively
 * into words, from LSB to MSB, with a miniblock possibly straddling two words. The  
 * procedure addressed miniblocks in this abstract model, but writes them at suitable 
 * positions inside $index$, which is an indexed string. 
 *
 * Remark: the procedure assumes that, when a miniblock straddles two words, the next word
 * is not a header?!?!?!?!?!? how is this enforced???????????
 */
static inline void DNA5_setMiniblock(unsigned int *index, const unsigned int miniblockID, const unsigned int value) {
	const unsigned int BIT_ID = miniblockID*BITS_PER_MINIBLOCK;
	const unsigned int WORD_ID = BIT_ID/DNA5_bits_per_word;
	const unsigned int OFFSET_IN_WORD = BIT_ID%DNA5_bits_per_word;
	const unsigned int BLOCK_ID = WORD_ID/DNA5_useful_words_per_block;
	unsigned int wordInIndex = BLOCK_ID*DNA5_words_per_block+DNA5_header_size_in_words+(WORD_ID%DNA5_useful_words_per_block);
	unsigned int tmpValue;

	tmpValue=(value&MINIBLOCK_MASK)<<OFFSET_IN_WORD;
	index[wordInIndex]&=~(MINIBLOCK_MASK<<OFFSET_IN_WORD);
	index[wordInIndex]|=tmpValue;
	if (OFFSET_IN_WORD>DNA5_bits_per_word-BITS_PER_MINIBLOCK) {
		wordInIndex++;
		tmpValue=value>>(DNA5_bits_per_word-OFFSET_IN_WORD);
		index[wordInIndex]&=0xFFFFFFFF<<(BITS_PER_MINIBLOCK-(DNA5_bits_per_word-OFFSET_IN_WORD));
		index[wordInIndex]|=tmpValue;
	}
}


/**
 * Returns the value of the $miniblockID$-th miniblock in the seven LSBs of the result.
 * 
 * Remark: the procedure assumes that, when a miniblock straddles two words, the next word
 * is not a header?!?!?!?!?!? how is this enforced???????????
 */
static inline unsigned int DNA5_getMiniblock(unsigned int *index, unsigned int miniblockID) {
	const unsigned int BIT_ID = miniblockID*BITS_PER_MINIBLOCK;
	const unsigned int WORD_ID = BIT_ID/DNA5_bits_per_word;
	const unsigned int OFFSET_IN_WORD = BIT_ID%DNA5_bits_per_word;
	const unsigned int BLOCK_ID = WORD_ID/DNA5_useful_words_per_block;
	const unsigned int wordInIndex = BLOCK_ID*DNA5_words_per_block+DNA5_header_size_in_words+(WORD_ID%DNA5_useful_words_per_block);
	unsigned int tmpValue;
	
	tmpValue=index[wordInIndex]>>OFFSET_IN_WORD;
	if (OFFSET_IN_WORD>DNA5_bits_per_word-BITS_PER_MINIBLOCK) tmpValue|=index[wordInIndex+1]<<(DNA5_bits_per_word-OFFSET_IN_WORD);
	return tmpValue&MINIBLOCK_MASK;
}


/**
 * @return $oldPointer$ moved forward by at least one page, and so that it coincides with 
 * the beginning of a page. The value of $oldPointer$ is stored immediately before the
 * pointer returned in output.
 */
static inline unsigned int *alignIndex(unsigned int *oldPointer) {
	unsigned int *newPointer = (unsigned int *)( (unsigned char *)oldPointer+(BYTES_PER_PAGE<<1)-((unsigned int)oldPointer)%BYTES_PER_PAGE );
	*(((unsigned int **)newPointer)-1)=oldPointer;
	return newPointer;
}


/**
 * Uses the pointer written by $alignIndex()$.
 */
void free_basic_DNA5_seq(unsigned int *index) {
	free( *(((unsigned int **)index)-1) );
}


/**
 * The procedure takes into account the partially-used extra space at the beginning of the
 * index, needed to align it to pages (see $alignIndex()$).
 *
 * @param textLength in characters;
 * @return the index size in bytes.
 */
inline unsigned int getIndexSize(const unsigned int textLength) {
	const unsigned int N_BLOCKS = textLength/DNA5_chars_per_block;
	const unsigned int REMAINING_CHARS = textLength-N_BLOCKS*DNA5_chars_per_block;
	const unsigned int REMAINING_MINIBLOCKS = DNA5_ceildiv(REMAINING_CHARS,DNA5_chars_per_miniblock);
	const unsigned int SIZE_IN_BITS = N_BLOCKS*DNA5_bits_per_block+DNA5_header_size_in_bits+REMAINING_MINIBLOCKS*BITS_PER_MINIBLOCK;
	const unsigned int SIZE_IN_PAGES = DNA5_ceildiv(SIZE_IN_BITS,BITS_PER_PAGE);
	
	return (SIZE_IN_PAGES+2)*BYTES_PER_PAGE+DNA5_bytes_per_block;
}


/**
 * Sets the $charID$-th character to the two LSBs in $value$.
 * The procedure assumes that the string is a sequence of 2-bit characters stored inside
 * consecutive miniblocks.
 */
extern inline void DNA5_set_char(unsigned int *index, unsigned int charID, unsigned char value) {
	unsigned int MINIBLOCK_ID = charID/DNA5_chars_per_miniblock;
	unsigned int OFFSET_IN_MINIBLOC = charID%DNA5_chars_per_miniblock;  // In chars
	unsigned int val = DNA5_getMiniblock(index,MINIBLOCK_ID);

	val+=DNA5_alpha_pows[OFFSET_IN_MINIBLOC]*value;  // why??????????????????????????
	DNA5_setMiniblock(index,MINIBLOCK_ID,val);
}


unsigned int *build_basic_DNA5_seq(unsigned char *text, unsigned int textLength, unsigned int *outputSize, unsigned int *characterCount) {
	unsigned int *pointer = NULL;
	unsigned int *index = NULL;
	unsigned char charID, miniblock;
	unsigned int i, j;
	unsigned int miniblockID, nAllocatedBytes;
	unsigned int cumulativeCounts[5];
	
	nAllocatedBytes=getIndexSize(textLength);
	pointer=(unsigned int *)calloc(1,nAllocatedBytes);
	if (pointer==NULL) {
		*outputSize=0;
		return NULL;
	}
	index=alignIndex(pointer);
	for (i=0; i<=4; i++) cumulativeCounts[i]=0;
	miniblockID=0; pointer=index;
	for (i=0; i<textLength; i+=DNA5_chars_per_miniblock) {
		// Block header
		if (miniblockID%DNA5_miniblocks_per_block==0) {
			for (j=0; j<=3; j++) pointer[j]=cumulativeCounts[j];
			pointer+=DNA5_words_per_block;
		}
		// Block payload
		miniblock=0;
		if (i+2<textLength) {
			charID=DNA_5_ascii2alphabet[text[i+2]];		
			cumulativeCounts[charID]++;
			miniblock+=charID;
			miniblock*=DNA5_alphabet_size;
			charID=DNA_5_ascii2alphabet[text[i+1]];
			cumulativeCounts[charID]++;
			miniblock+=charID;
			miniblock*=DNA5_alphabet_size;
		}
		else if (i+1<textLength) {
			charID=DNA_5_ascii2alphabet[text[i+1]];		
			cumulativeCounts[charID]++;
			miniblock+=charID;
			miniblock*=DNA5_alphabet_size;
		}
		charID=DNA_5_ascii2alphabet[text[i]];
		cumulativeCounts[charID]++;
		miniblock+=charID;
		DNA5_setMiniblock(index,miniblockID,miniblock);
		miniblockID++;
	}
	for (i=0; i<=3; i++) characterCount[i]=cumulativeCounts[i];
	*outputSize=nAllocatedBytes;
	return index;
}








#define subblock_size 32


/**
 * // The upcoming piece of code is not very elegant. It is quite dependend on the constants. 
// But it is crucial to make it run fast, so we sacrifice elegance for speed. 
 * 
 * @param output array of counts.
 */
void DNA5_get_char_pref_counts(unsigned int *count, unsigned int *index, unsigned int textPosition) {
	unsigned int i;
	const unsigned int BLOCK_ID = textPosition/DNA5_chars_per_block;
	const unsigned int CHAR_IN_BLOCK = textPosition%DNA5_chars_per_block;
	const unsigned int MINIBLOCK_ID = CHAR_IN_BLOCK/DNA5_chars_per_miniblock;
	const unsigned int CHAR_IN_MINIBLOCK = textPosition%DNA5_chars_per_miniblock;
	const unsigned int *block = &index[BLOCK_ID*DNA5_words_per_block];
	unsigned int word_idx = 0;
	unsigned int buffered_word;
	unsigned int tmp_counts;
	unsigned int val_7bits;
	unsigned int pos7bits_idx;
	
	
	// Counts before the block
	for (i=0; i<=3; i++) count[i]=block[i];
	
	// Counts before the miniblock
	word_idx=4; pos7bits_idx=0; tmp_counts=0;
	while (pos7bits_idx+subblock_size-1<=MINIBLOCK_ID) {
		buffered_word=block[word_idx];
		for (i=0; i<MINIBLOCKS_PER_WORD; i++) {
			val_7bits=buffered_word&MINIBLOCK_MASK;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			buffered_word>>=BITS_PER_MINIBLOCK;
		}
		word_idx++;
		buffered_word|=block[word_idx]<<4;
		for (i=0; i<MINIBLOCKS_PER_WORD; i++) {
			val_7bits=buffered_word&MINIBLOCK_MASK;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			buffered_word>>=BITS_PER_MINIBLOCK;
		}
		buffered_word|=(block[word_idx]>>28)<<4;
		word_idx++;
		buffered_word|=block[word_idx]<<8;
		for(i=0;i<4;i++)
		{
			val_7bits=buffered_word&MINIBLOCK_MASK;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			buffered_word>>=BITS_PER_MINIBLOCK;
		};
		buffered_word|=(block[word_idx]>>24)<<4;
		word_idx++;
		buffered_word|=block[word_idx]<<12;
		for(i=0;i<4;i++)
		{
			val_7bits=buffered_word&MINIBLOCK_MASK;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			buffered_word>>=BITS_PER_MINIBLOCK;
		};
		buffered_word|=(block[word_idx]>>20)<<4;
		word_idx++;
		buffered_word|=block[word_idx]<<16;
		for(i=0;i<4;i++)
		{
			val_7bits=buffered_word&MINIBLOCK_MASK;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			buffered_word>>=BITS_PER_MINIBLOCK;
		};
		buffered_word|=(block[word_idx]>>16)<<4;
		word_idx++;
		buffered_word|=block[word_idx]<<20;
		for(i=0;i<4;i++)
		{
			val_7bits=buffered_word&MINIBLOCK_MASK;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			buffered_word>>=BITS_PER_MINIBLOCK;
		};
		buffered_word|=(block[word_idx]>>12)<<4;
		word_idx++;
		buffered_word|=block[word_idx]<<24;
		for(i=0;i<4;i++)
		{
			val_7bits=buffered_word&MINIBLOCK_MASK;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			buffered_word>>=BITS_PER_MINIBLOCK;
		};
		buffered_word|=(block[word_idx]>>8)<<4;
		for(i=0;i<4;i++)
		{
			val_7bits=buffered_word&MINIBLOCK_MASK;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			buffered_word>>=BITS_PER_MINIBLOCK;
		};
		if(tmp_counts&0x80808080)
		{
			for(i=0;i<4;i++)
			{
				count[i]+=tmp_counts&0xff;
				tmp_counts>>=8;
			};
		};
		word_idx++;
		pos7bits_idx+=subblock_size;
	}
	
	
	
	if(MINIBLOCK_ID>=pos7bits_idx)
	{
		buffered_word=block[word_idx];
		i=0;
//		printf("the word number is %d\n",word_idx);
		do 
		{
			val_7bits=buffered_word&MINIBLOCK_MASK;
//			printf("val_7bits is %d\n",val_7bits);
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
//			printf("tmp_counts is %d\n",tmp_counts);
			if(pos7bits_idx==MINIBLOCK_ID)
				goto final;
			buffered_word>>=BITS_PER_MINIBLOCK;
			pos7bits_idx++;
			i++;
			if(i==4) 
				break;
		}while(1);		
		word_idx++;
		buffered_word|=block[word_idx]<<4;
		i=0;
		do 
		{
			val_7bits=buffered_word&MINIBLOCK_MASK;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			if(pos7bits_idx==MINIBLOCK_ID)
				goto final;
			buffered_word>>=BITS_PER_MINIBLOCK;
			pos7bits_idx++;
			i++;
			if(i==4) 
				break;
		}while(1);		
		buffered_word|=(block[word_idx]>>28)<<4;
		word_idx++;
		buffered_word|=block[word_idx]<<8;
		i=0;
		do 
		{
			val_7bits=buffered_word&MINIBLOCK_MASK;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			if(pos7bits_idx==MINIBLOCK_ID)
				goto final;
			buffered_word>>=BITS_PER_MINIBLOCK;
			pos7bits_idx++;
			i++;
			if(i==4) 
				break;

		}while(1);		
		buffered_word|=(block[word_idx]>>24)<<4;
		word_idx++;
		buffered_word|=block[word_idx]<<12;
		i=0;
		do 
		{
			val_7bits=buffered_word&MINIBLOCK_MASK;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			if(pos7bits_idx==MINIBLOCK_ID)
				goto final;
			buffered_word>>=BITS_PER_MINIBLOCK;
			pos7bits_idx++;
			i++;
			if(i==4) 
				break;
		}while(1);		
		buffered_word|=(block[word_idx]>>20)<<4;
		word_idx++;
		buffered_word|=block[word_idx]<<16;
		i=0;
		do 
		{
			val_7bits=buffered_word&MINIBLOCK_MASK;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			if(pos7bits_idx==MINIBLOCK_ID)
				goto final;
			buffered_word>>=BITS_PER_MINIBLOCK;
			pos7bits_idx++;
			i++;
			if(i==4) 
				break;
		}while(1);		
		buffered_word|=(block[word_idx]>>16)<<4;
		word_idx++;
		buffered_word|=block[word_idx]<<20;
		i=0;
		do 
		{
			val_7bits=buffered_word&MINIBLOCK_MASK;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			if(pos7bits_idx==MINIBLOCK_ID)
				goto final;
			buffered_word>>=BITS_PER_MINIBLOCK;
			pos7bits_idx++;
			i++;
			if(i==4) 
				break;
		}while(1);		
		buffered_word|=(block[word_idx]>>12)<<4;
		word_idx++;
		buffered_word|=block[word_idx]<<24;
		i=0;
		do 
		{
			val_7bits=buffered_word&MINIBLOCK_MASK;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			if(pos7bits_idx==MINIBLOCK_ID)
				goto final;
			buffered_word>>=BITS_PER_MINIBLOCK;
			pos7bits_idx++;
			i++;
			if(i==4) 
				break;

		}while(1);		
		buffered_word|=(block[word_idx]>>8)<<4;
		i=0;
		do 
		{
			val_7bits=buffered_word&MINIBLOCK_MASK;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			if(pos7bits_idx==MINIBLOCK_ID)
				goto final;
			i++;
			if(i==4) 
				break;
			buffered_word>>=BITS_PER_MINIBLOCK;
			pos7bits_idx++;
		}while(1);		
	};
	final:				
//	printf("count up to pos %d and val7 bits is %d\n",pos,val_7bits);
//	printf("The word is %d\n",block[word_idx]);
//	we add 0x02020202 to avoid underflows
	tmp_counts+=0x02020202;
	tmp_counts-=DNA_5_extract_suff_table[val_7bits+128*(2-CHAR_IN_MINIBLOCK)];

	i=0;
	do
	{
//		printf("increment count[%d] by %d\n",i,tmp_counts&255);
		count[i]+=(tmp_counts&0xff);
		count[i]-=2;
		i++;
		if(i==4)
			break;
		tmp_counts>>=8;
	} while(1);


};
// Compute multiple t rank queries given at given positsion 
// 
void DNA5_multipe_char_pref_counts(unsigned int * indexed_seq,
		unsigned int t,
		unsigned int * positions, 
		unsigned int * counts)
{
//	unsigned int * indexed_seq=alignIndex(indexed_seq0);
	if(t<=1)
	{
		if(t==0)
			return;
		DNA5_get_char_pref_counts(&counts[0],indexed_seq,positions[0]);
		return;
	};	
	unsigned int i,j;
	unsigned int block_pos0;
	unsigned int block_pos1;
	unsigned int char_pos_in_block0,char_pos_in_block1;
	unsigned int pos_7bits_in_block0,pos_7bits_in_block1;
	unsigned int charpos_in_7bits0,charpos_in_7bits1;
	unsigned int  * block;
	unsigned int word_idx;
	unsigned int buffered_word;
	unsigned int tmp_counts;
	unsigned int val_7bits=0;
	unsigned int bit_idx;
	unsigned int bit_idx_in_word;
	DNA5_get_char_pref_counts(&counts[0],indexed_seq,positions[0]);
	block_pos0=positions[0]/DNA5_chars_per_block;
	for(i=1;i<t;i++)
	{
		block_pos1=positions[i]/DNA5_chars_per_block;
		if(block_pos1!=block_pos0)
		{
			DNA5_get_char_pref_counts(&counts[i*4],indexed_seq,positions[i]);
			block_pos0=block_pos1;
			continue;
		};
		char_pos_in_block0=positions[i-1]-block_pos0*DNA5_chars_per_block;
		char_pos_in_block1=positions[i]-block_pos0*DNA5_chars_per_block;
		pos_7bits_in_block0=char_pos_in_block0/DNA5_chars_per_miniblock;
		pos_7bits_in_block1=char_pos_in_block1/DNA5_chars_per_miniblock;
		charpos_in_7bits0=positions[i-1]%DNA5_chars_per_miniblock;
		charpos_in_7bits1=positions[i]%DNA5_chars_per_miniblock;
		block=&indexed_seq[block_pos0*DNA5_words_per_block+
					DNA5_header_size_in_words];
		for(j=0;j<4;j++)
			counts[i*4+j]=counts[(i-1)*4+j];
		bit_idx=pos_7bits_in_block0*BITS_PER_MINIBLOCK;
		tmp_counts=0;
	//	printf("we count from pos %d until pos %d\n",pos0,pos1);
		if(charpos_in_7bits0<2)
		{
			word_idx=bit_idx/DNA5_bits_per_word;
			bit_idx_in_word=bit_idx%DNA5_bits_per_word;
			val_7bits=block[word_idx]>>bit_idx_in_word;
			if(bit_idx_in_word>DNA5_bits_per_word-BITS_PER_MINIBLOCK)
			{
				word_idx++;
				val_7bits|=block[word_idx]<<
					(DNA5_bits_per_word-bit_idx_in_word);
			};
			val_7bits&=MINIBLOCK_MASK;
	//		printf("val_7bits is %d\n",val_7bits);
	
			tmp_counts+=DNA_5_extract_suff_table[val_7bits+
					128*(2-charpos_in_7bits0)];
			if(pos_7bits_in_block0==pos_7bits_in_block1)
				goto final_final;
		};
		bit_idx+=BITS_PER_MINIBLOCK;
		while(bit_idx+21<=pos_7bits_in_block1*BITS_PER_MINIBLOCK)
		{
			word_idx=bit_idx/DNA5_bits_per_word;
			bit_idx_in_word=bit_idx%DNA5_bits_per_word;
			buffered_word=block[word_idx]>>bit_idx_in_word;
			if(bit_idx_in_word>DNA5_bits_per_word-28)
			{
				word_idx++;
				buffered_word|=block[word_idx]<<
					(DNA5_bits_per_word-bit_idx_in_word);
			};
			j=0;
			do 
			{
				val_7bits=buffered_word&MINIBLOCK_MASK;
				tmp_counts+=DNA5_char_counts_3gram[val_7bits];
				j++;
				if(j==4)
					break;
				buffered_word>>=BITS_PER_MINIBLOCK;
			}while(1);
			if(tmp_counts&0x80808080)
			{
				for(j=0;j<4;j++)
				{
					counts[i*4+j]+=tmp_counts&0xff;
					tmp_counts>>=8;
				};
			};
			bit_idx+=28;
		};
		while(bit_idx<=pos_7bits_in_block1*BITS_PER_MINIBLOCK)
		{
			word_idx=bit_idx/DNA5_bits_per_word;
			bit_idx_in_word=bit_idx%DNA5_bits_per_word;
			val_7bits=block[word_idx]>>bit_idx_in_word;
			if(bit_idx_in_word>DNA5_bits_per_word-BITS_PER_MINIBLOCK)
			{
				word_idx++;
				val_7bits|=block[word_idx]<<
					(DNA5_bits_per_word-bit_idx_in_word);
			};
			val_7bits&=MINIBLOCK_MASK;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			bit_idx+=BITS_PER_MINIBLOCK;
		};
		final_final:
		// We add 0x02020202 to avoid underflows
		tmp_counts+=0x02020202;
		tmp_counts-=DNA_5_extract_suff_table[val_7bits+128*(2-charpos_in_7bits1)];
		j=0;
		do
		{
			counts[4*i+j]+=tmp_counts&0xff;
			counts[4*i+j]-=2;
			j++;
			if(j==4)
				break;
			tmp_counts>>=8;
		} while(1);

	};
};




/* FC> unused?!
void complete_basic_DNA5_seq(unsigned int * indexed_seq,unsigned int seqlen)
{
	unsigned int i,j;
	unsigned int last_block=DNA5_floordiv(seqlen,DNA5_chars_per_block);
	unsigned int char_idx;
	unsigned int * block_ptr=indexed_seq;
	unsigned int counts[4];
	for(j=0;j<4;j++)
		counts[j]=0;
	i=0;
	char_idx=0;
	do
	{
		for(j=0;j<4;j++)
			block_ptr[j]=counts[j];
		if(i>=last_block)
			break;
		i++;
		char_idx+=DNA5_chars_per_block;
		block_ptr+=DNA5_words_per_block;
		DNA5_get_char_pref_counts(counts,indexed_seq,char_idx-1);
	}while(1);
};
*/

/* FC> unused?!
void DNA5_pack_indexed_seq_from_text(unsigned char * orig_text,
	unsigned int * indexed_seq,unsigned int textlen)
{
	unsigned int i;
	unsigned int packed_7bits;
	unsigned int word_pos;
	unsigned int val_write;
	unsigned int bit_pos=DNA5_header_size_in_bits;
	unsigned int bit_pos_in_word;
	unsigned int translated_chars[3];
	for(i=0;i+2<textlen;i+=3)
	{
		translated_chars[0]=DNA_5_ascii2alphabet[orig_text[i]];
		translated_chars[1]=DNA_5_ascii2alphabet[orig_text[i+1]];
		translated_chars[2]=DNA_5_ascii2alphabet[orig_text[i+2]];
		packed_7bits=translated_chars[0]+DNA5_alphabet_size*translated_chars[1]+
			(DNA5_alphabet_size*DNA5_alphabet_size)*translated_chars[2];
//		printf("write character in packed bwt %d\n",packed_7bits);
		bit_pos_in_word=bit_pos%DNA5_bits_per_word;
		word_pos=bit_pos/DNA5_bits_per_word;
		indexed_seq[word_pos]&=~(MINIBLOCK_MASK<<bit_pos_in_word);
		val_write=packed_7bits<<bit_pos_in_word;
		indexed_seq[word_pos]|=val_write;
		if(bit_pos_in_word>DNA5_bits_per_word-BITS_PER_MINIBLOCK)
		{
			word_pos++;
			val_write=packed_7bits>>(DNA5_bits_per_word-bit_pos_in_word);
			indexed_seq[word_pos]&=0xffffffff<<
				(BITS_PER_MINIBLOCK-(DNA5_bits_per_word-bit_pos_in_word));
			indexed_seq[word_pos]|=val_write;
		};
		bit_pos+=BITS_PER_MINIBLOCK;
		if((bit_pos%DNA5_bits_per_block)==0)
			bit_pos+=DNA5_header_size_in_bits;
	};
	if(i<textlen)
	{
		packed_7bits=DNA_5_ascii2alphabet[orig_text[i]];
		if(i+1<textlen)
			packed_7bits+=DNA5_alphabet_size*
				DNA_5_ascii2alphabet[orig_text[i+1]];
//		printf("write character in packed bwt %d\n",packed_7bits);
		bit_pos_in_word=bit_pos%DNA5_bits_per_word;
		word_pos=bit_pos/DNA5_bits_per_word;
		indexed_seq[word_pos]&=~(MINIBLOCK_MASK<<bit_pos_in_word);
		val_write=packed_7bits<<bit_pos_in_word;
		indexed_seq[word_pos]|=val_write;
		if(bit_pos_in_word>DNA5_bits_per_word-BITS_PER_MINIBLOCK)
		{
			word_pos++;
			val_write=packed_7bits>>(DNA5_bits_per_word-bit_pos_in_word);
			indexed_seq[word_pos]&=0xffffffff<<
				(BITS_PER_MINIBLOCK-(DNA5_bits_per_word-bit_pos_in_word));
			indexed_seq[word_pos]|=val_write;
		};
	};
};
*/

/* FC> UNUSED?!
inline void DNA5_set_triplet_at(unsigned int *indexed_seq, unsigned int pos, unsigned char *chars) {
	unsigned int val;
	
	val=DNA_5_ascii2alphabet[chars[0]]+
		DNA_5_ascii2alphabet[chars[1]]*5+
		DNA_5_ascii2alphabet[chars[2]]*25;
	DNA5_setMiniblock(indexed_seq,pos,val);
}
*/

/* FC> UNUSED?!
unsigned int DNA5_extract_char(unsigned int * indexed_seq,unsigned int charpos)
{
//	unsigned int * indexed_seq=alignIndex(indexed_seq0);
	unsigned int charpos_in_7bits=charpos%DNA5_chars_per_miniblock;
	unsigned int pos_of_7bits=charpos/DNA5_chars_per_miniblock;
	unsigned int read_val=DNA5_getMiniblock(indexed_seq,pos_of_7bits);
//	printf("extracted char %dat pos %d\n",
//	DNA_5_extract_table[read_val*DNA5_chars_per_miniblock+charpos_in_7bits],
//	charpos);
	return DNA_5_extract_table[read_val*DNA5_chars_per_miniblock+charpos_in_7bits];
};
*/


/*
unsigned int *new_basic_DNA5_seq(unsigned int textLength, unsigned int *_output_size) {
	unsigned int alloc_size = getIndexSize(textLength);
	unsigned int *indexed_seq0 = (unsigned int *)calloc(1,alloc_size);
	unsigned int *index = alignIndex(indexed_seq0);
	if (indexed_seq0==NULL) {
		(*_output_size)=0;
		return 0;
	}
	*(((unsigned int **)index)-1)=indexed_seq0;   // ??????????????????
	return index;
}
*/