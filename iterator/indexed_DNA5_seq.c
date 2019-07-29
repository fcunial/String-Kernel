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
#define BITS_IN_MINIBLOCKS_PER_WORD ((BITS_PER_MINIBLOCK)*(MINIBLOCKS_PER_WORD))

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

/**
 * A \emph{sub-block} is a group of 32 consecutive miniblocks, spanning seven 32-bit  
 * words, such that the 32-th miniblock ends at the end of the seventh word. 
 * Because of this periodicity, we use sub-blocks as units of computation in procedure
 * $DNA5_get_char_pref_counts()$.
 */
#define MINIBLOCKS_PER_SUBBLOCK 32
#define WORDS_PER_SUBBLOCK 7
#define CHARS_PER_SUBBLOCK ((MINIBLOCKS_PER_SUBBLOCK)*(DNA5_chars_per_miniblock))

/**
 * Tables of constants used throughout the code.
 */
extern unsigned char DNA_5_ascii2alphabet[256];
extern unsigned int DNA5_alpha_pows[3];
extern unsigned int DNA5_char_counts_3gram[128];
extern unsigned int DNA5_char_counts_3gram_fabio[128];
extern unsigned int DNA_5_extract_suff_table[128*3];
extern unsigned int DNA_5_extract_suff_table_fabio[128*4];
extern unsigned int DNA_5_miniblock_substring_table[128*4];



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
static inline void DNA5_setMiniblock(unsigned int *restrict index, const unsigned int miniblockID, const unsigned int value) {
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
static inline unsigned int DNA5_getMiniblock(unsigned int *restrict index, unsigned int miniblockID) {
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
extern inline void DNA5_set_char(unsigned int *restrict index, unsigned int charID, unsigned char value) {
	unsigned int MINIBLOCK_ID = charID/DNA5_chars_per_miniblock;
	unsigned int OFFSET_IN_MINIBLOC = charID%DNA5_chars_per_miniblock;  // In chars
	unsigned int val = DNA5_getMiniblock(index,MINIBLOCK_ID);

	val+=DNA5_alpha_pows[OFFSET_IN_MINIBLOC]*value;  // why??????????????????????????
	DNA5_setMiniblock(index,MINIBLOCK_ID,val);
}


/**
 * Every substring $T[i..i+2]$ of length 3 is transformed into a number 
 * $25*T[i+2] + 5*T[i+1] + 1*T[i]$. At the boundary, $T$ is assumed to be concatenated to 
 * three zeros.
 *
 */
unsigned int *build_basic_DNA5_seq(unsigned char *restrict text, unsigned int textLength, unsigned int *restrict outputSize, unsigned int *restrict characterCount) {
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


/**
 * Adds to $count$ the number of occurrences of all characters in A,C,G,T inside the 
 * interval that starts from the beginning of the $fromSubblock$-th sub-block of $block$, 
 * and that ends at $textPosition$ in the text, which is assumed to belong to $block$ as
 * well (but might belong to a different sub-block).
 *
 * Remark: the computation proceeds one sub-block at a time. The counts in a sub-block are
 * stored in a single word with binary representation $C_3 C_2 C_1 C_0$, where each $C_i$ 
 * takes 8 bits and is the number of times $i$ occurs inside a sub-block. Since a sub-
 * block contains 32 miniblocks, and each miniblock corresponds to 3 positions of the 
 * text, $C_i$ can be at most 96, so 7 bits suffice.
 *
 * @param block pointer to the first sub-block in the block, i.e. excluding the header of 
 * the block.
 */
static inline void countInBlock(const unsigned int *restrict block, const unsigned int fromSubblock, unsigned int toMiniblock, unsigned int charInToMiniblock, unsigned int *restrict count) {
	const unsigned char IS_LAST_MINIBLOCK_IN_SUBBLOCK = (toMiniblock+1)%MINIBLOCKS_PER_SUBBLOCK==0;
	unsigned int i;
	unsigned int wordID, miniblock, miniblockValue;
	register unsigned int tmpWord, tmpCounts;
	register unsigned int count0, count1, count2, count3;
	
	// Occurrences in all sub-blocks before the target miniblock.
	//
	// Remark: if $toMiniblock$ is the last one in its sub-block, all the characters in
	// the miniblock are cumulated to $tmpCounts$, rather than just the characters up to
	// $charInToMiniblock$.
	count0=0; count1=0; count2=0; count3=0; tmpCounts=0;
	miniblock=fromSubblock*MINIBLOCKS_PER_SUBBLOCK;
	wordID=fromSubblock*WORDS_PER_SUBBLOCK;
	while (miniblock+MINIBLOCKS_PER_SUBBLOCK-1<=toMiniblock) {
		tmpWord=block[wordID];
		for (i=0; i<MINIBLOCKS_PER_WORD; i++) {
			miniblockValue=tmpWord&MINIBLOCK_MASK;
			tmpCounts+=DNA5_char_counts_3gram_fabio[miniblockValue];
			tmpWord>>=BITS_PER_MINIBLOCK;
		}
		tmpWord|=block[wordID+1]<<4;
		for (i=0; i<MINIBLOCKS_PER_WORD; i++) {
			miniblockValue=tmpWord&MINIBLOCK_MASK;
			tmpCounts+=DNA5_char_counts_3gram_fabio[miniblockValue];
			tmpWord>>=BITS_PER_MINIBLOCK;
		}
		tmpWord=(block[wordID+1]>>24)|(block[wordID+2]<<8);
		for (i=0; i<MINIBLOCKS_PER_WORD; i++) {
			miniblockValue=tmpWord&MINIBLOCK_MASK;
			tmpCounts+=DNA5_char_counts_3gram_fabio[miniblockValue];
			tmpWord>>=BITS_PER_MINIBLOCK;
		}
		tmpWord=(block[wordID+2]>>20)|(block[wordID+3]<<12);
		for (i=0; i<MINIBLOCKS_PER_WORD; i++) {
			miniblockValue=tmpWord&MINIBLOCK_MASK;
			tmpCounts+=DNA5_char_counts_3gram_fabio[miniblockValue];
			tmpWord>>=BITS_PER_MINIBLOCK;
		}
		tmpWord=(block[wordID+3]>>16)|(block[wordID+4]<<16);
		for (i=0; i<MINIBLOCKS_PER_WORD; i++) {
			miniblockValue=tmpWord&MINIBLOCK_MASK;
			tmpCounts+=DNA5_char_counts_3gram_fabio[miniblockValue];
			tmpWord>>=BITS_PER_MINIBLOCK;
		}
		tmpWord=(block[wordID+4]>>12)|(block[wordID+5]<<20);
		for (i=0; i<MINIBLOCKS_PER_WORD; i++) {
			miniblockValue=tmpWord&MINIBLOCK_MASK;
			tmpCounts+=DNA5_char_counts_3gram_fabio[miniblockValue];
			tmpWord>>=BITS_PER_MINIBLOCK;
		}
		tmpWord=(block[wordID+5]>>8)|(block[wordID+6]<<24);
		for (i=0; i<MINIBLOCKS_PER_WORD; i++) {
			miniblockValue=tmpWord&MINIBLOCK_MASK;
			tmpCounts+=DNA5_char_counts_3gram_fabio[miniblockValue];
			tmpWord>>=BITS_PER_MINIBLOCK;
		}
		tmpWord=block[wordID+6]>>4;
		for (i=0; i<MINIBLOCKS_PER_WORD; i++) {
			miniblockValue=tmpWord&MINIBLOCK_MASK;
			tmpCounts+=DNA5_char_counts_3gram_fabio[miniblockValue];
			tmpWord>>=BITS_PER_MINIBLOCK;
		}
		count0+=tmpCounts&0xFF;
		tmpCounts>>=8;
		count1+=tmpCounts&0xFF;
		tmpCounts>>=8;
		count2+=tmpCounts&0xFF;
		tmpCounts>>=8;
		count3+=tmpCounts&0xFF;
		tmpCounts>>=8;
		// Now $tmpCounts$ equals zero
		wordID+=WORDS_PER_SUBBLOCK;
		miniblock+=MINIBLOCKS_PER_SUBBLOCK;
	}
	if (IS_LAST_MINIBLOCK_IN_SUBBLOCK) {
		// Removing from $count$ the extra counts inside $toMiniblock$.
		tmpCounts=DNA_5_extract_suff_table_fabio[(miniblockValue<<2)+charInToMiniblock];
		count[0]+=count0-(tmpCounts&0xFF);
		tmpCounts>>=8;
		count[1]+=count1-(tmpCounts&0xFF);
		tmpCounts>>=8;
		count[2]+=count2-(tmpCounts&0xFF);
		tmpCounts>>=8;
		count[3]+=count3-(tmpCounts&0xFF);
		return;
	}
	
	// Occurrences inside the sub-block to which the target miniblock belongs.
	//
	// Remark: all characters in $toMiniblock$ are cumulated to $tmpCounts$, rather
	// than just the characters up to $charInToMiniblock$.
	tmpWord=block[wordID];
	for (i=0; i<MINIBLOCKS_PER_WORD; i++) {
		miniblockValue=tmpWord&MINIBLOCK_MASK;
		tmpCounts+=DNA5_char_counts_3gram_fabio[miniblockValue];
		if (miniblock==toMiniblock) goto countInBlock_end;
		tmpWord>>=BITS_PER_MINIBLOCK;
		miniblock++;
	}
	tmpWord|=block[wordID+1]<<4;
	for (i=0; i<MINIBLOCKS_PER_WORD; i++) {
		miniblockValue=tmpWord&MINIBLOCK_MASK;
		tmpCounts+=DNA5_char_counts_3gram_fabio[miniblockValue];
		if (miniblock==toMiniblock) goto countInBlock_end;
		tmpWord>>=BITS_PER_MINIBLOCK;
		miniblock++;
	}
	tmpWord=(block[wordID+1]>>24)|(block[wordID+2]<<8);
	for (i=0; i<MINIBLOCKS_PER_WORD; i++) {
		miniblockValue=tmpWord&MINIBLOCK_MASK;
		tmpCounts+=DNA5_char_counts_3gram_fabio[miniblockValue];
		if (miniblock==toMiniblock) goto countInBlock_end;
		tmpWord>>=BITS_PER_MINIBLOCK;
		miniblock++;
	}
	tmpWord=(block[wordID+2]>>20)|(block[wordID+3]<<12);
	for (i=0; i<MINIBLOCKS_PER_WORD; i++) {
		miniblockValue=tmpWord&MINIBLOCK_MASK;
		tmpCounts+=DNA5_char_counts_3gram_fabio[miniblockValue];
		if (miniblock==toMiniblock) goto countInBlock_end;
		tmpWord>>=BITS_PER_MINIBLOCK;
		miniblock++;
	}
	tmpWord=(block[wordID+3]>>16)|(block[wordID+4]<<16);
	for (i=0; i<MINIBLOCKS_PER_WORD; i++) {
		miniblockValue=tmpWord&MINIBLOCK_MASK;
		tmpCounts+=DNA5_char_counts_3gram_fabio[miniblockValue];
		if (miniblock==toMiniblock) goto countInBlock_end;
		tmpWord>>=BITS_PER_MINIBLOCK;
		miniblock++;
	}
	tmpWord=(block[wordID+4]>>12)|(block[wordID+5]<<20);
	for (i=0; i<MINIBLOCKS_PER_WORD; i++) {	
		miniblockValue=tmpWord&MINIBLOCK_MASK;
		tmpCounts+=DNA5_char_counts_3gram_fabio[miniblockValue];
		if (miniblock==toMiniblock) goto countInBlock_end;
		tmpWord>>=BITS_PER_MINIBLOCK;
		miniblock++;
	}
	tmpWord=(block[wordID+5]>>8)|(block[wordID+6]<<24);
	for (i=0; i<MINIBLOCKS_PER_WORD; i++) {
		miniblockValue=tmpWord&MINIBLOCK_MASK;
		tmpCounts+=DNA5_char_counts_3gram_fabio[miniblockValue];
		if (miniblock==toMiniblock) goto countInBlock_end;
		tmpWord>>=BITS_PER_MINIBLOCK;
		miniblock++;
	}	
	tmpWord=block[wordID+6]>>4;
	for (i=0; i<MINIBLOCKS_PER_WORD; i++) {
		miniblockValue=tmpWord&MINIBLOCK_MASK;
		tmpCounts+=DNA5_char_counts_3gram_fabio[miniblockValue];
		if (miniblock==toMiniblock) goto countInBlock_end;
		tmpWord>>=BITS_PER_MINIBLOCK;
		miniblock++;
	}
	
	// Removing from $tmpCounts$ the extra counts inside $toMiniblock$.
countInBlock_end:
	tmpCounts-=DNA_5_extract_suff_table_fabio[(miniblockValue<<2)+charInToMiniblock];
	count[0]+=count0+(tmpCounts&0xFF);
	tmpCounts>>=8;
	count[1]+=count1+(tmpCounts&0xFF);
	tmpCounts>>=8;
	count[2]+=count2+(tmpCounts&0xFF);
	tmpCounts>>=8;
	count[3]+=count3+(tmpCounts&0xFF);
}


/**
 * Returns the number of occurrences of all characters in A,C,G,T inside the 
 * interval that starts from the beginning of the $fromMiniblock$-th miniblock of $block$,
 * and that ends at the $charInToMiniblock$-th character of the $toMiniblock$-th miniblock
 * of $block$ (included), where $fromMiniblock$ and $toMiniblock$ are assumed to belong to
 * the same sub-block.
 *
 * @param block pointer to the first sub-block in the block, i.e. excluding the header of 
 * the block;
 * @param toMiniblock can be equal to $fromMiniblock$;
 * @param charInToMiniblock in {0,1,2};
 * @return the counts, packed in a single integer as described in $countInBlock()$.
 */
static inline unsigned int countInSubblock(const unsigned int *restrict block, unsigned int fromMiniblock, unsigned int toMiniblock, unsigned int charInToMiniblock) {
	const unsigned int LAST_BIT = (toMiniblock+1)*BITS_PER_MINIBLOCK-1;
	unsigned int i;
	unsigned int bits, wordID, bitsInWord, miniblockValue;
	register unsigned int tmpWord, tmpCounts;
	
	// Occurrences in the following miniblocks, considered in chunks of 
	// $MINIBLOCKS_PER_WORD$ miniblocks.
	tmpCounts=0; bits=fromMiniblock*BITS_PER_MINIBLOCK;
	while (bits+BITS_IN_MINIBLOCKS_PER_WORD-1<=LAST_BIT) {
		wordID=bits/DNA5_bits_per_word;
		bitsInWord=bits%DNA5_bits_per_word;
		tmpWord=block[wordID]>>bitsInWord;
		if (bitsInWord>DNA5_bits_per_word-BITS_IN_MINIBLOCKS_PER_WORD) tmpWord|=block[wordID+1]<<(DNA5_bits_per_word-bitsInWord);
		for (i=0; i<MINIBLOCKS_PER_WORD; i++) {
			miniblockValue=tmpWord&MINIBLOCK_MASK;
			tmpCounts+=DNA5_char_counts_3gram_fabio[miniblockValue];
			tmpWord>>=BITS_PER_MINIBLOCK;
		}
		bits+=BITS_IN_MINIBLOCKS_PER_WORD;		
	}
	if ((toMiniblock-fromMiniblock+1)%MINIBLOCKS_PER_WORD==0) {
		// Removing the extra counts inside $toMiniblock$.
		return tmpCounts-DNA_5_extract_suff_table_fabio[(miniblockValue<<2)+charInToMiniblock];
	}
	
	// Occurrences fewer than a word away from the beginning of the last miniblock
	while (bits<LAST_BIT) {
		wordID=bits/DNA5_bits_per_word;
		bitsInWord=bits%DNA5_bits_per_word;
		miniblockValue=block[wordID]>>bitsInWord;
		if (bitsInWord>DNA5_bits_per_word-BITS_PER_MINIBLOCK) miniblockValue|=block[wordID+1]<<(DNA5_bits_per_word-bitsInWord);
		miniblockValue&=MINIBLOCK_MASK;
		tmpCounts+=DNA5_char_counts_3gram[miniblockValue];
		bits+=BITS_PER_MINIBLOCK;
	}
	
	// Removing from $tmpCounts$ the extra counts inside $toMiniblock$.
	return tmpCounts-DNA_5_extract_suff_table_fabio[(miniblockValue<<2)+charInToMiniblock];
}


/**
 * Answers all positions that lie in the same block, using a single scan of the block.
 */
void DNA5_multipe_char_pref_counts(unsigned int *index, unsigned int *restrict textPositions, unsigned int nTextPositions, unsigned int *restrict counts) {
	unsigned int i;
	unsigned int blockID, previousBlockID, miniblockID, previousMiniblockID;
	unsigned int charInBlock, previousCharInBlock, charInMiniblock, previousCharInMiniblock;
	unsigned int wordID, row, bits, bitsInWord, miniblockValue;
	unsigned int subBlockID, previousSubBlockID;
	register unsigned int tmpCounts;
	register unsigned int count0, count1, count2, count3;
	unsigned int *block;
	
	// First position
	previousBlockID=textPositions[0]/DNA5_chars_per_block;
	previousCharInBlock=textPositions[0]%DNA5_chars_per_block;
	previousMiniblockID=previousCharInBlock/DNA5_chars_per_miniblock;
	previousCharInMiniblock=textPositions[0]%DNA5_chars_per_miniblock;
	block=&index[previousBlockID*DNA5_words_per_block];
	counts[0]=block[0]; counts[1]=block[1]; counts[2]=block[2]; counts[3]=block[3];
	countInBlock(&block[DNA5_header_size_in_words],0,previousMiniblockID,previousCharInMiniblock,counts);
	if (nTextPositions==1) return;
	previousSubBlockID=previousCharInBlock/CHARS_PER_SUBBLOCK;
	
	// Other positions
	count0=counts[0]; count1=counts[1]; count2=counts[2]; count3=counts[3];
	for (i=1; i<nTextPositions; i++) {
		row=i<<2;
		blockID=textPositions[i]/DNA5_chars_per_block;
		charInBlock=textPositions[i]%DNA5_chars_per_block;
		subBlockID=charInBlock/CHARS_PER_SUBBLOCK;
		miniblockID=charInBlock/DNA5_chars_per_miniblock;
		charInMiniblock=textPositions[i]%DNA5_chars_per_miniblock;
		if (blockID!=previousBlockID) {
			// Counting just from the beginning of $blockID$.
			block=&index[blockID*DNA5_words_per_block];
			counts[row+0]=block[0]; counts[row+1]=block[1]; 
			counts[row+2]=block[2]; counts[row+3]=block[3];
			countInBlock(&block[DNA5_header_size_in_words],0,miniblockID,charInMiniblock,&counts[row]);
			goto DNA5_multipe_char_pref_counts_nextPosition;
		}
		
		// Positions $i$ and $i-1$ lie in the same block
		block=&index[blockID*DNA5_words_per_block+DNA5_header_size_in_words];
		
		// Occurrences inside the previous miniblock
		if (previousCharInMiniblock!=2) {
			bits=previousMiniblockID*BITS_PER_MINIBLOCK;
			wordID=bits/DNA5_bits_per_word;
			bitsInWord=bits%DNA5_bits_per_word;
			miniblockValue=block[wordID]>>bitsInWord;
			if (bitsInWord>DNA5_bits_per_word-BITS_PER_MINIBLOCK) miniblockValue|=block[wordID+1]<<(DNA5_bits_per_word-bitsInWord);
			miniblockValue&=MINIBLOCK_MASK;
			if (previousMiniblockID==miniblockID) {
				tmpCounts=DNA_5_miniblock_substring_table[(miniblockValue<<2)+(previousCharInMiniblock<<1)+charInMiniblock-1];
				counts[row+0]=count0+(tmpCounts&0xFF);
				tmpCounts>>=8;
				counts[row+1]=count1+(tmpCounts&0xFF);
				tmpCounts>>=8;
				counts[row+2]=count2+(tmpCounts&0xFF);
				tmpCounts>>=8;
				counts[row+3]=count3+(tmpCounts&0xFF);
				goto DNA5_multipe_char_pref_counts_nextPosition;
			}
			else tmpCounts=DNA_5_extract_suff_table_fabio[(miniblockValue<<2)+previousCharInMiniblock];
		}
		else tmpCounts=0;
		if (subBlockID==previousSubBlockID) {
			// Occurrences inside the common sub-block
			tmpCounts+=countInSubblock(block,previousMiniblockID+1,miniblockID,charInMiniblock);
			counts[row+0]=count0+(tmpCounts&0xFF);
			tmpCounts>>=8;
			counts[row+1]=count1+(tmpCounts&0xFF);
			tmpCounts>>=8;
			counts[row+2]=count2+(tmpCounts&0xFF);
			tmpCounts>>=8;
			counts[row+3]=count3+(tmpCounts&0xFF);
			goto DNA5_multipe_char_pref_counts_nextPosition;
		}
		if (((previousMiniblockID+1)%MINIBLOCKS_PER_SUBBLOCK)!=0) {
			// Occurrences inside the previous sub-block
			tmpCounts+=countInSubblock(block,previousMiniblockID+1,(previousSubBlockID+1)*MINIBLOCKS_PER_SUBBLOCK-1,2);
		}			
		counts[row+0]=count0+(tmpCounts&0xFF);
		tmpCounts>>=8;
		counts[row+1]=count1+(tmpCounts&0xFF);
		tmpCounts>>=8;
		counts[row+2]=count2+(tmpCounts&0xFF);
		tmpCounts>>=8;
		counts[row+3]=count3+(tmpCounts&0xFF);
		// Occurrences inside the following sub-blocks
		countInBlock(block,previousSubBlockID+1,miniblockID,charInMiniblock,&counts[row]);
		
		// Next iteration
DNA5_multipe_char_pref_counts_nextPosition:
		previousBlockID=blockID;
		previousCharInBlock=charInBlock;
		previousSubBlockID=subBlockID;
		previousMiniblockID=miniblockID;
		previousCharInMiniblock=charInMiniblock;
		count0=counts[row+0]; count1=counts[row+1];
		count2=counts[row+2]; count3=counts[row+3];
	}
}