#include<stdlib.h>
#include<stdio.h>
#include<stdint.h>
#include"indexed_DNA5_seq.h"

#define DNA5_chars_per_7bits ((3))
#define DNA5_bits_per_byte ((8))
#define DNA5_bytes_per_word ((sizeof(unsigned int)))
#define DNA5_bits_per_word (((DNA5_bytes_per_word)*(DNA5_bits_per_byte)))
#define DNA5_words_per_block ((32))
#define DNA5_bytes_per_block (((DNA5_words_per_block)*(DNA5_bytes_per_word)))
#define DNA5_bits_per_block (((DNA5_words_per_block)*(DNA5_bits_per_word)))
#define DNA5_header_size_in_words ((4))
#define DNA5_header_size_in_bytes (((DNA5_header_size_in_words)*(DNA5_bytes_per_word)))
#define DNA5_header_size_in_bits (((DNA5_header_size_in_words)*(DNA5_bits_per_word)))
#define DNA5_useful_words_per_block (((DNA5_words_per_block)-(DNA5_header_size_in_words)))
#define DNA5_useful_bytes_per_block ((DNA5_useful_words_per_block)*(DNA5_bytes_per_word))
#define DNA5_useful_bits_per_block (DNA5_useful_words_per_block*DNA5_bits_per_word)
#define DNA5_7bits_per_block (((DNA5_useful_bits_per_block)/7))
#define DNA5_chars_per_block (((DNA5_7bits_per_block)*(DNA5_chars_per_7bits)))
#define DNA5_alpha_size ((5))
#define DNA5_ceildiv(x,y) ((((x)+(y)-1)/(y)))
#define DNA5_ceilround(x,y) (((DNA5_ceildiv(x,y))*y))
#define DNA5_floordiv(x,y) ((x)/(y))
#define malloc_granularity 8 
#define bit_malloc_granularity (malloc_granularity*8)

static inline unsigned int * round_to_next_block_boundary(unsigned int * x)
{
        unsigned int pad_bytes=DNA5_bytes_per_block-
                ((uintptr_t)x-malloc_granularity)%DNA5_bytes_per_block;
        pad_bytes%=DNA5_bytes_per_block;
        pad_bytes+=malloc_granularity;
//      printf("pad bytes are %d and original add was %d mod 128\n",
//              pad_bytes,((unsigned int)x)%128);
        return (unsigned int *)((unsigned char *)x+pad_bytes);
//      return x;
};

// Return the allocation size in bytes, given the sequence 
// length in number of characters. The allocation size takes into
// account the padding needed to align to block boundaries. 
unsigned int get_DNA_index_seq_size(unsigned int seqlen)
{
        unsigned int nblocks=DNA5_floordiv(seqlen,DNA5_chars_per_block);
        unsigned int rem_chars=seqlen-DNA5_chars_per_block*nblocks;
        unsigned int rem_7_bits=DNA5_ceildiv(rem_chars,DNA5_chars_per_7bits);
        unsigned long long alloc_bits=DNA5_header_size_in_bits+
                rem_7_bits*7+(unsigned long long) nblocks*DNA5_bits_per_block;
        unsigned int malloc_units=DNA5_ceildiv(alloc_bits,bit_malloc_granularity);
//      printf("malloc units is %d\n",malloc_units);
//      printf("seqlen is %d and allocated size is %d\n",seqlen,(malloc_units-1)*malloc_granularity+DNA5_bytes_per_block);
//      printf("old allocated size was %d",DNA5_ceildiv(seqlen,DNA5_chars_per_block)*DNA5_bytes_per_block);
        return malloc_units*malloc_granularity+DNA5_bytes_per_block;
//      return DNA5_ceildiv(seqlen,DNA5_chars_per_block)*DNA5_bytes_per_block;
};

static inline void DNA5_set_7bits_at(unsigned int * indexed_seq,unsigned int pos,unsigned int val)
{
	unsigned int val_write;
	unsigned int bit_pos=7*pos;
	unsigned int word_pos=bit_pos/DNA5_bits_per_word;
	unsigned int block_pos=word_pos/DNA5_useful_words_per_block;
	unsigned int bit_pos_in_word=bit_pos%DNA5_bits_per_word;
	unsigned int real_word_pos=word_pos%DNA5_useful_words_per_block;
	real_word_pos+=DNA5_header_size_in_words+block_pos*DNA5_words_per_block;

	val&=127;
	val_write=val<<bit_pos_in_word;
	indexed_seq[real_word_pos]&=~(127<<bit_pos_in_word);
	indexed_seq[real_word_pos]|=val_write;
//	printf("or word number %d of with %d shifted by %d\n",
//		real_word_pos,val,bit_pos_in_word);
	if(bit_pos_in_word>DNA5_bits_per_word-7)
	{
		real_word_pos++;
		val_write=val>>(DNA5_bits_per_word-bit_pos_in_word);
		indexed_seq[real_word_pos]&=0xffffffff<<
			(7-(DNA5_bits_per_word-bit_pos_in_word));
		indexed_seq[real_word_pos]|=val_write;
	};
};
inline void DNA5_set_triplet_at(unsigned int * indexed_seq,unsigned int pos,unsigned char * chars)
{
	unsigned int val;
	val=DNA_5_alpha_trans_table[chars[0]]+
		DNA_5_alpha_trans_table[chars[1]]*5+
		DNA_5_alpha_trans_table[chars[2]]*25;
	DNA5_set_7bits_at(indexed_seq,pos,val);
};

static inline unsigned int DNA5_get_7bits_at(unsigned int * indexed_seq,unsigned int pos)
{
	unsigned int val_read;
	unsigned int bit_pos=7*pos;
	unsigned int word_pos=bit_pos/DNA5_bits_per_word;
	unsigned int block_pos=word_pos/DNA5_useful_words_per_block;
	unsigned int bit_pos_in_word=bit_pos%DNA5_bits_per_word;
	unsigned int real_word_pos=word_pos%DNA5_useful_words_per_block;
	real_word_pos+=DNA5_header_size_in_words+block_pos*DNA5_words_per_block;
	val_read=indexed_seq[real_word_pos]>>bit_pos_in_word;
//	printf("real_word_pos was %d,",real_word_pos);
//	printf("read_val was %d and bit_pos was %d,",val_read,bit_pos);
	if(bit_pos_in_word>DNA5_bits_per_word-7)
	{
		real_word_pos++;
		val_read|=indexed_seq[real_word_pos]<<(DNA5_bits_per_word-bit_pos_in_word);
	};
//	printf("read_val became %d\n",val_read);
	return val_read&127;
};

unsigned int DNA5_extract_char(unsigned int * indexed_seq,unsigned int charpos)
{
//	unsigned int * indexed_seq=round_to_next_block_boundary(indexed_seq0);
	unsigned int charpos_in_7bits=charpos%DNA5_chars_per_7bits;
	unsigned int pos_of_7bits=charpos/DNA5_chars_per_7bits;
	unsigned int read_val=DNA5_get_7bits_at(indexed_seq,pos_of_7bits);
//	printf("extracted char %dat pos %d\n",
//	DNA_5_extract_table[read_val*DNA5_chars_per_7bits+charpos_in_7bits],
//	charpos);
	return DNA_5_extract_table[read_val*DNA5_chars_per_7bits+charpos_in_7bits];
};

static unsigned int DNA5_alpha_pows[3]={1,5,25};
inline void DNA5_set_char(unsigned int * indexed_seq,
	unsigned int charpos,unsigned char char_val)
{
//	unsigned int * indexed_seq=round_to_next_block_boundary(indexed_seq0);
	unsigned int charpos_in_7bits=charpos%DNA5_chars_per_7bits;
	unsigned int pos_of_7bits=charpos/DNA5_chars_per_7bits;
	unsigned int val=DNA5_get_7bits_at(indexed_seq,pos_of_7bits);
//	printf("the 7 bits were %d,",val);
	val+=DNA5_alpha_pows[charpos_in_7bits]*char_val;
	DNA5_set_7bits_at(indexed_seq,pos_of_7bits,val);
//	printf("the 7 bits are set to %d,",val);
//	printf("val was incremented by %d, ",DNA5_alpha_pows[charpos_in_7bits]*char_val);
//	printf("the char was %d,",char_val);
//	printf("the pos of 7 bits was %d\n",pos_of_7bits);

};

unsigned int * new_basic_DNA5_seq(unsigned int seqlen,
	unsigned int *_output_size)
{
	unsigned int alloc_size=get_DNA_index_seq_size(seqlen);
	unsigned int * indexed_seq0=(unsigned int *)calloc(1,alloc_size);
	unsigned int * indexed_seq=round_to_next_block_boundary(indexed_seq0);
	if(indexed_seq0==NULL)
	{
		(*_output_size)=0;
		return 0;
	};
	*(((unsigned int **)indexed_seq)-1)=indexed_seq0;
//	printf("The new indexed_seq pointer is %d\n",indexed_seq);
	return indexed_seq;
};


unsigned int * build_basic_DNA5_seq(unsigned char * orig_seq,
	unsigned int seqlen,unsigned int *_output_size,unsigned int * count)
{
//	unsigned int nblocks=DNA5_ceildiv(seqlen,DNA5_chars_per_block);
//	unsigned int alloc_size=nblocks*DNA5_bytes_per_block;
	unsigned int alloc_size=get_DNA_index_seq_size(seqlen);
	unsigned int * indexed_seq0=(unsigned int *)calloc(1,alloc_size);
	unsigned int * indexed_seq=round_to_next_block_boundary(indexed_seq0);
	unsigned int i,j;
	unsigned char packed_7bits;
	unsigned char translated_char;	
	unsigned int counts[5];
	unsigned int * block_ptr=indexed_seq;
	unsigned int pos_7bits;
	*(((unsigned int **)indexed_seq)-1)=indexed_seq0;
//	printf("indexed_seq was %u and indexed_seq0 was %u\n",indexed_seq,indexed_seq0);
	
	if(indexed_seq0==NULL)
	{
		(*_output_size)=0;
		return 0;
	};
	pos_7bits=0;
	for(i=0;i<5;i++)
		counts[i]=0;
	for(i=0,pos_7bits=0;i<seqlen;i+=DNA5_chars_per_7bits,pos_7bits++)
	{
// First set the header counts
		if(pos_7bits%DNA5_7bits_per_block==0)
		{
//			printf("We are in block number %d\n",pos_7bits/DNA5_7bits_per_block);
			for(j=0;j<4;j++)
			{
				block_ptr[j]=counts[j];
//				printf("We set partial count of block of character number %d to %d\n",
//				j,counts[j]);
			};
			block_ptr+=DNA5_words_per_block;
		};
// Then pack the characters in the block
		packed_7bits=0;
		if(i+2<seqlen)
		{
			translated_char=DNA_5_alpha_trans_table[orig_seq[i+2]];		
			counts[translated_char]++;
			packed_7bits+=translated_char;
			packed_7bits*=DNA5_alpha_size;
			translated_char=DNA_5_alpha_trans_table[orig_seq[i+1]];		
			counts[translated_char]++;
			packed_7bits+=translated_char;
			packed_7bits*=DNA5_alpha_size;

		}
		else if(i+1<seqlen)
		{
			translated_char=DNA_5_alpha_trans_table[orig_seq[i+1]];		
			counts[translated_char]++;
			packed_7bits+=translated_char;
			packed_7bits*=DNA5_alpha_size;
		};
		translated_char=DNA_5_alpha_trans_table[orig_seq[i]];
		counts[translated_char]++;
		packed_7bits+=translated_char;
		DNA5_set_7bits_at(indexed_seq,pos_7bits,packed_7bits);
//		printf("packed 7 bits are set to %d\n",packed_7bits);
//		printf("original chars were %c,%c,%c\n",orig_seq[i],orig_seq[i+1],orig_seq[i+2]);
	};
	for(i=0;i<4;i++)
		count[i]=counts[i];
	(*_output_size)=alloc_size;
	return indexed_seq;
};

void free_basic_DNA5_seq(unsigned int * indexed_seq)
{
	unsigned int * real_pointer=*(((unsigned int **)indexed_seq)-1);
	free(real_pointer);
};

#define subblock_size 32

void DNA5_get_char_pref_counts(unsigned int * count,unsigned int * indexed_seq,unsigned int pos)
{
//	unsigned int * indexed_seq=round_to_next_block_boundary(indexed_seq0);
	unsigned int i;
	unsigned int block_pos=pos/DNA5_chars_per_block;
	unsigned int char_pos_in_block=pos%DNA5_chars_per_block;
	unsigned int pos_7bits_in_block=char_pos_in_block/DNA5_chars_per_7bits;
	unsigned int charpos_in_7bits=pos%DNA5_chars_per_7bits;
	unsigned int  * block=&indexed_seq[block_pos*DNA5_words_per_block];
	unsigned int word_idx=0;
	unsigned int buffered_word;
	unsigned int tmp_counts;
	unsigned int val_7bits;
	unsigned int pos7bits_idx;
//	unsigned char DNA5_char;
//	printf("char_pos_in_block=%d and pos_7bits_in_block=%d\n",char_pos_in_block,pos_7bits_in_block);
//	printf("charpos_in_7bits=%d and block_pos=%d\n",charpos_in_7bits,block_pos);
	for(i=0;i<4;i++)
	{
		count[i]=block[i];
//		printf("partial count[%d]=%d\n",i,count[i]);
	};
	word_idx=4;
// The upcoming piece of code is not very elegant. It is quite dependend on the constants. 
// But it is crucial to make it run fast, so we sacrifice elegance for speed. 
	pos7bits_idx=0;
	tmp_counts=0;
	while(pos_7bits_in_block>=pos7bits_idx+subblock_size-1)
	{
		buffered_word=block[word_idx];
		for(i=0;i<4;i++)
		{
			val_7bits=buffered_word&127;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			buffered_word>>=7;
		};		
		word_idx++;
		buffered_word|=block[word_idx]<<4;
		for(i=0;i<4;i++)
		{
			val_7bits=buffered_word&127;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			buffered_word>>=7;
		};
		buffered_word|=(block[word_idx]>>28)<<4;
		word_idx++;
		buffered_word|=block[word_idx]<<8;
		for(i=0;i<4;i++)
		{
			val_7bits=buffered_word&127;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			buffered_word>>=7;
		};
		buffered_word|=(block[word_idx]>>24)<<4;
		word_idx++;
		buffered_word|=block[word_idx]<<12;
		for(i=0;i<4;i++)
		{
			val_7bits=buffered_word&127;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			buffered_word>>=7;
		};
		buffered_word|=(block[word_idx]>>20)<<4;
		word_idx++;
		buffered_word|=block[word_idx]<<16;
		for(i=0;i<4;i++)
		{
			val_7bits=buffered_word&127;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			buffered_word>>=7;
		};
		buffered_word|=(block[word_idx]>>16)<<4;
		word_idx++;
		buffered_word|=block[word_idx]<<20;
		for(i=0;i<4;i++)
		{
			val_7bits=buffered_word&127;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			buffered_word>>=7;
		};
		buffered_word|=(block[word_idx]>>12)<<4;
		word_idx++;
		buffered_word|=block[word_idx]<<24;
		for(i=0;i<4;i++)
		{
			val_7bits=buffered_word&127;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			buffered_word>>=7;
		};
		buffered_word|=(block[word_idx]>>8)<<4;
		for(i=0;i<4;i++)
		{
			val_7bits=buffered_word&127;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			buffered_word>>=7;
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
	};
	if(pos_7bits_in_block>=pos7bits_idx)
	{
		buffered_word=block[word_idx];
		i=0;
//		printf("the word number is %d\n",word_idx);
		do 
		{
			val_7bits=buffered_word&127;
//			printf("val_7bits is %d\n",val_7bits);
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
//			printf("tmp_counts is %d\n",tmp_counts);
			if(pos7bits_idx==pos_7bits_in_block)
				goto final;
			buffered_word>>=7;
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
			val_7bits=buffered_word&127;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			if(pos7bits_idx==pos_7bits_in_block)
				goto final;
			buffered_word>>=7;
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
			val_7bits=buffered_word&127;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			if(pos7bits_idx==pos_7bits_in_block)
				goto final;
			buffered_word>>=7;
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
			val_7bits=buffered_word&127;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			if(pos7bits_idx==pos_7bits_in_block)
				goto final;
			buffered_word>>=7;
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
			val_7bits=buffered_word&127;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			if(pos7bits_idx==pos_7bits_in_block)
				goto final;
			buffered_word>>=7;
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
			val_7bits=buffered_word&127;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			if(pos7bits_idx==pos_7bits_in_block)
				goto final;
			buffered_word>>=7;
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
			val_7bits=buffered_word&127;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			if(pos7bits_idx==pos_7bits_in_block)
				goto final;
			buffered_word>>=7;
			pos7bits_idx++;
			i++;
			if(i==4) 
				break;

		}while(1);		
		buffered_word|=(block[word_idx]>>8)<<4;
		i=0;
		do 
		{
			val_7bits=buffered_word&127;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			if(pos7bits_idx==pos_7bits_in_block)
				goto final;
			i++;
			if(i==4) 
				break;
			buffered_word>>=7;
			pos7bits_idx++;
		}while(1);		
	};
	final:				
//	printf("count up to pos %d and val7 bits is %d\n",pos,val_7bits);
//	printf("The word is %d\n",block[word_idx]);
//	we add 0x02020202 to avoid underflows
	tmp_counts+=0x02020202;
	tmp_counts-=DNA_5_extract_suff_table[val_7bits+128*(2-charpos_in_7bits)];

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
//	unsigned int * indexed_seq=round_to_next_block_boundary(indexed_seq0);
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
		pos_7bits_in_block0=char_pos_in_block0/DNA5_chars_per_7bits;
		pos_7bits_in_block1=char_pos_in_block1/DNA5_chars_per_7bits;
		charpos_in_7bits0=positions[i-1]%DNA5_chars_per_7bits;
		charpos_in_7bits1=positions[i]%DNA5_chars_per_7bits;
		block=&indexed_seq[block_pos0*DNA5_words_per_block+
					DNA5_header_size_in_words];
		for(j=0;j<4;j++)
			counts[i*4+j]=counts[(i-1)*4+j];
		bit_idx=pos_7bits_in_block0*7;
		tmp_counts=0;
	//	printf("we count from pos %d until pos %d\n",pos0,pos1);
		if(charpos_in_7bits0<2)
		{
			word_idx=bit_idx/DNA5_bits_per_word;
			bit_idx_in_word=bit_idx%DNA5_bits_per_word;
			val_7bits=block[word_idx]>>bit_idx_in_word;
			if(bit_idx_in_word>DNA5_bits_per_word-7)
			{
				word_idx++;
				val_7bits|=block[word_idx]<<
					(DNA5_bits_per_word-bit_idx_in_word);
			};
			val_7bits&=127;
	//		printf("val_7bits is %d\n",val_7bits);
	
			tmp_counts+=DNA_5_extract_suff_table[val_7bits+
					128*(2-charpos_in_7bits0)];
			if(pos_7bits_in_block0==pos_7bits_in_block1)
				goto final_final;
		};
		bit_idx+=7;
		while(bit_idx+21<=pos_7bits_in_block1*7)
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
				val_7bits=buffered_word&127;
				tmp_counts+=DNA5_char_counts_3gram[val_7bits];
				j++;
				if(j==4)
					break;
				buffered_word>>=7;
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
		while(bit_idx<=pos_7bits_in_block1*7)
		{
			word_idx=bit_idx/DNA5_bits_per_word;
			bit_idx_in_word=bit_idx%DNA5_bits_per_word;
			val_7bits=block[word_idx]>>bit_idx_in_word;
			if(bit_idx_in_word>DNA5_bits_per_word-7)
			{
				word_idx++;
				val_7bits|=block[word_idx]<<
					(DNA5_bits_per_word-bit_idx_in_word);
			};
			val_7bits&=127;
			tmp_counts+=DNA5_char_counts_3gram[val_7bits];
			bit_idx+=7;
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
unsigned int append_DNA5_seq_to_opened_file(unsigned int * indexed_DNA5_seq,
		unsigned int length, FILE * file)
{
	unsigned int i;
	unsigned int written_bytes=0;
	unsigned int * block_ptr=indexed_DNA5_seq;
	unsigned int remwords;
	unsigned int rembits;
	unsigned int remchars;
	for(i=0;i+DNA5_chars_per_block<=length;i+=DNA5_chars_per_block)
	{
		fwrite(&block_ptr[4],DNA5_useful_words_per_block,DNA5_bytes_per_word,file);
		written_bytes+=DNA5_useful_bytes_per_block;
		block_ptr+=DNA5_words_per_block;
		printf("writing %d bytes with i=%d\n",DNA5_useful_bytes_per_block,i);
	};
	if(i<length)
	{
		remchars=length-i;
		rembits=((remchars+DNA5_chars_per_7bits-1)/DNA5_chars_per_7bits)*7;
		remwords=(rembits+DNA5_bits_per_word-1)/DNA5_bits_per_word;
		fwrite(&block_ptr[4],remwords,DNA5_bytes_per_word,file);
		written_bytes+=remwords*DNA5_bytes_per_word;
/*		printf("writing %d bytes with remchars=%d\n",
			remwords*DNA5_bytes_per_word,remchars);*/
	};
	return written_bytes;
};
unsigned int load_DNA5_seq_from_opened_file(unsigned int * indexed_DNA5_seq,
	unsigned int length, FILE * file)

{
	unsigned int i;
	unsigned int read_bytes=0;
	unsigned int * block_ptr=indexed_DNA5_seq;
	unsigned int remwords;
	unsigned int rembits;
	unsigned int remchars;
	for(i=0;i+DNA5_chars_per_block<=length;i+=DNA5_chars_per_block)
	{
		read_bytes+=fread(&block_ptr[4],DNA5_bytes_per_word,
			DNA5_useful_words_per_block,file)*DNA5_bytes_per_word;
		block_ptr+=DNA5_words_per_block;
	};
	if(i<length)
	{
		remchars=length-i;
		rembits=((remchars+DNA5_chars_per_7bits-1)/DNA5_chars_per_7bits)*7;
		remwords=(rembits+DNA5_bits_per_word-1)/DNA5_bits_per_word;
		read_bytes+=fread(&block_ptr[4],DNA5_bytes_per_word,remwords,file)*
			DNA5_bytes_per_word;
	};
	complete_basic_DNA5_seq(indexed_DNA5_seq,length);
	return read_bytes;
};

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
		translated_chars[0]=DNA_5_alpha_trans_table[orig_text[i]];
		translated_chars[1]=DNA_5_alpha_trans_table[orig_text[i+1]];
		translated_chars[2]=DNA_5_alpha_trans_table[orig_text[i+2]];
		packed_7bits=translated_chars[0]+DNA5_alpha_size*translated_chars[1]+
			(DNA5_alpha_size*DNA5_alpha_size)*translated_chars[2];
//		printf("write character in packed bwt %d\n",packed_7bits);
		bit_pos_in_word=bit_pos%DNA5_bits_per_word;
		word_pos=bit_pos/DNA5_bits_per_word;
		indexed_seq[word_pos]&=~(127<<bit_pos_in_word);
		val_write=packed_7bits<<bit_pos_in_word;
		indexed_seq[word_pos]|=val_write;
		if(bit_pos_in_word>DNA5_bits_per_word-7)
		{
			word_pos++;
			val_write=packed_7bits>>(DNA5_bits_per_word-bit_pos_in_word);
			indexed_seq[word_pos]&=0xffffffff<<
				(7-(DNA5_bits_per_word-bit_pos_in_word));
			indexed_seq[word_pos]|=val_write;
		};
		bit_pos+=7;
		if((bit_pos%DNA5_bits_per_block)==0)
			bit_pos+=DNA5_header_size_in_bits;
	};
	if(i<textlen)
	{
		packed_7bits=DNA_5_alpha_trans_table[orig_text[i]];
		if(i+1<textlen)
			packed_7bits+=DNA5_alpha_size*
				DNA_5_alpha_trans_table[orig_text[i+1]];
//		printf("write character in packed bwt %d\n",packed_7bits);
		bit_pos_in_word=bit_pos%DNA5_bits_per_word;
		word_pos=bit_pos/DNA5_bits_per_word;
		indexed_seq[word_pos]&=~(127<<bit_pos_in_word);
		val_write=packed_7bits<<bit_pos_in_word;
		indexed_seq[word_pos]|=val_write;
		if(bit_pos_in_word>DNA5_bits_per_word-7)
		{
			word_pos++;
			val_write=packed_7bits>>(DNA5_bits_per_word-bit_pos_in_word);
			indexed_seq[word_pos]&=0xffffffff<<
				(7-(DNA5_bits_per_word-bit_pos_in_word));
			indexed_seq[word_pos]|=val_write;
		};
	};
};
