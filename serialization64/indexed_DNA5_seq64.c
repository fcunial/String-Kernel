#include<stdlib.h>
#include<stdio.h>
#include"indexed_DNA5_seq64.h"

#ifndef DNA5_chars_per_7bits
#define DNA5_chars_per_7bits ((3))
#endif
#ifndef DNA5_bits_per_byte
#define DNA5_bits_per_byte ((8))
#endif
#ifndef DNA5_bytes_per_word
#define DNA5_bytes_per_word ((sizeof(unsigned int)))
#endif
#ifndef DNA5_bits_per_word
#define DNA5_bits_per_word (((DNA5_bytes_per_word)*(DNA5_bits_per_byte)))
#endif
#ifndef DNA5_words_per_block
#define DNA5_words_per_block ((32))
#endif
#ifndef DNA5_bytes_per_block
#define DNA5_bytes_per_block (((DNA5_words_per_block)*(DNA5_bytes_per_word)))
#endif
#ifndef DNA5_bits_per_block 
#define DNA5_bits_per_block (((DNA5_words_per_block)*(DNA5_bits_per_word)))
#endif
#ifndef DNA5_header_size_in_words
#define DNA5_header_size_in_words ((4))
#endif
#ifndef DNA5_header_size_in_bytes 
#define DNA5_header_size_in_bytes (((DNA5_header_size_in_words)*(DNA5_bytes_per_word)))
#endif
#ifndef DNA5_header_size_in_bits
#define DNA5_header_size_in_bits (((DNA5_header_size_in_words)*(DNA5_bits_per_word)))
#endif
#ifndef DNA5_useful_words_per_block
#define DNA5_useful_words_per_block (((DNA5_words_per_block)-(DNA5_header_size_in_words)))
#endif
#ifndef DNA5_useful_bytes_per_block
#define DNA5_useful_bytes_per_block ((DNA5_useful_words_per_block)*(DNA5_bytes_per_word))
#endif
#ifndef DNA5_useful_bits_per_block
#define DNA5_useful_bits_per_block (DNA5_useful_words_per_block*DNA5_bits_per_word)
#endif
#ifndef DNA5_7bits_per_block
#define DNA5_7bits_per_block (((DNA5_useful_bits_per_block)/7))
#endif
#ifndef DNA5_chars_per_block
#define DNA5_chars_per_block (((DNA5_7bits_per_block)*(DNA5_chars_per_7bits)))
#endif
#ifndef DNA5_alpha_size
#define DNA5_alpha_size ((5))
#endif
#ifndef DNA5_ceildiv
#define DNA5_ceildiv(x,y) ((((x)+(y)-1)/(y)))
#endif
#ifndef DNA5_ceilround
#define DNA5_ceilround(x,y) (((DNA5_ceildiv(x,y))*y))
#endif
#ifndef DNA5_floordiv
#define DNA5_floordiv(x,y) ((x)/(y))
#endif
#ifndef malloc_granularity
#define malloc_granularity 8 
#endif
#ifndef bit_malloc_granularity
#define bit_malloc_granularity (malloc_granularity*8)
#endif
#define DNA5_blocks_per_superblock ((1<<17))
#define DNA5_chars_per_superblock (((unsigned long long)DNA5_blocks_per_superblock*DNA5_chars_per_block))
#define DNA5_round_to_superblocks(x) (((x+DNA5_chars_per_superblock-1)/DNA5_chars_per_superblock))


static inline unsigned int * round_to_next_block_boundary64(unsigned int * x)
{
	unsigned long long pad_bytes=DNA5_bytes_per_block-
		((size_t)x-malloc_granularity)%DNA5_bytes_per_block;
	pad_bytes%=DNA5_bytes_per_block;
	pad_bytes+=malloc_granularity;
//	printf("pad bytes are %d and original add was %d mod 128\n",
//		pad_bytes,((unsigned int)x)%128);
	return (unsigned int *)((unsigned char *)x+pad_bytes);
//	return x;
};

// Return the allocation size in bytes, given the sequence 
// length in number of characters. The allocation size takes into
// account the padding needed to align to block boundaries. 
inline unsigned long long get_DNA_index_seq_size64(unsigned long long seqlen)
{
	unsigned long nblocks=DNA5_floordiv(seqlen,DNA5_chars_per_block);
	unsigned int rem_chars=seqlen-DNA5_chars_per_block*nblocks;
	unsigned int rem_7_bits=DNA5_ceildiv(rem_chars,DNA5_chars_per_7bits);
	unsigned long long alloc_bits=DNA5_header_size_in_bits+
		rem_7_bits*7+nblocks*DNA5_bits_per_block;
	unsigned long long malloc_units=DNA5_ceildiv(alloc_bits,bit_malloc_granularity);
//	printf("malloc units is %d\n",malloc_units);
//	printf("seqlen is %d and allocated size is %d\n",seqlen,(malloc_units-1)*malloc_granularity+DNA5_bytes_per_block);
//	printf("old allocated size was %d",DNA5_ceildiv(seqlen,DNA5_chars_per_block)*DNA5_bytes_per_block);
	return malloc_units*malloc_granularity+DNA5_bytes_per_block;
//	return DNA5_ceildiv(seqlen,DNA5_chars_per_block)*DNA5_bytes_per_block;
};

static inline void DNA5_set_7bits_at64(unsigned int * indexed_seq,
		unsigned int long long pos,unsigned int val)
{
	unsigned int val_write;
	unsigned long long bit_pos=7*pos;
	unsigned long long word_pos=bit_pos/DNA5_bits_per_word;
	unsigned long long block_pos=word_pos/DNA5_useful_words_per_block;
	unsigned long long bit_pos_in_word=bit_pos%DNA5_bits_per_word;
	unsigned long long real_word_pos=word_pos%DNA5_useful_words_per_block;
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
inline void _DNA5_set_triplet_at64(unsigned int * indexed_seq,unsigned long long pos,unsigned char * chars)
{
	unsigned int val;
	val=DNA_5_alpha_trans_table[chars[0]]+
		DNA_5_alpha_trans_table[chars[1]]*5+
		DNA_5_alpha_trans_table[chars[2]]*25;
	DNA5_set_7bits_at64(indexed_seq,pos,val);
};
inline void DNA5_set_triplet_at64(indexed_DNA5_seq64_t * indexed_seq64,
		unsigned long long pos,unsigned char * chars)
{
	_DNA5_set_triplet_at64(indexed_seq64->indexed_seq,pos,chars); 
};

static inline unsigned int DNA5_get_7bits_at64(unsigned int * indexed_seq,unsigned long long pos)
{
	unsigned int val_read;
	unsigned long long bit_pos=7*pos;
	unsigned long long word_pos=bit_pos/DNA5_bits_per_word;
	unsigned long long block_pos=word_pos/DNA5_useful_words_per_block;
	unsigned long long bit_pos_in_word=bit_pos%DNA5_bits_per_word;
	unsigned long long real_word_pos=word_pos%DNA5_useful_words_per_block;
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

unsigned int DNA5_extract_char64(indexed_DNA5_seq64_t * indexed_DNA5_seq64,
			unsigned long long charpos)
{
//	unsigned int * indexed_seq=round_to_next_block_boundary(indexed_seq0);
	unsigned int charpos_in_7bits=charpos%DNA5_chars_per_7bits;
	unsigned long long pos_of_7bits=charpos/DNA5_chars_per_7bits;
	unsigned int read_val=DNA5_get_7bits_at64(indexed_DNA5_seq64->indexed_seq,pos_of_7bits);
//	printf("extracted char %dat pos %d\n",
//	DNA_5_extract_table[read_val*DNA5_chars_per_7bits+charpos_in_7bits],
//	charpos);
	return DNA_5_extract_table[read_val*DNA5_chars_per_7bits+charpos_in_7bits];
};

static unsigned int DNA5_alpha_pows[3]={1,5,25};
inline void _DNA5_set_char64(unsigned int * indexed_seq,
	unsigned long long charpos,unsigned char char_val)
{
//	unsigned int * indexed_seq=round_to_next_block_boundary(indexed_seq0);
	unsigned long long charpos_in_7bits=charpos%DNA5_chars_per_7bits;
	unsigned long long pos_of_7bits=charpos/DNA5_chars_per_7bits;
	unsigned long long val=DNA5_get_7bits_at64(indexed_seq,pos_of_7bits);
//	printf("the 7 bits were %d,",val);
	val+=DNA5_alpha_pows[charpos_in_7bits]*char_val;
	DNA5_set_7bits_at64(indexed_seq,pos_of_7bits,val);
//	printf("the 7 bits are set to %d,",val);
//	printf("val was incremented by %d, ",DNA5_alpha_pows[charpos_in_7bits]*char_val);
//	printf("the char was %d,",char_val);
//	printf("the pos of 7 bits was %d\n",pos_of_7bits);

};
inline void DNA5_set_char64(indexed_DNA5_seq64_t * indexed_seq64,
	unsigned long long charpos,unsigned char char_val)
{
	_DNA5_set_char64(indexed_seq64->indexed_seq,charpos,char_val); 
};

indexed_DNA5_seq64_t * new_basic_DNA5_seq64(unsigned int seqlen)
{
	unsigned int alloc_size=get_DNA_index_seq_size64(seqlen);
	unsigned int * indexed_seq0;
	unsigned int nsuperblocks=DNA5_round_to_superblocks(seqlen);
	indexed_DNA5_seq64_t * indexed_DNA5_seq64=(indexed_DNA5_seq64_t*)
				malloc(sizeof(indexed_DNA5_seq64_t));
	if(indexed_DNA5_seq64==0)
	{
//		(*_output_size)=0;
		return 0;
	};
	indexed_DNA5_seq64->superblock_counts=
		(unsigned long long *)calloc(nsuperblocks*4,
					sizeof(unsigned long long));
	if(indexed_DNA5_seq64->superblock_counts==0)
	{
		free(indexed_DNA5_seq64);
		return 0;

	};
	indexed_seq0=(unsigned int *)calloc(1,alloc_size);
	if(indexed_seq0==NULL)
	{
		free(indexed_DNA5_seq64);
		free(indexed_DNA5_seq64->superblock_counts);
//		(*_output_size)=0;
		return 0;
	};

	indexed_DNA5_seq64->indexed_seq=round_to_next_block_boundary64(indexed_seq0);
	*(((unsigned int **)indexed_DNA5_seq64->indexed_seq)-1)=indexed_seq0;
	indexed_DNA5_seq64->length=seqlen;
	return indexed_DNA5_seq64;
};


indexed_DNA5_seq64_t * build_basic_DNA5_seq64(unsigned char * orig_seq,
	unsigned long long seqlen,unsigned long long * counts64)
{
	unsigned long long i;
	unsigned int j,k;
	unsigned char packed_7bits;
	unsigned char translated_char;	
	unsigned int counts[5];
	unsigned int * block_ptr;
	unsigned long pos_7bits;
	unsigned long alloc_size=get_DNA_index_seq_size64(seqlen);
	unsigned int * indexed_seq0;
	unsigned int nsuperblocks=DNA5_round_to_superblocks(seqlen);
	unsigned int superblock_idx;
	unsigned int curr_superblock_len;
	indexed_DNA5_seq64_t * indexed_DNA5_seq64=(indexed_DNA5_seq64_t*)
				malloc(sizeof(indexed_DNA5_seq64_t));
//	printf("the number of superblocks is %d\n",nsuperblocks);
//	printf("block size is %u superblock size is %llu and blocks per superblock is %u\n",
//		DNA5_chars_per_block,DNA5_chars_per_superblock,DNA5_blocks_per_superblock);
	if(indexed_DNA5_seq64==0)
	{
//		(*_output_size)=0;
		return 0;
	};
	indexed_DNA5_seq64->superblock_counts=
		(unsigned long long *)calloc(nsuperblocks*4,
			sizeof(unsigned long long));
	if(indexed_DNA5_seq64->superblock_counts==0)
	{
		free(indexed_DNA5_seq64);
		return 0;

	};
	indexed_seq0=(unsigned int *)calloc(1,alloc_size);
	if(indexed_seq0==NULL)
	{
		free(indexed_DNA5_seq64);
		free(indexed_DNA5_seq64->superblock_counts);
//		(*_output_size)=0;
		return 0;
	};
	indexed_DNA5_seq64->indexed_seq=round_to_next_block_boundary64(indexed_seq0);
	*(((unsigned int **)indexed_DNA5_seq64->indexed_seq)-1)=indexed_seq0;
	indexed_DNA5_seq64->length=seqlen;
	block_ptr=indexed_DNA5_seq64->indexed_seq;
	for(j=0;j<4;j++)
		indexed_DNA5_seq64->superblock_counts[j]=0;
	i=0;
	pos_7bits=0;
	for(superblock_idx=0;superblock_idx<nsuperblocks;superblock_idx++)
	{
		if(superblock_idx==nsuperblocks-1)
			curr_superblock_len=seqlen%DNA5_chars_per_superblock;
		else
			curr_superblock_len=DNA5_chars_per_superblock;
		printf("current super block length is %d\n",
			curr_superblock_len);
		for(j=0;j<5;j++)
			counts[j]=0;
		for(k=0;k<curr_superblock_len;
			i+=DNA5_chars_per_7bits,k+=DNA5_chars_per_7bits,pos_7bits++)
		{
// First set the header counts
			if(pos_7bits%DNA5_7bits_per_block==0)
			{
//				printf("We are in block number %d\n",pos_7bits/DNA5_7bits_per_block);
				for(j=0;j<4;j++)
				{
					block_ptr[j]=counts[j];
//					printf("We set partial count of block of character number %d to %d\n",
//					j,counts[j]);
				};
				block_ptr+=DNA5_words_per_block;
			};
// Then	pack the characters in the block
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
//			printf("value of i is %d and char is %c ",i,orig_seq[i]);
			translated_char=DNA_5_alpha_trans_table[orig_seq[i]];
			counts[translated_char]++;
			packed_7bits+=translated_char;
			DNA5_set_7bits_at64(indexed_DNA5_seq64->indexed_seq,
				pos_7bits,packed_7bits);
//			printf("packed 7 bits are set to %d\n",packed_7bits);
//			printf("original chars were %c,%c,%c\n",orig_seq[i],orig_seq[i+1],orig_seq[i+2]);
		};	
		if(superblock_idx<nsuperblocks-1)
			for(j=0;j<4;j++)
				indexed_DNA5_seq64->
					superblock_counts[(superblock_idx+1)*4+j]=
					indexed_DNA5_seq64->superblock_counts
						[superblock_idx*4+j]+counts[j];
	}
	for(j=0;j<4;j++)
		counts64[j]=indexed_DNA5_seq64->superblock_counts[
				(nsuperblocks-1)*4+j]+counts[j];
//	(*_output_size)=alloc_size;
	return indexed_DNA5_seq64;
};	

void free_basic_DNA5_seq64(indexed_DNA5_seq64_t * indexed_DNA5_seq64)
{
	unsigned int * real_indexed_seq_pointer=
		*(((unsigned int **)indexed_DNA5_seq64->indexed_seq)-1);
	free(real_indexed_seq_pointer);
	free(indexed_DNA5_seq64->superblock_counts);
	free(indexed_DNA5_seq64);
};

#define subblock_size 32

void DNA5_get_char_pref_counts64(unsigned long long * count64,
	indexed_DNA5_seq64_t * indexed_DNA5_seq64,unsigned long long pos)
{
//	unsigned int * indexed_seq=round_to_next_block_boundary(indexed_seq0);
	unsigned int i;
	unsigned int * indexed_seq=indexed_DNA5_seq64->indexed_seq;
	unsigned long long block_pos=pos/DNA5_chars_per_block;
	unsigned long long word_idx=0;
	unsigned long long superblock_pos=block_pos/DNA5_blocks_per_superblock;
	unsigned int char_pos_in_block=pos%DNA5_chars_per_block;
	unsigned int pos_7bits_in_block=char_pos_in_block/DNA5_chars_per_7bits;
	unsigned int charpos_in_7bits=pos%DNA5_chars_per_7bits;
	unsigned int  * block=&indexed_seq[block_pos*DNA5_words_per_block];
	unsigned int buffered_word;
	unsigned int tmp_counts;
	unsigned int val_7bits;
	unsigned int pos7bits_idx;
	unsigned int count[4];
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
	for(i=0;i<4;i++)
		count64[i]=count[i]+indexed_DNA5_seq64->
			superblock_counts[i+4*superblock_pos];
};
// Compute multiple t rank queries given at given posisions 

void DNA5_multipe_char_pref_counts64(indexed_DNA5_seq64_t * indexed_DNA5_seq64,
		unsigned int t,
		unsigned long long * positions, 
		unsigned long long * counts64)
{
//	unsigned int * indexed_seq=round_to_next_block_boundary(indexed_seq0);
	if(t<=1)
	{
		if(t==0)
			return;
		DNA5_get_char_pref_counts64(counts64,indexed_DNA5_seq64,positions[0]);
//		printf("Do a single counting with pos %d\n",positions[0]);
//		printf("number of c is %d\n",counts64[1]);
		return;
	};
	unsigned int i,j;
	unsigned int * indexed_seq=indexed_DNA5_seq64->indexed_seq;
	unsigned long long word_idx;
	unsigned long long block_pos0;
	unsigned long long block_pos1;
	unsigned int char_pos_in_block0,char_pos_in_block1;
	unsigned int pos_7bits_in_block0,pos_7bits_in_block1;
	unsigned int charpos_in_7bits0,charpos_in_7bits1;
	unsigned int  * block;
	unsigned int buffered_word;
	unsigned int tmp_counts;
	unsigned int val_7bits=0;
	unsigned int bit_idx;
	unsigned int bit_idx_in_word;
	unsigned int counts[4];

	DNA5_get_char_pref_counts64(&counts64[0],indexed_DNA5_seq64,positions[0]);
	block_pos0=positions[0]/DNA5_chars_per_block;
	for(i=1;i<t;i++)
	{
		block_pos1=positions[i]/DNA5_chars_per_block;
		if(block_pos1!=block_pos0)
		{
			DNA5_get_char_pref_counts64(&counts64[i*4],indexed_DNA5_seq64,positions[i]);
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
			counts[j]=0;
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
					counts[j]+=tmp_counts&0xff;
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
			counts[j]+=tmp_counts&0xff;
			counts[j]-=2;
			j++;
			if(j==4)
				break;
			tmp_counts>>=8;
		} while(1);
		for(j=0;j<4;j++)
			counts64[4*i+j]=counts64[4*(i-1)+j]+counts[j];
	};
};

void complete_basic_DNA5_seq64(indexed_DNA5_seq64_t * indexed_DNA5_seq64)
{
	unsigned int i,j;
	unsigned int superblock_idx;
	unsigned long long seqlen=indexed_DNA5_seq64->length;
	unsigned int char_idx;
	unsigned int * block_ptr=indexed_DNA5_seq64->indexed_seq;
	unsigned long long counts64[4];
	unsigned long long old_counts64[4];
	unsigned int counts[4];
	unsigned int nsuperblocks=DNA5_round_to_superblocks(seqlen);
	unsigned int curr_superblock_size;
	for(j=0;j<4;j++)
		indexed_DNA5_seq64->superblock_counts[j]=0;
	for(j=0;j<4;j++)
		counts64[j]=0;
	char_idx=0;
	for(superblock_idx=0;superblock_idx<nsuperblocks;superblock_idx++)
	{
		for(j=0;j<4;j++)
			indexed_DNA5_seq64->superblock_counts[j+4*superblock_idx]=counts64[j];
		for(j=0;j<4;j++)
		{
			counts[j]=0;
			old_counts64[j]=counts64[j];
		};
		if(superblock_idx<nsuperblocks-1)
			curr_superblock_size=DNA5_chars_per_superblock;
		else
			curr_superblock_size=(indexed_DNA5_seq64->length+
					DNA5_chars_per_superblock-1)%
					DNA5_chars_per_superblock+1;
		for(i=0;i<curr_superblock_size;i+=DNA5_chars_per_block)
		{
		
			for(j=0;j<4;j++)
				block_ptr[j]=counts[j];
			char_idx+=DNA5_chars_per_block;
			block_ptr+=DNA5_words_per_block;
			DNA5_get_char_pref_counts64(counts64,indexed_DNA5_seq64,char_idx-1);
			for(j=0;j<4;j++)
				counts[j]=counts64[j]-old_counts64[j];
		}
	};
}
unsigned long long append_DNA5_seq64_to_opened_file(indexed_DNA5_seq64_t * indexed_DNA5_seq, FILE * file)
{
	unsigned long long i;
	unsigned long long written_bytes=0;
	unsigned int * block_ptr=indexed_DNA5_seq->indexed_seq;
	unsigned int remwords;
	unsigned int rembits;
	unsigned int remchars;
	for(i=0;i+DNA5_chars_per_block<=indexed_DNA5_seq->length;i+=DNA5_chars_per_block)
	{
		written_bytes+=fwrite(&block_ptr[4],DNA5_bytes_per_word,
			DNA5_useful_words_per_block,file)*DNA5_bytes_per_word;
		block_ptr+=DNA5_words_per_block;
	};
	if(i<indexed_DNA5_seq->length)
	{
		remchars=indexed_DNA5_seq->length-i;
		rembits=((remchars+DNA5_chars_per_7bits-1)/DNA5_chars_per_7bits)*7;
		remwords=(rembits+DNA5_bits_per_word-1)/DNA5_bits_per_word;
		written_bytes+=fwrite(&block_ptr[4],DNA5_bytes_per_word,remwords,file)*
				DNA5_bytes_per_word;
	};
	return written_bytes;
};
unsigned long long load_DNA5_seq64_from_opened_file(indexed_DNA5_seq64_t * indexed_DNA5_seq, FILE * file)
{
	unsigned long long i;
	unsigned long long read_bytes=0;
	unsigned int * block_ptr=indexed_DNA5_seq->indexed_seq;
	unsigned int remwords;
	unsigned int rembits;
	unsigned int remchars;
	for(i=0;i+DNA5_chars_per_block<=indexed_DNA5_seq->length;i+=DNA5_chars_per_block)
	{
		read_bytes+=fread(&block_ptr[4],DNA5_bytes_per_word,
			DNA5_useful_words_per_block,file)*DNA5_bytes_per_word;
		block_ptr+=DNA5_words_per_block;
	};
	if(i<indexed_DNA5_seq->length)
	{
		remchars=indexed_DNA5_seq->length-i;
		rembits=((remchars+DNA5_chars_per_7bits-1)/DNA5_chars_per_7bits)*7;
		remwords=(rembits+DNA5_bits_per_word-1)/DNA5_bits_per_word;
		read_bytes+=fread(&block_ptr[4],DNA5_bytes_per_word,remwords,file)*
			DNA5_bytes_per_word;
	};
	complete_basic_DNA5_seq64(indexed_DNA5_seq);
	return read_bytes;
};

/*
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
};*/

