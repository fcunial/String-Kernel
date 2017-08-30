#ifndef indexed_DNA5_seq64_h
#define indexed_DNA5_seq64_h
#include<stdio.h>
typedef struct
{
	unsigned long long * superblock_counts;
	unsigned int * indexed_seq;
	unsigned long long length;
} indexed_DNA5_seq64_t;
extern unsigned int DNA5_char_counts_3gram[128];
extern unsigned char DNA_5_extract_table[128*3];
extern unsigned char DNA_5_alpha_trans_table[256];
extern unsigned int DNA_5_extract_suff_table[128*3];
unsigned int DNA5_extract_char64(indexed_DNA5_seq64_t * indexed_DNA5_seq64,unsigned long long charpos);
indexed_DNA5_seq64_t * new_basic_DNA5_seq64(unsigned int seqlen);
unsigned long long append_DNA5_seq64_to_opened_file(indexed_DNA5_seq64_t * indexed_DNA5_seq, FILE * file);
unsigned long long load_DNA5_seq64_from_opened_file(indexed_DNA5_seq64_t * indexed_DNA5_seq, FILE * file);
void complete_basic_DNA5_seq64(indexed_DNA5_seq64_t * indexed_DNA5_seq64);
indexed_DNA5_seq64_t * build_basic_DNA5_seq64(unsigned char * orig_seq,
	unsigned long long seqlen,unsigned long long * counts64);
void free_basic_DNA5_seq64(indexed_DNA5_seq64_t * indexed_DNA5_seq64);
inline void DNA5_set_char64(indexed_DNA5_seq64_t * indexed_DNA5_seq64,
	unsigned long long charpos,unsigned char char_val);
void DNA5_get_char_pref_counts64(unsigned long long * count,
	indexed_DNA5_seq64_t * indexed_DNA5_seq64,unsigned long long pos);
inline unsigned long long get_DNA_index_seq_size64(unsigned long long seqlen);
inline void DNA5_set_triplet_at64(indexed_DNA5_seq64_t * indexed_DNA5_seq64,
		unsigned long long pos,unsigned char * chars);
void DNA5_pack_indexed_seq_from_text64(unsigned char * orig_text,
	indexed_DNA5_seq64_t * indexed_DNA5_seq64,unsigned long long textlen);
void DNA5_joint_get_char_pref_counts64(unsigned long long * count0,unsigned long long * count1,
	indexed_DNA5_seq64_t * indexed_DNA5_seq64,unsigned long long pos0,unsigned long long pos1);

void DNA5_multipe_char_pref_counts64( 
		indexed_DNA5_seq64_t * indexed_DNA5_seq64,
		unsigned int t,
		unsigned long long * positions, 
		unsigned long long * counts);

#endif
