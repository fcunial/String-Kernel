#ifndef indexed_DNA5_seq_h
#define indexed_DNA5_seq_h
extern unsigned int DNA5_char_counts_3gram[128];
extern unsigned char DNA_5_extract_table[128*3];
extern unsigned char DNA_5_alpha_trans_table[256];
extern unsigned int DNA_5_extract_suff_table[128*3];
extern unsigned int DNA5_extract_char(unsigned int * indexed_seq,unsigned int charpos);
unsigned int * new_basic_DNA5_seq(unsigned int seqlen,
	unsigned int *_output_size);
unsigned int append_DNA5_seq_to_opened_file(unsigned int * indexed_DNA5_seq,
		unsigned int length, FILE * file);
unsigned int load_DNA5_seq_from_opened_file(unsigned int * indexed_DNA5_seq,
	unsigned int length, FILE * file);
void complete_basic_DNA5_seq(unsigned int * indexed_seq,unsigned int seqlen);
unsigned int * build_basic_DNA5_seq(unsigned char * orig_seq,
	unsigned int seqlen,unsigned int *_output_size,unsigned int * count);
void free_basic_DNA5_seq(unsigned int * indexed_seq);
inline void DNA5_set_char(unsigned int * indexed_seq,
	unsigned int charpos,unsigned char char_val);
void DNA5_get_char_pref_counts(unsigned int * count,unsigned int * indexed_seq,unsigned int pos);
inline unsigned int get_DNA_index_seq_size(unsigned int seqlen);
inline void DNA5_set_triplet_at(unsigned int * indexed_seq,unsigned int pos,unsigned char * chars);
void DNA5_pack_indexed_seq_from_text(unsigned char * orig_text,
	unsigned int * indexed_seq,unsigned int textlen);
void DNA5_joint_get_char_pref_counts(unsigned int * count0,unsigned int * count1,
	unsigned int * indexed_seq,unsigned int pos0,unsigned int pos1);

void DNA5_multipe_char_pref_counts(unsigned int * indexed_seq,
		unsigned int t,
		unsigned int * positions, 
		unsigned int * counts);

#endif
