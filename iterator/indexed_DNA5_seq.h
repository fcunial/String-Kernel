/**
 * @author Djamal Belazzougui, Fabio Cunial
 */
#ifndef indexed_DNA5_seq_h
#define indexed_DNA5_seq_h


/**
 * 
 */
unsigned int *new_basic_DNA5_seq(unsigned int textLength, unsigned int *_output_size);


/** 
 * Builds the index on string $text$ of length $textLength$. 
 *
 * @param characterCount output array: the procedure stores here the total number of 
 * occurrences of each character (0=A, 1=C, 2=G, 3=T);
 * @param outputSize output value: the procedure stores here the size of the data 
 * structure, in bytes;
 * @return NULL if construction failed.
 */
unsigned int *build_basic_DNA5_seq(unsigned char *text, unsigned int textLength, unsigned int *outputSize, unsigned int *characterCount);


void free_basic_DNA5_seq(unsigned int *index);




void complete_basic_DNA5_seq(unsigned int *indexed_seq, unsigned int seqlen);


void DNA5_set_char(unsigned int *indexed_seq, unsigned int charpos, unsigned char char_val);

void DNA5_get_char_pref_counts(unsigned int *count, unsigned int *indexed_seq, unsigned int pos);

inline unsigned int get_DNA_index_seq_size(unsigned int textLength);

void DNA5_set_triplet_at(unsigned int *indexed_seq, unsigned int pos, unsigned char *chars);

void DNA5_pack_indexed_seq_from_text(unsigned char *orig_text, unsigned int *indexed_seq, unsigned int textlen);

void DNA5_joint_get_char_pref_counts(unsigned int *count0, unsigned int *count1, unsigned int *indexed_seq, unsigned int pos0, unsigned int pos1);

void DNA5_multipe_char_pref_counts(unsigned int *indexed_seq, unsigned int t, unsigned int *positions, unsigned int *counts);


extern unsigned char DNA_5_ascii2alphabet[256];

extern unsigned int DNA5_alpha_pows[3];

extern unsigned int DNA5_char_counts_3gram[128];

extern unsigned char DNA_5_extract_table[128*3];



extern unsigned int DNA_5_extract_suff_table[128*3];

extern unsigned int DNA5_extract_char(unsigned int *indexed_seq, unsigned int charpos);


#endif