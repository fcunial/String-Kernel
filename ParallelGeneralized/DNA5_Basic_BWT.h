/** 
 * Interface of the BWT
 */
#ifndef DNA5_Basic_BWT_h
#define DNA5_Basic_BWT_h


#define Basic_bwt_no_free_text 0
#define Basic_bwt_free_text 0


#include <stdlib.h>
#include "indexed_DNA5_seq.h"


typedef struct {
	unsigned int *indexed_BWT;
	unsigned int char_base[5];
	unsigned int size;
	unsigned int primary_idx;
	unsigned int textlen;
} Basic_BWT_t;


static inline unsigned char DNA5_BWT_get_prev_char(Basic_BWT_t *Basic_BWT, unsigned int suff_idx) {
	suff_idx++;
	if (suff_idx==Basic_BWT->primary_idx) return 255;
	return DNA5_extract_char(Basic_BWT->indexed_BWT,suff_idx);
}


Basic_BWT_t *new_Basic_BWT();


void free_Basic_BWT(Basic_BWT_t *Basic_BWT);


int build_sequence_index ( unsigned char *raw_seq,
						   unsigned int seqlen, 
						   unsigned int **_indexed_seq,
						   unsigned int *_alloc_size,
						   unsigned int *char_base );

						   
Basic_BWT_t *Build_BWT_index_from_text(unsigned char *text, unsigned int textlen, unsigned int options);


int patt_count(unsigned char *P, unsigned int m, Basic_BWT_t *Basic_BWT, unsigned int _SA_interval[2]);


int Backward_step( unsigned int *in_interval, 
				   unsigned int *out_interval, 
				   unsigned char c, 
				   Basic_BWT_t *Basic_BWT );


unsigned int LF_map(unsigned int in_pos, Basic_BWT_t *Basic_BWT);


int Basic_BWT_batch_extract( unsigned int *bitvector, 
 						     unsigned int nelements, 
							 unsigned int *output_vector,
							 Basic_BWT_t *Basic_BWT );


#endif