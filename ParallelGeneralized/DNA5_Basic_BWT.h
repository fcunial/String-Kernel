#ifndef DNA5_Basic_BWT_h
#define DNA5_Basic_BWT_h
#include"indexed_DNA5_seq.h"
#include<stdlib.h>
#define Basic_bwt_no_free_text 0
#define Basic_bwt_free_text 1


typedef struct 
{
	unsigned int * indexed_BWT;
	unsigned int char_base[5];
	unsigned int size;
	unsigned int primary_idx;
	unsigned int textlen;
} Basic_BWT_t;
static inline unsigned char DNA5_BWT_get_prev_char(Basic_BWT_t * Basic_BWT,unsigned int suff_idx)
{
	unsigned char c;
	suff_idx++;	
	c=DNA5_extract_char(Basic_BWT->indexed_BWT,suff_idx);
	return suff_idx==Basic_BWT->primary_idx?255:c;
};
// Create a new BBWT structure 
Basic_BWT_t * new_Basic_BWT();
// Load a BBWT from a file in compacted format
int load_Basic_BWT_from_file(Basic_BWT_t ** BBWT,char * file_name);
// Save BBWT to a file in compacted format
int save_Basic_BWT_to_file(Basic_BWT_t * BBWT,char * file_name);
// Compare two BBWT structures
int cmp_BBWTs(Basic_BWT_t * BBWT1, Basic_BWT_t * BBWT2);
// The following call back is used to load the BWT of a text 
// Into a BBWT structure. The BWT is loaded into blocks.  
// Each block is stored into a buffer and the pointer to the buffer 
// put into _buffer. The length of the loaded block is returned by the callback. 
// The callback signals the end of the loading by returning a zero length 
// buffer (setting _buffer to zero or returning zero as the length of the block).
typedef unsigned int (*BWT_load_callback_t)(unsigned char ** _buffer, void * intern_state);
int load_Basic_BWT_from_callback(Basic_BWT_t ** _BBWT, unsigned int BWT_length,
	unsigned int primary_idx, BWT_load_callback_t load_callback, void * intern_state);
// Free a BBWT structure. 
void free_Basic_BWT(Basic_BWT_t * Basic_BWT);
int build_sequence_index (unsigned char * raw_seq,
		unsigned int seqlen, unsigned int ** _indexed_seq,
		unsigned int * _alloc_size,unsigned int * char_base);
// Build a BBWT index from a given text. 
Basic_BWT_t * Build_BWT_index_from_text(unsigned char * text,
	unsigned int textlen,unsigned int options);
// Counts the number of occurrences of pattern in a text indexed in BBWT. 
int patt_count(unsigned char *P,unsigned int m,
	Basic_BWT_t * Basic_BWT,unsigned int _SA_interval[2]);
// Do a bacward step on a BBWT, given an input interval and a character. 
int Backward_step(unsigned int * in_interval,unsigned int * out_interval,
		unsigned char c,Basic_BWT_t * Basic_BWT);
// Does an LF mapping on the BBWT index. 
unsigned int LF_map(unsigned int in_pos,Basic_BWT_t * Basic_BWT);
// Does a batch extraction on the BBWT where the suffix 
// array positions to be extracted are given in bitvector. 
int Basic_BWT_batch_extract(unsigned int * bitvector,
		unsigned int nelements,unsigned int * output_vector,
		Basic_BWT_t * Basic_BWT);

#endif
