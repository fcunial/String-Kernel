#ifndef DNA5_Basic_BWT64_h
#define DNA5_Basic_BWT64_h
#include"indexed_DNA5_seq64.h"
#include<stdlib.h>
#ifndef Basic_bwt_no_free_text
#define Basic_bwt_no_free_text 0
#endif
#ifndef Basic_bwt_free_text
#define Basic_bwt_free_text 1
#endif

typedef struct 
{
	indexed_DNA5_seq64_t * indexed_BWT;
	unsigned long long char_base[5];
//	unsigned long long size;
	unsigned long long primary_idx;
	unsigned long long textlen;
} Basic_BWT64_t;
static inline unsigned char DNA5_BWT_get_prev_char(Basic_BWT64_t * Basic_BWT,unsigned long long suff_idx)
{
	unsigned char c;
	suff_idx++;	
	c=DNA5_extract_char64(Basic_BWT->indexed_BWT,suff_idx);
	return suff_idx==Basic_BWT->primary_idx?255:c;
};
// Create a new BBWT structure 
Basic_BWT64_t * new_Basic_BWT64();
// Load a BBWT from a file in compacted format
long long load_Basic_BWT64_from_file(Basic_BWT64_t ** _BBWT,char * file_name);
// Save BBWT to a file in compacted format
long long save_Basic_BWT64_to_file(Basic_BWT64_t * BBWT,char * file_name);
// Compare two BBWT structures
int cmp_BBWT64s(Basic_BWT64_t * BBWT1, Basic_BWT64_t * BBWT2);
// The following call back is used to load the BWT of a text 
// Into a BBWT structure. The BWT is loaded into blocks.  
// Each block is stored into a buffer and the pointer to the buffer 
// put into _buffer. The length of the loaded block is returned by the callback. 
// The callback signals the end of the loading by returning a zero length 
// buffer (setting _buffer to zero or returning zero as the length of the block).

typedef unsigned long long (*BWT64_load_callback_t)(unsigned char ** _buffer, void * intern_state);

int load_Basic_BWT64_from_callback(Basic_BWT64_t ** _BBWT, unsigned long long BWT_length,
	unsigned long long primary_idx, BWT64_load_callback_t load_callback, void * intern_state);
// Free a BBWT structure. 
void free_Basic_BWT64(Basic_BWT64_t * Basic_BWT);


int build_sequence_index64(unsigned char * raw_seq,
		unsigned long long seqlen,indexed_DNA5_seq64_t ** _indexed_DNA5_seq64,
		unsigned long long * _alloc_size,unsigned int * char_base);
// Build a BBWT index from a given text. 
Basic_BWT64_t * Build_BWT_index_from_text64(unsigned char * text,
	unsigned long long textlen,unsigned int options);
// Counts the number of occurrences of pattern in a text indexed in BBWT. 
unsigned long long patt_count64(unsigned char *P,unsigned long long m,
	Basic_BWT64_t * Basic_BWT,unsigned long long _SA_interval[2]);
// Do a backward step on a BBWT, given an input interval and a character. 
int Backward_step64(unsigned long long * in_interval,unsigned long long * out_interval,
		unsigned char c,Basic_BWT64_t * Basic_BWT);
// Does an LF mapping on the BBWT index. 
unsigned long long LF_map64(unsigned long long in_pos,Basic_BWT64_t * Basic_BWT);
// Does a batch extraction on the BBWT where the suffix 
// array positions to be extracted are given in bitvector. 
int Basic_BWT_batch_extract64(unsigned int * bitvector,
		unsigned long long nelements,unsigned long long * output_vector,
		Basic_BWT64_t * Basic_BWT);
#endif
