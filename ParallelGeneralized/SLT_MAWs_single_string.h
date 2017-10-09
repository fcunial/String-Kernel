#include "SLT_single_string.h"

#ifndef SLT_MAWs_h
#define SLT_MAWs_h
#endif


/** 
 * @param mem 0 iff MAWs should not be written to the output; otherwise, MAWs are written
 * to output file "maws.txt";
 * @return the number of MAWs found.
 */
unsigned int SLT_find_MAWs_single_string( Basic_BWT_t *BBWT1, 
  										  unsigned int minlen,
										  unsigned int mem );
