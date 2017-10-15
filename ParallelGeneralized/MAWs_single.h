#ifndef MAWs_single_h
#define MAWs_single_h


#include "SLT_single_string.h"


/** 
 * @param mem 0 iff MAWs should not be written to the output; otherwise, MAWs are written
 * to output file $filePath$;
 * @return the number of MAWs found.
 */
unsigned int find_MAWs_single(Basic_BWT_t *BBWT, unsigned int minlen, unsigned int mem, unsigned char *filePath);


#endif