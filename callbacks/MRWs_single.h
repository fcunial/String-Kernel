#ifndef MRWs_single_h
#define MRWs_single_h


#include "../io/io.h"
#include "../iterator/SLT_single_string.h"


/** 
 * Detects minimal rare words $W$ such that $minFreq \leq f(W) < maxFreq$ and  
 * $f(V) \geq maxFreq$ for every substring $V$ of $W$.
 *
 * @param minLength considers only MRWs of length at least $minLength$;
 * @param writeMRWs 0 iff MRWs should not be written to the output; otherwise, MRWs are  
 * appended to file $filePath$;
 * @return the number of MRWs found.
 */
unsigned int find_MRWs_single(Basic_BWT_t *BBWT, unsigned int textLength, unsigned int minLength, unsigned int minFreq, unsigned int maxFreq, unsigned char writeMRWs, unsigned char computeScore, char *filePath);

#endif