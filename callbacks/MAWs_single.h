#ifndef MAWs_single_h
#define MAWs_single_h


#include "../io/io.h"
#include "../iterator/SLT_single_string.h"


/**  
 * @param minLength considers only MAWs of length at least $minLength$;
 * @param writeMAWs 0 iff MAWs should not be written to the output; otherwise, MAWs are  
 * appended to file $filePath$;
 * @return the number of MAWs found.
 */
unsigned int find_MAWs_single(Basic_BWT_t *BBWT, unsigned int minLength, unsigned char writeMAWs, char *filePath);


#endif