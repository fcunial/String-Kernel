#ifndef DNA5_Basic_BWT_h
#define DNA5_Basic_BWT_h

#include <stdint.h>
#include <stdlib.h>
#include "../iterator/indexed_DNA5_seq.h"
#include "../io/io.h"


/**
 * Creates the index.
 * 
 * @param options for dbwt.
 */
uint8_t *useDivsufsort(char *text, uint64_t length, uint64_t *sharpPosition);

#endif
