#ifndef io_h
#define io_h
#include "alphabet.h"


/**
 * Memory representation of the concatenation of all DNA sequences stored in a  
 * multi-FASTA file. Characters not in $DNA_ALPHABET$ are not loaded. The concatenation 
 * is terminated with the special character $SEPARATOR$, but sequences in the FASTA file
 * are not separated. At the end of the concatenation, its reverse-complement might be
 * appended as well, terminated again by $SEPARATOR$.
 */
typedef struct {
	unsigned char *buffer;
	unsigned long long length;  // Number of characters in memory, including RC.
	unsigned long long inputLength;  // Number of non-header characters in the file
	unsigned char hasRC;  // Reverse-complement present (1) or absent (0).
} Concatenation;


Concatenation *loadFASTA(unsigned char *inputFilePath, unsigned char appendRC);


#endif