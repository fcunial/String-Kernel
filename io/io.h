/**
 * Basic input/output procedures.
 *
 * @author Fabio Cunial
 */
#ifndef io_h
#define io_h

#include <stdint.h>


// IO constants used throughout the code
#ifndef DNA5_alphabet_size
#define DNA5_alphabet_size 5  // Includes #
#endif
#ifndef CONCATENATION_SEPARATOR
#define CONCATENATION_SEPARATOR 'z'
#endif
#ifndef OUTPUT_SEPARATOR_1
#define OUTPUT_SEPARATOR_1 ','
#endif
#ifndef OUTPUT_SEPARATOR_2
#define OUTPUT_SEPARATOR_2 '\n'
#endif
#ifndef BUFFER_CHUNK
#define BUFFER_CHUNK 1024  // Size of a buffer chunk, in bytes.
#endif
#ifndef ALLOC_GROWTH_NUM
#define ALLOC_GROWTH_NUM 4  // Reallocation rate
#endif
#ifndef ALLOC_GROWTH_DENOM
#define ALLOC_GROWTH_DENOM 3  // Reallocation rate
#endif
extern char *DNA_ALPHABET;  // Characters of the alphabet
extern double DNA_ALPHABET_PROBABILITIES[4];  // Empirical probability of each character
extern double LOG_DNA_ALPHABET_PROBABILITIES[4];  // log_e of the above


/**
 * In-memory representation of the concatenation of all DNA sequences stored in a multi-
 * FASTA file. Characters not in $DNA_ALPHABET$ are replaced by $CONCATENATION_SEPARATOR$. 
 * If the concatenation contains more than one string, each string except the last one is 
 * terminated by the special character $CONCATENATION_SEPARATOR$ not in the alphabet. 
 * At the end of the concatenation, its reverse-complement might be appended as well, 
 * terminated again by $CONCATENATION_SEPARATOR$.
 */
typedef struct {
	char *buffer;
	uint64_t length;  // Number of characters in memory, including RC.
	uint64_t inputLength;  // Number of non-header characters in the input file
	uint8_t hasRC;  // Reverse-complement present (1/0).
} Concatenation;


Concatenation loadFASTA(char *inputFilePath, uint8_t appendRC);


/**
 * In microseconds.
 */
double getTime();


#endif