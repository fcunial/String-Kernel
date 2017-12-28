#ifndef io_h
#define io_h

#ifndef CONCATENATION_SEPARATOR
#define CONCATENATION_SEPARATOR 'z'
#endif
#ifndef OUTPUT_SEPARATOR
#define OUTPUT_SEPARATOR '\n'
#endif
#ifndef BUFFER_CHUNK
#define BUFFER_CHUNK 1048576  // In bytes. Default=2^20.
#endif
#ifndef ALLOC_GROWTH_NUM  // Stack reallocation rate
#define ALLOC_GROWTH_NUM 4
#endif
#ifndef ALLOC_GROWTH_DENOM
#define ALLOC_GROWTH_DENOM 3
#endif

extern const char *DNA_ALPHABET;
extern const unsigned char SEPARATOR;


/**
 * In-memory representation of the concatenation of all DNA sequences stored in a multi-
 * FASTA file. Characters not in $DNA_ALPHABET$ are replaced by $CONCATENATION_SEPARATOR$. 
 * If the concatenation contains more than one string, each string except the last one is 
 * terminated by the special character $CONCATENATION_SEPARATOR$ not in the alphabet. 
 * At the end of the concatenation, its reverse-complement might be appended as well, 
 * terminated again by $CONCATENATION_SEPARATOR$.
 */
typedef struct {
	unsigned char *buffer;
	unsigned long length;  // Number of characters in memory, including RC.
	unsigned long inputLength;  // Number of non-header characters in the input file
	unsigned char hasRC;  // Reverse-complement present (1/0).
} Concatenation;


Concatenation loadFASTA(char *inputFilePath, unsigned char appendRC);


/**
 * In microseconds
 */
double getTime();


#endif