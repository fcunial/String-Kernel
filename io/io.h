#ifndef io_h
#define io_h

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
#define BUFFER_CHUNK 1024  // In bytes. Default=2^10.
#endif
#ifndef ALLOC_GROWTH_NUM  // Stack reallocation rate
#define ALLOC_GROWTH_NUM 4
#endif
#ifndef ALLOC_GROWTH_DENOM
#define ALLOC_GROWTH_DENOM 3
#endif
#ifndef MY_CEIL  // ceil(N/D) where N and D are integers.
#define MY_CEIL(N,D) (1+((N)-1)/(D))
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