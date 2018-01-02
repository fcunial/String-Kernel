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
#define BUFFER_CHUNK 1048576  // In bytes. Default=2^20.
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
#ifndef BIT_MASK
#define BIT_MASK 1L  // 1-bit selector
#endif
#ifndef TWO_BIT_MASK
#define TWO_BIT_MASK 3L  // 2-bit selector
#endif


extern const char *DNA_ALPHABET;  // Characters of the alphabet
extern double DNA_ALPHABET_PROBABILITIES[4];  // Empirical probability of each character
extern double LOG_DNA_ALPHABET_PROBABILITIES[4];  // log_e of the above
extern const unsigned char BITS_PER_LONG, INITIAL_REM;
extern const unsigned long INITIAL_MASK;


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


/**
 * For debugging only
 */
void printLong(unsigned long number);


/**
 * Read the $i$-th pair of bits from $buffer$.
 */
char readTwoBits(unsigned long *buffer, unsigned int i);


/**
 * Writes $value$ in the $i$-th pair of bits from $buffer$.
 * $value$ is assumed to use just the two LSBs.
 */
void writeTwoBits(unsigned long *buffer, unsigned int i, unsigned char value);


char readBit(unsigned long *buffer, unsigned int i);


/** 
 * @param value 1/0.
 */
void writeBit(unsigned long *buffer, unsigned int i, unsigned char value);


/**
 * Appends to $out$ bits $[0..lastBit]$ of $bitmap$, as characters.
 *
 * @param outSize the number of characters that have already been written to $out$;
 * @return the new value of $outSize$.
 */
unsigned int printBits(unsigned long *bitmap, unsigned int lastBit, char *out, unsigned int outSize);


/**
 * Let $array$ be an array of 2-bit numbers. The procedure appends to $out$ all numbers
 * in $array[0..lastElement]$, in reverse order, interpreting each number as a position in
 * $alphabet$.
 *
 * @param outSize the number of characters that have already been written to $out$;
 * @return the new value of $outSize$.
 */
unsigned int printTwoBitsReverse(unsigned long *array, unsigned int lastElement, char *out, unsigned int outSize, const char *alphabet);

#endif