#ifndef io_h
#define io_h

#ifndef CONCATENATION_SEPARATOR
#define CONCATENATION_SEPARATOR 'Z'
#endif
#ifndef OUTPUT_SEPARATOR
#define OUTPUT_SEPARATOR '\n'
#endif

extern const char *DNA_ALPHABET;
extern const unsigned char SEPARATOR;


/**
 * In-memory representation of the concatenation of all DNA sequences stored in a  
 * multi-FASTA file. Characters not in $DNA_ALPHABET$ are not loaded. The concatenation 
 * is terminated with the special character $CONCATENATION_SEPARATOR$, but sequences in 
 * the FASTA file are not separated. At the end of the concatenation, its reverse- 
 * complement might be appended as well, terminated again by $CONCATENATION_SEPARATOR$.
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