#ifndef bits_h
#define bits_h


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


#endif