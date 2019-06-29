#ifndef bits_h
#define bits_h


/**
 * (For debugging only)
 * Prints the bits in $number$ from MSB to LSB.
 */
void printLong(unsigned long number);


/**
 * Read the $i$-th pair of bits from $buffer$.
 *
 * Remark: bits inside each long of $buffer$ are assumed to be stored from LSB to MSB.
 */
char readTwoBits(unsigned long *buffer, unsigned int i);


/**
 * Writes $value$ in the $i$-th pair of bits from $buffer$. $value$ is assumed to use just
 * the two LSBs.
 *
 * Remark: bits inside each long of $buffer$ are assumed to be stored from LSB to MSB.
 */
void writeTwoBits(unsigned long *buffer, unsigned int i, unsigned char value);


/**
 * Remark: bits inside each long of $buffer$ are assumed to be stored from LSB to MSB.
 */
char readBit(unsigned long *buffer, unsigned int i);


/** 
 * @param value 1/0.
 *
 * Remark: bits inside each long of $buffer$ are assumed to be stored from LSB to MSB.
 */
void writeBit(unsigned long *buffer, unsigned int i, unsigned char value);


/**
 * @return 1 iff $bitvector[0..lastBit]$ (coordinates in bits) contains a one-bit.
 */
char hasOneBit(unsigned long *bitvector, unsigned int lastBit);


#endif