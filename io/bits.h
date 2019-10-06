/**
 * Basic operations on bitvectors.
 *
 * @author Fabio Cunial
 */
#ifndef bits_h
#define bits_h

#include <stdint.h>


// Bit constants used throughout the code
#ifndef BYTES_PER_CHAR
#define BYTES_PER_CHAR (sizeof(char))
#endif
#ifndef BYTES_PER_WORD
#define BYTES_PER_WORD ((sizeof(uint32_t)))
#endif
#ifndef BYTES_PER_LONG
#define BYTES_PER_LONG (sizeof(uint64_t))
#endif
#ifndef BITS_PER_BYTE
#define BITS_PER_BYTE 8
#endif
#ifndef BITS_PER_WORD
#define BITS_PER_WORD (((BYTES_PER_WORD)*(BITS_PER_BYTE)))
#endif
#ifndef BITS_PER_LONG
#define BITS_PER_LONG ((BYTES_PER_LONG)<<3)
#endif
#ifndef BIT_MASK
#define BIT_MASK (1L)  // 1-bit long selector
#endif
#ifndef TWO_BIT_MASK
#define TWO_BIT_MASK (3L)  // 2-bit long selector
#endif
#ifndef ALL_ONES_8
#define ALL_ONES_8 (0xFF)
#endif
#ifndef ALL_ONES_32
#define ALL_ONES_32 (0xFFFFFFFF)
#endif
#ifndef ALL_ONES_64
#define ALL_ONES_64 (0xFFFFFFFFFFFFFFFF)
#endif
#ifndef MY_CEIL
#define MY_CEIL(N,D) (1+((N)-1)/(D))  // ceil(N/D) where N and D are integers.
#endif


/**
 * Read the $i$-th pair of bits from $buffer$.
 *
 * Remark: bits inside each long of $buffer$ are assumed to be stored from LSB to MSB.
 */
uint8_t readTwoBits(uint64_t *buffer, uint64_t i);


/**
 * Writes $value$ in the $i$-th pair of bits from $buffer$. $value$ is assumed to use just
 * the two LSBs.
 *
 * Remark: bits inside each long of $buffer$ are assumed to be stored from LSB to MSB.
 */
void writeTwoBits(uint64_t *buffer, uint64_t i, uint8_t value);


/**
 * Remark: bits inside each long of $buffer$ are assumed to be stored from LSB to MSB.
 */
uint8_t readBit(uint64_t *buffer, uint64_t i);


/** 
 * @param value 1/0.
 *
 * Remark: bits inside each long of $buffer$ are assumed to be stored from LSB to MSB.
 */
void writeBit(uint64_t *buffer, uint64_t i, uint8_t value);


/**
 * @return 1 iff $bitvector[0..lastBit]$ (coordinates in bits) contains a one-bit.
 */
uint8_t hasOneBit(uint64_t *bitvector, uint64_t lastBit);


/**
 * (For debugging only)
 * Prints the bits in $number$ from MSB to LSB.
 */
void printLong(uint64_t number);


#endif