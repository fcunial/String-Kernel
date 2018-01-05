#ifndef buffered_file_writer_h
#define buffered_file_writer_h


#include <stdio.h>


typedef struct {
	char *buffer;
	unsigned int capacity;  // Maximum number of characters in the buffer
	unsigned int size;  // Number of characters currently in the buffer
	FILE *file;
} buffered_file_writer_t;


void initializeBufferedFileWriter(buffered_file_writer_t *file, char *path);


void finalizeBufferedFileWriter(buffered_file_writer_t *file);


/**
 * Writes character $c$ to $to->file$.
 */
void writeChar(char c, buffered_file_writer_t *to);


/**
 * Writes to $to->file$ all characters in $from[0..last]$.
 */
void writeChars(char *from, unsigned int last, buffered_file_writer_t *to);


/**
 * Let $from$ be an array of bits. The procedure appends to $to$ all bits in 
 * $from[0..lastBit]$, as charactersb.
 */
void writeBits(unsigned long *from, unsigned int lastBit, buffered_file_writer_t *to);


/**
 * Let $from$ be an array of 2-bit numbers. The procedure appends to $to$ all numbers in 
 * $from[0..last]$, in reverse order, interpreting each number as a position in 
 * $alphabet$.
 */
void writeTwoBitsReversed(unsigned long *from, unsigned int last, buffered_file_writer_t *to, char *alphabet);


#endif