/**
 * 
 *
 * @author Fabio Cunial
 */
#ifndef buffered_file_writer_h
#define buffered_file_writer_h


#include <stdio.h>


typedef struct {
	char *buffer;
	unsigned int capacity;  // Maximum number of characters in the buffer
	unsigned int size;  // Number of characters currently in the buffer
	FILE *file;
} buffered_file_writer_t;


/** 
 * Remark: $file->file$ is opened in append mode, so its content is preserved.
 */
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
 * $from[0..lastBit]$ (coordinates refer to bits), as characters.
 *
 * Remark: bits inside each long of $from$ are assumed to be stored from LSB to MSB.
 */
void writeBits(unsigned long *from, unsigned int lastBit, buffered_file_writer_t *to);


/**
 * Let $from$ be an array of 2-bit numbers. The procedure appends to $to$ all numbers in 
 * $from[0..last]$ (coordinates refer to numbers), in reverse order, interpreting each 
 * number as a position in $alphabet$.
 *
 * Remark: bits inside each long of $from$ are assumed to be stored from LSB to MSB.
 */
void writeTwoBitsReversed(unsigned long *from, unsigned int last, buffered_file_writer_t *to, char *alphabet);


#endif