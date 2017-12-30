#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/time.h>
#include "io.h"


const char *DNA_ALPHABET = "acgt";
double DNA_ALPHABET_PROBABILITIES[4];


Concatenation loadFASTA(char *inputFilePath, unsigned char appendRC) {
	int i;
	unsigned int j;
	int c;
	unsigned int lineLength;  // Length of a line of the file
	unsigned int stringLength;  // Length of a FASTA sequence
	unsigned int bufferLength;
	unsigned int inputLength;  // Total length of the input, including non-DNA characters.
	unsigned int outputLength;  // Total length of the output, including non-DNA characters.
	unsigned int outputLengthDNA;  // Number of DNA characters in the output
	unsigned int outputLengthPrime;  // Temporary
	char *pointer;
	unsigned char *buffer;
	FILE *file;
	Concatenation out;
	
	file=fopen(inputFilePath,"r");
	if (file==NULL) {
		fprintf(stderr,"ERROR: cannot open input file \n");
		exit(EXIT_FAILURE);
	}
	fclose(file);
	
	// Loading the multi-FASTA input file
	for (i=0; i<4; i++) DNA_ALPHABET_PROBABILITIES[i]=0.0;
	file=fopen(inputFilePath,"r");
	if (file==NULL) {
		fprintf(stderr,"ERROR: cannot open input file \n");
		exit(EXIT_FAILURE);
	}
	buffer=(unsigned char *)malloc(BUFFER_CHUNK);
	bufferLength=BUFFER_CHUNK; inputLength=0; outputLength=0; outputLengthDNA=0;
	c=fgetc(file);
	do {
		if (c!='>') {
			fprintf(stderr,"ERROR: input file not in FASTA format \n");
			exit(EXIT_FAILURE);
		}
		// Header
		c=fgetc(file);
		while (c!='\n' && c!=EOF) c=fgetc(file);
		if (c==EOF) {
			fprintf(stderr,"Omitting empty string \n");
			break;
		}
		// String
		stringLength=0; lineLength=0;
		c=fgetc(file);
		while (c!=EOF && c!='>') {
			if (c=='\n') {
				c=fgetc(file);
				if (lineLength==0) fprintf(stderr,"Omitting empty line \n");
				lineLength=0;
				continue;
			}
			lineLength++; stringLength++; inputLength++; c=tolower(c);
			pointer=strchr(DNA_ALPHABET,c);
			if (pointer==NULL) c=CONCATENATION_SEPARATOR;
			else {
				outputLengthDNA++;
				DNA_ALPHABET_PROBABILITIES[pointer-DNA_ALPHABET]+=1.0;
			}
			if (outputLength==bufferLength) {
				bufferLength+=BUFFER_CHUNK;
				buffer=(unsigned char *)realloc(buffer,bufferLength*sizeof(unsigned char));
			}
			buffer[outputLength++]=c;
			c=fgetc(file);
		}
		if (c!=EOF) {
			if (stringLength==0) fprintf(stderr,"Omitting empty string \n");
			else {
				if (outputLength==bufferLength) {
					bufferLength+=BUFFER_CHUNK;
					buffer=(unsigned char *)realloc(buffer,bufferLength*sizeof(unsigned char));
				}
				buffer[outputLength++]=CONCATENATION_SEPARATOR;
			}
		}
	} while (c!=EOF);
	fclose(file);
	for (i=0; i<4; i++) DNA_ALPHABET_PROBABILITIES[i]/=outputLengthDNA;
	
	// Appending reverse-complement, if needed.
	if (appendRC==1) {
		outputLengthPrime=(outputLength<<1)+1;
		if (bufferLength<outputLengthPrime) {
			bufferLength=outputLengthPrime;
			buffer=(unsigned char *)realloc(buffer,bufferLength*sizeof(unsigned char));
		}
		buffer[outputLength]=CONCATENATION_SEPARATOR;
		i=outputLength-1; j=outputLength+1;
		while (i>=0) {
			switch (buffer[i]) {
				case 'a': buffer[j]='t'; break;
				case 'c': buffer[j]='g'; break;
				case 'g': buffer[j]='c'; break;
				case 't': buffer[j]='a'; break;
				case CONCATENATION_SEPARATOR: buffer[j]=CONCATENATION_SEPARATOR; break;
			}
			i--; j++;
		}
	}
	else outputLengthPrime=outputLength;
	
	out.buffer=buffer;
	out.length=outputLengthPrime;
	out.inputLength=inputLength;
	out.hasRC=appendRC;
	return out;
}


double getTime() {
	struct timeval ttime;
	gettimeofday(&ttime,0);
	return ttime.tv_sec+ttime.tv_usec*0.000001;
}
