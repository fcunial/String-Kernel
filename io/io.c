#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "io.h"

#define BUFFER_CHUNK 1048576  // In bytes. Default=2^20.


const char *DNA_ALPHABET = "ACGT";


Concatenation loadFASTA(char *inputFilePath, unsigned char appendRC) {
	long i;
	unsigned long j;
	int c;
	unsigned long lineLength;
	unsigned long bufferLength, stringLength, inputLength, outputLength, outputLengthPrime;
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
	file=fopen(inputFilePath,"r");
	if (file==NULL) {
		fprintf(stderr,"ERROR: cannot open input file \n");
		exit(EXIT_FAILURE);
	}
	buffer=(unsigned char *)malloc(BUFFER_CHUNK);
	bufferLength=BUFFER_CHUNK; inputLength=0; outputLength=0;
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
			lineLength++; inputLength++;
			if (strchr(DNA_ALPHABET,c)!=NULL) {	
				if (outputLength==bufferLength) {
					bufferLength+=BUFFER_CHUNK;
					buffer=(unsigned char *)realloc(buffer,bufferLength*sizeof(unsigned char));
				}
				buffer[outputLength++]=c;
			}
			c=fgetc(file);
		}
	} while (c!=EOF);
	fclose(file);
	
	// Appending reverse-complement, if needed.
	if (appendRC==1) {
		outputLengthPrime=(outputLength<<1)+2;
		if (bufferLength<outputLengthPrime) {
			bufferLength=outputLengthPrime;
			buffer=(unsigned char *)realloc(buffer,bufferLength*sizeof(unsigned char));
		}
		buffer[outputLength]=CONCATENATION_SEPARATOR;
		buffer[outputLengthPrime-1]=CONCATENATION_SEPARATOR;
		i=outputLength-1; j=outputLength+1;
		while (i>=0) {
			switch (buffer[i]) {
				case 'A': buffer[j]='T'; break;
				case 'C': buffer[j]='G'; break;
				case 'G': buffer[j]='C'; break;
				case 'T': buffer[j]='A'; break;
			}
			i--; j++;
		}
	}
	else {
		outputLengthPrime=outputLength+1;
		if (bufferLength<outputLengthPrime) {
			bufferLength=outputLengthPrime;
			buffer=(unsigned char *)realloc(buffer,bufferLength*sizeof(unsigned char));
		}
		buffer[outputLength]=CONCATENATION_SEPARATOR;
	}
	
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