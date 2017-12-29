#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/time.h>
#include "io.h"


const char *DNA_ALPHABET = "acgt";


Concatenation loadFASTA(char *inputFilePath, unsigned char appendRC) {
	long i;
	unsigned long j;
	int c;
	unsigned long lineLength, stringLength;
	unsigned long bufferLength, inputLength, outputLength, outputLengthPrime;
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
			lineLength++; stringLength++; inputLength++; c=tolower(c);
			if (strchr(DNA_ALPHABET,c)==NULL) c=CONCATENATION_SEPARATOR;
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
