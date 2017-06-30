#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include <string.h>
#include"mt19937ar.h"
#include"DNA5_Basic_BWT.h"
#include"SLT_MAWs.h"
#include<time.h>
#include "../malloc_count-master/malloc_count.h"
#include "naive_MAWs.h"

double gettime( void )
{
	struct timeval ttime;
	gettimeofday( &ttime , 0 );
	return ttime.tv_sec + ttime.tv_usec * 0.000001;
};
int str_compar(const void *_str1, const void *_str2)
{
	return strcmp(*(char **)_str1,*(char **)_str2);
};


#define min_MAW_len 2
#define ALLOC_SIZE 1048576
#define DNA                     "ACGT"                         //DNA alphabet

int main(int argc, char **argv)
{
	unsigned char * text1=NULL;
	unsigned char * text2=NULL;
	unsigned int textlen1=0;
	unsigned int textlen2=0;
	unsigned int nMAWs1;
	unsigned int nMAWs2;
	unsigned int nMAWs;
	unsigned int memory=atoi(argv[3]);
	unsigned int RC=atoi(argv[4]);
	unsigned int cores= atoi(argv[5]);
	double LW=0;
	Basic_BWT_t * BBWT1;
	Basic_BWT_t * BBWT2;

	unsigned int i;
	unsigned int num=66;
	char files[num][20];
	strcpy(files[0],"../data/BA.fa");
	strcpy(files[1],"../data/BS.fa");
	strcpy(files[2],"../data/EC.fa");
	strcpy(files[3],"../data/HI.fa");
	strcpy(files[4],"../data/HP.fa");
	strcpy(files[5],"../data/LC.fa");
	strcpy(files[6],"../data/LL.fa");
	strcpy(files[7],"../data/MG.fa");
	strcpy(files[8],"../data/SA.fa");
	strcpy(files[9],"../data/SP.fa");
	strcpy(files[10],"../data/XC.fa");
	strcpy(files[11],"../data/AT1.fasta");
	strcpy(files[12],"../data/AT2.fasta");
	strcpy(files[13],"../data/AT3.fasta");
	strcpy(files[14],"../data/AT4.fasta");
	strcpy(files[15],"../data/AT5.fasta");
	strcpy(files[16],"../data/DM2L.fasta");
	strcpy(files[17],"../data/DM2R.fasta");
	strcpy(files[18],"../data/DM3L.fasta");
	strcpy(files[19],"../data/DM3R.fasta");
	strcpy(files[20],"../data/DMX.fasta");
	for(i=21;i<45;i++) {
		char t[20];
		snprintf(t, 20, "../data/HS%d.fasta", i-20);
		strcpy(files[i],t);
	}
	for(i=45;i<66;i++) {
		char t[20];
		snprintf(t, 20, "../data/MM%d.fasta", i-44);
		strcpy(files[i],t);
	}

	// clearing output file
	FILE *x=fopen("output.txt", "w");
	fclose(x);

	// Reading File 1
	FILE *f1;
	/* Read the (Multi)FASTA file in memory */
	if ( ! ( f1 = fopen ( files[atoi(argv[1])], "r") ) )
	{
		fprintf ( stderr, " Error: Cannot open file!" );
		return ( 1 );
	}

	char c;
	c = fgetc( f1 );
	unsigned int max_alloc_seq = 0;
	do
	{
		if ( c != '>' )
		{
			fprintf ( stderr, " Error: input file is not in FASTA format!\n" );
			return ( 1 );
		}
		else
		{
			while ( ( c = fgetc( f1 ) ) != EOF && c != '\n' )
			{
			}
		}
		while ( ( c = fgetc( f1 ) ) != EOF && c != '>' )
		{
			if( textlen1 == 0 && c == '\n' )
			{
				fprintf ( stderr, " Omitting empty sequence in file!\n" );
				c = fgetc( f1 );
				break;
			}
			if( c == '\n' ) continue;

			if ( textlen1 >= max_alloc_seq )
			{
				text1 = ( unsigned char * ) realloc ( text1,   ( max_alloc_seq + ALLOC_SIZE ) * sizeof ( unsigned char ) );
				max_alloc_seq += ALLOC_SIZE;
			}
			if( strchr ( DNA, c ) )
			{
				text1[ textlen1++ ] = c;
			}
		}
	} while( c != EOF );
	unsigned long int l=textlen1;
			if(RC) {
				text1 = ( unsigned char * ) realloc ( text1,   (textlen1*2+2) * sizeof ( unsigned char ) );
				unsigned long int r=textlen1-1;
				text1[l++]='Z';
				while ( l<textlen1*2+1)
				{
					switch ( text1[r--] )
					{
					case 'A':
						text1[l++] = 'T';
						break;
					case 'C':
						text1[l++] = 'G';
						break;
					case 'G':
						text1[l++] = 'C';
						break;
					case 'T':
						text1[l++] = 'A';
						break;
					default:
						return ( 0 );
					}
				}
				text1[l++]='Z';
			}
	fclose(f1);
	// Reading File 2
	if ( ! ( f1 = fopen ( files[atoi(argv[2])], "r") ) )
	{
		fprintf ( stderr, " Error: Cannot open file!" );
		return ( 1 );
	}

	c = fgetc( f1 );
	max_alloc_seq = 0;
	do
	{
		if ( c != '>' )
		{
			fprintf ( stderr, " Error: input file is not in FASTA format!\n" );
			return ( 1 );
		}
		else
		{
			while ( ( c = fgetc( f1 ) ) != EOF && c != '\n' )
			{
			}

		}
		while ( ( c = fgetc( f1 ) ) != EOF && c != '>' )
		{
			if( textlen2 == 0 && c == '\n' )
			{
				fprintf ( stderr, " Omitting empty sequence in file!\n" );
				c = fgetc( f1 );
				break;
			}
			if( c == '\n' ) continue;

			if ( textlen2 >= max_alloc_seq )
			{
				text2 = ( unsigned char * ) realloc ( text2,   ( max_alloc_seq + ALLOC_SIZE ) * sizeof ( unsigned char ) );
				max_alloc_seq += ALLOC_SIZE;
			}
			if( strchr ( DNA, c ) )
			{
				text2[ textlen2++ ] = c;
			}
		}
	} while( c != EOF );
	l=textlen2;
			if(RC) {
				text2 = ( unsigned char * ) realloc ( text2,   (textlen2*2+2) * sizeof ( unsigned char ) );
				unsigned long int r=textlen2-1;
				text2[l++]='Z';
				while ( l<textlen2*2+1)
				{
					switch ( text2[r--] )
					{
					case 'A':
						text2[l++] = 'T';
						break;
					case 'C':
						text2[l++] = 'G';
						break;
					case 'G':
						text2[l++] = 'C';
						break;
					case 'T':
						text2[l++] = 'A';
						break;
					default:
						return ( 0 );
					}
				}
				text2[l++]='Z';
			}
	fclose(f1);

	printf("Text len is: %d - %d\n",textlen1, textlen2);

	double t1= gettime();
	// Build a BWT index on the text.
	BBWT1=Build_BWT_index_from_text(text1,textlen1,Basic_bwt_free_text);
	BBWT2=Build_BWT_index_from_text(text2,textlen2,Basic_bwt_free_text);
	// Launch the SLT based algorithm
	double t2= gettime();
	nMAWs=SLT_find_MAWs(BBWT1,BBWT2,min_MAW_len,&nMAWs1,&nMAWs2,&LW, memory, cores);
	naive_find_MAWs(text1, textlen1, 2);
	double t3=gettime();
	double jaccard= (double) nMAWs/(nMAWs1+nMAWs2-nMAWs);
	FILE *results= fopen("../results_gen", "a");
	fprintf(results,"Computing %s and %s; Common MAWs is %d, Maws1: %d, Maws2: %d; Jaccard: %f; LW: %f\n", files[atoi(argv[1])], files[atoi(argv[2])], nMAWs, nMAWs1, nMAWs2, jaccard,LW);
	//fprintf(results, "Time BWT: %f; Time MAWs: %f; Our peak memory allocation: %lld; Number of cores: %d\n",t2-t1, t3-t2,(long long)malloc_count_peak(),cores);
	free_Basic_BWT(BBWT1);
	free_Basic_BWT(BBWT2);
	return 0; 
}
