#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include"naive_MAWs.h"

#define alloc_growth_num 4
#define alloc_growth_denom 3

int contains_N(unsigned char* str,unsigned int str_len)
{
	unsigned int i;
	for(i=0;i<str_len;i++)
		if(str[i]=='Z')
			return 1;
	return 0;
}
double g(unsigned int y, unsigned int length) {
	return (double) (length-y+2)/(length-y+1)*(length-y+2)/(length-y+3);
}

static unsigned char alpha4_to_ACGT[4]={'A','C','G','T'};

unsigned int naive_find_MAWs(unsigned char * text, unsigned int textlen,
		unsigned int minlen)
{
	unsigned int k,h,fw, faw, fwb;
	unsigned int i,j;
	unsigned char cl,cr;
	unsigned char * mark_vector=(unsigned char *) malloc(textlen);
	unsigned char left_exists[5];
	unsigned char right_exists[5];
	unsigned char left_right_exists[5][5];
	unsigned int nMAWs=0;
	unsigned char inv_alphabet[256];
	unsigned char right_char;
	unsigned char left_char;
	double D=0;
	for(i=0;i<256;i++)
	{
		switch(i)
		{
		case 'A': inv_alphabet[i]=0;break;
		case 'C': inv_alphabet[i]=1;break;
		case 'G': inv_alphabet[i]=2;break;
		case 'T': inv_alphabet[i]=3;break;
		default: inv_alphabet[i]=4;
		};
	};
	if(minlen<=2)
	{
		for(cl=0;cl<5;cl++)
		{
			left_exists[cl]=0;
			right_exists[cl]=0;
			for(cr=0;cr<5;cr++)
				left_right_exists[cl][cr]=0;
		};
		for(i=0;i<textlen-1;i++)
		{
			left_char=inv_alphabet[text[i]];
			right_char=inv_alphabet[text[i+1]];
			left_exists[left_char]+=1;
			right_exists[right_char]+=1;
			left_right_exists[left_char][right_char]+=1;
		};
		left_char=inv_alphabet[text[i]];
		//first string
		left_exists[4]+=1;
		right_exists[inv_alphabet[text[0]]]+=1;
		left_right_exists[4][inv_alphabet[text[0]]]+=1;
		//last string
		left_exists[left_char]+=1;
		right_exists[4]+=1;
		left_right_exists[left_char][4]+=1;

		for(cl=0;cl<5;cl++)
		{
			if(!left_exists[cl])
				continue;
			for(cr=0;cr<5;cr++)
				if(right_exists[cr]) {
					if(left_right_exists[cl][cr]==0) {
						if(cl!=4 && cr!=4) {
							// We have a minimal absent word. We output it.
							nMAWs++;
							D++;
						}
					}
					else {
						fw= 0;
						faw= 0;
						fwb= 0;
						for(h=0; h<5; h++) {
							faw+= left_right_exists[cl][h];
							fwb+= left_right_exists[h][cr];
							for(k=0; k<5; k++)
								fw+= left_right_exists[h][k];
						}

						D+=(g(2,textlen+2)*left_right_exists[cl][cr]*fw/(faw*fwb)-1)*
								(g(2,textlen+2)*left_right_exists[cl][cr]*fw/(faw*fwb)-1);
					}
				}
		};
	};
	for(k=1;k<=textlen;k++)
	{
		for(i=0;i<textlen;i++)
			mark_vector[i]=0;
		for(i=0;i<textlen-k+1;i++)
		{
			if(mark_vector[i])
				continue;
			for(cl=0;cl<5;cl++)
			{
				left_exists[cl]=0;
				right_exists[cl]=0;
			};	
			for(cl=0;cl<5;cl++)
				for(cr=0;cr<5;cr++)
					left_right_exists[cl][cr]=0;
			if(i>0)
				left_char=inv_alphabet[text[i-1]];
			else
				left_char=4;
			if(i+k==textlen)
				right_char=4;
			else
				right_char=inv_alphabet[text[i+k]];
			left_exists[left_char]+=1;
			right_exists[right_char]+=1;
			left_right_exists[left_char][right_char]+=1;
			for(j=i+1;j<textlen-k+1;j++)
				if(memcmp(&text[i],&text[j],k)==0)
				{
					mark_vector[j]=1;
					left_char=inv_alphabet[text[j-1]];
					if(j<textlen-k)
						right_char=inv_alphabet[text[j+k]];
					else
						right_char=4;
					left_exists[left_char]+=1;
					right_exists[right_char]+=1;
					left_right_exists[left_char][right_char]+=1;
				};
			for(cl=0;cl<5;cl++)
			{
				if(left_exists[cl]==0)
					continue;
				for(cr=0;cr<5;cr++) {
					if(right_exists[cr]) {
						if(left_right_exists[cl][cr]==0) {
							if(cl!=4 && cr!=4) {
								// We have a minimal absent word. We output it.
								nMAWs++;
								D++;
							}
						}
						else {
							unsigned int z;
							fw= 0;
							faw=0;
							fwb=0;
							for(h=0; h<5; h++) {
								faw+= left_right_exists[cl][h];
								fwb+= left_right_exists[h][cr];
								for(z=0; z<5; z++)
									fw+= left_right_exists[h][z];
							}

							D+=(g(k+2,textlen+2)*left_right_exists[cl][cr]*fw/(faw*fwb)-1)*
									(g(k+2,textlen+2)*left_right_exists[cl][cr]*fw/(faw*fwb)-1);
						}
					}
				}
			};
		};
	};
	free(mark_vector);
	printf("Naive Markovian kernel: %f MAWs: %d\n",D, nMAWs);
	return nMAWs;
};

