#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include"naive_MAWs.h"

#define alloc_growth_num 4
#define alloc_growth_denom 3

static inline unsigned int contains_N(unsigned char* str,unsigned int str_len)
{
	unsigned int i;
	for(i=0;i<str_len;i++)
		if(str[i]=='Z')
			return 1;
	return 0;
};
static unsigned char alpha4_to_ACGT[4]={'A','C','G','T'};

unsigned int naive_find_MAWs(unsigned char * text, unsigned int textlen,
		unsigned int minlen,unsigned char *** _MAW_ptr,
		unsigned char ** _MAW_buffer)
{
	unsigned int k;
	unsigned int i,j;
	unsigned char cl,cr;
	unsigned char * mark_vector=(unsigned char *) malloc(textlen);
	unsigned char left_exists[5];
	unsigned char right_exists[5];
	unsigned char left_right_exists[5][5];
	unsigned int nMAWs=0;
	unsigned char ** MAW_ptr;
	unsigned int * MAW_idx=0;
	unsigned int nMAW_capacity=0;
	unsigned char * MAW_buffer=0;
	unsigned int MAW_buffer_idx=0;
	unsigned int MAW_buffer_size=0;
	unsigned int MAW_length;
	unsigned int is_repeat;
	unsigned char inv_alphabet[256];
	unsigned char right_char;
	unsigned char left_char;
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
		for(cl=0;cl<4;cl++)
		{
			left_exists[cl]=0;
			right_exists[cl]=0;
		};	
		for(cl=0;cl<4;cl++)
			for(cr=0;cr<4;cr++)
				left_right_exists[cl][cr]=0;
		for(i=0;i<textlen-1;i++)
		{
			left_char=inv_alphabet[text[i]];
			right_char=inv_alphabet[text[i+1]];
			left_exists[left_char]=1;
			right_exists[left_char]=1;
			left_right_exists[left_char][right_char]=1;
		};
		left_char=inv_alphabet[text[i]];
		left_exists[left_char]=1;
		right_exists[left_char]=1;
		for(cl=0;cl<4;cl++)
		{
			if(!left_exists[cl])
				continue;
			for(cr=0;cr<4;cr++)
				if(right_exists[cr] && left_right_exists[cl][cr]==0)
				{
				// We have a minimal absent word. We output it. 					
					nMAWs++;
					if(nMAWs>nMAW_capacity)
					{
						nMAW_capacity=(
							nMAW_capacity+1)
							*alloc_growth_num/alloc_growth_denom;
						MAW_idx=(unsigned int*)realloc(
								MAW_idx,nMAW_capacity*
							sizeof(unsigned int ));
					};
					MAW_idx[nMAWs-1]=MAW_buffer_idx;
					MAW_length=3;
					if(MAW_buffer_idx+MAW_length>MAW_buffer_size)
					{
						for(;;)						
						{
							MAW_buffer_size=(MAW_buffer_size+1)
								*alloc_growth_num/alloc_growth_denom;
							if(MAW_buffer_idx+MAW_length<=
								MAW_buffer_size)
								break;
						};
						MAW_buffer=(unsigned char *)realloc(
							MAW_buffer,MAW_buffer_size);				
					};
					MAW_buffer[MAW_buffer_idx++]=alpha4_to_ACGT[cl];
					MAW_buffer[MAW_buffer_idx++]=alpha4_to_ACGT[cr];
					MAW_buffer[MAW_buffer_idx++]=0;
				};
			
		};			

	};
	for(k=minlen>2?minlen-2:1;k<textlen;k++)
	{
//		printf("finding MAWs of length %d\n",k+2);
		for(i=0;i<textlen;i++)
			mark_vector[i]=0;
		for(i=0;i<textlen-k;i++)
		{
			if(mark_vector[i] || contains_N(&text[i],k))
				continue;
			for(cl=0;cl<4;cl++)
			{
				left_exists[cl]=0;
				right_exists[cl]=0;
			};	
			for(cl=0;cl<4;cl++)
				for(cr=0;cr<4;cr++)
					left_right_exists[cl][cr]=0;
			if(i>0)
				left_char=inv_alphabet[text[i-1]];
			else
				left_char=4;
			right_char=inv_alphabet[text[i+k]];
			left_exists[left_char]=1;
			right_exists[right_char]=1;
			left_right_exists[left_char][right_char]=1;
			is_repeat=0;
			for(j=i+1;j<textlen-k+1;j++)
				if(memcmp(&text[i],&text[j],k)==0)
				{
					mark_vector[j]=1;
					is_repeat=1;
					left_char=inv_alphabet[text[j-1]];
					if(j<textlen-k)
						right_char=inv_alphabet[text[j+k]];
					else
						right_char=4;
					left_exists[left_char]=1;
					right_exists[right_char]=1;
					left_right_exists[left_char][right_char]=1;
				};
			if(!is_repeat)
				continue;
			for(cl=0;cl<4;cl++)
			{
				if(!left_exists[cl])
					continue;
				for(cr=0;cr<4;cr++)
					if(right_exists[cr] && left_right_exists[cl][cr]==0)
					{
					// We have a minimal absent word. We output it. 					
						nMAWs++;
						if(nMAWs>nMAW_capacity)
						{
							nMAW_capacity=(
								nMAW_capacity+1)
								*alloc_growth_num/alloc_growth_denom;
							MAW_idx=(unsigned int*)realloc(
								MAW_idx,nMAW_capacity*
								sizeof(unsigned int ));
						};
						MAW_idx[nMAWs-1]=MAW_buffer_idx;
						MAW_length=k+3;
						if(MAW_buffer_idx+MAW_length>MAW_buffer_size)
						{
							for(;;)						
							{
								MAW_buffer_size=(MAW_buffer_size+1)
									*alloc_growth_num/alloc_growth_denom;
								if(MAW_buffer_idx+MAW_length<=
									MAW_buffer_size)
									break;
							};
							MAW_buffer=(unsigned char *)realloc(
								MAW_buffer,MAW_buffer_size);				
						};
						MAW_buffer[MAW_buffer_idx++]=alpha4_to_ACGT[cl];
						memcpy(&MAW_buffer[MAW_buffer_idx], 
							&text[i],k);
						MAW_buffer_idx+=k;
						MAW_buffer[MAW_buffer_idx++]=alpha4_to_ACGT[cr];
						MAW_buffer[MAW_buffer_idx++]=0;
							
					};
				
			};			
			
		};


	};
	free(mark_vector);
	MAW_buffer=(unsigned char *)realloc(
		MAW_buffer,MAW_buffer_idx);
	(*_MAW_buffer)=MAW_buffer;
	MAW_idx=(unsigned int *)realloc(
		MAW_idx,nMAWs*
		sizeof(unsigned int));
	MAW_ptr=(unsigned char **) malloc(nMAWs*sizeof(unsigned int));
	for(i=0;i<nMAWs;i++)
		MAW_ptr[i]=&MAW_buffer[MAW_idx[i]];
	free(MAW_idx);
	(*_MAW_ptr)=MAW_ptr;
	return nMAWs;
};

