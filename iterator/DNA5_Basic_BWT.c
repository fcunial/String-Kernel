#include <stdio.h>
#include <stdlib.h>
#include "basic_bitvec.h"
#include "DNA5_Basic_BWT.h"

#define use_dbwt 
#ifdef use_dbwt
#include "../dbwt/dbwt.h"
//#define free_text_in_bwt_build 0
//#define free_text_in_bwt_build 1
#else
#include<divsufsort.h>
#endif


Basic_BWT_t * new_Basic_BWT()
{
	return (Basic_BWT_t *) calloc(1,sizeof(Basic_BWT_t));
};
void free_Basic_BWT(Basic_BWT_t * Basic_BWT)
{
	if(Basic_BWT)
	{
		if(Basic_BWT->indexed_BWT)
		{
			free_basic_DNA5_seq(Basic_BWT->indexed_BWT);
			Basic_BWT->indexed_BWT=NULL;
		};
		free(Basic_BWT);
	};
};


unsigned int LF_map(unsigned int in_pos,Basic_BWT_t * Basic_BWT)
{
	unsigned int counts0[4];
	unsigned int c;
	in_pos++;
	if(in_pos==Basic_BWT->primary_idx)
		return 0xffffffff;
	c=DNA5_extract_char(Basic_BWT->indexed_BWT,in_pos);
	DNA5_get_char_pref_counts(counts0,Basic_BWT->indexed_BWT,in_pos-1);
	if(c==0 && in_pos>Basic_BWT->primary_idx)
		counts0[0]--;
	return Basic_BWT->char_base[c]+counts0[c];
};

int Basic_BWT_batch_extract(unsigned int * bitvector,
		unsigned int nelements,unsigned int * output_vector,
		Basic_BWT_t * Basic_BWT)
{
	unsigned char c;
	unsigned int out_idx=0;
	unsigned int SA_pos;
	unsigned int txt_pos=Basic_BWT->textlen;
	unsigned int counts[5];
	c=DNA5_extract_char(Basic_BWT->indexed_BWT,0);
	counts[c]=0;
	do
	{
		txt_pos--;
//		printf("text pos is %d\n",txt_pos);
//		printf("character is %d\n",c);
		SA_pos=counts[c]+Basic_BWT->char_base[c];
//		printf("SA pos is %d\n",SA_pos);
		if(ismarkedbit(SA_pos,bitvector))
		{
			output_vector[out_idx++]=txt_pos;	
			output_vector[out_idx++]=SA_pos;
			if(out_idx==2*nelements)
				break;
		};
//		printf("SA[%d]=%d\n",SA_pos,txt_pos);
		c=DNA5_extract_char(Basic_BWT->indexed_BWT,SA_pos+1);
		DNA5_get_char_pref_counts(counts,Basic_BWT->indexed_BWT,SA_pos);
		if(c==0 && SA_pos>=Basic_BWT->primary_idx)
			counts[0]--;
		else if(c==4)
		{
			counts[4]=SA_pos+1-(counts[0]+counts[1]+counts[2]+counts[3]);
			//printf("counts[4]=%d for position %d\n",counts[4],SA_pos);
		};
	} while(txt_pos>0);
	return out_idx/2;
};

int Backward_step(unsigned int * in_interval,unsigned int * out_interval,
		unsigned char c,Basic_BWT_t * Basic_BWT)
{
	unsigned int counts0[4],counts1[4];
	unsigned int _in_interval[2];
	_in_interval[0]=in_interval[0]+1;
	_in_interval[1]=in_interval[1]+1;
/*	printf("backward step called with in_interval[0]=%d and in_interval[1]=%d \n",
		in_interval[0],in_interval[1]);*/
	DNA5_get_char_pref_counts(counts0,Basic_BWT->indexed_BWT,_in_interval[0]-1);
	DNA5_get_char_pref_counts(counts1,Basic_BWT->indexed_BWT,_in_interval[1]);
	if(c==0)
	{
		if(_in_interval[0]>Basic_BWT->primary_idx)
		{
			counts0[0]--;
			counts1[0]--;
		} else
		{
			if(_in_interval[1]>=Basic_BWT->primary_idx)
				counts1[0]--;
		};
	};
	if(counts0[c]==counts1[c])
		return -1;
	out_interval[0]=Basic_BWT->char_base[c]+counts0[c];
	out_interval[1]=Basic_BWT->char_base[c]+counts1[c]-1;
	return 0;
};

int patt_count(unsigned char *P,unsigned int m,
		Basic_BWT_t * Basic_BWT,unsigned int _SA_interval[2])
{
	int i;
	unsigned int op_res;
	unsigned char c;
	unsigned int SA_interval[2];
	if(m==0)
		return -1;
	c=DNA_5_alpha_trans_table[P[m-1]];
	if(c>=4)
		return -1;
	SA_interval[0]=Basic_BWT->char_base[c];
	if(SA_interval[0]==Basic_BWT->char_base[c+1])
		return -1;
	SA_interval[1]=Basic_BWT->char_base[c+1]-1;
/*	printf("Initial SA interval for character %d (%c) is (%d,%d)\n",
		c,P[m-1],SA_interval[0],SA_interval[1]);*/
	if(m>1)
		for(i=m-2;i>=0;i--)
		{
			c=DNA_5_alpha_trans_table[P[i]];
			if(c>=4)
				return -1;
			op_res=Backward_step(SA_interval,SA_interval,c,Basic_BWT);
			if(op_res!=0)
				return -1;
		};
	_SA_interval[0]=SA_interval[0];
	_SA_interval[1]=SA_interval[1];
	return 0;
};


Basic_BWT_t * Build_BWT_index_from_text(unsigned char * text,
	unsigned int textlen,unsigned int options)
{
	Basic_BWT_t * Basic_BWT=new_Basic_BWT();
	unsigned int i;
	unsigned int char_count[4];
	unsigned char * temp_BWT=0;
#ifdef use_dbwt
	temp_BWT=dbwt_bwt(text,textlen,&Basic_BWT->primary_idx,
		options);
	temp_BWT[Basic_BWT->primary_idx]='A';
//	printf("The computed bwt is : ");
//	for(i=0;i<=textlen;i++)
//		printf("%c",temp_BWT[i]);
//	printf("\n");
#else
	int build_res;
	int * SA_array;
	unsigned int SA_val;
	temp_BWT=(unsigned char *)malloc(textlen+1);
	SA_array=(int *)malloc(textlen*sizeof(int));
	build_res=divsufsort(text,SA_array,textlen);	
	if(build_res!=0)
	{
		free_Basic_BWT(Basic_BWT);
		Basic_BWT=NULL;
		goto return_point;
	}
	temp_BWT[0]=text[textlen-1];
	for(i=0;i<textlen;i++)
	{
		SA_val=SA_array[i];
		if(SA_val==0)
		{
			Basic_BWT->primary_idx=i+1;		
			temp_BWT[i+1]='A';
		}
		else
			temp_BWT[i+1]=text[SA_val-1];
	};
#endif
	Basic_BWT->indexed_BWT=build_basic_DNA5_seq(temp_BWT,textlen+1,
		&Basic_BWT->size,char_count);
	if(Basic_BWT->indexed_BWT==NULL)
	{
		free_Basic_BWT(Basic_BWT);
		Basic_BWT=NULL;
		goto return_point;
	}
	Basic_BWT->char_base[0]=0;
	Basic_BWT->char_base[1]=char_count[0]-1;
	for(i=2;i<5;i++)
		Basic_BWT->char_base[i]=Basic_BWT->char_base[i-1]+char_count[i-1];
//	printf("char bases are %d,%d,%d,%d,%d\n",Basic_BWT->char_base[0],
//		Basic_BWT->char_base[1],Basic_BWT->char_base[2],
//		Basic_BWT->char_base[3],Basic_BWT->char_base[4]);
	Basic_BWT->textlen=textlen;
return_point:
#ifndef use_dbwt
	free(SA_array);
#endif
	free(temp_BWT);
	return Basic_BWT;
};

