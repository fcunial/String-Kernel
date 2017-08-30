#include<stdio.h>
#include<stdlib.h>
#include"vbyte.h"
#include"basic_bitvec64.h"
#include"DNA5_Basic_BWT64.h"

#define use_dbwt 
#ifdef use_dbwt
#include"dbwt.h"
//#define free_text_in_bwt_build 0
//#define free_text_in_bwt_build 1
#else
#include<divsufsort.h>
#endif


Basic_BWT64_t * new_Basic_BWT64()
{
	return (Basic_BWT64_t *) calloc(1,sizeof(Basic_BWT64_t));
};
long long load_Basic_BWT64_from_file(Basic_BWT64_t ** _BBWT,char * file_name)
{
	Basic_BWT64_t * BBWT =new_Basic_BWT64();
	unsigned char vbyte_buffer[60];
	unsigned int i;
	unsigned long long read_size;
	FILE * file=fopen(file_name,"r");
	if(file==0)
	{
		free(BBWT);
		(*_BBWT)=0;
		return -1;
	};
	for(i=0;i<60;i++)
		vbyte_buffer[i]=0;
	read_size=fread(vbyte_buffer,1,60,file);// Just to shut a compiler warning
	read_size=vbyte_decode(&BBWT->textlen,vbyte_buffer);
	read_size+=vbyte_decode(&BBWT->primary_idx,&vbyte_buffer[read_size]);
	BBWT->char_base[0]=0;
	for(i=1;i<5;i++)
		read_size+=vbyte_decode(&BBWT->char_base[i],&vbyte_buffer[read_size]);
	if(BBWT->textlen==0)
	{
		fclose(file);
		free(BBWT);
		(*_BBWT)=0;
		return -2;
	};
	fseek(file,read_size, SEEK_SET);
	BBWT->indexed_BWT=new_basic_DNA5_seq64(BBWT->textlen+1);
	read_size+=load_DNA5_seq64_from_opened_file(BBWT->indexed_BWT,file);
	fclose(file);
	(*_BBWT)=BBWT;
	return read_size;
};
long long save_Basic_BWT64_to_file(Basic_BWT64_t * BBWT,char * file_name)
{
	unsigned char vbyte_buffer[60];
	unsigned int i;
	unsigned int write_size;
	FILE * file=fopen(file_name,"w");
	if(file==0)
		return -1;
	for(i=0;i<60;i++)
		vbyte_buffer[i]=0;
	write_size=vbyte_encode(BBWT->textlen,vbyte_buffer);
	write_size+=vbyte_encode(BBWT->primary_idx,&vbyte_buffer[write_size]);
	for(i=1;i<5;i++)
		write_size+=vbyte_encode(BBWT->char_base[i],&vbyte_buffer[write_size]);
	fwrite(vbyte_buffer,1,write_size,file);
	BBWT->indexed_BWT->length=BBWT->textlen;
	write_size+=append_DNA5_seq64_to_opened_file(BBWT->indexed_BWT,file);
	fclose(file);
	return write_size;
};
int cmp_BBWT64s(Basic_BWT64_t * BBWT1, Basic_BWT64_t * BBWT2)
{
	unsigned long long i;
	unsigned int j;
	unsigned long long counts1[4],counts2[4];
	if(BBWT1->textlen!=BBWT2->textlen)
		return -5;
	if(BBWT1->primary_idx!=BBWT2->primary_idx)
		return -4;
	for(i=0;i<5;i++)
		if(BBWT1->char_base[i]!=BBWT2->char_base[i])
			return -3;
	for(i=0;i<BBWT1->textlen+1;i++)
	{
		if(DNA5_extract_char64(BBWT1->indexed_BWT,i)!=
			DNA5_extract_char64(BBWT2->indexed_BWT,i))
			{
				printf("character %llu was different\n",i);
				printf("It was %d in BBWT1 and %d in BBWT2\n",
				DNA5_extract_char64(BBWT1->indexed_BWT,i),
				DNA5_extract_char64(BBWT2->indexed_BWT,i));
				return -2;
			};
		DNA5_get_char_pref_counts64(counts1,BBWT1->indexed_BWT,i);
		DNA5_get_char_pref_counts64(counts2,BBWT2->indexed_BWT,i);
		for(j=0;j<4;j++)
			if(counts1[j]!=counts2[j])
				return -1;
	};
	return 0;
};

int load_Basic_BWT64_from_callback(Basic_BWT64_t ** _BBWT, unsigned long long BWT_length,
	unsigned long long primary_idx, BWT64_load_callback_t load_callback, 
 	void * intern_state)
{
	unsigned char * buffer;
	unsigned char * triplet;
	unsigned long long buffer_len;
	unsigned long long curr_pos;
	unsigned long long curr_triplet_pos;
	unsigned long long i;
	unsigned int j;
	unsigned long long char_count[5];
	unsigned char curr_char;
	Basic_BWT64_t * BBWT=(Basic_BWT64_t *) malloc(sizeof(Basic_BWT64_t));
	BBWT->primary_idx=primary_idx;
	BBWT->textlen=BWT_length;
	BBWT->indexed_BWT=new_basic_DNA5_seq64(BBWT->textlen+1);
	for(i=0;i<5;i++)
		char_count[i]=0;
	curr_pos=0;
	while(1) 
	{
		buffer_len=load_callback(&buffer,intern_state);
		if(buffer_len==0 || buffer==0)
			break;
		i=0;
		while((curr_pos%3)>0)
		{	
			curr_char=DNA_5_alpha_trans_table[buffer[i]];
			char_count[curr_char]++;
			DNA5_set_char64(BBWT->indexed_BWT,curr_pos,curr_char);
			curr_pos++;
			i++;
			if(i==buffer_len || curr_pos==BWT_length+1)
				break;
		};
		if(curr_pos==BWT_length+1)
			break;
		if(i==buffer_len)
			continue;
		curr_triplet_pos=curr_pos/3;
		while(i+2<buffer_len && curr_pos+2<BWT_length+1)
		{
			triplet=&buffer[i];
			DNA5_set_triplet_at64(BBWT->indexed_BWT,curr_triplet_pos,triplet);
			for(j=0;j<3;j++,i++)
			{
				curr_char=DNA_5_alpha_trans_table[buffer[i]];
				char_count[curr_char]++;
			};
			curr_triplet_pos++;
			curr_pos+=3;
		};
		while(i<buffer_len && curr_pos<BWT_length+1)
		{	
			curr_char=DNA_5_alpha_trans_table[buffer[i]];
			char_count[curr_char]++;
			DNA5_set_char64(BBWT->indexed_BWT,curr_pos,curr_char);
			curr_pos++;
			i++;
		};
		if(curr_pos==BWT_length+1)
			break;
	};
	BBWT->char_base[0]=0;
	BBWT->char_base[1]=char_count[0]-1;
	for(i=2;i<5;i++)
		BBWT->char_base[i]=BBWT->char_base[i-1]+char_count[i-1];
	complete_basic_DNA5_seq64(BBWT->indexed_BWT);
	(*_BBWT)=BBWT;
	return 0;
};

void free_Basic_BWT64(Basic_BWT64_t * Basic_BWT)
{
	if(Basic_BWT)
	{
		if(Basic_BWT->indexed_BWT)
		{
			free_basic_DNA5_seq64(Basic_BWT->indexed_BWT);
			Basic_BWT->indexed_BWT=NULL;
		};
		free(Basic_BWT);
	};
};


unsigned long long LF_map64(unsigned long long in_pos,Basic_BWT64_t * Basic_BWT)
{
	unsigned long long counts0[4];
	unsigned int c;
	in_pos++;
	if(in_pos==Basic_BWT->primary_idx)
		return 0xffffffff;
	c=DNA5_extract_char64(Basic_BWT->indexed_BWT,in_pos);
	DNA5_get_char_pref_counts64(counts0,Basic_BWT->indexed_BWT,in_pos-1);
	if(c==0 && in_pos>Basic_BWT->primary_idx)
		counts0[0]--;
	return Basic_BWT->char_base[c]+counts0[c];
};

int Basic_BWT_batch_extract64(unsigned int * bitvector,
		unsigned long long nelements,unsigned long long * output_vector,
		Basic_BWT64_t * Basic_BWT)
{
	unsigned char c;
	unsigned long long out_idx=0;
	unsigned long long SA_pos;
	unsigned long txt_pos=Basic_BWT->textlen;
	unsigned long long counts[5];
	c=DNA5_extract_char64(Basic_BWT->indexed_BWT,0);
	counts[c]=0;
	do
	{
		txt_pos--;
//		printf("text pos is %d\n",txt_pos);
//		printf("character is %d\n",c);
		SA_pos=counts[c]+Basic_BWT->char_base[c];
//		printf("SA pos is %d\n",SA_pos);
		if(ismarkedbit64(SA_pos,bitvector))
		{
			output_vector[out_idx++]=txt_pos;	
			output_vector[out_idx++]=SA_pos;
			if(out_idx==2*nelements)
				break;
		};
//		printf("SA[%d]=%d\n",SA_pos,txt_pos);
		c=DNA5_extract_char64(Basic_BWT->indexed_BWT,SA_pos+1);
		DNA5_get_char_pref_counts64(counts,Basic_BWT->indexed_BWT,SA_pos);
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

int Backward_step64(unsigned long long * in_interval,unsigned long long * out_interval,
		unsigned char c,Basic_BWT64_t * Basic_BWT)
{
	unsigned long long counts0[4],counts1[4];
	unsigned long long _in_interval[2];
	_in_interval[0]=in_interval[0]+1;
	_in_interval[1]=in_interval[1]+1;
/*	printf("backward step called with in_interval[0]=%d and in_interval[1]=%d \n",
		in_interval[0],in_interval[1]);*/
	DNA5_get_char_pref_counts64(counts0,Basic_BWT->indexed_BWT,_in_interval[0]-1);
	DNA5_get_char_pref_counts64(counts1,Basic_BWT->indexed_BWT,_in_interval[1]);
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

unsigned long long patt_count64(unsigned char *P,unsigned long long m,
		Basic_BWT64_t * Basic_BWT,unsigned long long _SA_interval[2])
{
	int i;
	unsigned long long op_res;
	unsigned char c;
	unsigned long long SA_interval[2];
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
			op_res=Backward_step64(SA_interval,SA_interval,c,Basic_BWT);
			if(op_res!=0)
				return -1;
		};
	_SA_interval[0]=SA_interval[0];
	_SA_interval[1]=SA_interval[1];
	return 0;
};
// The bwt should be the bwt of a text from alphabet {a,c,g,t,z} terminated 
// '#' which a character lexicographically smaller than 'a'. The primary_idx 
// should be the position of an 'a' that represents the position of character #. 
Basic_BWT64_t * Build_BWT_index_from_bwt64(unsigned char * bwt,
	unsigned long long bwtlen,unsigned long long primary_idx)
{
	unsigned int i;
	Basic_BWT64_t * Basic_BWT=new_Basic_BWT64();	
	unsigned long long char_count[4];
	Basic_BWT->indexed_BWT=build_basic_DNA5_seq64(bwt,bwtlen,char_count);
	if(Basic_BWT->indexed_BWT==NULL)
	{
		free_Basic_BWT64(Basic_BWT);
		return 0;
	};
	Basic_BWT->textlen=bwtlen-1;
	Basic_BWT->primary_idx=primary_idx;
	Basic_BWT->char_base[0]=0;
	Basic_BWT->char_base[1]=char_count[0]-1;
	for(i=2;i<5;i++)
		Basic_BWT->char_base[i]=Basic_BWT->char_base[i-1]+char_count[i-1];
	return Basic_BWT;


};

Basic_BWT64_t * Build_BWT_index_from_text64(unsigned char * text,
	unsigned long long textlen,unsigned int options)
{
	Basic_BWT64_t * Basic_BWT=new_Basic_BWT64();
	unsigned int i;
	unsigned long long char_count[4];
	unsigned char * temp_BWT=0;
	unsigned int short_primary_idx;
//	unsigned long long primary_idx;
#ifdef use_dbwt
	temp_BWT=dbwt_bwt(text,textlen,&short_primary_idx,
		options);
	Basic_BWT->primary_idx=short_primary_idx;
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
	Basic_BWT->indexed_BWT=build_basic_DNA5_seq64(temp_BWT,textlen+1,char_count);
	if(Basic_BWT->indexed_BWT==NULL)
	{
		free_Basic_BWT64(Basic_BWT);
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

