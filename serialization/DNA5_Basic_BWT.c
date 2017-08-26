#include<stdio.h>
#include<stdlib.h>
#include"basic_bitvec.h"
#include"DNA5_Basic_BWT.h"
#include"vbyte.h"

#define use_dbwt 
#ifdef use_dbwt
#include"dbwt.h"
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
			//printf("The indexed_seq pointer to be freed is %d\n",Basic_BWT->indexed_BWT);
			free_basic_DNA5_seq(Basic_BWT->indexed_BWT);
			Basic_BWT->indexed_BWT=NULL;
		};
		free(Basic_BWT);
	};
};
int load_Basic_BWT_from_file(Basic_BWT_t ** _BBWT,char * file_name)
{
	Basic_BWT_t * BBWT= (Basic_BWT_t *) malloc(sizeof(Basic_BWT_t));
	unsigned char vbyte_buffer[35];
	unsigned int i;
	unsigned int read_size;
	unsigned long long tmp_int;
	unsigned int output_size;
	FILE * file=fopen(file_name,"r");
	if(file==0)
		return -1;
	for(i=0;i<35;i++)
		vbyte_buffer[i]=0;
	tmp_int=fread(vbyte_buffer,1,35,file);
	read_size=vbyte_decode(&tmp_int,vbyte_buffer);
//	printf("Text length is %llu\n",tmp_int);
	if(tmp_int>(3ul*(1<<30)) && tmp_int==0)
	{
		fclose(file);
		free(BBWT);
		(*_BBWT)=0;
		return -2;
	};
	BBWT->textlen=tmp_int;
	read_size+=vbyte_decode(&tmp_int,&vbyte_buffer[read_size]);
//	printf("primary idx is %llu\n",tmp_int);
	if(tmp_int>BBWT->textlen)
	{
		fclose(file);
		free(BBWT);
		(*_BBWT)=0;
		return -2;
	};
	BBWT->primary_idx=tmp_int;
	BBWT->char_base[0]=0;
	for(i=1;i<5;i++)
	{
		read_size+=vbyte_decode(&tmp_int,&vbyte_buffer[read_size]);
//		printf("char_base[i]=%llu\n",tmp_int);
		if(tmp_int>BBWT->textlen)
		{
			fclose(file);
			free(BBWT);
			(*_BBWT)=0;
			return -2;
		};
		BBWT->char_base[i]=tmp_int;
	};
//	printf("success to read integer components of BBWT\n");
	fseek(file,read_size, SEEK_SET);
	BBWT->indexed_BWT=new_basic_DNA5_seq(BBWT->textlen+1,&output_size);
	read_size+=load_DNA5_seq_from_opened_file(BBWT->indexed_BWT,BBWT->textlen+1,file);
	fclose(file);
	(*_BBWT)=BBWT;
	return read_size;
};
int save_Basic_BWT_to_file(Basic_BWT_t * BBWT,char * file_name)
{
	unsigned char vbyte_buffer[30];
	unsigned int i;
	unsigned int write_size;
	FILE * file=fopen(file_name,"w");
	if(file==0)
		return -1;
	for(i=0;i<30;i++)
		vbyte_buffer[i]=0;
	write_size=vbyte_encode(BBWT->textlen,vbyte_buffer);
	write_size+=vbyte_encode(BBWT->primary_idx,&vbyte_buffer[write_size]);
	for(i=1;i<5;i++)
		write_size+=vbyte_encode(BBWT->char_base[i],&vbyte_buffer[write_size]);
	fwrite(vbyte_buffer,write_size,1,file);
	write_size+=append_DNA5_seq_to_opened_file(BBWT->indexed_BWT,BBWT->textlen+1,file);
	fclose(file);
	return write_size;
};

int cmp_BBWTs(Basic_BWT_t * BBWT1, Basic_BWT_t * BBWT2)
{
	unsigned int i;
	unsigned int j;
	unsigned int counts1[4],counts2[4];
	if(BBWT1->textlen!=BBWT2->textlen)
		return -4;
	if(BBWT1->primary_idx!=BBWT2->primary_idx)
		return -3;
	for(i=0;i<5;i++)
		if(BBWT1->char_base[i]!=BBWT2->char_base[i])
			return -2;
	for(i=0;i<BBWT1->textlen+1;i++)
	{
		if(DNA5_extract_char(BBWT1->indexed_BWT,i)!=
			DNA5_extract_char(BBWT2->indexed_BWT,i))
				return -1;
		DNA5_get_char_pref_counts(counts1,BBWT1->indexed_BWT,i);
		DNA5_get_char_pref_counts(counts2,BBWT2->indexed_BWT,i);
		for(j=0;j<4;j++)
			if(counts1[j]!=counts2[j])
				return -1;
	};
	return 0;

};
int load_Basic_BWT_from_callback(Basic_BWT_t * BBWT, unsigned int BWT_length,
	 unsigned int primary_idx, BWT_load_callback_t load_callback)
{
	unsigned char * buffer;
	unsigned char triplet[3];
	unsigned int buffer_len;
	unsigned int curr_pos;
	unsigned int curr_triplet_pos;
	unsigned int i;
	unsigned int j;
	unsigned int char_count[5];
	unsigned char curr_char;
	unsigned int output_size;
	BBWT->primary_idx=primary_idx;
	BBWT->textlen=BWT_length;
	BBWT->indexed_BWT=new_basic_DNA5_seq(BBWT->textlen+1,&output_size);
	for(i=0;i<5;i++)
		char_count[i]=0;
	curr_pos=0;
	while(1) 
	{
		buffer_len=load_callback(&buffer);
		if(buffer_len==0)
			break;
		i=0;
		while((curr_pos%3)>0)
		{	
			curr_char=DNA_5_alpha_trans_table[buffer[i]];
			char_count[curr_char]++;
			DNA5_set_char(BBWT->indexed_BWT,curr_pos,curr_char);
			curr_pos++;
			i++;
			if(i==buffer_len || curr_pos==BWT_length)
				break;
		};
		if(curr_pos==BWT_length)
			break;
		if(i==buffer_len)
			continue;
		curr_triplet_pos=curr_pos/3;
		for(;i<buffer_len && curr_pos<BWT_length;)
		{
			curr_triplet_pos++;
			curr_pos+=3;
			for(j=0;j<3;j++,i++)
			{
				curr_char=DNA_5_alpha_trans_table[buffer[i]];
				triplet[j]=curr_char;
				char_count[curr_char]++;
			};
			DNA5_set_triplet_at(BBWT->indexed_BWT,curr_triplet_pos,triplet);
		};
		curr_pos=curr_triplet_pos*3;
		while(i<buffer_len && curr_pos<BWT_length)
		{	
			curr_char=DNA_5_alpha_trans_table[buffer[i]];
			char_count[curr_char]++;
			DNA5_set_char(BBWT->indexed_BWT,curr_pos,curr_char);
			curr_pos++;
			i++;
		};
		if(curr_pos==BWT_length)
			break;
	};
	BBWT->char_base[0]=0;
	BBWT->char_base[1]=char_count[0]-1;
	for(i=2;i<5;i++)
		BBWT->char_base[i]=BBWT->char_base[i-1]+char_count[i-1];
	complete_basic_DNA5_seq(BBWT->indexed_BWT,BWT_length);
	return 0;
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

