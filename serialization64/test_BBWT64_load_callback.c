#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include <string.h>
#include"mt19937ar.h"
#include"DNA5_Basic_BWT64.h"
#include"dbwt.h"
#define textsize ((1<<20))
//#define textsize (1024)
static void naive_str_print(unsigned char * str1, unsigned int str_len)
{
	unsigned int i;
	for(i=0;i<str_len;i++)
		printf("%c",str1[i]);

};

static unsigned char local_DNA_trans_table[5]={'A','C','G','T','Z'};
typedef struct 
{
	FILE * BWT_file;
	unsigned char * buffer;
	unsigned int block_size;
	unsigned int end_flag;
} callback_state_t;


unsigned long long BWT_load_callback (unsigned char ** _buffer, void * intern_state)
{
	unsigned int len;
	callback_state_t * state=(callback_state_t*) intern_state;
	if(state->end_flag)
	{
		(*_buffer)=0;
		return 0;
	};
	len=fread(state->buffer,1,state->block_size,state->BWT_file);
	if(len<state->block_size)
		state->end_flag=1;
	(*_buffer)=state->buffer;
	return len;
};
int test_build_BBWT_from_text(unsigned char * text, unsigned int len, 
	Basic_BWT64_t ** _BBWT)
{
	FILE * out_file;	
	unsigned int primary_idx;
	unsigned char * temp_BWT=0;
	unsigned long long written_bytes;
	const unsigned int block_size=(1<<20);
	int ret_code;
	callback_state_t state;
	int print_string=0;
	out_file=fopen("tmp_file.bwt","w");
	if(out_file==0)
		return -1;
	temp_BWT=dbwt_bwt(text,len,&primary_idx,Basic_bwt_no_free_text);
	if(temp_BWT==0)
	{
		fclose(out_file);
		return -2; 
	};
	temp_BWT[primary_idx]='A';
	if(print_string)
	{
		printf("dbwt string is : ");
		naive_str_print(temp_BWT,len+1);
		printf("\n");
	};
	written_bytes=fwrite(temp_BWT,1,len+1,out_file);
	free(temp_BWT);
	if(written_bytes<len+1)
	{
		printf("wrote %llu bytes instead of %d\n",written_bytes,len+1);
		fclose(out_file);
		return -3;
	};
	fclose(out_file);
	state.buffer=malloc(block_size);
	state.block_size=block_size;
	state.end_flag=0;
	state.BWT_file=fopen("tmp_file.bwt","r");
	if(state.BWT_file==0 || state.buffer==0)
	{
		if(state.buffer!=0)
			free(state.buffer);
		if(state.BWT_file!=0)
			fclose(state.BWT_file);
		return -4;
	};
//	printf("We are going to call the callback\n");
	ret_code=load_Basic_BWT64_from_callback(_BBWT,len,primary_idx,BWT_load_callback,&state);
	free(state.buffer);
	fclose(state.BWT_file);
	if(ret_code!=0)
	{
		printf("error in call back building function with code %d\n",ret_code);
		return -5;
	};
//	printf("No apparent error in callback\n");
	return 0;
};


int main()
{
	unsigned char * text;
	unsigned int textlen;
	unsigned int i;
	unsigned char tmp_char;
	unsigned int print_string=0;
	int err_code;
	Basic_BWT64_t * BBWT1;
	Basic_BWT64_t * BBWT2;
	textlen=textsize;
	text=(unsigned char *)malloc(textlen+1);
	text[textlen]='0';
// Generate a random text.
	init_genrand(time(0));
	for(i=0;i<textlen;i++)
	{
		tmp_char=genrand_int32()%4;
		text[i]=local_DNA_trans_table[tmp_char];
	};
	printf("Text len is %d\n",textlen);
	if(print_string)
	{
		printf("Input string is : ");
		naive_str_print(text,textlen);
		printf("\n");
	};
// Build a BWT index on the text. 
	BBWT1=Build_BWT_index_from_text64(text,textlen,Basic_bwt_no_free_text);
	err_code=test_build_BBWT_from_text(text,textlen,&BBWT2);
	if(err_code!=0)
	{
		printf("Failure in building the second BBWT based on callback\n");
		printf("Error code was %d\n",err_code);
		goto end_label;
	};
	if((err_code=cmp_BBWT64s(BBWT1,BBWT2)))
		printf("Error in callback bwt load mechanism with error code %d\n",err_code);
	else
		printf("Success in callback bwt load mechanism\n");
	free_Basic_BWT64(BBWT2);
end_label:
	free_Basic_BWT64(BBWT1);
	free(text);
	return 0; 
};
