#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include"SLT_joint64.h"

typedef struct 
{
	unsigned long long string_depth;
	unsigned long long interval_start1;
	unsigned long long interval_start2;
	unsigned long long child_freqs1[6];
	unsigned long long child_freqs2[6];
//	unsigned int revbwt_start1;
//	unsigned int revbwt_start2;
	unsigned long long interval_size1;
	unsigned long long interval_size2;
	unsigned char WL_char;
} SLT64_stack_item_t;



static inline void swap2_stack_items(SLT64_stack_item_t * SLT_stack_item1,
			SLT64_stack_item_t * SLT_stack_item2)
{
	unsigned int i;
	SLT_stack_item1->string_depth^=SLT_stack_item2->string_depth;
	SLT_stack_item1->interval_start1^=SLT_stack_item2->interval_start1;
	SLT_stack_item1->interval_start2^=SLT_stack_item2->interval_start2;
	SLT_stack_item1->interval_size1^=SLT_stack_item2->interval_size1;
	SLT_stack_item1->interval_size2^=SLT_stack_item2->interval_size2;
	SLT_stack_item1->WL_char^=SLT_stack_item2->WL_char;

	SLT_stack_item2->string_depth^=SLT_stack_item1->string_depth;
	SLT_stack_item2->interval_start1^=SLT_stack_item1->interval_start1;
	SLT_stack_item2->interval_start2^=SLT_stack_item1->interval_start2;
	SLT_stack_item2->interval_size1^=SLT_stack_item1->interval_size1;
	SLT_stack_item2->interval_size2^=SLT_stack_item1->interval_size2;
	SLT_stack_item2->WL_char^=SLT_stack_item1->WL_char;

	SLT_stack_item1->string_depth^=SLT_stack_item2->string_depth;
	SLT_stack_item1->interval_start1^=SLT_stack_item2->interval_start1;
	SLT_stack_item1->interval_start2^=SLT_stack_item2->interval_start2;
	SLT_stack_item1->interval_size1^=SLT_stack_item2->interval_size1;
	SLT_stack_item1->interval_size2^=SLT_stack_item2->interval_size2;
	SLT_stack_item1->WL_char^=SLT_stack_item2->WL_char;

	for(i=0;i<3;i++)
	{
		SLT_stack_item1->child_freqs1[i*2]^=SLT_stack_item2->child_freqs1[i*2];
		SLT_stack_item1->child_freqs1[i*2+1]^=SLT_stack_item2->child_freqs1[i*2+1];
		SLT_stack_item2->child_freqs1[i*2]^=SLT_stack_item1->child_freqs1[i*2];
		SLT_stack_item2->child_freqs1[i*2+1]^=SLT_stack_item1->child_freqs1[i*2+1];
		SLT_stack_item1->child_freqs1[i*2]^=SLT_stack_item2->child_freqs1[i*2];
		SLT_stack_item1->child_freqs1[i*2+1]^=SLT_stack_item2->child_freqs1[i*2+1];

		SLT_stack_item1->child_freqs2[i*2]^=SLT_stack_item2->child_freqs2[i*2];
		SLT_stack_item1->child_freqs2[i*2+1]^=SLT_stack_item2->child_freqs2[i*2+1];
		SLT_stack_item2->child_freqs2[i*2]^=SLT_stack_item1->child_freqs2[i*2];
		SLT_stack_item2->child_freqs2[i*2+1]^=SLT_stack_item1->child_freqs2[i*2+1];
		SLT_stack_item1->child_freqs2[i*2]^=SLT_stack_item2->child_freqs2[i*2];
		SLT_stack_item1->child_freqs2[i*2+1]^=SLT_stack_item2->child_freqs2[i*2+1];

	};
};

#define min_SLT_stack_size 4

void SLT64_joint_execute_iterator(SLT64_joint_iterator_t * SLT_iterator)
{
	unsigned long long curr_stack_size=min_SLT_stack_size;
	unsigned long long curr_stack_idx=0;
	SLT64_joint_params_t SLT_params;
	unsigned long long string_depth;
	unsigned long long char_pref_counts1[28];
	unsigned long long char_pref_counts2[28];
	unsigned long long last_char_pref_counts1[7];
	unsigned long long pref_count_query_points1[7];
	unsigned long long last_char_pref_counts2[7];
	unsigned long long pref_count_query_points2[7];
	unsigned long long last_char_pref_count1;
	unsigned long long last_char_pref_count2;
	unsigned int i,j,k;
	unsigned int extension_exists1;
	unsigned int extension_exists2;
	unsigned int nchildren1;
	unsigned int nchildren2;
	unsigned int nchildren;
	unsigned int includes_EOT_char1;
	unsigned int includes_EOT_char2;
	unsigned long long last_char_freq1;
	unsigned long long last_char_freq2;
	unsigned int npref_query_points1;
	unsigned int npref_query_points2;
	unsigned long long interval_size1;
	unsigned long long interval_size2;
	unsigned long long sum_interval_size;
	unsigned long long max_sum_interval_size;
	unsigned int nexplicit_WL;
	unsigned int max_interval_idx;
	unsigned long long ntraversed_nodes=0;
//	unsigned int revbwt_start;
	Basic_BWT64_t * BBWT1=SLT_iterator->BBWT1;
	Basic_BWT64_t * BBWT2=SLT_iterator->BBWT2;
	unsigned int options=SLT_iterator->options;
	unsigned int j1,j2;
// Allocate the stack
//	printf("length of the text given to SLT program is %d\n",
//		BBWT->textlen);

	SLT64_stack_item_t * SLT_stack=(SLT64_stack_item_t *)
		malloc((min_SLT_stack_size+1)*sizeof(SLT64_stack_item_t));
// Push the root node on the stack
	SLT_stack[curr_stack_idx].WL_char=0;
	SLT_stack[curr_stack_idx].string_depth=0;
	SLT_stack[curr_stack_idx].interval_start1=0;
	SLT_stack[curr_stack_idx].interval_start2=0;
	SLT_stack[curr_stack_idx].child_freqs1[0]=1;
	SLT_stack[curr_stack_idx].child_freqs2[0]=1;
	SLT_stack[curr_stack_idx].interval_size1=BBWT1->textlen;
	SLT_stack[curr_stack_idx].interval_size2=BBWT1->textlen;
	for(i=1;i<5;i++)
	{
		SLT_stack[curr_stack_idx].child_freqs1[i]=
			BBWT1->char_base[i]-BBWT1->char_base[i-1];
		SLT_stack[curr_stack_idx].child_freqs2[i]=
			BBWT2->char_base[i]-BBWT2->char_base[i-1];
	}
	SLT_stack[curr_stack_idx].child_freqs1[5]=BBWT1->textlen-BBWT1->char_base[4];
	SLT_stack[curr_stack_idx].child_freqs2[5]=BBWT2->textlen-BBWT2->char_base[4];
//	for(i=0;i<=5;i++)
//		printf("freq %d is %d\n",i,
//			SLT_stack[curr_stack_idx].child_freqs[i]);
	curr_stack_idx++;
// Enter the main loop
	do
	{
		ntraversed_nodes++;
// Pop a node from the stacj 
		curr_stack_idx--;
// Set the first rank query points
//		printf("start generating query points\n");
		pref_count_query_points1[0]=
			SLT_stack[curr_stack_idx].interval_start1-1;
		if((pref_count_query_points1[0]+1)==0)
			for(i=0;i<4;i++)
				char_pref_counts1[i]=0;
		pref_count_query_points2[0]=
			SLT_stack[curr_stack_idx].interval_start2-1;
		if((pref_count_query_points2[0]+1)==0)
			for(i=0;i<4;i++)
				char_pref_counts2[i]=0;

//		printf("first query point is %d\n",pref_count_query_points[0]);
// Set the data related to the top node to be given as parameter to the call back
// function. Also set the remaining rank query points. 

		SLT_params.WL_char=SLT_stack[curr_stack_idx].WL_char;
		SLT_params.string_depth=SLT_stack[curr_stack_idx].string_depth;
		SLT_params.interval_size1=SLT_stack[curr_stack_idx].interval_size1;
		SLT_params.interval_size2=SLT_stack[curr_stack_idx].interval_size2;
		SLT_params.right_extension_bitmap1=0;
		SLT_params.right_extension_bitmap2=0;
		for(i=1,j=1;i<7;i++)
		{
			if(SLT_stack[curr_stack_idx].child_freqs1[i-1])
			{
				SLT_params.right_extension_bitmap1|=(1<<(i-1));
				pref_count_query_points1[j]=
					pref_count_query_points1[j-1]+
					SLT_stack[curr_stack_idx].child_freqs1[i-1];
				j++;
			};
		};
		npref_query_points1=j;
		SLT_params.nright_extensions1=npref_query_points1-1;
		for(i=1,j=1;i<7;i++)
		{
			if(SLT_stack[curr_stack_idx].child_freqs2[i-1])
			{
				SLT_params.right_extension_bitmap2|=(1<<(i-1));
				pref_count_query_points2[j]=
					pref_count_query_points2[j-1]+
					SLT_stack[curr_stack_idx].child_freqs2[i-1];
				j++;
			};
		};
		npref_query_points2=j;
		SLT_params.nright_extensions2=npref_query_points2-1;

		if(SLT_params.nright_extensions1)
		{
			if((pref_count_query_points1[0]+1)==0)
			{
				npref_query_points1--;
				DNA5_multipe_char_pref_counts64(BBWT1->indexed_BWT,npref_query_points1,
					&pref_count_query_points1[1],&char_pref_counts1[4]);
				npref_query_points1++;
			}
			else
				DNA5_multipe_char_pref_counts64(BBWT1->indexed_BWT,npref_query_points1,
					pref_count_query_points1,char_pref_counts1);
			includes_EOT_char1=((BBWT1->primary_idx>=(pref_count_query_points1[0]+1))&&
			(BBWT1->primary_idx<=pref_count_query_points1[npref_query_points1-1]));
			SLT_params.nleft_extensions1=includes_EOT_char1;
			SLT_params.left_extension_bitmap1=includes_EOT_char1;
// Set pref counts for the last character
			for(i=0;i<npref_query_points1;i++)
			{
				last_char_pref_count1=pref_count_query_points1[i]+1;
				for(j=0;j<4;j++)
					last_char_pref_count1-=char_pref_counts1[j+i*4];
				last_char_pref_counts1[i]=last_char_pref_count1;
			};
			last_char_freq1=0;

		}
		if(SLT_params.nright_extensions2)
		{
			if((pref_count_query_points2[0]+1)==0)
			{
				npref_query_points2--;
				DNA5_multipe_char_pref_counts64(BBWT2->indexed_BWT,npref_query_points2,
					&pref_count_query_points2[1],&char_pref_counts2[4]);
				npref_query_points2++;
			}
			else
				DNA5_multipe_char_pref_counts64(BBWT2->indexed_BWT,npref_query_points2,
					pref_count_query_points2,char_pref_counts2);
			includes_EOT_char2=((BBWT2->primary_idx>=(pref_count_query_points2[0]+1))&&
				(BBWT2->primary_idx<=pref_count_query_points2[npref_query_points2-1]));
			SLT_params.nleft_extensions2=includes_EOT_char2;
			SLT_params.left_extension_bitmap2=includes_EOT_char2;
// Set pref counts for the last character
			for(i=0;i<npref_query_points2;i++)
			{
				last_char_pref_count2=pref_count_query_points2[i]+1;
				for(j=0;j<4;j++)
					last_char_pref_count2-=char_pref_counts2[j+i*4];
				last_char_pref_counts2[i]=last_char_pref_count2;
			};
			last_char_freq2=0;

		}
// Compute the frequencies of all combinations of left and right extensions
		memset(SLT_params.left_right_extension_freqs1,0, 
			sizeof(SLT_params.left_right_extension_freqs1));
		memset(SLT_params.left_right_extension_freqs2,0, 
			sizeof(SLT_params.left_right_extension_freqs2));
		for(i=1,j1=1,j2=1;i<7;i++)
		{
			if(SLT_params.right_extension_bitmap1&(1<<(i-1)))
			{	
				includes_EOT_char1=((BBWT1->primary_idx>=(pref_count_query_points1[j1-1]+1))&&
						(BBWT1->primary_idx<=pref_count_query_points1[j1]));
				SLT_params.left_right_extension_freqs1[0][i-1]=includes_EOT_char1;
				SLT_params.left_right_extension_freqs1[1][i-1]=char_pref_counts1[j1*4]-
						char_pref_counts1[(j1-1)*4]-
						includes_EOT_char1;
				for(k=2;k<5;k++)
					SLT_params.left_right_extension_freqs1[k][i-1]=
						char_pref_counts1[j1*4+k-1]-
						char_pref_counts1[(j1-1)*4+k-1];
				SLT_params.left_right_extension_freqs1[5][i-1]=
					last_char_pref_counts1[j1]-
					last_char_pref_counts1[j1-1];
				last_char_freq1+=SLT_params.left_right_extension_freqs1[5][i-1];
				j1++;
			}
			if(SLT_params.right_extension_bitmap2&(1<<(i-1)))
			{	
				includes_EOT_char2=((BBWT2->primary_idx>=(pref_count_query_points2[j2-1]+1))&&
						(BBWT2->primary_idx<=pref_count_query_points2[j2]));
				SLT_params.left_right_extension_freqs2[0][i-1]=includes_EOT_char2;
				SLT_params.left_right_extension_freqs2[1][i-1]=char_pref_counts2[j2*4]-
						char_pref_counts2[(j2-1)*4]-
						includes_EOT_char2;
				for(k=2;k<5;k++)
					SLT_params.left_right_extension_freqs2[k][i-1]=
						char_pref_counts2[j2*4+k-1]-
						char_pref_counts2[(j2-1)*4+k-1];
				SLT_params.left_right_extension_freqs2[5][i-1]=
					last_char_pref_counts2[j2]-
					last_char_pref_counts2[j2-1];
				last_char_freq2+=SLT_params.left_right_extension_freqs2[5][i-1];
				j2++;
			}

		};
		extension_exists1=(last_char_freq1>0);
		extension_exists2=(last_char_freq2>0);
		SLT_params.nleft_extensions1+=extension_exists1;
		SLT_params.nleft_extensions2+=extension_exists2;
		SLT_params.left_extension_bitmap1|=(extension_exists1<<5);
		SLT_params.left_extension_bitmap2|=(extension_exists2<<5);
		string_depth=SLT_params.string_depth+1;
// Now generate the elements to be put in the stack and complete the 
// param structure to be passed to the callback function. 
		max_interval_idx=0;
		max_sum_interval_size=2;
		nexplicit_WL=0;
// First push the node labelled with character 1 if it exists
		nchildren=0;
		nchildren1=0;
		nchildren2=0;
		interval_size1=0;
		interval_size2=0;
		sum_interval_size=0;
		for(j=0;j<6;j++)
		{
			nchildren1+=(SLT_params.left_right_extension_freqs1[1][j]!=0);
			nchildren2+=(SLT_params.left_right_extension_freqs2[1][j]!=0);
			nchildren+=(SLT_params.left_right_extension_freqs1[1][j]+
				SLT_params.left_right_extension_freqs2[1][j]>0); 
// We speculatively write into the stack
			SLT_stack[curr_stack_idx].child_freqs1[j]=
				SLT_params.left_right_extension_freqs1[1][j];
			SLT_stack[curr_stack_idx].child_freqs2[j]=
				SLT_params.left_right_extension_freqs2[1][j];
			sum_interval_size+=SLT_params.left_right_extension_freqs1[1][j]+
					SLT_params.left_right_extension_freqs2[1][j];
			interval_size1+=SLT_params.left_right_extension_freqs1[1][j];		
			interval_size2+=SLT_params.left_right_extension_freqs2[1][j];	
		};
		extension_exists1=(nchildren1>0);
		extension_exists2=(nchildren2>0);
		SLT_params.nleft_extensions1+=extension_exists1;
		SLT_params.nleft_extensions2+=extension_exists2;
		SLT_params.left_extension_bitmap1|=(extension_exists1<<1);
		SLT_params.left_extension_bitmap2|=(extension_exists2<<1);
		if(((options & SLT_joint_or_enum)==0 && nchildren1>0 && nchildren2>0 &&( nchildren>1 || 
				SLT_params.left_right_extension_freqs1[1][0]+
				SLT_params.left_right_extension_freqs1[1][5]+
				SLT_params.left_right_extension_freqs2[1][0]+
				SLT_params.left_right_extension_freqs2[1][5]>=2))
		|| ((options & SLT_joint_or_enum) && (nchildren1>1 || nchildren2>1 || nchildren>1 ||
         			SLT_params.left_right_extension_freqs1[1][0]+
				SLT_params.left_right_extension_freqs1[1][5]+
				SLT_params.left_right_extension_freqs2[1][0]+
				SLT_params.left_right_extension_freqs2[1][5]>=2)))
		{
			// Push a new node in the stack.				
//			printf("We push a node with character 1 and string depth %d\n",
//				string_depth);
			if(curr_stack_size<=curr_stack_idx)
			{
				curr_stack_size*=2;
				SLT_stack=(SLT64_stack_item_t*)realloc(SLT_stack, 
					sizeof(SLT64_stack_item_t)*(curr_stack_size+1));
			};
			SLT_stack[curr_stack_idx].WL_char=1;
			SLT_stack[curr_stack_idx].string_depth=string_depth;
			includes_EOT_char1=(BBWT1->primary_idx<(pref_count_query_points1[0]+1));
			includes_EOT_char2=(BBWT2->primary_idx<(pref_count_query_points2[0]+1));
			SLT_stack[curr_stack_idx].interval_start1=
				char_pref_counts1[0]+BBWT1->char_base[0]+
				1-includes_EOT_char1;
			SLT_stack[curr_stack_idx].interval_start2=
				char_pref_counts2[0]+BBWT2->char_base[0]+
				1-includes_EOT_char2;
			SLT_stack[curr_stack_idx].interval_size1=interval_size1;
			SLT_stack[curr_stack_idx].interval_size2=interval_size2;
			sum_interval_size=interval_size1+interval_size2;
			max_sum_interval_size=sum_interval_size;
			nexplicit_WL++;
			curr_stack_idx++;
		}
// Then push nodes labelled with other characters
		for(i=1;i<4;i++)
		{
			nchildren=0;
			nchildren1=0;
			nchildren2=0;
			interval_size1=0;
			interval_size2=0;
			for(j=0;j<6;j++)
			{
				nchildren1+=(SLT_params.left_right_extension_freqs1[i+1][j]!=0);
				nchildren2+=(SLT_params.left_right_extension_freqs2[i+1][j]!=0);
				nchildren+=(SLT_params.left_right_extension_freqs1[i+1][j]+
					SLT_params.left_right_extension_freqs2[i+1][j]>0); 
// We speculatively write into the stack
				SLT_stack[curr_stack_idx].child_freqs1[j]=
					SLT_params.left_right_extension_freqs1[i+1][j];
				SLT_stack[curr_stack_idx].child_freqs2[j]=
					SLT_params.left_right_extension_freqs2[i+1][j];
				interval_size1+=SLT_params.left_right_extension_freqs1[i+1][j];		
				interval_size2+=SLT_params.left_right_extension_freqs2[i+1][j];		
			};
			extension_exists1=(nchildren1>0);
			extension_exists2=(nchildren2>0);
			SLT_params.nleft_extensions1+=extension_exists1;
			SLT_params.nleft_extensions2+=extension_exists2;
			SLT_params.left_extension_bitmap1|=(extension_exists1<<(i+1));
			SLT_params.left_extension_bitmap2|=(extension_exists2<<(i+1));
			if(((options & SLT_joint_or_enum)==0 && nchildren1>0 && nchildren2>0 &&( nchildren>1 || 
					SLT_params.left_right_extension_freqs1[i+1][0]+
					SLT_params.left_right_extension_freqs1[i+1][5]+
					SLT_params.left_right_extension_freqs2[i+1][0]+
					SLT_params.left_right_extension_freqs2[i+1][5]>=2))
			|| ((options & SLT_joint_or_enum) && (nchildren1>1 || nchildren2>1 || nchildren>1 ||
	         			SLT_params.left_right_extension_freqs1[i+1][0]+
					SLT_params.left_right_extension_freqs1[i+1][5]+
					SLT_params.left_right_extension_freqs2[i+1][0]+
					SLT_params.left_right_extension_freqs2[i+1][5]>=2)))
			{
				// Push a new node in the stack.				
//				printf("We push a node with character %d and string depth %d\n",i+1,
//					string_depth);
				if(curr_stack_size<=curr_stack_idx)
				{
					curr_stack_size*=2;
					SLT_stack=(SLT64_stack_item_t*)realloc(SLT_stack, 
						sizeof(SLT64_stack_item_t)*(curr_stack_size+1));
				};
				SLT_stack[curr_stack_idx].WL_char=i+1;
				SLT_stack[curr_stack_idx].string_depth=string_depth;
				SLT_stack[curr_stack_idx].interval_start1=
					char_pref_counts1[i]+BBWT1->char_base[i]+1;
				SLT_stack[curr_stack_idx].interval_start2=
					char_pref_counts2[i]+BBWT2->char_base[i]+1;

				SLT_stack[curr_stack_idx].interval_size1=interval_size1;
				SLT_stack[curr_stack_idx].interval_size2=interval_size2;
				sum_interval_size=interval_size1+interval_size2;
				if((options&SLT_stack_trick) && sum_interval_size>max_sum_interval_size)
				{
					max_sum_interval_size=sum_interval_size;
					max_interval_idx=nexplicit_WL;
				};
				nexplicit_WL++;
				curr_stack_idx++;
			};
		};

		if((options&SLT_stack_trick) && max_interval_idx)
			swap2_stack_items(&SLT_stack[curr_stack_idx-nexplicit_WL],
				&SLT_stack[curr_stack_idx-nexplicit_WL+max_interval_idx]);
		if((options&SLT_lex_order))
		{
			for(j=0;j<nexplicit_WL/2;j++)
				swap2_stack_items(&SLT_stack[curr_stack_idx-nexplicit_WL+j],
					&SLT_stack[curr_stack_idx-j-1]);
		};
		SLT_iterator->SLT_callback(&SLT_params,SLT_iterator->intern_state);
	}while(curr_stack_idx);
	printf("The number of traversed suffix tree nodes is %llu\n",ntraversed_nodes);
	free(SLT_stack);
};
