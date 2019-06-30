/**
 * 
 *
 * @author Djamal Belazzougui, Fabio Cunial
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SLT_single_string.h"

#define MIN_SLT_STACK_SIZE 4


SLT_iterator_t_single_string new_SLT_iterator(SLT_callback_t SLT_callback, void *intern_state, Basic_BWT_t *BBWT, unsigned short options) {
	SLT_iterator_t_single_string SLT_iterator;
	SLT_iterator.SLT_callback=SLT_callback;
	SLT_iterator.intern_state=intern_state;
	SLT_iterator.BBWT=BBWT;
	SLT_iterator.options=options;
	return SLT_iterator;
}


/**
 * A frame of the stack
 */
typedef struct {
	unsigned int string_depth;
	unsigned int interval_start;
	unsigned int child_freqs[6];
	unsigned int interval_size;
	unsigned char WL_char;
} SLT_stack_item_t;


static inline void swap2_stack_items(SLT_stack_item_t *SLT_stack_item1, SLT_stack_item_t *SLT_stack_item2) {
	SLT_stack_item1->string_depth^=SLT_stack_item2->string_depth;
	SLT_stack_item1->interval_start^=SLT_stack_item2->interval_start;
	SLT_stack_item1->interval_size^=SLT_stack_item2->interval_size;
	SLT_stack_item1->WL_char^=SLT_stack_item2->WL_char;

	SLT_stack_item2->string_depth^=SLT_stack_item1->string_depth;
	SLT_stack_item2->interval_start^=SLT_stack_item1->interval_start;
	SLT_stack_item2->interval_size^=SLT_stack_item1->interval_size;
	SLT_stack_item2->WL_char^=SLT_stack_item1->WL_char;

	SLT_stack_item1->string_depth^=SLT_stack_item2->string_depth;
	SLT_stack_item1->interval_start^=SLT_stack_item2->interval_start;
	SLT_stack_item1->interval_size^=SLT_stack_item2->interval_size;
	SLT_stack_item1->WL_char^=SLT_stack_item2->WL_char;

	SLT_stack_item1->child_freqs[0]^=SLT_stack_item2->child_freqs[0];
	SLT_stack_item1->child_freqs[1]^=SLT_stack_item2->child_freqs[1];
	SLT_stack_item1->child_freqs[2]^=SLT_stack_item2->child_freqs[2];
	SLT_stack_item1->child_freqs[3]^=SLT_stack_item2->child_freqs[3];
	SLT_stack_item1->child_freqs[4]^=SLT_stack_item2->child_freqs[4];
	SLT_stack_item1->child_freqs[5]^=SLT_stack_item2->child_freqs[5];

	SLT_stack_item2->child_freqs[0]^=SLT_stack_item1->child_freqs[0];
	SLT_stack_item2->child_freqs[1]^=SLT_stack_item1->child_freqs[1];
	SLT_stack_item2->child_freqs[2]^=SLT_stack_item1->child_freqs[2];
	SLT_stack_item2->child_freqs[3]^=SLT_stack_item1->child_freqs[3];
	SLT_stack_item2->child_freqs[4]^=SLT_stack_item1->child_freqs[4];
	SLT_stack_item2->child_freqs[5]^=SLT_stack_item1->child_freqs[5];

	SLT_stack_item1->child_freqs[0]^=SLT_stack_item2->child_freqs[0];
	SLT_stack_item1->child_freqs[1]^=SLT_stack_item2->child_freqs[1];
	SLT_stack_item1->child_freqs[2]^=SLT_stack_item2->child_freqs[2];
	SLT_stack_item1->child_freqs[3]^=SLT_stack_item2->child_freqs[3];
	SLT_stack_item1->child_freqs[4]^=SLT_stack_item2->child_freqs[4];
	SLT_stack_item1->child_freqs[5]^=SLT_stack_item2->child_freqs[5];
}





void SLT_execute_iterator(SLT_iterator_t_single_string *SLT_iterator) {
	const Basic_BWT_t *bwt = SLT_iterator->BBWT;
	unsigned char containsSharp, tmpBitmap, npref_query_points;
	unsigned int i, j, k;
	unsigned long count, stackSize, stackPointer, nTraversedNodes;
	SLT_stack_item_t *stack;
	SLT_params_t SLT_params;
	
	// Let $[i..j]$ be the BWT interval of a right-maximal string $W$. The array contains
	// the sorted list of positions $i-1,e_1,e_2,...,e_k$, where $e_p$ is the last
	// position of every sub-interval of $[i..j]$ induced by a right-extension of $W$.
	// Note that there can be at most 7 elements in the array, but there might be fewer.
	unsigned int pref_count_query_points[7];  // Should be long
	
	// 7 blocks of 4 elements each. The $j$-th element of block $i$ contains the number of
	// occurrences of character $j$ (0=A, 1=C, 2=G, 3=T) up to the $i$-th element of
	// $pref_count_query_points$ (included).
	unsigned int char_pref_counts[28];  // Should be long
	
	// Position $i$ contains the number of occurrences of character N, up to the $i$-th
	// element of $pref_count_query_points$ (included).
	unsigned int last_char_pref_counts[7];  // Should be long
	
	
	
	
	unsigned int extension_exists;
	unsigned int nchildren;
	unsigned int last_char_freq;
	unsigned int interval_size;
	unsigned int max_interval_size;
	unsigned int nexplicit_WL;
	unsigned int max_interval_idx;
	unsigned int options = SLT_iterator->options;
	unsigned int string_depth;
	
	
	stackSize=MIN_SLT_STACK_SIZE;  // In frames
	stack=(SLT_stack_item_t *)malloc((1+MIN_SLT_STACK_SIZE)*sizeof(SLT_stack_item_t));
	stack[0].WL_char=0;
	stack[0].string_depth=0;
	stack[0].interval_start=0;
	stack[0].child_freqs[0]=1;
	for (i=1; i<5; i++) stack[0].child_freqs[i]=bwt->char_base[i]-bwt->char_base[i-1];
	stack[0].child_freqs[5]=bwt->textlen-bwt->char_base[4];
	stack[0].interval_size=bwt->textlen+1;
	stackPointer=1L; nTraversedNodes=0L;
	do {
		stackPointer--;
		nTraversedNodes++;

		// Computing all query points, and all their ranks in a single call.
		tmpBitmap=0;
		j=0;
		pref_count_query_points[j]=stack[stackPointer].interval_start-1;
		for (i=0; i<=5; i++) {
			count=stack[stackPointer].child_freqs[i];
			if (count>0) {
				tmpBitmap|=1<<i;
				j++;
				pref_count_query_points[j]=pref_count_query_points[j-1]+count;
			}
		}
		containsSharp=(bwt->primary_idx>=pref_count_query_points[0]+1)&&(bwt->primary_idx<=pref_count_query_points[j]);
		npref_query_points=j+1;
		if (pref_count_query_points[0]+1==0) {
			for (i=0; i<4; i++) char_pref_counts[i]=0;
			DNA5_multipe_char_pref_counts(bwt->indexed_BWT,npref_query_points-1,&pref_count_query_points[1],&char_pref_counts[4]);
		}
		else DNA5_multipe_char_pref_counts(bwt->indexed_BWT,npref_query_points,pref_count_query_points,char_pref_counts);
		for (i=0; i<npref_query_points; i++) {
			count=pref_count_query_points[i]+1;
			for (j=0; j<4; j++) count-=char_pref_counts[(i<<2)+j];
			last_char_pref_counts[i]=count;
		}

		// Set the data related to the top node to be given as parameter to the call back
		// function.
		SLT_params.WL_char=stack[stackPointer].WL_char;
		SLT_params.string_depth=stack[stackPointer].string_depth;
		SLT_params.bwt_start=stack[stackPointer].interval_start;
		SLT_params.interval_size=stack[stackPointer].interval_size;
		SLT_params.right_extension_bitmap=tmpBitmap;
		SLT_params.nright_extensions=npref_query_points-1;
		SLT_params.nleft_extensions=containsSharp;
		SLT_params.left_extension_bitmap=containsSharp;
		
		// Set the starting point of left children in bwt in SLT_params
		for (i=0; i<=3; i++) SLT_params.left_ext_bwt_start[i]=bwt->char_base[i]+char_pref_counts[i]+1;
		if (bwt->primary_idx<=pref_count_query_points[0]) {
			// We subtract because to position $primary_idx$ it was assigned character A
			// in the BWT, rather than the extra real character #.
			SLT_params.left_ext_bwt_start[0]--;
		}
		SLT_params.left_ext_bwt_start[4]=bwt->char_base[4]+last_char_pref_counts[0]+1;
		
		
		
		
// Compute the frequencies of all combinations of left and right extensions
		last_char_freq=0;
		memset(SLT_params.left_right_extension_freqs,0, 
			sizeof(SLT_params.left_right_extension_freqs));
		for(i=1,j=1;i<7;i++)
		{
			if((SLT_params.right_extension_bitmap&(1<<(i-1)))==0)
				continue;
			containsSharp=((bwt->primary_idx>=(pref_count_query_points[j-1]+1))&&
					(bwt->primary_idx<=pref_count_query_points[j]));
			SLT_params.left_right_extension_freqs[0][i-1]=containsSharp;
			SLT_params.left_right_extension_freqs[1][i-1]=char_pref_counts[j*4]-
					char_pref_counts[(j-1)*4]-
					containsSharp;
//			printf("containsSharp=%d\n",containsSharp);
//			printf("freq of char 1 is =%d\n",
//				SLT_params.left_right_extension_freqs[1][i-1]);
//			printf("pref count of char 1 up to %d is %d and to %d is %d\n",j,
//				char_pref_counts[j*4],j-1,char_pref_counts[(j-1)*4]);
			for(k=2;k<5;k++)
				SLT_params.left_right_extension_freqs[k][i-1]=
					char_pref_counts[j*4+k-1]-
					char_pref_counts[(j-1)*4+k-1];
			SLT_params.left_right_extension_freqs[5][i-1]=
				last_char_pref_counts[j]-
				last_char_pref_counts[j-1];
			last_char_freq+=SLT_params.left_right_extension_freqs[5][i-1];
			j++;
		};
		extension_exists=(last_char_freq>0);
		SLT_params.nleft_extensions+=extension_exists;
		SLT_params.left_extension_bitmap|=(extension_exists<<5);
//		containsSharp=(bwt->primary_idx<=char_pref_counts[0]);
//		char_pref_counts[0]-=containsSharp;
		string_depth=SLT_params.string_depth+1;
// Now generate the elements to be put in the stack and complete the 
// param structure to be passed to the callback function. 
		max_interval_idx=0;
		max_interval_size=2;
		nexplicit_WL=0;
// First push the node labelled with character 1 if it exists
		nchildren=0;
		interval_size=0;
		for(j=0;j<6;j++)
		{
			nchildren+=(SLT_params.left_right_extension_freqs[1][j]!=0);
// We speculatively write into the stack
			stack[stackPointer].child_freqs[j]=
				SLT_params.left_right_extension_freqs[1][j];
			interval_size+=SLT_params.left_right_extension_freqs[1][j];
		};
		extension_exists=(nchildren>0);
		SLT_params.nleft_extensions+=extension_exists;
		SLT_params.left_extension_bitmap|=(extension_exists<<1);
		if(nchildren>1 || SLT_params.left_right_extension_freqs[1][5]>=2)
		{
			// Push a new node in the stack.				
//			printf("We push a node with character 1 and string depth %d\n",
//				string_depth);
			if(stackSize<=stackPointer)
			{
				stackSize*=2;
				stack=(SLT_stack_item_t*)realloc(stack, 
					sizeof(SLT_stack_item_t)*(stackSize+1));
			};
			stack[stackPointer].WL_char=1;
			stack[stackPointer].string_depth=string_depth;
			containsSharp=(bwt->primary_idx<(pref_count_query_points[0]+1));
			stack[stackPointer].interval_start=
				char_pref_counts[0]+bwt->char_base[0]+
				1-containsSharp;
			stack[stackPointer].interval_size=interval_size;
			max_interval_size=interval_size;
			nexplicit_WL++;
			stackPointer++;
		}
// Then push nodes labelled with other characters
		for(i=1;i<4;i++)
		{
			nchildren=0;
			interval_size=0;
			for(j=0;j<6;j++)
			{
				nchildren+=(SLT_params.left_right_extension_freqs[i+1][j]!=0);
// We speculatively write into the stack
				stack[stackPointer].child_freqs[j]=
					SLT_params.left_right_extension_freqs[i+1][j];
				interval_size+=SLT_params.left_right_extension_freqs[i+1][j];		
			};
			extension_exists=(nchildren>0);
			SLT_params.nleft_extensions+=extension_exists;
			SLT_params.left_extension_bitmap|=(extension_exists<<(i+1));
			if(nchildren>1 || SLT_params.left_right_extension_freqs[i+1][5]>=2)
			{
				// Push a new node in the stack.				
//				printf("We push a node with character %d and string depth %d\n",i+1,
//					string_depth);
				if(stackSize<=stackPointer)
				{
					stackSize*=2;
					stack=(SLT_stack_item_t*)realloc(stack, 
						sizeof(SLT_stack_item_t)*(stackSize+1));
				};
				stack[stackPointer].WL_char=i+1;
				stack[stackPointer].string_depth=string_depth;
				stack[stackPointer].interval_start=
					char_pref_counts[i]+bwt->char_base[i]+1;
				stack[stackPointer].interval_size=interval_size;
				if(options==SLT_stack_trick && interval_size>max_interval_size)
				{
					max_interval_size=interval_size;
					max_interval_idx=nexplicit_WL;
				};
				nexplicit_WL++;
				stackPointer++;
			};
		};

		if(options==SLT_stack_trick && max_interval_idx)
			swap2_stack_items(&stack[stackPointer-nexplicit_WL],
				&stack[stackPointer-nexplicit_WL+max_interval_idx]);
		if(options==SLT_lex_order)
		{
			for(j=0;j<nexplicit_WL/2;j++)
				swap2_stack_items(&stack[stackPointer-nexplicit_WL+j],
					&stack[stackPointer-j-1]);
		};
		SLT_iterator->SLT_callback(SLT_params,SLT_iterator->intern_state);
	}while(stackPointer);
	printf("The number of traversed suffix tree nodes is %ld \n",nTraversedNodes);
	free(stack);
};
