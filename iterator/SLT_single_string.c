/**
 * 
 *
 * @author Djamal Belazzougui, Fabio Cunial
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SLT_single_string.h"


/**
 * Initial size of the iterator stack (in stack frames).
 */
static const unsigned char MIN_SLT_STACK_SIZE = 4;



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


static inline void swapStackFrames(SLT_stack_item_t *SLT_stack_item1, SLT_stack_item_t *SLT_stack_item2) {
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


/**
 * Computes all distinct right-extensions $Wa$ of the string $W$ encoded in $stackFrame$,
 * as well as all their ranks. The results are written in the data structures given in 
 * input.
 *  
 * @param rightExtensionBitmap Output value. The $i$-th LSB is set to one iff character 
 * $i$ (0=#, 1=A, 2=C, 3=G, 4=T, 5=N) is a right-extension of $W$.
 *
 * @param pref_count_query_points Output array of length at least 7. Let $[i..j]$ be the 
 * BWT interval of $W$. The array contains the sorted list of positions $i-1,e_1,e_2,...,
 * e_k$, where $e_p$ is the last position of every sub-interval of $[i..j]$ induced by a 
 * right-extension of $W$, and $k<=7$ is returned inside $npref_query_points$.
 *
 * @param char_pref_counts Output array of length at least 28. It consists of at most 7 
 * blocks of 4 elements each. The $j$-th element of block $i$ contains the number of 
 * occurrences of character $j$ (0=A, 1=C, 2=G, 3=T) up to the $i$-th element of 
 * $pref_count_query_points$ (included).
 *
 * @param last_char_pref_counts Output array of length at least 7. Position $i$ contains 
 * the number of occurrences of character N, up to the $i$-th element of 
 * $pref_count_query_points$ (included).
 *
 * @param containsSharp Output value. True iff the BWT interval of $W$ contains the sharp.
 */
static inline void getRanksOfRightExtensions(const SLT_stack_item_t *stackFrame, const Basic_BWT_t *bwt, unsigned char *rightExtensionBitmap, unsigned int *pref_count_query_points, unsigned char *npref_query_points, unsigned int *char_pref_counts, unsigned int *last_char_pref_counts, unsigned char *containsSharp) {
	unsigned int i, j;
	unsigned long count;
	
	*rightExtensionBitmap=0;
	j=0;
	pref_count_query_points[j]=stackFrame->interval_start-1;
	for (i=0; i<=5; i++) {
		count=stackFrame->child_freqs[i];
		if (count>0) {
			*rightExtensionBitmap|=1<<i;
			j++;
			pref_count_query_points[j]=pref_count_query_points[j-1]+count;
		}
	}
	*containsSharp=(bwt->primary_idx>=pref_count_query_points[0]+1)&&(bwt->primary_idx<=pref_count_query_points[j]);
	*npref_query_points=j+1;
	if (pref_count_query_points[0]+1==0) {
		for (i=0; i<4; i++) char_pref_counts[i]=0;
		DNA5_multipe_char_pref_counts(bwt->indexed_BWT,*npref_query_points-1,&pref_count_query_points[1],&char_pref_counts[4]);
	}
	else DNA5_multipe_char_pref_counts(bwt->indexed_BWT,*npref_query_points,pref_count_query_points,char_pref_counts);
	for (i=0; i<*npref_query_points; i++) {
		count=pref_count_query_points[i]+1;
		for (j=0; j<4; j++) count-=char_pref_counts[(i<<2)+j];
		last_char_pref_counts[i]=count;
	}
}


/**
 * 
 */
static inline void buildCallbackState(SLT_params_t *SLT_params, const SLT_stack_item_t *stackFrame, const Basic_BWT_t *bwt, const unsigned char rightExtensionBitmap, const unsigned int *pref_count_query_points, const unsigned char npref_query_points, const unsigned int *char_pref_counts, const unsigned int *last_char_pref_counts, const unsigned char containsSharp) {
	unsigned char containsSharpTmp, leftExtensionBitmap;
	unsigned int i, j, k;
	
	SLT_params->string_depth=stackFrame->string_depth;
	SLT_params->bwt_start=stackFrame->interval_start;
	SLT_params->interval_size=stackFrame->interval_size;
	SLT_params->WL_char=stackFrame->WL_char;
	SLT_params->nright_extensions=npref_query_points-1;
	SLT_params->right_extension_bitmap=rightExtensionBitmap;
	for (i=0; i<=3; i++) SLT_params->left_ext_bwt_start[i]=bwt->char_base[i]+char_pref_counts[i]+1;
	if (bwt->primary_idx<=pref_count_query_points[0]) {
		// We subtract one because character A, and not the actual sharp, is assigned
		// to position $primary_idx$ in the BWT.
		SLT_params->left_ext_bwt_start[0]--;
	}
	SLT_params->left_ext_bwt_start[4]=bwt->char_base[4]+last_char_pref_counts[0]+1;
	
	// Computing the frequencies of all combinations of left and right extensions
	memset(SLT_params->left_right_extension_freqs,0,sizeof(SLT_params->left_right_extension_freqs));
	j=0; leftExtensionBitmap=0;
	for (i=0; i<=5; i++) {
		if ((rightExtensionBitmap&(1<<i))==0) continue;
		j++;
		containsSharpTmp=((bwt->primary_idx>=(pref_count_query_points[j-1]+1))&&(bwt->primary_idx<=pref_count_query_points[j]));
		SLT_params->left_right_extension_freqs[0][i]=containsSharpTmp;
		leftExtensionBitmap|=containsSharpTmp;
		SLT_params->left_right_extension_freqs[1][i]=char_pref_counts[j<<2]-char_pref_counts[(j-1)<<2]-containsSharpTmp;  // We subtract $containsSharpTmp$ because character A, and not the actual sharp, is assigned to position $primary_idx$ in the BWT.
		leftExtensionBitmap|=(!!SLT_params->left_right_extension_freqs[1][i])<<1;
		for (k=1; k<=3; k++) {
			SLT_params->left_right_extension_freqs[k+1][i]=char_pref_counts[(j<<2)+k]-char_pref_counts[((j-1)<<2)+k];
			leftExtensionBitmap|=(!!SLT_params->left_right_extension_freqs[k+1][i])<<(k+1);
		}
		SLT_params->left_right_extension_freqs[5][i]=last_char_pref_counts[j]-last_char_pref_counts[j-1];
		leftExtensionBitmap|=(!!SLT_params->left_right_extension_freqs[5][i])<<5;
	}
	SLT_params->left_extension_bitmap=leftExtensionBitmap;
	SLT_params->nleft_extensions=0;
	for (i=0; i<=5; i++) SLT_params->nleft_extensions+=(leftExtensionBitmap&(1<<i))!=0;
}




void SLT_execute_iterator(SLT_iterator_t_single_string *SLT_iterator) {
	const Basic_BWT_t *bwt = SLT_iterator->BBWT;
	unsigned char containsSharp, rightExtensionBitmap, npref_query_points;
	unsigned int i, j;
	unsigned long stackSize, stackPointer, nTraversedNodes;
	SLT_stack_item_t *stack;
	SLT_params_t SLT_params = {0};
	
	unsigned int pref_count_query_points[7];  // Should be long
	unsigned int char_pref_counts[28];  // Should be long
	unsigned int last_char_pref_counts[7];  // Should be long
	
	
	
	
	unsigned int extension_exists;
	unsigned int nchildren;
	unsigned int interval_size;
	unsigned int max_interval_size;
	unsigned int nexplicit_WL;
	unsigned int max_interval_idx;
	unsigned int options = SLT_iterator->options;
	unsigned int string_depth;
	
	
	stackSize=MIN_SLT_STACK_SIZE;
	stack=(SLT_stack_item_t *)malloc((1+MIN_SLT_STACK_SIZE)*sizeof(SLT_stack_item_t));
	stack[0].WL_char=0;
	stack[0].string_depth=0;
	stack[0].interval_start=0;
	stack[0].child_freqs[0]=1;
	for (i=1; i<=4; i++) stack[0].child_freqs[i]=bwt->char_base[i]-bwt->char_base[i-1];
	stack[0].child_freqs[5]=bwt->textlen-bwt->char_base[4];
	stack[0].interval_size=bwt->textlen+1;
	stackPointer=1L; nTraversedNodes=0L;
	do {
		stackPointer--;
		nTraversedNodes++;

		getRanksOfRightExtensions(&stack[stackPointer],bwt,&rightExtensionBitmap,pref_count_query_points,&npref_query_points,char_pref_counts,last_char_pref_counts,&containsSharp);
		buildCallbackState(&SLT_params,&stack[stackPointer],bwt,rightExtensionBitmap,pref_count_query_points,npref_query_points,char_pref_counts,last_char_pref_counts,containsSharp);


		

		string_depth=SLT_params.string_depth+1;
// Now generate the elements to be put in the stack and complete the 
// param structure to be passed to the callback function. 
		max_interval_idx=0;
		max_interval_size=2;
		nexplicit_WL=0;
// First push the node labelled with character 1 if it exists
		nchildren=0;
		interval_size=0;
		for (j=0; j<=5; j++) {
			nchildren+=(SLT_params.left_right_extension_freqs[1][j]!=0);
// We speculatively write into the stack
			stack[stackPointer].child_freqs[j]=SLT_params.left_right_extension_freqs[1][j];
			interval_size+=SLT_params.left_right_extension_freqs[1][j];
		}
		extension_exists=(nchildren>0);
		if (nchildren>1 || SLT_params.left_right_extension_freqs[1][5]>=2) {
			// Push a new node in the stack.
			if (stackPointer>=stackSize) {
				stackSize<<=1;
				stack=(SLT_stack_item_t*)realloc(stack,sizeof(SLT_stack_item_t)*stackSize);
			}
			stack[stackPointer].WL_char=1;
			stack[stackPointer].string_depth=string_depth;
			containsSharp=(bwt->primary_idx<(pref_count_query_points[0]+1));
			stack[stackPointer].interval_start=bwt->char_base[0]+char_pref_counts[0]+1-containsSharp;
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
			swapStackFrames(&stack[stackPointer-nexplicit_WL],
				&stack[stackPointer-nexplicit_WL+max_interval_idx]);
		if(options==SLT_lex_order)
		{
			for(j=0;j<nexplicit_WL/2;j++)
				swapStackFrames(&stack[stackPointer-nexplicit_WL+j],
					&stack[stackPointer-j-1]);
		};
		SLT_iterator->SLT_callback(SLT_params,SLT_iterator->intern_state);
	}while(stackPointer);
	printf("The number of traversed suffix tree nodes is %ld \n",nTraversedNodes);
	free(stack);
};
