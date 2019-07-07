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
static void getRanksOfRightExtensions(const SLT_stack_item_t *stackFrame, const Basic_BWT_t *bwt, unsigned char *rightExtensionBitmap, unsigned int *pref_count_query_points, unsigned char *npref_query_points, unsigned int *char_pref_counts, unsigned int *last_char_pref_counts, unsigned char *containsSharp) {
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
 * Sets all fields of $SLT_params$ based on $stackFrame$.
 * See function $getRanksOfRightExtensions()$ for details on the input parameters.
 *
 * Remark: the procedure assumes that $SLT_params->left_right_extension_freqs$ contains
 * only zeros.
 *
 * @param nRightExtensions Cell $a \in [0..5]$ contains the number of (at most 6) distinct 
 * right-extensions of string $aW$. The array is assumed to be initialized to all zeros.
 *
 * @param intervalSize Cell $a \in [0..5]$ contains the size of the BWT interval of $aW$.
 * The array is assumed to be initialized to all zeros.
 */
static void buildCallbackState(SLT_params_t *SLT_params, const SLT_stack_item_t *stackFrame, const Basic_BWT_t *bwt, const unsigned char rightExtensionBitmap, const unsigned int *pref_count_query_points, const unsigned char npref_query_points, const unsigned int *char_pref_counts, const unsigned int *last_char_pref_counts, const unsigned char containsSharp, unsigned char *nRightExtensions, unsigned int *intervalSize) {
	unsigned char containsSharpTmp, extensionExists, leftExtensionBitmap;
	unsigned int i, j, k;
	
	SLT_params->string_depth=stackFrame->string_depth;
	SLT_params->bwt_start=stackFrame->interval_start;
	SLT_params->interval_size=stackFrame->interval_size;
	SLT_params->WL_char=stackFrame->WL_char;
	SLT_params->nright_extensions=npref_query_points-1;
	SLT_params->right_extension_bitmap=rightExtensionBitmap;
	for (i=0; i<=3; i++) SLT_params->left_ext_bwt_start[i]=bwt->char_base[i]+char_pref_counts[i]+1;
	if (bwt->primary_idx<(pref_count_query_points[0]+1)) {
		// We subtract one because character A, and not the actual sharp, is assigned
		// to position $primary_idx$ in the BWT.
		SLT_params->left_ext_bwt_start[0]--;
	}
	SLT_params->left_ext_bwt_start[4]=bwt->char_base[4]+last_char_pref_counts[0]+1;
	
	// Computing the frequencies of all combinations of left and right extensions
	j=0; leftExtensionBitmap=0; nRightExtensions[0]=1; intervalSize[0]=1;
	for (i=0; i<=5; i++) {  // For every right-extension
		if ((rightExtensionBitmap&(1<<i))==0) continue;
		j++;
		// Left-extension by #
		containsSharpTmp=((bwt->primary_idx>=(pref_count_query_points[j-1]+1))&&(bwt->primary_idx<=pref_count_query_points[j]));
		SLT_params->left_right_extension_freqs[0][i]=containsSharpTmp;
		leftExtensionBitmap|=containsSharpTmp;
		// Left-extension by A
		SLT_params->left_right_extension_freqs[1][i]=char_pref_counts[j<<2]-char_pref_counts[(j-1)<<2]-containsSharpTmp;  // We subtract $containsSharpTmp$ because character A, and not the actual sharp, is assigned to position $primary_idx$ in the BWT.
		extensionExists=!!SLT_params->left_right_extension_freqs[1][i];
		leftExtensionBitmap|=extensionExists<<1;
		nRightExtensions[1]+=extensionExists;
		intervalSize[1]+=SLT_params->left_right_extension_freqs[1][i];
		// Left-extension by C,G,T.
		for (k=1; k<=3; k++) {
			SLT_params->left_right_extension_freqs[k+1][i]=char_pref_counts[(j<<2)+k]-char_pref_counts[((j-1)<<2)+k];
			extensionExists=!!SLT_params->left_right_extension_freqs[k+1][i];
			leftExtensionBitmap|=extensionExists<<(k+1);
			nRightExtensions[k+1]+=extensionExists;
			intervalSize[k+1]+=SLT_params->left_right_extension_freqs[k+1][i];
		}
		// Left-extension by N
		SLT_params->left_right_extension_freqs[5][i]=last_char_pref_counts[j]-last_char_pref_counts[j-1];
		extensionExists=!!SLT_params->left_right_extension_freqs[5][i];
		leftExtensionBitmap|=extensionExists<<5;
		nRightExtensions[5]+=extensionExists;
		intervalSize[5]+=SLT_params->left_right_extension_freqs[5][i];
	}
	SLT_params->left_extension_bitmap=leftExtensionBitmap;
	SLT_params->nleft_extensions=0;
	for (i=0; i<=5; i++) SLT_params->nleft_extensions+=(leftExtensionBitmap&(1<<i))!=0;
}


/**
 * Tries to push $AW$ onto $stack$ (where A has character ID equal to one).
 *
 * //////////////////// TODO: Remark: $stackPointer$ is incremented before pushing.
 * 
 * @return 0 if $AW$ was not pushed on the stack; otherwise, the size of the BWT interval
 * of $AW$.
 */
static inline unsigned int pushA(const SLT_params_t *SLT_params, const Basic_BWT_t *bwt, SLT_stack_item_t **stack, unsigned int *stackSize, unsigned int *stackPointer, const unsigned int stringDepth, const unsigned int *pref_count_query_points, const unsigned int *char_pref_counts, const unsigned char *nRightExtensionsOfLeft, const unsigned int *intervalSizeOfLeft) {
	unsigned char containsSharp;
	unsigned int i;
	
	if (nRightExtensionsOfLeft[1]<2 && SLT_params->left_right_extension_freqs[1][5]<2) return 0;
	if (*stackPointer>=*stackSize) {
		*stackSize=(*stackSize)<<=1;
		*stack=(SLT_stack_item_t *)realloc(*stack,sizeof(SLT_stack_item_t)*(*stackSize));
	}
	(*stack)[*stackPointer].WL_char=1;
	(*stack)[*stackPointer].string_depth=stringDepth;
	containsSharp=bwt->primary_idx<(pref_count_query_points[0]+1);
	(*stack)[*stackPointer].interval_start=bwt->char_base[0]+char_pref_counts[0]+1-containsSharp;
	(*stack)[*stackPointer].interval_size=intervalSizeOfLeft[1];
	for (i=0; i<=5; i++) (*stack)[*stackPointer].child_freqs[i]=SLT_params->left_right_extension_freqs[1][i];
	*stackPointer=*stackPointer+1;
	return intervalSizeOfLeft[1];
}


/**
 * @param characterID >=2.
 * @return 0 if $bW$ was not pushed on the stack; otherwise, the size of the BWT interval
 * of $bW$, where $b=characterID$.
 */
static inline unsigned int pushNonA(int characterID, const SLT_params_t *SLT_params, const Basic_BWT_t *bwt, SLT_stack_item_t **stack, unsigned int *stackSize, unsigned int *stackPointer, const unsigned int stringDepth, const unsigned int *pref_count_query_points, const unsigned int *char_pref_counts, const unsigned char *nRightExtensionsOfLeft, const unsigned int *intervalSizeOfLeft) {
	unsigned int i;
	
	if (nRightExtensionsOfLeft[characterID]<2 && SLT_params->left_right_extension_freqs[characterID][5]<2) return 0;
	if (*stackPointer>=*stackSize) {
		*stackSize=(*stackSize)<<=1;
		*stack=(SLT_stack_item_t *)realloc(*stack,sizeof(SLT_stack_item_t)*(*stackSize));
	}
	(*stack)[*stackPointer].WL_char=characterID;
	(*stack)[*stackPointer].string_depth=stringDepth;
	(*stack)[*stackPointer].interval_start=bwt->char_base[characterID-1]+char_pref_counts[characterID-1]+1;
	(*stack)[*stackPointer].interval_size=intervalSizeOfLeft[characterID];
	for (i=0; i<=5; i++) (*stack)[*stackPointer].child_freqs[i]=SLT_params->left_right_extension_freqs[characterID][i];
	*stackPointer=*stackPointer+1;	
	return intervalSizeOfLeft[characterID];
}



//// TODO: maybe it's possible to simplify the code above and remove all these pointers by declaring global variables, very simply. in this way we would also remove all those input parameters to the functions.


// Above: add a destroyer that cancels new_SLT_iterator and the globals. 


/**
 * 
 */
void SLT_execute_iterator(SLT_iterator_t_single_string *SLT_iterator) {
	const Basic_BWT_t *BWT = SLT_iterator->BBWT;
	unsigned char i;
	unsigned char maxIntervalID, nExplicitWL, containsSharp, rightExtensionBitmap, npref_query_points;
	unsigned int stackSize, stackPointer, stringDepth, nTraversedNodes, intervalSize, maxIntervalSize;  // Should be 64-bit
	SLT_params_t SLT_params = {0};
	SLT_stack_item_t *stack;
	unsigned int pref_count_query_points[7];  // Should be 64-bit
	unsigned int char_pref_counts[28];  // Should be 64-bit
	unsigned int last_char_pref_counts[7];  // Should be 64-bit
	unsigned char nRightExtensionsOfLeft[6];
	unsigned int intervalSizeOfLeft[6];  // Should be 64-bit
	
	stackSize=MIN_SLT_STACK_SIZE;
	stack=(SLT_stack_item_t *)malloc((1+MIN_SLT_STACK_SIZE)*sizeof(SLT_stack_item_t));
	stack[0].WL_char=0;
	stack[0].string_depth=0;
	stack[0].interval_start=0;
	stack[0].child_freqs[0]=1;
	for (i=1; i<=4; i++) stack[0].child_freqs[i]=BWT->char_base[i]-BWT->char_base[i-1];
	stack[0].child_freqs[5]=BWT->textlen-BWT->char_base[4];
	stack[0].interval_size=BWT->textlen+1;
	nTraversedNodes=0; stackPointer=1;
	do {
		stackPointer--;
		nTraversedNodes++;
		
		// Computing ranks
		getRanksOfRightExtensions(&stack[stackPointer],BWT,&rightExtensionBitmap,pref_count_query_points,&npref_query_points,char_pref_counts,last_char_pref_counts,&containsSharp);
		
		// Issuing the callback function
		memset(SLT_params.left_right_extension_freqs,0,sizeof(SLT_params.left_right_extension_freqs));
		memset(nRightExtensionsOfLeft,0,sizeof(nRightExtensionsOfLeft));
		memset(intervalSizeOfLeft,0,sizeof(intervalSizeOfLeft));
		buildCallbackState(&SLT_params,&stack[stackPointer],BWT,rightExtensionBitmap,pref_count_query_points,npref_query_points,char_pref_counts,last_char_pref_counts,containsSharp,nRightExtensionsOfLeft,intervalSizeOfLeft);
		SLT_iterator->SLT_callback(SLT_params,SLT_iterator->intern_state);
		
		// Pushing $aW$, for $a \in {A,C,G,T}$ only, if it exists and it is right-maximal.
		stringDepth=SLT_params.string_depth+1;
		maxIntervalSize=pushA(&SLT_params,BWT,&stack,&stackSize,&stackPointer,stringDepth,pref_count_query_points,char_pref_counts,nRightExtensionsOfLeft,intervalSizeOfLeft);
		maxIntervalID=0;
		nExplicitWL=!!maxIntervalSize;
		for (i=2; i<=4; i++) {
			intervalSize=pushNonA(i,&SLT_params,BWT,&stack,&stackSize,&stackPointer,stringDepth,pref_count_query_points,char_pref_counts,nRightExtensionsOfLeft,intervalSizeOfLeft);
			if (!intervalSize) continue;
			if (intervalSize>maxIntervalSize) {
				maxIntervalSize=intervalSize;
				maxIntervalID=nExplicitWL;
			}
			nExplicitWL++;	
		}
		if (maxIntervalID && (SLT_iterator->options&SLT_stack_trick)) swapStackFrames(&stack[stackPointer-nExplicitWL],&stack[stackPointer-nExplicitWL+maxIntervalID]);
		if (SLT_iterator->options&SLT_lex_order) {
			for (i=0; i<nExplicitWL>>1; i++) swapStackFrames(&stack[stackPointer-nExplicitWL+i],&stack[stackPointer-1-i]);
		}
	} while(stackPointer);
	free(stack);
	
	printf("The number of traversed suffix tree nodes is %d \n",nTraversedNodes);
}