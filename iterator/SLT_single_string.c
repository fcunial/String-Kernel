/**
 * @author Djamal Belazzougui, Fabio Cunial
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SLT_single_string.h"


/**
 * Initial size of the iterator stack (in stack frames).
 */
static const unsigned char MIN_SLT_STACK_SIZE = 16;


UnaryIterator_t newIterator(SLT_callback_t SLT_callback, void *applicationData, Basic_BWT_t *BBWT, unsigned int maxLength) {
	UnaryIterator_t iterator;
	iterator.SLT_callback=SLT_callback;
	iterator.applicationData=applicationData;
	iterator.BBWT=BBWT;
	iterator.maxLength=maxLength;
	return iterator;
}


/**
 * A stack frame
 */
typedef struct {
	unsigned int length;
	unsigned int bwtStart;
	unsigned int frequency;
	unsigned char firstCharacter;
	unsigned int frequency_right[6];  // 0=#, 1=A, 2=C, 3=G, 4=T, 5=N.
} StackFrame_t;


static inline void swapStackFrames(StackFrame_t *SLT_stack_item1, StackFrame_t *SLT_stack_item2) {
	SLT_stack_item1->length^=SLT_stack_item2->length;
	SLT_stack_item1->bwtStart^=SLT_stack_item2->bwtStart;
	SLT_stack_item1->frequency^=SLT_stack_item2->frequency;
	SLT_stack_item1->firstCharacter^=SLT_stack_item2->firstCharacter;

	SLT_stack_item2->length^=SLT_stack_item1->length;
	SLT_stack_item2->bwtStart^=SLT_stack_item1->bwtStart;
	SLT_stack_item2->frequency^=SLT_stack_item1->frequency;
	SLT_stack_item2->firstCharacter^=SLT_stack_item1->firstCharacter;

	SLT_stack_item1->length^=SLT_stack_item2->length;
	SLT_stack_item1->bwtStart^=SLT_stack_item2->bwtStart;
	SLT_stack_item1->frequency^=SLT_stack_item2->frequency;
	SLT_stack_item1->firstCharacter^=SLT_stack_item2->firstCharacter;

	SLT_stack_item1->frequency_right[0]^=SLT_stack_item2->frequency_right[0];
	SLT_stack_item1->frequency_right[1]^=SLT_stack_item2->frequency_right[1];
	SLT_stack_item1->frequency_right[2]^=SLT_stack_item2->frequency_right[2];
	SLT_stack_item1->frequency_right[3]^=SLT_stack_item2->frequency_right[3];
	SLT_stack_item1->frequency_right[4]^=SLT_stack_item2->frequency_right[4];
	SLT_stack_item1->frequency_right[5]^=SLT_stack_item2->frequency_right[5];

	SLT_stack_item2->frequency_right[0]^=SLT_stack_item1->frequency_right[0];
	SLT_stack_item2->frequency_right[1]^=SLT_stack_item1->frequency_right[1];
	SLT_stack_item2->frequency_right[2]^=SLT_stack_item1->frequency_right[2];
	SLT_stack_item2->frequency_right[3]^=SLT_stack_item1->frequency_right[3];
	SLT_stack_item2->frequency_right[4]^=SLT_stack_item1->frequency_right[4];
	SLT_stack_item2->frequency_right[5]^=SLT_stack_item1->frequency_right[5];

	SLT_stack_item1->frequency_right[0]^=SLT_stack_item2->frequency_right[0];
	SLT_stack_item1->frequency_right[1]^=SLT_stack_item2->frequency_right[1];
	SLT_stack_item1->frequency_right[2]^=SLT_stack_item2->frequency_right[2];
	SLT_stack_item1->frequency_right[3]^=SLT_stack_item2->frequency_right[3];
	SLT_stack_item1->frequency_right[4]^=SLT_stack_item2->frequency_right[4];
	SLT_stack_item1->frequency_right[5]^=SLT_stack_item2->frequency_right[5];
}


/**
 * Computes all distinct right-extensions $Wa$ of the string $W$ encoded in $stackFrame$,
 * as well as all their ranks. The results are written in the data structures given in 
 * input.
 *  
 * @param rightExtensionBitmap Output value. The $i$-th LSB is set to one iff character 
 * $i$ (0=#, 1=A, 2=C, 3=G, 4=T, 5=N) is a right-extension of $W$.
 *
 * @param rankPoints Output array of length at least 7. Let $[i..j]$ be the 
 * BWT interval of $W$. The array contains the sorted list of positions $i-1,e_1,e_2,...,
 * e_k$, where $e_p$ is the last position of every sub-interval of $[i..j]$ induced by a 
 * right-extension of $W$, and $k<=7$ is returned inside $npref_query_points$.
 *
 * @param rankValues Output array of length at least 28. It consists of at most 7 
 * blocks of 4 elements each. The $j$-th element of block $i$ contains the number of 
 * occurrences of character $j$ (0=A, 1=C, 2=G, 3=T) up to the $i$-th element of 
 * $rankPoints$ (included).
 *
 * @param rankValuesN Output array of length at least 7. Position $i$ contains 
 * the number of occurrences of character N, up to the $i$-th element of 
 * $rankPoints$ (included).
 *
 * @param containsSharp Output value. True iff the BWT interval of $W$ contains the sharp.
 */
static void getRanksOfRightExtensions(const StackFrame_t *stackFrame, const Basic_BWT_t *bwt, unsigned char *rightExtensionBitmap, unsigned int *rankPoints, unsigned char *npref_query_points, unsigned int *rankValues, unsigned int *rankValuesN, unsigned char *containsSharp) {
	unsigned int i, j;
	unsigned long count;
	
	*rightExtensionBitmap=0;
	j=0;
	rankPoints[j]=stackFrame->bwtStart-1;
	for (i=0; i<=5; i++) {
		count=stackFrame->frequency_right[i];
		if (count>0) {
			*rightExtensionBitmap|=1<<i;
			j++;
			rankPoints[j]=rankPoints[j-1]+count;
		}
	}
	*containsSharp=(bwt->primary_idx>=rankPoints[0]+1)&&(bwt->primary_idx<=rankPoints[j]);
	*npref_query_points=j+1;
	if (rankPoints[0]+1==0) {
		for (i=0; i<4; i++) rankValues[i]=0;
		DNA5_multipe_char_pref_counts(bwt->indexed_BWT,*npref_query_points-1,&rankPoints[1],&rankValues[4]);
	}
	else DNA5_multipe_char_pref_counts(bwt->indexed_BWT,*npref_query_points,rankPoints,rankValues);
	for (i=0; i<*npref_query_points; i++) {
		count=rankPoints[i]+1;
		for (j=0; j<4; j++) count-=rankValues[(i<<2)+j];
		rankValuesN[i]=count;
	}
}


/**
 * Sets all fields of $rightMaximalString$ based on $stackFrame$.
 * See function $getRanksOfRightExtensions()$ for details on the input parameters.
 *
 * Remark: the procedure assumes that $rightMaximalString->frequency_leftRight$ contains
 * only zeros.
 *
 * @param nRightExtensions cell $a \in [0..5]$ contains the number of (at most 6) distinct 
 * right-extensions of string $aW$; the array is assumed to be initialized to all zeros.
 *
 * @param intervalSize cell $a \in [0..5]$ contains the size of the BWT interval of $aW$;
 * the array is assumed to be initialized to all zeros.
 */
static void buildCallbackState(RightMaximalString_t *rightMaximalString, const StackFrame_t *stackFrame, const Basic_BWT_t *bwt, const unsigned char rightExtensionBitmap, const unsigned int *rankPoints, const unsigned char npref_query_points, const unsigned int *rankValues, const unsigned int *rankValuesN, const unsigned char containsSharp, unsigned char *nRightExtensions, unsigned int *intervalSize) {
	unsigned char containsSharpTmp, extensionExists, leftExtensionBitmap;
	unsigned int i, j, k;
	
	rightMaximalString->length=stackFrame->length;
	rightMaximalString->bwtStart=stackFrame->bwtStart;
	rightMaximalString->frequency=stackFrame->frequency;
	rightMaximalString->firstCharacter=stackFrame->firstCharacter;
	rightMaximalString->nRightExtensions=npref_query_points-1;
	rightMaximalString->rightExtensionBitmap=rightExtensionBitmap;
	for (i=0; i<=3; i++) rightMaximalString->bwtStart_left[i]=bwt->char_base[i]+rankValues[i]+1;
	if (bwt->primary_idx<(rankPoints[0]+1)) {
		// We subtract one because character A, and not the actual sharp, is assigned
		// to position $primary_idx$ in the BWT.
		rightMaximalString->bwtStart_left[0]--;
	}
	rightMaximalString->bwtStart_left[4]=bwt->char_base[4]+rankValuesN[0]+1;
	
	// Computing the frequencies of all combinations of left and right extensions
	j=0; leftExtensionBitmap=0; nRightExtensions[0]=1; intervalSize[0]=1;
	for (i=0; i<=5; i++) {  // For every right-extension
		if ((rightExtensionBitmap&(1<<i))==0) continue;
		j++;
		// Left-extension by #
		containsSharpTmp=((bwt->primary_idx>=(rankPoints[j-1]+1))&&(bwt->primary_idx<=rankPoints[j]));
		rightMaximalString->frequency_leftRight[0][i]=containsSharpTmp;
		leftExtensionBitmap|=containsSharpTmp;
		// Left-extension by A
		rightMaximalString->frequency_leftRight[1][i]=rankValues[j<<2]-rankValues[(j-1)<<2]-containsSharpTmp;  // We subtract $containsSharpTmp$ because character A, and not the actual sharp, is assigned to position $primary_idx$ in the BWT.
		extensionExists=!!rightMaximalString->frequency_leftRight[1][i];
		leftExtensionBitmap|=extensionExists<<1;
		nRightExtensions[1]+=extensionExists;
		intervalSize[1]+=rightMaximalString->frequency_leftRight[1][i];
		// Left-extension by C,G,T.
		for (k=1; k<=3; k++) {
			rightMaximalString->frequency_leftRight[k+1][i]=rankValues[(j<<2)+k]-rankValues[((j-1)<<2)+k];
			extensionExists=!!rightMaximalString->frequency_leftRight[k+1][i];
			leftExtensionBitmap|=extensionExists<<(k+1);
			nRightExtensions[k+1]+=extensionExists;
			intervalSize[k+1]+=rightMaximalString->frequency_leftRight[k+1][i];
		}
		// Left-extension by N
		rightMaximalString->frequency_leftRight[5][i]=rankValuesN[j]-rankValuesN[j-1];
		extensionExists=!!rightMaximalString->frequency_leftRight[5][i];
		leftExtensionBitmap|=extensionExists<<5;
		nRightExtensions[5]+=extensionExists;
		intervalSize[5]+=rightMaximalString->frequency_leftRight[5][i];
	}
	rightMaximalString->leftExtensionBitmap=leftExtensionBitmap;
	rightMaximalString->nLeftExtensions=0;
	for (i=0; i<=5; i++) rightMaximalString->nLeftExtensions+=(leftExtensionBitmap&(1<<i))!=0;
}


/**
 * @return 1 if the left-extension of $RightMaximalString_t$ by character $b$ is 
 * right-maximal by the current definition, zero otherwise.
 */
static inline unsigned char isLeftExtensionRightMaximal(unsigned char b, const RightMaximalString_t *rightMaximalString, const unsigned char *nRightExtensionsOfLeft) {
	unsigned char nRightExtensions;
	unsigned int i;
	
	if (TRAVERSAL_MAXIMALITY==0) {
		if (nRightExtensionsOfLeft[b]<2) return 0;
	}
	else if (TRAVERSAL_MAXIMALITY==1) {
		if (nRightExtensionsOfLeft[b]<2 && rightMaximalString->frequency_leftRight[b][5]<2) return 0;
	}
	else if (TRAVERSAL_MAXIMALITY==2) {
		nRightExtensions=0;
		for (i=1; i<=4; i++) nRightExtensions+=!!rightMaximalString->frequency_leftRight[b][i];
		if (nRightExtensions<2) return 0;
	}
	return 1;
}


/**
 * Tries to push $aW$ onto $stack$, where a=A has character ID equal to one.
 * 
 * @param stackPointer pointer to the first free frame in the stack; the procedure 
 * increments $stackPointer$ at the end;
 * @return 0 if $AW$ was not pushed on the stack; otherwise, the size of the BWT interval
 * of $AW$.
 */
static inline unsigned int pushA(const RightMaximalString_t *rightMaximalString, const Basic_BWT_t *bwt, StackFrame_t **stack, unsigned int *stackSize, unsigned int *stackPointer, const unsigned int length, const unsigned int *rankPoints, const unsigned int *rankValues, const unsigned char *nRightExtensionsOfLeft, const unsigned int *intervalSizeOfLeft) {
	unsigned char containsSharp;
	unsigned int i;
	
	if (!isLeftExtensionRightMaximal(1,rightMaximalString,nRightExtensionsOfLeft)) return 0;
	if (*stackPointer>=*stackSize) {
		*stackSize=(*stackSize)<<=1;
		*stack=(StackFrame_t *)realloc(*stack,sizeof(StackFrame_t)*(*stackSize));
	}
	(*stack)[*stackPointer].firstCharacter=1;
	(*stack)[*stackPointer].length=length;
	containsSharp=bwt->primary_idx<(rankPoints[0]+1);
	(*stack)[*stackPointer].bwtStart=bwt->char_base[0]+rankValues[0]+1-containsSharp;
	(*stack)[*stackPointer].frequency=intervalSizeOfLeft[1];
	for (i=0; i<=5; i++) (*stack)[*stackPointer].frequency_right[i]=rightMaximalString->frequency_leftRight[1][i];
	*stackPointer=*stackPointer+1;
	return intervalSizeOfLeft[1];
}


/**
 * Tries to push $bW$ onto $stack$, where $b \in {C,G,T,N}$.
 *
 * @param b the character ID >=2 of the character to push;
 * @param stackPointer pointer to the first free frame in the stack; the procedure 
 * increments $stackPointer$ at the end;
 * @return 0 if $bW$ was not pushed on the stack; otherwise, the size of the BWT interval
 * of $bW$.
 */
static inline unsigned int pushNonA(unsigned char b, const RightMaximalString_t *rightMaximalString, const Basic_BWT_t *bwt, StackFrame_t **stack, unsigned int *stackSize, unsigned int *stackPointer, const unsigned int length, const unsigned int *rankPoints, const unsigned int *rankValues, const unsigned char *nRightExtensionsOfLeft, const unsigned int *intervalSizeOfLeft) {
	unsigned int i;
	
	if (!isLeftExtensionRightMaximal(b,rightMaximalString,nRightExtensionsOfLeft)) return 0;
	if (*stackPointer>=*stackSize) {
		*stackSize=(*stackSize)<<=1;
		*stack=(StackFrame_t *)realloc(*stack,sizeof(StackFrame_t)*(*stackSize));
	}
	(*stack)[*stackPointer].firstCharacter=b;
	(*stack)[*stackPointer].length=length;
	(*stack)[*stackPointer].bwtStart=bwt->char_base[b-1]+rankValues[b-1]+1;
	(*stack)[*stackPointer].frequency=intervalSizeOfLeft[b];
	for (i=0; i<=5; i++) (*stack)[*stackPointer].frequency_right[i]=rightMaximalString->frequency_leftRight[b][i];
	*stackPointer=*stackPointer+1;
	return intervalSizeOfLeft[b];
}


void run(UnaryIterator_t *iterator) {
	const Basic_BWT_t *BWT = iterator->BBWT;
	const unsigned int MAX_LENGTH = iterator->maxLength;
	unsigned char i;
	unsigned char maxIntervalID, nExplicitWL, containsSharp, rightExtensionBitmap, npref_query_points;
	unsigned int length, nTraversedNodes, intervalSize, maxIntervalSize;  // Should be 64-bit
	unsigned int stackPointer;  // Pointer to the first free frame
	unsigned int stackSize;  // In frames
	RightMaximalString_t rightMaximalString = {0};
	StackFrame_t *stack;
	unsigned int rankPoints[7];  // Should be 64-bit
	unsigned int rankValues[28];  // Should be 64-bit
	unsigned int rankValuesN[7];  // Should be 64-bit
	unsigned char nRightExtensionsOfLeft[6];
	unsigned int intervalSizeOfLeft[6];  // Should be 64-bit
	
	// Initializing the stack
	stackSize=MIN_SLT_STACK_SIZE;
	stack=(StackFrame_t *)malloc((1+MIN_SLT_STACK_SIZE)*sizeof(StackFrame_t));
	stack[0].firstCharacter=0;
	stack[0].length=0;
	stack[0].bwtStart=0;
	stack[0].frequency_right[0]=1;
	for (i=1; i<=4; i++) stack[0].frequency_right[i]=BWT->char_base[i]-BWT->char_base[i-1];
	stack[0].frequency_right[5]=BWT->textlen-BWT->char_base[4];
	stack[0].frequency=BWT->textlen+1;
	
	// Traversing the suffix-link tree
	nTraversedNodes=0; stackPointer=1;
	do {
		nTraversedNodes++;
		stackPointer--;
		
		// Computing ranks
		getRanksOfRightExtensions(&stack[stackPointer],BWT,&rightExtensionBitmap,rankPoints,&npref_query_points,rankValues,rankValuesN,&containsSharp);
		
		// Issuing the callback function on the top of the stack
		memset(rightMaximalString.frequency_leftRight,0,sizeof(rightMaximalString.frequency_leftRight));
		memset(nRightExtensionsOfLeft,0,sizeof(nRightExtensionsOfLeft));
		memset(intervalSizeOfLeft,0,sizeof(intervalSizeOfLeft));
		buildCallbackState(&rightMaximalString,&stack[stackPointer],BWT,rightExtensionBitmap,rankPoints,npref_query_points,rankValues,rankValuesN,containsSharp,nRightExtensionsOfLeft,intervalSizeOfLeft);
		iterator->SLT_callback(rightMaximalString,iterator->applicationData);
		
		// Pushing $aW$ for $a \in {A,C,G,T}$ only, if it exists and it is right-maximal.
		length=rightMaximalString.length+1;
		if (length>MAX_LENGTH) continue;
		maxIntervalSize=pushA(&rightMaximalString,BWT,&stack,&stackSize,&stackPointer,length,rankPoints,rankValues,nRightExtensionsOfLeft,intervalSizeOfLeft);
		maxIntervalID=0;
		nExplicitWL=!!maxIntervalSize;
		for (i=2; i<=4; i++) {
			intervalSize=pushNonA(i,&rightMaximalString,BWT,&stack,&stackSize,&stackPointer,length,rankPoints,rankValues,nRightExtensionsOfLeft,intervalSizeOfLeft);
			if (!intervalSize) continue;
			if (intervalSize>maxIntervalSize) {
				maxIntervalSize=intervalSize;
				maxIntervalID=nExplicitWL;
			}
			nExplicitWL++;	
		}
		if (!nExplicitWL) continue;
		
		// Sorting the new left-extensions, if required.
		if (TRAVERSAL_ORDER==1) {
			if (maxIntervalID) swapStackFrames(&stack[stackPointer-nExplicitWL],&stack[stackPointer-nExplicitWL+maxIntervalID]);
		}
		else if (TRAVERSAL_ORDER==2) {
			for (i=0; i<nExplicitWL>>1; i++) swapStackFrames(&stack[stackPointer-nExplicitWL+i],&stack[stackPointer-1-i]);
		}
	} while (stackPointer);
	free(stack);
	
	printf("The number of traversed suffix tree nodes is %d \n",nTraversedNodes);
}