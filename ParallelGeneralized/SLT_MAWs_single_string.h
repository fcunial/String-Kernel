#ifndef SLT_MAWs_h
#define SLT_MAWs_h
#include"SLT_single_string.h"

void SLT_MAWs_callback(const SLT_params_t * SLT_params,void * intern_state, unsigned int mem);
unsigned int SLT_find_MAWs_single_string(Basic_BWT_t * BBWT1, unsigned int minlen, unsigned int * _nMAWs1,
		double * _output_result, unsigned int mem);
void convert_MAWs_to_ACGT(unsigned char ** MAW_ptr,unsigned int nMAWs);


#endif
