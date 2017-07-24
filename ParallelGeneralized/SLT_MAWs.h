#ifndef SLT_MAWs_h
#define SLT_MAWs_h
#include"SLT.h"

void SLT_callback_MAWs(const SLT_joint_params_t * SLT_params,void * intern_state, unsigned int memory);
void SLT_callback_MAWs_present(const SLT_joint_params_t * SLT_params,void * intern_state, unsigned int memory);
void SLT_callback_MAWs_kernel(const SLT_joint_params_t * SLT_params,void * intern_state, unsigned int memory);
void* SLT_cloner(void* p, unsigned int t);
void SLT_combiner(void** intern_state, void* state, unsigned int t, unsigned int mem);
void SLT_free(void* intern_state, unsigned int mem);
unsigned int SLT_find_MAWs(Basic_BWT_t * BBWT1,Basic_BWT_t * BBWT2,
		unsigned int minlen, unsigned int * nMAWs1,unsigned int * nMAWs2, double * LW,unsigned int mem, unsigned int cores, unsigned int result);
void convert_MAWs_to_ACGT(unsigned char ** MAW_ptr,unsigned int nMAWs);
double g1(int y);
double g2(int y);

#endif
