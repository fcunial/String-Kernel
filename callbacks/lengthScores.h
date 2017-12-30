#ifndef lengthScores_h
#define lengthScores_h


#include <math.h>


double LENGTH_SCORE_EPSILON = 0.1;


/**
 * Length score used in \cite{crochemore2016linear}.     
 */
double lengthScore1(unsigned int length) {
	return 1.0/(length*length);
}


/**
 * Length score used in \cite{smola2003fast}.     
 */
double lengthScore2(unsigned int length) {
	return pow(LENGTH_SCORE_EPSILON,length);
}


#endif