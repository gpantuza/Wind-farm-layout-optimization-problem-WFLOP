
#include <cstdio>
#include <set>
#include <list>
#include <vector>
#include <cmath>
#include <sys/times.h>
//#include "gurobi_c.h"
//#include "gurobi_c++.h"
//#include "wflop.hpp"
//#include "wflopUTIL.h"


//imprime relatorio
 void writeLog(char *fname, int k, double lvalue, double pcost, double tempo);

void writeLogCal(char *fname, int k, double lvalue, double pcost, double tempo, int iter, int green, int yellow, int red, int aint, double lmd, double alphai, double alpham, double alphaf, double bestCost);

/* Gera inteiro aleatorio entre min e max */
int randint(int mini, int maxi);

int *cria_vetor(int tam);

void iniciaVetorDouble(double *vetor, int tam, double v);
void iniciaVetor(int *vetor, int tam, int v);
double *criaVetorDouble(int tam);
int *criaVetor(int tam);
//
