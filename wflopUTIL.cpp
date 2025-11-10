#include <cstdio>
#include <set>
#include <list>
#include <vector>
#include <cmath>
#include <sys/times.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
//#include <time.h>
//#include "gurobi_c.h"
//#include "gurobi_c++.h"
#include "wflopUTIL.h"


#define MAX(x,y) ((x)<(y) ? (y) : (x))
/////////////////

void writeLog(char *fname, int k, double lvalue, double pcost, double tempo)
{
   
  char pname[100];
  sprintf(pname, "/home/gpantuza/wflop/wflopVOL/wflopDATA/Moseti95Case%s_BestSol.dat",fname);     
  
  FILE *f = NULL;

  f = fopen(pname, "r");
  FILE *file = fopen(pname, "r");
  if (!file) {
     printf("Failure to open WFLOP datafile: %s\n ", pname);
     abort();
  }

  double bestSol = 0;
  int pare = 0;
  while(pare == 0){
     //fscanf(f,"%d	%lf",&bestSol);
     fscanf(f,"%d %lf",&pare, &bestSol);
     if(pare == k) break;
     else pare = 0;
  }
  fclose(f);
  //printf("\nk  %d  Pare %d   bs %f",k,pare,bestSol);

  //char status[10];
  double gapBS = -111;
  double gapL = -111;
  double gapVOL = -111;
  
  if(bestSol > 0) gapL = -100*(bestSol - lvalue)/bestSol;
  if(bestSol > 0) gapBS = 100*(pcost - bestSol)/bestSol;
  if (pcost == lvalue){
      gapVOL = 0.0;
      //status = "Optimal";
  }
  else if (pcost > 0) gapVOL = 100*(pcost - lvalue)/pcost;



  FILE *arq = NULL;
  arq = fopen("wflopVol.log","a+");
  if (!arq) printf("O Arquivo rel %s nao pode ser aberto.\n","Relatorio.log");
  fprintf(arq,"vol; %s; %d; %f; %f; %1.0f; %2.2f; %2.2f; %2.2f\n", 
                  fname, k, pcost, lvalue, tempo, gapVOL, gapBS, gapL);
  fclose(arq);

}

void writeLogCal(char *fname, int k, double lvalue, double pcost, double tempo, int iter, int green, int yellow, int red, int aint, double lmd, double alphai, double alpham, double alphaf, double bestCost)
{

  char pname[100];
  sprintf(pname, "/home/gpantuza/wflop/wflopVOL/wflopDATA/Moseti95Case%s_BestSol.dat",fname);     
 
  FILE *f = NULL;

  f = fopen(pname, "r");
  //FILE *file = fopen(pname, "r");
  if (!f) {
     printf("Failure to open WFLOP datafile: %s\n ", pname);
     abort();
  }

  double bestSol = 0;
  int pare = 0;
  while(pare == 0){
     //fscanf(f,"%d	%lf",&bestSol);
     fscanf(f,"%d %lf",&pare, &bestSol);
     if(pare == k) break;
     else pare = 0;
  }
  fclose(f);
  //printf("\nk  %d  Pare %d   bs %f",k,pare,bestSol);

  double gapBS = 100*(bestSol - pcost)/bestSol;
  double gapL = 100*(bestSol - lvalue)/bestSol;
  double gapVOL = -1;
  if (pcost > 0) gapVOL = 100*(pcost - lvalue)/pcost;
  
  
  FILE *ar = NULL;
  ar = fopen("/home/gpantuza/wflop/wflopVOL/wflopDATA/mTurner_LP.dat", "r");
  if (!ar) {
     printf("Failure to open WFLOP datafile: %s\n ", "mTurner_LP.dat");
     abort();
  }
  
  double lpSol = -1;
  pare = 0;      //Moseti95Case1
  int cse = -1;
  int cs = atoi(fname) ;
  while(pare == 0){
     fscanf(ar,"%d %d %lf",&cse, &pare, &lpSol);
     //printf("\n cse %d  cs %d; pare %d   k %d,  lpspol %f",cse, cs, pare, k, lpSol); getchar();
     if(pare == k && cse == cs)  break;
     else pare = 0;
  }
  fclose(ar);

  double gapLP = 100*(lpSol - lvalue)/lpSol;

  /**/


  FILE *arq = NULL;
  arq = fopen("wflopVolCal.log","a+");
  if (!arq) printf("O Arquivo rel %s nao pode ser aberto.\n","wflopVolCal.log");
  fprintf(arq,"vol; %s; %d; %f; %f; %f; %1.0f; %2.2f; %2.2f; %2.2f; %2.2f; %d; %d; %d; %d; %d; %1.4f; %1.4f; %1.4f; %1.4f \n", fname, k, bestCost, pcost, lvalue, tempo, gapBS, gapL,gapLP, gapVOL, iter, green, yellow, red, aint, lmd, alphai, alpham, alphaf);
  fclose(arq);


}
/*-------------------------------------------------------------------------------------------*/
/* Gera inteiro aleatorio entre min e max */
int randint(int mini, int maxi)
{
  int r;
  r = mini +(rand()% maxi);

  return r;
}





////////////////////////////////////////////////////////////////////////////////////////////
/* cria memoria para um vetor de tam posicoes */
int *cria_vetor(int tam)
{
  int *vetor;

  vetor = (int *) malloc(tam*sizeof(int));
  if (!vetor){
  	printf("Falta memoria para alocar o vetor de ponteiros");
    exit(1);
  }
  return vetor;
}
/*-------------------------------------------------------------------------------------------*/
/* Cria vetor inteiro*/
int *criaVetor(int tam)
{
        int *vetor;
        vetor = (int *) malloc(tam*sizeof(int));
        if (!vetor) printf("Falta memoria para alocar o vetor de ponteiros");
        return vetor;
}
/*-------------------------------------------------------------------------------------------*/
/* Cria vetor double*/
double *criaVetorDouble(int tam)
{
        double *vetor;
        vetor = (double*) malloc(tam*sizeof(double));
        if (!vetor) printf("Falta memoria para alocar o vetor de ponteiros");
        return vetor;
}
/*-------------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------*/
/* Inicializa o vetor int com o valor v*/
void iniciaVetor(int *vetor, int tam, int v)
{
    int j;
    for (j=0; j<tam; j++) vetor[j] = v;
}

/*-------------------------------------------------------------------------------------------*/
/* Inicializa o vetor double com o valor v*/
void iniciaVetorDouble(double *vetor, int tam, double v)
{
    int j;
    for (j=0; j<tam; j++) vetor[j] = v;
}
/*-------------------------------------------------------------------------------------------*/


/////////////////////////////////////////////////////
///////////////////////////////////////////////////
//////////////////////////////////////////////////
/////////////////////////////////////////////////
///////////////////////////////////////////////
