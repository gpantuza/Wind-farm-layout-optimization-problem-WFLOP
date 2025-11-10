// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).
// Code based on Francisco Barahona implementation

// This is an implementation of the Volume algorithm for WFLOP (wind farm layout optimization problem) 
// Guido Pantuza - April/2024 

///////////////////////////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <set>
#include <list>
#include <vector>
#include <cmath>
#include <sys/times.h>
//#include <time.h>

//#include "gurobi_c.h"
//#include "gurobi_c++.h"
#include "wflop.hpp"
//#include "wflopUTIL.h"

//--------------------------------------------------------------------------------------------------

//Load data
int WFLOPreadData(const char* fname, WFLOP& data);

//--------------------------------------------------------------------------------------------------


int main(int argc, char* argv[]) {

  printf("\n\n WFLOP - VOLUME \n");
  srand(time(NULL));
  int erro = 0;

  WFLOP  wflop_data;
  wflop_data.k = strtol(argv[2], NULL, 10); //le K
  wflop_data.fname = argv[1];
  char pname[200];
  sprintf(pname, "/home/gpantuza/wflop/wflopVOL/wflopDATA/Moseti95Case%s.dat",wflop_data.fname);  

  erro = WFLOPreadData(wflop_data.fname, wflop_data);
  const int n = wflop_data.n;
  //wflop_data.hash.allocate(100*n);
  //wflop_data.hash = -1;
 
  VOL_problem volp;

//inici parametrosdo alg
  volp.psize = n*n+n ;      //# primal var's
//  volp.dsize = n*n;   //# dual var's (relaxed constraints) Formulação compacta
  volp.dsize = 3*n*n;   //# dual var's (relaxed constraints)
  bool ifdual = false; 
  wflop_data.iter = 0;
  wflop_data.icost = 9999999;
  //wflop_data.bestCost = 9999999;
  //wflop_data.bSol.allocate(n);
  //for(int i = 0 ; i < n; ++i)   wflop_data.bSol[i] = 0;

  volp.parm.greentestinvl= 1;
  volp.parm.yellowtestinvl= 400; //200, 300, 400, 500
  volp.parm.redtestinvl= 50; //20, 50 //cal: 10, 20, 35, 50, 70
  volp.parm.alphaint= 50;  //50, 100, 150, 200, 300, 400 ,500
  volp.parm.lambdainit= 0.1 ;
  volp.parm.alphainit= 0.1;
  volp.parm.alphamin= 0.0001;
  volp.parm.alphafactor= 0.99; 
 
  ///VOL_PARAMETERS: print
  volp.parm.printflag = 0; //log level
  //VOL_PARAMETERS: heuristic
  volp.parm.heurinvl= 100; // n: number of times that the primal heuristic is run 
  //VOL_PARAMETERS: stop criteria
  volp.parm.maxsgriters=99999;
  volp.parm.primal_abs_precision=0.01;
  volp.parm.gap_abs_precision=0.01;
  volp.parm.gap_rel_precision=0.01;
  volp.parm.granularity=0.01;

   // allocates memo: index matriz
   volp.dual_lb.allocate(volp.dsize);
   //volp.dual_lb = -1e31;
   volp.dual_ub.allocate(volp.dsize);
   volp.dual_ub = 1e31;
   for(int i = 0 ; i < volp.dsize; ++i)   volp.dual_lb[i] = 0;

  
  //tempo de execução;
  clock_t inicio, fim;
  struct tms timearr; 
  clock_t tres;
  inicio = clock();
  tres = times(&timearr); 
  wflop_data.t0 = timearr.tms_utime;   

  // invoke volume algorithm
  printf(" invoke volume algorithm ....\n");   
  if (volp.solve(wflop_data, ifdual) < 0) 
  //getchar();
  
  //printf("\n Best integer solution value: %f\n", wflop_data.icost);
  //printf(" Lower bound: %f\n", volp.value);

  // end time measurement
  fim = clock();
  double tempo_decorrido = ((double)(fim - inicio) / CLOCKS_PER_SEC);
  tres = times(&timearr);
  double t =  (timearr.tms_utime-wflop_data.t0)/100.;
  printf("\nTotal Time: %f secs\n", t);
  printf("\nTotal Time: %f secs\n", tempo_decorrido);
  
  writeLog(wflop_data.fname, wflop_data.k, volp.value, wflop_data.icost, t);

   printf("---FIM-----------\n");

   return 0;
}
//############################################################################
//Load data
int WFLOPreadData(const char* fname, WFLOP& data)
{

   int erro;
   FILE *f = NULL;
   char name[100];
   sprintf(name, "/home/gpantuza/wflop/wflopVOL/wflopDATA/Moseti95Case%s.dat",fname); 

   VOL_dvector& w = data.w;
   VOL_dvector& ws = data.ws; 
   VOL_dvector& p = data.p;
   int& n = data.n;
   int& k = data.k;
   int& wd = data.wd;
   int& caso = data.caso;
   //char& *fname = ;
   double& alpha = data.alpha;
   double& r0 = data.r0;
   double& a = data.a;
   double& thetainit = data.thetainit;
   int& nws = data.nws;
   
   f = fopen(name, "r");
   FILE *file = fopen(name, "r");
   if (!file) {
      printf("Failure to open WFLOP datafile: %s\n ", name);
      return 1;
   }

   erro = fscanf(f,"%d\n",&caso);
   erro = fscanf(f,"%d\n",&wd);
   erro = fscanf(f,"%d\n",&n);
   erro = fscanf(f,"%lf %lf %lf %lf\n",&alpha, &r0, &a, &thetainit);
   erro = fscanf(f,"%d",&nws);

   w.allocate(n*n);
   ws.allocate(nws);
   p.allocate(wd);

   for(int i = 0; i < nws; ++i) erro = fscanf(f,"%lf",&ws[i]);
   for(int i = 0; i < wd; ++i)  erro = fscanf(f,"%lf\n",&p[i]);
   for(int i = 0; i < n*n; ++i) erro = fscanf(f,"%lf",&w[i]);

   fclose(f);

   return 0;
   
}
//############################################################################
//###### USER HOOKS
// compute reduced costs
int WFLOP::compute_rc(const VOL_dvector& u, VOL_dvector& rc)
{
   //printf("\n rc1 ...");
   
   for(int i = 0; i < n*n+n; ++i) rc[i] = 0;

   for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) if(i > j)
      rc[i*n+j] = w[i*n+j] - u[i*n+j] + u[n*n + i*n+j] + u[2*n*n + i*n+j];

   for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j){ 
      if(i > j)  rc[n*n+i] += (u[i*n+j] - u[n*n + i*n+j]); 
      else if(i < j)  rc[n*n+i] += (u[j*n +i] - u[2*n*n + i*n+j]);
   }

/* Formulação compacta
   for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) if(i > j)
      rc[i*n+j] = w[i*n+j] - u[i*n+j];

   for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j){ 
      if(i > j)  rc[n*n+i] += u[i*n+j] ; 
      else if(i < j)  rc[n*n+i] += u[j*n +i] ;
   }
*/

   //printf("... rc2\n");
   return 0;
}

//############################################################################
// IN: dual vector u
// OUT: primal solution to the Lagrangian subproblem (x)
//      optimal value of Lagrangian subproblem (lcost)
//      v = difference between the rhs and lhs when substituting
//                  x into the relaxed constraints (v)
//      objective value of x substituted into the original problem (pcost)
//      xrc
// return value: -1 (volume should quit) 0 ow

int WFLOP::solve_subproblem(const VOL_dvector& u, const VOL_dvector& rc, double& lcost, 
                          VOL_dvector& x, VOL_dvector& v, double& pcost)
{
   
   //printf("... solve_subproblem\n");
   lcost = 0.0;

   for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) 
       if(i > j) lcost -= u[i*n+j]; //soma parte constante da foD
 
   for (int i = 0; i < n*n+n; ++i) x[i] = 0; //atribui zero para todass as vars primais

//solve sub problem X
   for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) 
       if(rc[i*n+j] < 0)  x[i*n+j] = 1;
    
//solve sub problem Y
  for(int j = 0; j < k; ++j){
      int pos = 0;
      double menor = 999999;
      for(int i = 0; i < n; ++i) if(x[n*n+i] == 0) if(rc[n*n+i] < menor ){
         menor = rc[n*n+i];
         pos = i;
      }
      x[n*n+pos] = 1;
   }

//Calcula o valor dual
   for (int i = 0; i < n; ++i){
      if(x[n*n+i] == 1){ 
         lcost+= (rc[n*n+i])*(x[n*n+i]);
         //printf("\n Y %d : rc %f, y %f\n",i,rc[n*n+i],x[n*n+i]);
      }
      for (int j = 0; j < n; ++j) if(i > j){
         if(x[i*n+j] == 1){ 
            lcost+= (rc[i*n+j])*(x[i*n+j]);
         }
      }
   }

///compute primal solution
///Se o x[ij] = 0 p/ todo ij, então solucao primal também igual a zero
   pcost = 0;
   for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) 
       if(i > j) pcost += w[i*n+j]*x[i*n+j];
 
//Calcula vetor v = rhs - lhs (relaxed constraints)
   v = 0;
   double rhs = 0;
   double lhs = 0;
   for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) if(i > j){
      rhs = -1;
      lhs =  x[i*n+j] - x[n*n+i] - x[n*n+j];
      v[i*n+j] = rhs - lhs;
      
      rhs = 0;
      lhs =  x[n*n+i] - x[i*n+j]; 
      v[n*n + i*n+j] = rhs - lhs;
      
      rhs = 0;
      lhs =  x[n*n+j] - x[i*n+j]; 
      v[2*n*n + i*n+j] = rhs - lhs;

   }

//compute v = rhs - lhs (relaxed constraints)
//Formulação Compacta
/*   v = 0;
   double rhs = 0;
   double lhs = 0;
   for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) if(i > j){
      rhs = -1;
      lhs =  x[i*n+j] - x[n*n+i] - x[n*n+j];
      if( lhs < rhs) v[i*n+j] = rhs - lhs;  //Calcula V apenas se a restrição for violada
      //v[i*n+j] = rhs - lhs;                //Calcula V p/ toda restrição
   }
 /**/

  
/**/
   //Busca local VND
   double fo = 0;
   VOL_ivector sol;
   sol.allocate(n*n+n);
   sol = 0;
   VOL_ivector s;
   s.allocate(n); 
   s = 0;
   for(int i = 0; i < n; ++i) s[i] = x[n*n + i];
   VOL_ivector sl;
   sl.allocate(n); 
   sl = 0;
   double fol = 99999;  
   
   int aux = 0;
   for (int i = 0; i < n; ++i) aux += s[i];

   if(aux == k){
   //if(aux == k && (iter%2) == 0){
       
      fo = 0;
      for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) 
          if(i > j) if(s[i] == 1 && s[j] == 1) fo += w[i*n+j];
      
      fol = vndV2(n, s, fo, sl);
            
      aux = 0;
      for(int i = 0; i < n; ++i){ 
          sol[n*n + i] = sl[i];
          aux +=sl[i];
      }
      
      for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j)
   	  if(i > j) if(sol[n*n+i] == 1 && sol[n*n+j] == 1) sol[i*n+j] = 1;
     
      if (fol < icost && aux == k) {
          icost = fol;
          ix = sol;
      }
   }
/**/ 
   ++iter;
   //printf("\n iter %d",iter);
   return 0;
   
   
}

//############################################################################
// IN:  fractional primal solution (x),
//      best feasible soln value so far (icost)
// OUT: integral primal solution (ix) if better than best so far
//      and primal value (icost)
// returns -1 if Volume should stop, 0/1 if feasible solution wasn't/was
//         found.
int WFLOP::heuristics(const VOL_problem& p, const VOL_dvector& x, double& new_icost)
{
//constroi solução inicial inteira a partir da solução relaxada x, apenas p/ Y
   //printf("\n heuristica"); 
   double value = 0;
   double fo = 0;
   VOL_ivector sol;
   sol.allocate(n*n+n);
   for(int i = 0; i < n*n+n; ++i) sol[i] = 0;
   int ct = 0;
   VOL_ivector sl;
   sl.allocate(n); 
   for(int i = 0; i < n; ++i) sl[i] = 0;
   VOL_ivector s;
   s.allocate(n); 
   for(int i = 0; i < n; ++i) s[i] = 0;

//Verifica os valores iguais a 1 primeiro (apenas para var Yi)
   ct = 0;
   for(int i = 0; i < n; ++i) if(x[n*n+i] == 1){
      sol[n*n+i] = 1; 
      ++ct;
   }
  
//Verifica e seleciona os k restantes de acordo com o maior valor de x
   int pos = -1;
   for(int j = 0; j < (k - ct); ++j){
      double maior = 0;
      for(int i = 0; i < n; ++i) if(x[n*n+i] != 1 && sol[n*n+i] != 1) if(x[n*n+i] > maior){
         maior = x[n*n+i];
         pos = i;
      }
      sol[n*n+pos] = 1;
   } 

//Calcula o valor da solução inicial     
   value = 0;                
   for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j){ 
   	if(i > j) if(sol[n*n+i] == 1 && sol[n*n+j] == 1){
   	    value += w[i*n+j];
            sol[i*n+j] = 1;	    
   	}
   }	
   int aux = 0;
   for (int i = 0; i < n; ++i) aux += sol[n*n+i];
   //printf("\nSOL: aux %d  k %d",aux,k);//getchar();

/**/
//Atualiza a melhor solução até o momento
   new_icost = value;
   if (value < icost && aux == k) {
      icost = value;
      ix = sol;  
   }
   
/**/
//Busca local VND
   for(int i = 0; i < n; ++i) s[i] = sol[n*n + i];
   aux = 0;
   for (int i = 0; i < n; ++i) aux += s[i];
   if(aux == k){
       
       //fo = vnd(n, s, value, sl);
       fo = vndV2(n, s, value, sl);

       aux = 0;
       for (int i = 0; i < n; ++i) aux += sl[i];
   
       for(int i = 0; i < (n*n+n); ++i) sol[i] = 0;
       for(int i = 0; i < n; ++i) sol[n*n + i] = sl[i];
       for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j)
   	   if(i > j) if(sol[n*n+i] == 1 && sol[n*n+j] == 1) sol[i*n+j] = 1;
 
       new_icost = fo;
       if (fo < icost) {
          icost = fo;
          ix = sol;
       }
   } //else goto QUIT;
/**/  

/**/
///Agora o VNS
   int it = 0;
   int MelhorIter = 0;
   int ki = 1;
   int VNS_max = 100; //numero máximo de iterações sem melhora (a calibrar)
   int r = 6;        // número de pertubações (a calibrar)
   double fo2l = 9999;
   double fol = 9999;
   int nvezes = 0;
   ct = 0;
   int ctmax = int (VNS_max/2);
   VOL_ivector s2l;
   s2l.allocate(n);
   s2l = 0;
   for(int i = 0; i < n; ++i) s[i] = sol[n*n + i];
   
   aux = 0;
   for (int i = 0; i < n; ++i) aux += s[i];
   //printf("\n aux %d  k %d",aux,k);getchar();
   if(aux != k) goto QUIT; 
   
   while (it - MelhorIter < VNS_max){
    
     //printf("\n it %d  ki %d  k %d",it,ki,k);getchar(); 
     it++;
     ki = 1;
     
     while (ki <= r){
        //printf("\n it %d  ki %d <= %d r",it,ki, r); 
        
        for(int i = 0; i < n; ++i) sl[i] = s[i]; 
        fol = gera_k_vizinhos_aleatorios(n, sl, ki); // Escolher viz qualquer na k-vizinhanca 
        aux = 0;
        for (int i = 0; i < n; ++i) aux += sl[i];
        if(aux != k) goto QUIT;
   
        //fo2l = vnd(n, sl, fol, s2l);    
        fo2l = vndV2(n, sl, fol, s2l);
        aux = 0;
        for (int i = 0; i < n; ++i) aux += s2l[i];
        if(aux != k) goto QUIT;

        if (fo2l < fo){
           fo = fo2l;
           for(int i = 0; i < n; ++i) s[i] = s2l[i]; 
           MelhorIter = it;
           //printf("Iter VNS = %3d \t fo_star = %10.4f \t Vizinhanca = %3d\n", iter, fo, ki);
           ki = 1; //Basic-VNS
           nvezes = 1;
           
        } else {
           if(nvezes > 10) {
               ++ki;
               nvezes = 1;
           } else{ 
               ++nvezes;
           }
           ++ct;
        }     
     }
     
     //Adicionar: se várias iterações sem melhora, recomeçar vns da melhor solução ix 
     if(ct == ctmax){
         for(int i = 0; i < n; ++i) s[i] = ix[n*n+i];
         fo = icost;
         MelhorIter = it;
         ct = 0;
     }  
   }

   //for(int i = 0; i < n; ++i) if(s[i]) printf("\nFIM VNS s[%d] %d",i,s[i]);
   for(int i = 0; i < (n*n+n); ++i) sol[i] = 0;
   for(int i = 0; i < n; ++i) sol[n*n + i] = s[i];
   for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j)
   	if(i > j) if(sol[n*n+i] == 1 && sol[n*n+j] == 1) sol[i*n+j] = 1;

   new_icost = fo;
   if (fo < icost) {
       icost = fo;
       ix = sol;
   }
/**/
   //printf("\nFIM heuristica");
   //if( icost <= dcost) return -1;
   return 1;
   
   QUIT:
   //printf("\nFIM heuristica");
   return 0;
  
}

/***************************************************************************/
//############################################################################
//############################################################################
/***************************************************************************/
//Busca local VND
//Entrada: solução corrente através do vetor sol, e a melhor fo até o momento (icost), 
//Retorna: fo se houver melhora (ou icost ) e o melhor vizinho (se houver melhora) através de sl
//         se não houver melhora, retorna a solução corrente sol através de sl 

double WFLOP::vndV2(int n, VOL_ivector s, double fo, VOL_ivector& sl)
{
    VOL_ivector s2l;
    s2l.allocate(n);
    s2l = 0;
    for (int i = 0; i < n; i++) sl[i] = s[i];
    double fos2l = 9999;
    int km = 1;
    double fo_aux = fo;

    while(km < 3) { 
       
       for (int i = 0; i < n; ++i) s2l[i] = sl[i];
       
       if(km == 1) fos2l = mov1Grid(s2l, n);
       if(km == 2) fos2l = mov2Grid(s2l, n);

       if (fos2l < fo_aux) {
          for (int i = 0; i < n; i++) sl[i] = s2l[i];
          fo_aux = fos2l;
          km = 1;
       }  else km++;
       
    }//while 

    //printf("-------------------------------------------------\n");
    //printf("Melhor solucao encontrada pelo VND:\n");
    //printf("Funcao objetivo: %f\n", fo_aux);
    //printf("-------------------------------------------------\n");
    //for(int i = 0; i < n; ++i) if(sl[i]) printf("\nDentro: sl[%d] %d",i, sl[i]); //getchar(); 

    return fo_aux;
    
}//vnd
/***************************************************************************/
//############################################################################
//############################################################################
/***************************************************************************/
double WFLOP::gera_um_vizinho_aleatorio(int n, VOL_ivector& s)
{
  int i,j,aux;
  double fo_viz;
  
  do{
    i = randint(0,99);
    j = randint(0,99);
    //printf(" loop...  ");
  } while(s[i] == 0 || s[j] == 1);

  //printf("\n s[%d] = %d (0); s[%d] = %d (1)",i,s[i],j,s[j]); // getchar();
  s[i] = 0;
  s[j] = 1;
  
  fo_viz = 0;
  for(i = 0; i < n; ++i) for(j = 0; j < n; ++j) if(i > j) 
     if(s[i] == 1 && s[j] == 1) fo_viz += w[i*n+j];
  
  return fo_viz;
}// fim gera_um_vizinho_aleatorio
/***************************************************************************/
//############################################################################
//############################################################################
/***************************************************************************/
double WFLOP::gera_k_vizinhos_aleatorios(int n, VOL_ivector& s, int ki)
{
  double fo_viz;
  for (int i=0; i<ki; i++) fo_viz = gera_um_vizinho_aleatorio(n,s);

  return fo_viz;

}// fim gera_k_vizinhos_aleatorios
/***************************************************************************/
//############################################################################
//############################################################################
/***************************************************************************/
//Analisa todos os 8 vizinhos e retorna o melhor viz
double  WFLOP::mov1Grid(VOL_ivector& s2l, int n)
{
   double fo = 999;
   int i, j ;
   VOL_ivector stemp;
   stemp.allocate(n);
   stemp = 0;
   VOL_dvector foV;
   foV.allocate(9);
   foV = 9999;
   int mviz = 8;
   //for(int ii = 0; ii < 8; ++ii) foV[ii] = 9999;
   for(int ii = 0; ii < n; ++ii)  stemp[ii] = s2l[ii] ;
   
   do{
       i = randint(0,9); //mudar para (0,9) quando todos estiverem prontos
       j = randint(0,9);
   }while ( stemp[i*10+j] == 0 ); //enquanto for zero ou i menor que j	
   
   stemp[(i)*10+(j)] = 0;     //Atribui zero para o escolhido
   
   if( 0 < i && i < 9 && 0 < j && j < 9) {
      //faz o movimento (8 vizinhos) 
      if( stemp[(i-1)*10+(j-1)] == 0){   
          stemp[(i-1)*10+(j-1)] = 1; //faz o movimento 1° vizinho
          foV[0] = 0;
          for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
              if(stemp[u] == 1 && stemp[v] == 1) foV[0] += w[u*n+v];
          stemp[(i-1)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho
      }  
         
      if( stemp[(i-1)*10+(j)] == 0){     
	      stemp[(i-1)*10+(j)] = 1; //faz o movimento 2° vizinho
	      foV[1] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[1] += w[u*n+v];
	      stemp[(i-1)*10+(j)] = 0;
      }
      
      if( stemp[(i-1)*10+(j+1)] == 0){
	      stemp[(i-1)*10+(j+1)] = 1; //faz o movimento 3° vizinho
	      foV[2] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[2] += w[u*n+v];
	      stemp[(i-1)*10+(j+1)] = 0; 
      }      
      
      if( stemp[(i)*10+(j-1)] == 0){
	      stemp[(i)*10+(j-1)] = 1; //faz o movimento 2° vizinho
	      foV[3] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[3] += w[u*n+v];
	      stemp[(i)*10+(j-1)] = 0;
      }
      
      if( stemp[(i)*10+(j+1)] == 0){      
		  stemp[(i)*10+(j+1)] = 1; //faz o movimento 2° vizinho
		  foV[4] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[4] += w[u*n+v];
		  stemp[(i)*10+(j+1)] = 0;      
      }
      if( stemp[(i+1)*10+(j-1)] == 0){
		  stemp[(i+1)*10+(j-1)] = 1; //faz o movimento 2° vizinho
		  foV[5] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[5] += w[u*n+v];
		  stemp[(i+1)*10+(j-1)] = 0;      
      }
      if( stemp[(i+1)*10+(j)] == 0){
		  stemp[(i+1)*10+(j)] = 1; //faz o movimento 2° vizinho
		  foV[6] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[6] += w[u*n+v];
		  stemp[(i+1)*10+(j)] = 0;
      }
      if( stemp[(i+1)*10+(j+1)] == 0){
		  stemp[(i+1)*10+(j+1)] = 1; //faz o movimento 2° vizinho
		  foV[7] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[7] += w[u*n+v];
		  stemp[(i+1)*10+(j+1)] = 0;
      }        
   } else if(i == 0 && j != 0 && i != 9){
   
      if( stemp[(i)*10+(j-1)] == 0){
		  stemp[(i)*10+(j-1)] = 1; //faz o movimento 2° vizinho
		  foV[3] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[3] += w[u*n+v];
		  stemp[(i)*10+(j-1)] = 0;
      }
      if( stemp[(i)*10+(j+1)] == 0){
		  stemp[(i)*10+(j+1)] = 1; //faz o movimento 2° vizinho
		  foV[4] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[4] += w[u*n+v];
		  stemp[(i)*10+(j+1)] = 0;      
      }
      if( stemp[(i+1)*10+(j-1)] == 0){
		  stemp[(i+1)*10+(j-1)] = 1; //faz o movimento 2° vizinho
		  foV[5] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[5] += w[u*n+v];
		  stemp[(i+1)*10+(j-1)] = 0;      
      }
      if( stemp[(i+1)*10+(j)] == 0){
		  stemp[(i+1)*10+(j)] = 1; //faz o movimento 2° vizinho
		  foV[6] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[6] += w[u*n+v];
		  stemp[(i+1)*10+(j)] = 0;
      }
      if( stemp[(i+1)*10+(j+1)] == 0){
		  stemp[(i+1)*10+(j+1)] = 1; //faz o movimento 2° vizinho
		  foV[7] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[7] += w[u*n+v];
		  stemp[(i+1)*10+(j+1)] = 0;
      }
   } else if(i == 0 && j == 0 ){
   
      if( stemp[(i)*10+(j+1)] == 0){
		  stemp[(i)*10+(j+1)] = 1; //faz o movimento 2° vizinho
		  foV[4] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[4] += w[u*n+v];
		  stemp[(i)*10+(j+1)] = 0;      
      }
      if( stemp[(i+1)*10+(j)] == 0){
		  stemp[(i+1)*10+(j)] = 1; //faz o movimento 2° vizinho
		  foV[6] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[6] += w[u*n+v];
		  stemp[(i+1)*10+(j)] = 0;
      }
      if( stemp[(i+1)*10+(j+1)] == 0){
		  stemp[(i+1)*10+(j+1)] = 1; //faz o movimento 2° vizinho
		  foV[7] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[7] += w[u*n+v];
		  stemp[(i+1)*10+(j+1)] = 0;
      }
   } else if(i == 0 && j == 9 ){
   
      if( stemp[(i)*10+(j-1)] == 0){
		  stemp[(i)*10+(j-1)] = 1; //faz o movimento 2° vizinho
		  foV[3] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[3] += w[u*n+v];
		  stemp[(i)*10+(j-1)] = 0;
      }
      if( stemp[(i+1)*10+(j-1)] == 0){
          stemp[(i+1)*10+(j-1)] = 1; //faz o movimento 2° vizinho
	      foV[5] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[5] += w[u*n+v];
	      stemp[(i+1)*10+(j-1)] = 0;      
      }
      if( stemp[(i+1)*10+(j)] == 0){
	      stemp[(i+1)*10+(j)] = 1; //faz o movimento 2° vizinho
	      foV[6] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[6] += w[u*n+v];
	      stemp[(i+1)*10+(j)] = 0;
      }
   } else if(i != 0 && i != 9 && j == 0){
      if( stemp[(i-1)*10+(j)] == 0){
	      stemp[(i-1)*10+(j)] = 1; //faz o movimento 2° vizinho
	      foV[1] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[1] += w[u*n+v];
	      stemp[(i-1)*10+(j)] = 0;
      }
      if( stemp[(i-1)*10+(j+1)] == 0){
	      stemp[(i-1)*10+(j+1)] = 1; //faz o movimento 3° vizinho
	      foV[2] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[2] += w[u*n+v];
	      stemp[(i-1)*10+(j+1)] = 0; 
      }
      if( stemp[(i)*10+(j+1)] == 0){
	      stemp[(i)*10+(j+1)] = 1; //faz o movimento 2° vizinho
	      foV[4] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[4] += w[u*n+v];
	      stemp[(i)*10+(j+1)] = 0;      
      }      
      if( stemp[(i+1)*10+(j)] == 0){
	      stemp[(i+1)*10+(j)] = 1; //faz o movimento 2° vizinho
	      foV[6] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[6] += w[u*n+v];
	      stemp[(i+1)*10+(j)] = 0;
      }
      if( stemp[(i+1)*10+(j+1)] == 0){
	      stemp[(i+1)*10+(j+1)] = 1; //faz o movimento 2° vizinho
	      foV[7] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[7] += w[u*n+v];
	      stemp[(i+1)*10+(j+1)] = 0;
      }
   } else if(i != 0 && i != 9 && j == 9){
      if( stemp[(i-1)*10+(j-1)] == 0){
	      stemp[(i-1)*10+(j-1)] = 1; //faz o movimento 1° vizinho
	      foV[0] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[0] += w[u*n+v];
	      stemp[(i-1)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho
      }        
      if( stemp[(i-1)*10+(j)] == 0){
	      stemp[(i-1)*10+(j)] = 1; //faz o movimento 2° vizinho
	      foV[1] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[1] += w[u*n+v];
	      stemp[(i-1)*10+(j)] = 0;
      }   
      if( stemp[(i)*10+(j-1)] == 0){
	      stemp[(i)*10+(j-1)] = 1; //faz o movimento 2° vizinho
	      foV[3] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[3] += w[u*n+v];
	      stemp[(i)*10+(j-1)] = 0;
      }
      if( stemp[(i+1)*10+(j-1)] == 0){
	      stemp[(i+1)*10+(j-1)] = 1; //faz o movimento 2° vizinho
	      foV[5] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[5] += w[u*n+v];
	      stemp[(i+1)*10+(j-1)] = 0;      
      }
      if( stemp[(i+1)*10+(j)] == 0){
	      stemp[(i+1)*10+(j)] = 1; //faz o movimento 2° vizinho
	      foV[6] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[6] += w[u*n+v];
	      stemp[(i+1)*10+(j)] = 0;
      }
    } else if(i == 9 && j == 0){
  
      if( stemp[(i-1)*10+(j)] == 0){
	      stemp[(i-1)*10+(j)] = 1; //faz o movimento 2° vizinho
	      foV[1] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[1] += w[u*n+v];
	      stemp[(i-1)*10+(j)] = 0;
      }
      if( stemp[(i-1)*10+(j+1)] == 0){
	      stemp[(i-1)*10+(j+1)] = 1; //faz o movimento 3° vizinho
	      foV[2] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[2] += w[u*n+v];
	      stemp[(i-1)*10+(j+1)] = 0; 
      }
      if( stemp[(i)*10+(j+1)]  == 0){
	      stemp[(i)*10+(j+1)] = 1; //faz o movimento 2° vizinho
	      foV[4] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[4] += w[u*n+v];
	      stemp[(i)*10+(j+1)] = 0;      
      }
    } else if(i == 9 && j == 9){  
      if( stemp[(i-1)*10+(j-1)] == 0){
	      stemp[(i-1)*10+(j-1)] = 1; //faz o movimento 1° vizinho
	      foV[0] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[0] += w[u*n+v];
	      stemp[(i-1)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho
      }         
      if( stemp[(i-1)*10+(j)]  == 0){
	      stemp[(i-1)*10+(j)] = 1; //faz o movimento 2° vizinho
	      foV[1] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[1] += w[u*n+v];
	      stemp[(i-1)*10+(j)] = 0;
      }
      
      if( stemp[(i)*10+(j-1)] == 0){
	      stemp[(i)*10+(j-1)] = 1; //faz o movimento 2° vizinho
	      foV[3] = 0;
	      for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		  if(stemp[u] == 1 && stemp[v] == 1) foV[3] += w[u*n+v];
	      stemp[(i)*10+(j-1)] = 0;
      }
   }
 /**/
   
   //verifica o melhor vizinho
   for(int u = 0; u < 8; ++u){
       if(foV[u] < fo) {
           fo = foV[u];
           mviz = u;
       }
   }


 
   //Retorna s2l com o melhor vizinho
   s2l[i*10+j] = 0;
   if(mviz == 0){
      s2l[(i-1)*10+(j-1)] = 1; 
   } else if(mviz == 1) {
      s2l[(i-1)*10+(j)] = 1;
   } else if(mviz == 2) {
      s2l[(i-1)*10+(j+1)] = 1;
   } else if(mviz == 3) {
      s2l[(i)*10+(j-1)] = 1;
   } else if(mviz == 4) {
      s2l[(i)*10+(j+1)] = 1;
   } else if(mviz == 5) {
      s2l[(i+1)*10+(j-1)] = 1;
   } else if(mviz == 6) {
      s2l[(i+1)*10+(j)] = 1;
   } else if(mviz == 7) {
      s2l[(i+1)*10+(j+1)] = 1;
   } else if(mviz == 8) {
      s2l[i*10+j] = 1;
   }  
  /**/   
   return fo;
   
}
/***************************************************************************/
//############################################################################
//############################################################################
/***************************************************************************/
//Incompleto
//Analisa todos os 16 vizinhos e retorna o melhor viz
double  WFLOP::mov2Grid(VOL_ivector& s2l, int n)
{
   int i, j ;
   VOL_ivector stemp;
   stemp.allocate(n);
   VOL_dvector foV;
   foV.allocate(17);
   int mviz = 16;
   double fo = 9999;
   for(int ii = 0; ii < 16; ++ii) foV[ii] = 99999;   
   for(int ii = 0; ii < n; ++ii)  stemp[ii] = s2l[ii] ;
   
   do{
       i = randint(0,9); 
       j = randint(0,9);
   }while ( stemp[i*10+j] == 0 ); //enquanto for zero
   
   stemp[(i)*10+(j)] = 0;     //Atribui zero para o escolhido

//Calcula a Fo dos vizinhos
   if( (j > 1 && j < 8) && (i > 1 && i < 8) ){
       
      if(stemp[(i-2)*10+(j-2)]  == 0){
		  stemp[(i-2)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[0] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[0] += w[u*n+v];
		  stemp[(i-2)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }
       
      if(stemp[(i-2)*10+(j-1)] == 0){  
          stemp[(i-2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[1] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[1] += w[u*n+v];
		  stemp[(i-2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
      }
      
      if(stemp[(i-2)*10+(j)] == 0){   
		  stemp[(i-2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[2] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[2] += w[u*n+v];
		  stemp[(i-2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }
      
      if(stemp[(i-2)*10+(j+1)] == 0){ 
		  stemp[(i-2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[3] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[3] += w[u*n+v];
		  stemp[(i-2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
      }
      
      if(stemp[(i-2)*10+(j+2)] == 0){    
		  stemp[(i-2)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[4] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[4] += w[u*n+v];
		  stemp[(i-2)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }
         
      if(stemp[(i-1)*10+(j-2)] == 0){    
		  stemp[(i-1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[5] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[5] += w[u*n+v];
		  stemp[(i-1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }
	     
      if(stemp[(i-1)*10+(j+2)] == 0){   
		  stemp[(i-1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[6] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[6] += w[u*n+v];
		  stemp[(i-1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }
      
      if(stemp[(i)*10+(j-2)] == 0){    
		  stemp[(i)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[7] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[7] += w[u*n+v];
		  stemp[(i)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }
      
      if(stemp[(i)*10+(j+2)] == 0){ 
		  stemp[(i)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[8] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[8] += w[u*n+v];
		  stemp[(i)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }
      
      if(stemp[(i+1)*10+(j-2)] == 0){
		  stemp[(i+1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[9] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[9] += w[u*n+v];
		  stemp[(i+1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }
         
      if(stemp[(i+1)*10+(j+2)] == 0){
		  stemp[(i+1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[10] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[10] += w[u*n+v];
		  stemp[(i+1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }
      
      if(stemp[(i+2)*10+(j-2)]  == 0){      
		  stemp[(i+2)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[11] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[11] += w[u*n+v];
		  stemp[(i+2)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }
      
      if( stemp[(i+2)*10+(j-1)] == 0){ 
		  stemp[(i+2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[12] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[12] += w[u*n+v];
		  stemp[(i+2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
      }
       
      if( stemp[(i+2)*10+(j)] == 0){
		  stemp[(i+2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[13] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[13] += w[u*n+v];
		  stemp[(i+2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }
       
      if(stemp[(i+2)*10+(j+1)] == 0){
		  stemp[(i+2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[14] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[14] += w[u*n+v];
		  stemp[(i+2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
      }
      
      if(stemp[(i+2)*10+(j+2)]  == 0){
		  stemp[(i+2)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[15] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[15] += w[u*n+v];
		  stemp[(i+2)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }
      
   } else if( i == 0 && j == 0) {
        
      if(stemp[(i)*10+(j+2)] == 0){     
		  stemp[(i)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[8] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[8] += w[u*n+v];
		  stemp[(i)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho
      } 
      if(stemp[(i+1)*10+(j+2)] == 0){
		  stemp[(i+1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[10] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[10] += w[u*n+v];
		  stemp[(i+1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho
	  }    
      if(stemp[(i+2)*10+(j)] == 0){     
		  stemp[(i+2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[13] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[13] += w[u*n+v];
		  stemp[(i+2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }  
      if(stemp[(i+2)*10+(j+1)] == 0){
		  stemp[(i+2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[14] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[14] += w[u*n+v];
		  stemp[(i+2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
      } 
      if(stemp[(i+2)*10+(j+2)] == 0){
		  stemp[(i+2)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[15] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[15] += w[u*n+v];
		  stemp[(i+2)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho      
	   }
   } else if(i == 0 && j == 1){
       
      if(stemp[(i)*10+(j+2)] == 0){   
		  stemp[(i)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[8] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[8] += w[u*n+v];
		  stemp[(i)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho
      } 
      if(stemp[(i+1)*10+(j+2)] == 0){
		  stemp[(i+1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[10] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[10] += w[u*n+v];
		  stemp[(i+1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho
      }   
      if(stemp[(i+2)*10+(j-1)] == 0){    
		  stemp[(i+2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[12] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[12] += w[u*n+v];
		  stemp[(i+2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho
      } 
      if( stemp[(i+2)*10+(j)] == 0){
		  stemp[(i+2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[13] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[13] += w[u*n+v];
		  stemp[(i+2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
       } 
      if(stemp[(i+2)*10+(j+1)] == 0){
		  stemp[(i+2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[14] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[14] += w[u*n+v];
		  stemp[(i+2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
      } 
      if(stemp[(i+2)*10+(j+2)] == 0){
		  stemp[(i+2)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[15] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[15] += w[u*n+v];
		  stemp[(i+2)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho      
      }        
   } else if(i == 0 && j == 8) {
       
      if(stemp[(i)*10+(j-2)] == 0){   
		  stemp[(i)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[7] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[7] += w[u*n+v];
		  stemp[(i)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      } 
      if(stemp[(i+1)*10+(j-2)] == 0){     
		  stemp[(i+1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[9] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[9] += w[u*n+v];
		  stemp[(i+1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho      
	  }   
      if(stemp[(i+2)*10+(j-2)] == 0){
		  stemp[(i+2)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[11] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[11] += w[u*n+v];
		  stemp[(i+2)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      } 
      if(stemp[(i+2)*10+(j-1)] == 0){ 
		  stemp[(i+2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[12] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[12] += w[u*n+v];
		  stemp[(i+2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
	  }   
      if(stemp[(i+2)*10+(j)] == 0){ 
		  stemp[(i+2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[13] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[13] += w[u*n+v];
		  stemp[(i+2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }  
      if(stemp[(i+2)*10+(j+1)] == 0){
		  stemp[(i+2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[14] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[14] += w[u*n+v];
		  stemp[(i+2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
      }
   } else if(i == 0 && j == 9) {
       
      if( stemp[(i)*10+(j-2)] == 0){   
		  stemp[(i)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[7] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[7] += w[u*n+v];
		  stemp[(i)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
	  }    
      if( stemp[(i+1)*10+(j-2)] == 0){  
		  stemp[(i+1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[9] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[9] += w[u*n+v];
		  stemp[(i+1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho      
      } 
      if( stemp[(i+2)*10+(j-2)] == 0){
		  stemp[(i+2)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[11] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[11] += w[u*n+v];
		  stemp[(i+2)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
      if( stemp[(i+2)*10+(j-1)] == 0){
		  stemp[(i+2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[12] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[12] += w[u*n+v];
		  stemp[(i+2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
      } 
      if(  stemp[(i+2)*10+(j)] == 0){ 
		  stemp[(i+2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[13] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[13] += w[u*n+v];
		  stemp[(i+2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
     }
   } else if(i == 1 && j == 0) {
       
      if( stemp[(i-1)*10+(j+2)] == 0){   
		  stemp[(i-1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[6] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[6] += w[u*n+v];
		  stemp[(i-1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
	  }	   
      if(stemp[(i)*10+(j+2)] == 0){
		  stemp[(i)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[8] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[8] += w[u*n+v];
		  stemp[(i)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho
	  }   
      if( stemp[(i+1)*10+(j+2)] == 0){  
		  stemp[(i+1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[10] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[10] += w[u*n+v];
		  stemp[(i+1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
      if( stemp[(i+2)*10+(j)] == 0){  
		  stemp[(i+2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[13] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[13] += w[u*n+v];
		  stemp[(i+2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }  
      if( stemp[(i+2)*10+(j+1)] == 0){
		  stemp[(i+2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[14] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[14] += w[u*n+v];
		  stemp[(i+2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
      } 
      if( stemp[(i+2)*10+(j+2)] == 0){
		  stemp[(i+2)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[15] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[15] += w[u*n+v];
		  stemp[(i+2)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho      
	  } 
   } else if(i == 1 && j == 1) {
       
      if( stemp[(i-1)*10+(j+2)] == 0){   
		  stemp[(i-1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[6] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[6] += w[u*n+v];
		  stemp[(i-1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      } 
      if( stemp[(i)*10+(j+2)] == 0){
		  stemp[(i)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[8] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[8] += w[u*n+v];
		  stemp[(i)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho
      } 
      if( stemp[(i+1)*10+(j+2)] == 0){ 
		  stemp[(i+1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[10] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[10] += w[u*n+v];
		  stemp[(i+1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
	  }   
      if( stemp[(i+2)*10+(j-1)] == 0){
		  stemp[(i+2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[12] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[12] += w[u*n+v];
		  stemp[(i+2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
	  }   
      if( stemp[(i+2)*10+(j)] == 0){ 
		  stemp[(i+2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[13] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[13] += w[u*n+v];
		  stemp[(i+2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
	  }	   
      if( stemp[(i+2)*10+(j+1)] == 0){ 
		  stemp[(i+2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[14] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[14] += w[u*n+v];
		  stemp[(i+2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
	  }	   
      if( stemp[(i+2)*10+(j+2)] == 0){
		  stemp[(i+2)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[15] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[15] += w[u*n+v];
		  stemp[(i+2)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho    
	  }
   } else if(i == 1 && j == 8) {
       
      if(stemp[(i-1)*10+(j-2)] == 0){   
		  stemp[(i-1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[5] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[5] += w[u*n+v];
		  stemp[(i-1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      } 
      if( stemp[(i)*10+(j-2)] == 0){    
		  stemp[(i)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[7] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[7] += w[u*n+v];
		  stemp[(i)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      } 
      if( stemp[(i+1)*10+(j-2)] == 0){
		  stemp[(i+1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[9] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[9] += w[u*n+v];
		  stemp[(i+1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      } 
      if( stemp[(i+2)*10+(j-2)] == 0){
		  stemp[(i+2)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[11] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[11] += w[u*n+v];
		  stemp[(i+2)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
      if( stemp[(i+2)*10+(j-1)] == 0){
		  stemp[(i+2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[12] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[12] += w[u*n+v];
		  stemp[(i+2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
      if( stemp[(i+2)*10+(j)] == 0){
		  stemp[(i+2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[13] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[13] += w[u*n+v];
		  stemp[(i+2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }  
      if( stemp[(i+2)*10+(j+1)] == 0){
		  stemp[(i+2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[14] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[14] += w[u*n+v];
		  stemp[(i+2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
      }
   } else if(i == 1 && j == 9) {
        
      if(stemp[(i-1)*10+(j-2)] == 0){  
		  stemp[(i-1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[5] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[5] += w[u*n+v];
		  stemp[(i-1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      } 
      if(stemp[(i)*10+(j-2)] == 0){    
		  stemp[(i)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[7] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[7] += w[u*n+v];
		  stemp[(i)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      } 
      if( stemp[(i+1)*10+(j-2)] == 0){
		  stemp[(i+1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[9] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[9] += w[u*n+v];
		  stemp[(i+1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      } 
      if( stemp[(i+2)*10+(j-2)] == 0){
		  stemp[(i+2)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[11] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[11] += w[u*n+v];
		  stemp[(i+2)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
      if( stemp[(i+2)*10+(j-1)] == 0){
		  stemp[(i+2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[12] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[12] += w[u*n+v];
		  stemp[(i+2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
      if( stemp[(i+2)*10+(j)] == 0){
		  stemp[(i+2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[13] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[13] += w[u*n+v];
		  stemp[(i+2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }
   } else if(i > 1 && i < 8 && j == 0) {
       
      if(stemp[(i-2)*10+(j)] == 0){   
		  stemp[(i-2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[2] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[2] += w[u*n+v];
		  stemp[(i-2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }
      		   
      if(stemp[(i-2)*10+(j+1)] == 0){ 
		  stemp[(i-2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[3] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[3] += w[u*n+v];
		  stemp[(i-2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(stemp[(i-2)*10+(j+2)] == 0){   
		  stemp[(i-2)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[4] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[4] += w[u*n+v];
		  stemp[(i-2)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-1)*10+(j+2)] == 0){    
		  stemp[(i-1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[6] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[6] += w[u*n+v];
		  stemp[(i-1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i)*10+(j+2)] == 0){   
		  stemp[(i)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[8] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[8] += w[u*n+v];
		  stemp[(i)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho
       }  
        
      if( stemp[(i+1)*10+(j+2)] == 0){    
		  stemp[(i+1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[10] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[10] += w[u*n+v];
		  stemp[(i+1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+2)*10+(j)] == 0){   
		  stemp[(i+2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[13] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[13] += w[u*n+v];
		  stemp[(i+2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+2)*10+(j+1)] == 0){    
		  stemp[(i+2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[14] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[14] += w[u*n+v];
		  stemp[(i+2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+2)*10+(j+2)] == 0){   
		  stemp[(i+2)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[15] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[15] += w[u*n+v];
		  stemp[(i+2)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho    
      }
   } else if(i > 1 && i < 8 && j == 1) { 
        
      if( stemp[(i-2)*10+(j-1)] == 0){     
		  stemp[(i-2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[1] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[1] += w[u*n+v];
		  stemp[(i-2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho
      }  
        
      if( stemp[(i-2)*10+(j)] == 0){         
		  stemp[(i-2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[2] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[2] += w[u*n+v];
		  stemp[(i-2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-2)*10+(j+1)] == 0){    
		  stemp[(i-2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[3] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[3] += w[u*n+v];
		  stemp[(i-2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
       }  
        
      if( stemp[(i-2)*10+(j+2)] == 0){      
		  stemp[(i-2)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[4] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[4] += w[u*n+v];
		  stemp[(i-2)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-1)*10+(j+2)] == 0){    
		  stemp[(i-1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[6] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[6] += w[u*n+v];
		  stemp[(i-1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i)*10+(j+2)] == 0){   
		  stemp[(i)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[8] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[8] += w[u*n+v];
		  stemp[(i)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho
      }  
        
      if( stemp[(i+1)*10+(j+2)] == 0){     
		  stemp[(i+1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[10] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[10] += w[u*n+v];
		  stemp[(i+1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+2)*10+(j-1)] == 0){   
		  stemp[(i+2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[12] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[12] += w[u*n+v];
		  stemp[(i+2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+2)*10+(j)] == 0){    
          stemp[(i+2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[13] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[13] += w[u*n+v];
		  stemp[(i+2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(stemp[(i+2)*10+(j+1)] == 0){    
		  stemp[(i+2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[14] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[14] += w[u*n+v];
		  stemp[(i+2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(stemp[(i+2)*10+(j+2)] == 0){   
		  stemp[(i+2)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[15] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[15] += w[u*n+v];
		  stemp[(i+2)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinh
	  }
   
   } else if(i > 1 && i < 8 && j == 8) {
        
      if( stemp[(i-2)*10+(j-2)] == 0){   
          stemp[(i-2)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[0] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[0] += w[u*n+v];
		  stemp[(i-2)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(   stemp[(i-2)*10+(j-1)] == 0){      
		  stemp[(i-2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[1] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[1] += w[u*n+v];
		  stemp[(i-2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-2)*10+(j)] == 0){      
		  stemp[(i-2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[2] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[2] += w[u*n+v];
		  stemp[(i-2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(stemp[(i-2)*10+(j+1)] == 0){    
		  stemp[(i-2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[3] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[3] += w[u*n+v];
		  stemp[(i-2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-1)*10+(j-2)] == 0){   
		  stemp[(i-1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[5] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[5] += w[u*n+v];
		  stemp[(i-1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i)*10+(j-2)] == 0){       
		  stemp[(i)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[7] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[7] += w[u*n+v];
		  stemp[(i)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(  stemp[(i+1)*10+(j-2)] == 0){   
		  stemp[(i+1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[9] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[9] += w[u*n+v];
		  stemp[(i+1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+2)*10+(j-2)]  == 0){   
		  stemp[(i+2)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[11] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[11] += w[u*n+v];
		  stemp[(i+2)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
       }  
        
	   if( stemp[(i+2)*10+(j-1)] == 0){   
		  stemp[(i+2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[12] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[12] += w[u*n+v];
		  stemp[(i+2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
	  if(stemp[(i+2)*10+(j)] == 0){    
		  stemp[(i+2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[13] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[13] += w[u*n+v];
		  stemp[(i+2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
       }  
        
	   if(stemp[(i+2)*10+(j+1)] == 0){   
		  stemp[(i+2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[14] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[14] += w[u*n+v];
		  stemp[(i+2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
       }
   } else if(i > 1 && i < 8 && j == 9) {
        
	  if(stemp[(i-2)*10+(j-2)] == 0){      
		  stemp[(i-2)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[0] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[0] += w[u*n+v];
		  stemp[(i-2)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
       }  
        
      if( stemp[(i-2)*10+(j-1)] == 0){     
		  stemp[(i-2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[1] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[1] += w[u*n+v];
		  stemp[(i-2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-2)*10+(j)] == 0){      
		  stemp[(i-2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[2] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[2] += w[u*n+v];
		  stemp[(i-2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-1)*10+(j-2)] == 0){   
		  stemp[(i-1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[5] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[5] += w[u*n+v];
		  stemp[(i-1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i)*10+(j-2)] == 0){       
		  stemp[(i)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[7] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[7] += w[u*n+v];
		  stemp[(i)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+1)*10+(j-2)] == 0){   
		  stemp[(i+1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[9] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[9] += w[u*n+v];
		  stemp[(i+1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+2)*10+(j-2)] == 0){   
		  stemp[(i+2)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[11] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[11] += w[u*n+v];
		  stemp[(i+2)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
       }  
        
      if( stemp[(i+2)*10+(j-1)] == 0){   
		  stemp[(i+2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[12] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[12] += w[u*n+v];
		  stemp[(i+2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+2)*10+(j)] == 0){    
		  stemp[(i+2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[13] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[13] += w[u*n+v];
		  stemp[(i+2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }   
   } else if(i == 8 && j == 0) {  
        
      if( stemp[(i-2)*10+(j)] == 0){      
		  stemp[(i-2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[2] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[2] += w[u*n+v];
		  stemp[(i-2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(  stemp[(i-2)*10+(j+1)] == 0){    
		  stemp[(i-2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[3] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[3] += w[u*n+v];
		  stemp[(i-2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(  stemp[(i-2)*10+(j+2)] == 0){       
		  stemp[(i-2)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[4] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[4] += w[u*n+v];
		  stemp[(i-2)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-1)*10+(j+2)] == 0){   
		  stemp[(i-1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[6] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[6] += w[u*n+v];
		  stemp[(i-1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i)*10+(j+2)] == 0){   
		  stemp[(i)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[8] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[8] += w[u*n+v];
		  stemp[(i)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho
      }  
        
      if( stemp[(i+1)*10+(j+2)] == 0){     
		  stemp[(i+1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[10] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[10] += w[u*n+v];
		  stemp[(i+1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }
   } else if(i == 8 && j == 1) {
        
      if(stemp[(i-2)*10+(j-1)] == 0){     
		  stemp[(i-2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[1] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[1] += w[u*n+v];
		  stemp[(i-2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
       }  
        
      if(stemp[(i-2)*10+(j)] == 0){     
		  stemp[(i-2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[2] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[2] += w[u*n+v];
		  stemp[(i-2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
       }  
        
      if(stemp[(i-2)*10+(j+1)] == 0){    
		  stemp[(i-2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[3] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[3] += w[u*n+v];
		  stemp[(i-2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(stemp[(i-2)*10+(j+2)] == 0){       
		  stemp[(i-2)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[4] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[4] += w[u*n+v];
		  stemp[(i-2)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(stemp[(i-1)*10+(j+2)] == 0){      
		  stemp[(i-1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[6] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[6] += w[u*n+v];
		  stemp[(i-1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i)*10+(j+2)] == 0){   
		  stemp[(i)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[8] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[8] += w[u*n+v];
		  stemp[(i)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho
      }  
        
      if( stemp[(i+1)*10+(j+2)]== 0){     
		  stemp[(i+1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[10] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[10] += w[u*n+v];
		  stemp[(i+1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }
   } else if(i == 8 && j == 8) {   
        
      if( stemp[(i-2)*10+(j-2)] == 0){      
		  stemp[(i-2)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[0] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[0] += w[u*n+v];
		  stemp[(i-2)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
       }  
        
      if( stemp[(i-2)*10+(j-1)] == 0){     
		  stemp[(i-2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[1] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[1] += w[u*n+v];
		  stemp[(i-2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
       }  
        
      if( stemp[(i-2)*10+(j)] == 0){     
		  stemp[(i-2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[2] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[2] += w[u*n+v];
		  stemp[(i-2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
       }  
        
      if( stemp[(i-2)*10+(j+1)] == 0){    
		  stemp[(i-2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[3] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[3] += w[u*n+v];
		  stemp[(i-2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-1)*10+(j-2)] == 0){     
		  stemp[(i-1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[5] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[5] += w[u*n+v];
		  stemp[(i-1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i)*10+(j-2)] == 0){       
          stemp[(i)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[7] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[7] += w[u*n+v];
		  stemp[(i)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+1)*10+(j-2)] == 0){   
		  stemp[(i+1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[9] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[9] += w[u*n+v];
		  stemp[(i+1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho     
      }  
           
   } else if(i == 8 && j == 9) { 
        
      if( stemp[(i-2)*10+(j-2)] == 0){      
		  stemp[(i-2)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[0] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[0] += w[u*n+v];
		  stemp[(i-2)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
       }  
        
      if( stemp[(i-2)*10+(j-1)] == 0){     
		  stemp[(i-2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[1] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[1] += w[u*n+v];
		  stemp[(i-2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-2)*10+(j)] == 0){      
		  stemp[(i-2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[2] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[2] += w[u*n+v];
		  stemp[(i-2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
       }  
        
      if( stemp[(i-1)*10+(j-2)] == 0){     
		  stemp[(i-1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[5] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[5] += w[u*n+v];
		  stemp[(i-1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i)*10+(j-2)] == 0){       
		  stemp[(i)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[7] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[7] += w[u*n+v];
		  stemp[(i)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+1)*10+(j-2)] == 0){   
		  stemp[(i+1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[9] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[9] += w[u*n+v];
		  stemp[(i+1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho     
      }  
         
   } else if(i == 9 && j == 0) {
        
      if( stemp[(i-2)*10+(j)] == 0){   
		  stemp[(i-2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[2] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[2] += w[u*n+v];
		  stemp[(i-2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-2)*10+(j+1)] == 0){     
      	  stemp[(i-2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[3] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[3] += w[u*n+v];
		  stemp[(i-2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
       }  
        
      if( stemp[(i-2)*10+(j+2)] == 0){    
		  stemp[(i-2)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[4] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[4] += w[u*n+v];
		  stemp[(i-2)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-1)*10+(j+2)] == 0){   
		  stemp[(i-1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[6] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[6] += w[u*n+v];
		  stemp[(i-1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i)*10+(j+2)] == 0){   
		  stemp[(i)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[8] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[8] += w[u*n+v];
		  stemp[(i)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho
      }  
        
   
   } else if(i == 9 && j == 1) {
        
      if(stemp[(i-2)*10+(j-1)] == 0){    
		  stemp[(i-2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[1] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[1] += w[u*n+v];
		  stemp[(i-2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-2)*10+(j)] == 0){   
		  stemp[(i-2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[2] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[2] += w[u*n+v];
		  stemp[(i-2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-2)*10+(j+1)] == 0){     
		  stemp[(i-2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[3] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[3] += w[u*n+v];
		  stemp[(i-2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-2)*10+(j+2)] == 0){     
		  stemp[(i-2)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[4] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[4] += w[u*n+v];
		  stemp[(i-2)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-1)*10+(j+2)] == 0){   
		  stemp[(i-1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[6] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[6] += w[u*n+v];
		  stemp[(i-1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i)*10+(j+2)]  == 0){   
		  stemp[(i)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[8] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[8] += w[u*n+v];
		  stemp[(i)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho
      }  
          
   } else if(i == 9 && j == 8) { 
        
      if(stemp[(i-2)*10+(j-2)] == 0){    
		  stemp[(i-2)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[0] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[0] += w[u*n+v];
		  stemp[(i-2)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(  stemp[(i-2)*10+(j-1)] == 0){   
		  stemp[(i-2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[1] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[1] += w[u*n+v];
		  stemp[(i-2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(stemp[(i-2)*10+(j)] == 0){   
		  stemp[(i-2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[2] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[2] += w[u*n+v];
		  stemp[(i-2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-2)*10+(j+1)]== 0){     
		  stemp[(i-2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[3] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[3] += w[u*n+v];
		  stemp[(i-2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
       }  
        
      if( stemp[(i-1)*10+(j-2)] == 0){   
		  stemp[(i-1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[5] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[5] += w[u*n+v];
		  stemp[(i-1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i)*10+(j-2)] == 0){       
		  stemp[(i)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[7] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[7] += w[u*n+v];
		  stemp[(i)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho   
      }  
										     
   } else if(i == 9 && j == 9) {				
        
      if(stemp[(i-2)*10+(j-2)] == 0){   
		  stemp[(i-2)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[0] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[0] += w[u*n+v];
		  stemp[(i-2)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-2)*10+(j-1)] == 0){   
		  stemp[(i-2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[1] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[1] += w[u*n+v];
		  stemp[(i-2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-2)*10+(j)] == 0){   
		  stemp[(i-2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[2] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[2] += w[u*n+v];
		  stemp[(i-2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-1)*10+(j-2)] == 0){    
		  stemp[(i-1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[5] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[5] += w[u*n+v];
		  stemp[(i-1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i)*10+(j-2)] == 0){       
		  stemp[(i)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[7] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[7] += w[u*n+v];
		  stemp[(i)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho   
      }  
           
   } else if(i == 0 && j > 1 && j < 8) { 
        
      if(stemp[(i)*10+(j-2)] == 0){   
		  stemp[(i)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[7] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[7] += w[u*n+v];
		  stemp[(i)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(stemp[(i)*10+(j+2)] == 0){    
		  stemp[(i)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[8] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[8] += w[u*n+v];
		  stemp[(i)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+1)*10+(j-2)] == 0){   
      	  stemp[(i+1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[9] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[9] += w[u*n+v];
		  stemp[(i+1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+1)*10+(j+2)] == 0){   
		  stemp[(i+1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[10] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[10] += w[u*n+v];
		  stemp[(i+1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+2)*10+(j-2)] == 0){         
		  stemp[(i+2)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[11] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[11] += w[u*n+v];
		  stemp[(i+2)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+2)*10+(j-1)] == 0){    
		  stemp[(i+2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[12] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[12] += w[u*n+v];
		  stemp[(i+2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
       }  
        
      if( stemp[(i+2)*10+(j)] == 0){   
		  stemp[(i+2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[13] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[13] += w[u*n+v];
		  stemp[(i+2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+2)*10+(j+1)] == 0){    
		  stemp[(i+2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[14] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[14] += w[u*n+v];
		  stemp[(i+2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+2)*10+(j+2)] == 0){   
		  stemp[(i+2)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[15] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[15] += w[u*n+v];
		  stemp[(i+2)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
          
   } else if(i == 1 && j > 1 && j < 8) {
        
      if(stemp[(i-1)*10+(j-2)] == 0){   
		  stemp[(i-1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[5] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[5] += w[u*n+v];
		  stemp[(i-1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(stemp[(i-1)*10+(j+2)] == 0){   
		  stemp[(i-1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[6] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[6] += w[u*n+v];
		  stemp[(i-1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(stemp[(i)*10+(j-2)] == 0){   
		  stemp[(i)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[7] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[7] += w[u*n+v];
		  stemp[(i)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(stemp[(i)*10+(j+2)] == 0){    
		  stemp[(i)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[8] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[8] += w[u*n+v];
		  stemp[(i)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(stemp[(i+1)*10+(j-2)] == 0){   
		  stemp[(i+1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[9] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[9] += w[u*n+v];
		  stemp[(i+1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(stemp[(i+1)*10+(j+2)] == 0){   
		  stemp[(i+1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[10] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[10] += w[u*n+v];
		  stemp[(i+1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(stemp[(i+2)*10+(j-2)] == 0){        
		  stemp[(i+2)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[11] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[11] += w[u*n+v];
		  stemp[(i+2)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
       }  
        
      if(stemp[(i+2)*10+(j-1)] == 0){   
		  stemp[(i+2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[12] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[12] += w[u*n+v];
		  stemp[(i+2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+2)*10+(j)] == 0){    
		  stemp[(i+2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[13] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[13] += w[u*n+v];
		  stemp[(i+2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+2)*10+(j+1)] == 0){    
		  stemp[(i+2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[14] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[14] += w[u*n+v];
		  stemp[(i+2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(stemp[(i+2)*10+(j+2)] == 0){   
		  stemp[(i+2)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[15] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[15] += w[u*n+v];
		  stemp[(i+2)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
           
   } else if(i == 8 && j > 1 && j < 8) { 
        
      if(stemp[(i-2)*10+(j-2)] == 0){    
		  stemp[(i-2)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[0] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[0] += w[u*n+v];
		  stemp[(i-2)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(stemp[(i-2)*10+(j-1)] == 0){      
		  stemp[(i-2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[1] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[1] += w[u*n+v];
		  stemp[(i-2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
       }  
        
      if( stemp[(i-2)*10+(j)] == 0){     
		  stemp[(i-2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[2] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[2] += w[u*n+v];
		  stemp[(i-2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-2)*10+(j+1)] == 0){    
		  stemp[(i-2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[3] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[3] += w[u*n+v];
		  stemp[(i-2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-2)*10+(j+2)] == 0){       
		  stemp[(i-2)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[4] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[4] += w[u*n+v];
		  stemp[(i-2)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(stemp[(i-1)*10+(j-2)] == 0){       
		  stemp[(i-1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[5] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[5] += w[u*n+v];
		  stemp[(i-1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
       }  
        
      if( stemp[(i-1)*10+(j+2)] == 0){      
      stemp[(i-1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[6] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[6] += w[u*n+v];
		  stemp[(i-1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if(stemp[(i)*10+(j-2)] == 0){      
		  stemp[(i)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[7] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[7] += w[u*n+v];
		  stemp[(i)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i)*10+(j+2)] == 0){    
		  stemp[(i)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[8] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[8] += w[u*n+v];
		  stemp[(i)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+1)*10+(j-2)] == 0){   
		  stemp[(i+1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[9] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[9] += w[u*n+v];
		  stemp[(i+1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i+1)*10+(j+2)] == 0){   
		  stemp[(i+1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[10] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[10] += w[u*n+v];
		  stemp[(i+1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        	     
   } else if(i == 9 && j > 1 && j < 8) {			
        
      if( stemp[(i-2)*10+(j-2)] == 0){   
		  stemp[(i-2)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[0] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[0] += w[u*n+v];
		  stemp[(i-2)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-2)*10+(j-1)] == 0){      
		  stemp[(i-2)*10+(j-1)] = 1; //faz o movimento 1° vizinho
		  foV[1] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[1] += w[u*n+v];
		  stemp[(i-2)*10+(j-1)] = 0; //Desfaz o movimento 1° vizinho  
       }  
        
      if( stemp[(i-2)*10+(j)] == 0){     
		  stemp[(i-2)*10+(j)] = 1; //faz o movimento 1° vizinho
		  foV[2] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[2] += w[u*n+v];
		  stemp[(i-2)*10+(j)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-2)*10+(j+1)] == 0){    
		  stemp[(i-2)*10+(j+1)] = 1; //faz o movimento 1° vizinho
		  foV[3] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[3] += w[u*n+v];
		  stemp[(i-2)*10+(j+1)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-2)*10+(j+2)] == 0){       
		  stemp[(i-2)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[4] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[4] += w[u*n+v];
		  stemp[(i-2)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-1)*10+(j-2)] == 0){       
		  stemp[(i-1)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[5] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[5] += w[u*n+v];
		  stemp[(i-1)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i-1)*10+(j+2)] == 0){       
		  stemp[(i-1)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[6] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[6] += w[u*n+v];
		  stemp[(i-1)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i)*10+(j-2)] == 0){       
		  stemp[(i)*10+(j-2)] = 1; //faz o movimento 1° vizinho
		  foV[7] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[7] += w[u*n+v];
		  stemp[(i)*10+(j-2)] = 0; //Desfaz o movimento 1° vizinho  
      }  
        
      if( stemp[(i)*10+(j+2)] == 0){    
		  stemp[(i)*10+(j+2)] = 1; //faz o movimento 1° vizinho
		  foV[8] = 0;
		  for(int u = 0; u < n; ++u) for(int v = 0; v < n; ++v) if(u > v)
		      if(stemp[u] == 1 && stemp[v] == 1) foV[8] += w[u*n+v];
		  stemp[(i)*10+(j+2)] = 0; //Desfaz o movimento 1° vizinho
      }   
   }

//verifica o melhor vizinho
   for(int u = 0; u < 16; ++u){
       if(foV[u] < fo) {
           fo = foV[u];
           mviz = u;
       }
   }   
   
   /**/

//Retorna s2l com o melhor vizinho
   s2l[(i)*10+(j)] = 0;
   if(mviz == 0){
      s2l[(i-2)*10+(j-2)] = 1;
   } else if(mviz == 1){
      s2l[(i-2)*10+(j-1)] = 1;    
   } else if(mviz == 2){
      s2l[(i-2)*10+(j)] = 1;   
   } else if(mviz == 3){
      s2l[(i-2)*10+(j+1)] = 1;   
   } else if(mviz == 4){
      s2l[(i-2)*10+(j+2)] = 1;   
   } else if(mviz == 5){
      s2l[(i-1)*10+(j-2)] = 1;   
   } else if(mviz == 6){
      s2l[(i-1)*10+(j+2)] = 1;  
   } else if(mviz == 7){
      s2l[(i)*10+(j-2)] = 1;   
   } else if(mviz == 8){
      s2l[(i)*10+(j+2)] = 1;   
   } else if(mviz == 9){
      s2l[(i+1)*10+(j-2)] = 1;   
   } else if(mviz == 10){
      s2l[(i+1)*10+(j+2)] = 1;   
   } else if(mviz == 11){
      s2l[(i+2)*10+(j-2)] = 1;   
   } else if(mviz == 12){
      s2l[(i+2)*10+(j-1)] = 1;  
   } else if(mviz == 13){
      s2l[(i+2)*10+(j)] = 1;  
   } else if(mviz == 14){
      s2l[(i+2)*10+(j+1)] = 1;  
   } else if(mviz == 15){
      s2l[(i+2)*10+(j+2)] = 1;     
   } else if(mviz == 16){
      s2l[(i)*10+(j)] = 1; //caso não encontre um vizinho     
   }
/**/

   
   return fo;
}

//############################################################################
/***************************************************************************/
double  WFLOP::calcFO(VOL_ivector s)
{
   double fo = 0;
   for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) 
       if(i > j) if(s[i] == 1 && s[j] == 1) fo += w[i*n+j];  

   return fo;
}
/***************************************************************************/
//############################################################################
//############################################################################
/***************************************************************************/
double WFLOP::vnd(int n, VOL_ivector s, double fo, VOL_ivector& sl)
{
    VOL_ivector s2l(n);
    
    for (int i = 0; i < n; i++) sl[i] = s[i];
    //for(int i = 0; i < n; ++i) if(sl[i]) printf("\nsl[%d] %d",i, sl[i]); getchar();
    double fos2l = fo;
    int p = 0;
    int km = 1;

    while(km <= 3) {
       for (int i = 0; i < n; ++i) s2l[i] = sl[i];
       for (int j = 0; j < km; ++j) fos2l = realocaTurbinas(s2l, n, fo);
       if (fos2l < fo) {
          for (int i = 0; i < n; i++) sl[i] = s2l[i];
          fo = fos2l;
          km = 1;
       }  else km++;
    }//while 1

    //printf("-------------------------------------------------\n");
    //printf("Melhor solucao encontrada pelo VND:\n");
    //printf("Funcao objetivo: %d\n", fo);
    //printf("-------------------------------------------------\n");

    return fo;
}//vnd

/***************************************************************************/
//############################################################################
//############################################################################
/***************************************************************************/
double  WFLOP::realocaTurbinas(VOL_ivector& s2l, int n, double fo)
{
   double fos2l ;
   int p1, p ;
   
   do{
       p = randint(0,99);
       //printf("\ns2l[%d] %d ",p,s2l[p]); 
   }while ( s2l[p] == 0 );
   //printf("\n Fim do-while"); getchar();
   
   if(p == 99){
        s2l[p] = 0;
   	s2l[p-11] = 1;
   	fos2l = 0;
   	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) 
   	    if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
   	if(fos2l < fo) return fos2l;
   
	s2l[p-11] = 0;
	s2l[p-10] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) 
	    if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;
	       
	s2l[p-10] = 0;
	s2l[p-9] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;

	s2l[p-9] = 0;
	s2l[p-1] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;
	
   } else if(p >= 91 && p <= 98){
        s2l[p] = 0;
   	s2l[p-11] = 1;
   	fos2l = 0;
   	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
   	if(fos2l < fo) return fos2l;
   
	s2l[p-11] = 0;
	s2l[p-10] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;
	       
	s2l[p-10] = 0;
	s2l[p-9] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;

	s2l[p-9] = 0;
	s2l[p-1] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;

	s2l[p-1] = 0;
	s2l[p+1] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;
   
   } else if(p == 90){
        s2l[p-10] = 1;
	s2l[p] = 0;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;

	s2l[p-10] = 0;
	s2l[p-9] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;
	
	s2l[p-9] = 0;
	s2l[p+1] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;
	
   }else if(p == 10){
        s2l[p-10] = 1;
	s2l[p] = 0;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;

	s2l[p-10] = 0;
	s2l[p-9] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;
	
	s2l[p-9] = 0;
	s2l[p+1] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;
	
        s2l[p+1] = 0;
	s2l[p+10] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;
	
	s2l[p+10] = 0;
	s2l[p+11] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;
	
   } else if(p == 9){
        s2l[p-1] = 1;
	s2l[p] = 0;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;

	s2l[p-1] = 0;
	s2l[p+9] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;

	s2l[p+9] = 0;
	s2l[p+10] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;
   
   }else if(p >= 1 && p <= 8 ){
        s2l[p-1] = 1;
	s2l[p] = 0;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;
	
        s2l[p-1] = 0;
	s2l[p+1] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;

	s2l[p+1] = 0;
	s2l[p+9] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;

	s2l[p+9] = 0;
	s2l[p+10] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;
	   
	s2l[p+10] = 0;
	s2l[p+11] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;
	
   } else if(p == 0){
        s2l[p] = 0;
	s2l[p+1] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;
	
	s2l[p+1] = 0;
	s2l[p+9] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;

	s2l[p+9] = 0;
	s2l[p+10] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;
	   
	s2l[p+10] = 0;
	s2l[p+11] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;
	
   } else if(p >= 11 && p <= 89){
   	s2l[p] = 0;
   	s2l[p-11] = 1;
   	fos2l = 0;
   	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
   	if(fos2l < fo) return fos2l;
   
	s2l[p-11] = 0;
	s2l[p-10] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;
	       
	s2l[p-10] = 0;
	s2l[p-9] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;

	s2l[p-9] = 0;
	s2l[p-1] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;

	s2l[p-1] = 0;
	s2l[p+1] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;

	s2l[p+1] = 0;
	s2l[p+9] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;

	s2l[p+9] = 0;
	s2l[p+10] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;
	   
	s2l[p+10] = 0;
	s2l[p+11] = 1;
	fos2l = 0;
	for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(s2l[i] == 1 && s2l[j] == 1) fos2l += w[i*n+j];
	if(fos2l < fo) return fos2l;
   }
   
   s2l[p] = 1;
   
   return fo;
}
/***************************************************************************/
//############################################################################
//############################################################################
/***************************************************************************/
/*
   //printf("\n heuristica"); 
   double value=0;
   VOL_ivector sol(n);
   sol = 0;
   int ct = 0;
   VOL_ivector sl(n);

//constroi solução inteira a partir da solução relaxada 
   ct = 0;
   for(int i = 0; i < n; ++i) if(x[n*n+i] == 1){
      sol[i] = 1; 
      ++ct;
   }
   
   int pos = -1;
   for(int j = 0; j < (k - ct); ++j){
      double maior = 0;
      for(int i = 0; i < n; ++i) if(x[n*n+i] != 1 && sol[i] != 1) if(x[n*n+i] > maior){
         maior = x[n*n+i];
         pos = i;
      }
      sol[pos] = 1;
   }
   //for(int i = 0; i < n; ++i) printf("\n sol1[%d] = %d ",i,sol[i]); getchar();
                        
   value = 0;
   for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(i > j) if(sol[n*n+i] == 1 && sol[n*n+j] == 1) value += w[i*n+j];

   new_icost = value;
   if (value < icost) {
      icost = value;
      ix = sol;
   }
   if(value< bestCost){
       bestCost = value;
       for(int i = 0; i < n; ++i) bSol[i] = sol[i]  ;
   }

///Agora o ils
   //param ILS a calibrar
   int ILSmax    = 10; //# iter sem melhora ILS
   int vezesnivel = 3; // nro de vezes sem melhora em um dado nivel
   //
   int p1, p2;
   int nivel, iter, MelhorIter;
   float fo, fol;

   fo = vnd(n, sol,  value, sl);

   iter = MelhorIter = 0;
   nivel = 1;
   while (iter - MelhorIter < ILSmax){
       iter++;
       // s_2l <- s
       for(int i = 0; i < n; ++i) sl[i] = sol[i];
       int vezes = 0;
       while (vezes < vezesnivel){
           int ntrocasmax = nivel + 1;
           int ntrocas = 0;
           for(int i = 0; i < n; ++i) sl[i] = sol[i];
           //for(int i = 0; i < n; ++i) printf("\n s_2l[%d] = %d ",i,s_2l[i]); getchar();
           //for(int i = 0; i < n; ++i) printf("\n sol[%d] = %d ",i,sol[i]); getchar();
           fol = value;
           while (ntrocas < ntrocasmax){
               ntrocas++;
               p1 = randint(0,n-1);
               do{
                   p2 = randint(0,n-1);
               }while ( (p1 == p2) || (sl[p1]-sl[p2]==0) );
               //printf("\n yi= %d = %d || yj= %d = %d",p1,s_2l[p1],p2,s_2l[p2]); getchar();
               int aux;
               aux = sl[p1];
               sl[p1] = sl[p2];
               sl[p2] = aux;
           } // fim while ntrocasmax
           fol =  vnd(n, sol,  value, sl);
           if (fol < value){
               value = fol;
               for(int i = 0; i < n; ++i) sol[i] = sl[i]  ;
               vezes = 0;
               nivel = 1;
               MelhorIter = iter;
               //printf("fo = %12.4f \n", value);
               if(value< bestCost){
                   bestCost = value;
                   for(int i = 0; i < n; ++i) bSol[i] = sol[i]  ;
               }
           }
           vezes++;
       } // fim while vezesnivel
       nivel++;
       //printf("Aumentando o nivel perturbacao para %2d \n",nivel);
   } // fim while

   //printf("int sol %f\n", new_icost);


/**/
//####################################################################
//####################################################################
 
