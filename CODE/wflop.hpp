// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef __WFLOP_HPP__
#define __WFLOP_HPP__

#include <cstdio>
#include <cfloat>
#include <string>
#include <vector>
#include <cstdio>
#include <iomanip>
#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <sstream>
#include <limits>
#include <list>

#include "wflopUTIL.h"
#include <VolVolume.hpp>

//#include <cplex.h>

//#include <gurobi_c++.h>


//extern "C" {
//	#include <concorde.h>
//        #include <gurobi_c.h>
//}




// parameters controlled by the user
class WFLOP_parms {
public:
  std::string fdata; // file with the data
  std::string dualfile; // file with an initial dual solution
  std::string dual_savefile; // file to save final dual solution
  std::string int_savefile; // file to save primal integer solution
  int h_iter; // number of times that the primal heuristic will be run
  // after termination of the volume algorithm
   
  WFLOP_parms(const char* filename);
  ~WFLOP_parms() {}
};



class WFLOP : public VOL_user_hooks {

public:
  // for all hooks: return value of -1 means that volume should quit
  // compute reduced costs
   int compute_rc(const VOL_dvector& u, VOL_dvector& rc);
  // solve relaxed problem
   int solve_subproblem(const VOL_dvector& u, const VOL_dvector& rc,
			double& lcost, VOL_dvector& x, VOL_dvector&v,
			double& pcost);
  // primal heuristic
  // return DBL_MAX in heur_val if feas sol wasn't/was found 
   int heuristics(const VOL_problem& p, const VOL_dvector& x, double& heur_val);
   double realocaTurbinas(VOL_ivector& s2l, int n, double fo);
   double vnd(int n, VOL_ivector s, double fo, VOL_ivector& sl);
   double vndV2(int n, VOL_ivector s, double fo, VOL_ivector& sl);
   double mov1Grid(VOL_ivector& s2l, int n);
   double mov2Grid(VOL_ivector& s2l, int n);
   double calcFO(VOL_ivector s);
   double gera_k_vizinhos_aleatorios(int n, VOL_ivector& s, int ki);
   double gera_um_vizinho_aleatorio(int n, VOL_ivector& s);

public: 
  VOL_dvector w;     // wake efect
  VOL_dvector ws;     
  VOL_dvector p;     
  VOL_ivector ix;   // best integer feasible solution so far
  //VOL_ivector hash;


  int n;             // number of nodes
  int k;
  int wd;
  int caso;
  char *fname;
  double alpha;
  double r0;
  double a;
  double thetainit;
  int nws;
  //cpu
  double t0;
  //double t;
  double      icost;  // value of best integer feasible solution 
  double      dcost;  // value of the current dual solution 
  int iter;

 
public:
  WFLOP() : icost(DBL_MAX) {}
  virtual ~WFLOP() {}  
};

//#############################################################################
//########  Member functions          #########################################
//#############################################################################

//****** WFLOP_parms
// reading parameters specific to facility location
WFLOP_parms::WFLOP_parms(const char *filename) :
   fdata(""), 
   h_iter(10)
{

   char s[500];
   FILE * file = fopen(filename, "r");
   if (!file) { 
      printf("Failure to open WFLOP datafile: %s\n ", filename);
      abort();
   }
   
   while (fgets(s, 500, file)) {
      const int len = strlen(s) - 1;
      if (s[len] == '\n')
	 s[len] = 0;
      std::string ss;
      ss = s;
      
      if (ss.find("fdata") == 0) {
	 int j = ss.find("=");
	 int j1 = ss.length() - j + 1;
	 fdata = ss.substr(j+1, j1);
	 
      } else if (ss.find("dualfile") == 0) {
	 int j = ss.find("=");
	 int j1 = ss.length() - j + 1;
	 dualfile = ss.substr(j+1, j1);
	 
      } else if (ss.find("dual_savefile") == 0) {
	 int j = ss.find("=");
	 int j1 = ss.length() - j + 1;
	 dual_savefile = ss.substr(j+1, j1);

      } else if (ss.find("int_savefile") == 0) {
	 int j = ss.find("=");
	 int j1 = ss.length() - j + 1;
	 int_savefile = ss.substr(j+1, j1);
	 
      } else if (ss.find("h_iter") == 0) {
	 int i = ss.find("=");  
	 h_iter = atoi(&s[i+1]);
      }	
   }
   fclose(file);
   
}

#endif


 


