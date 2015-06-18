#ifndef __RATGRPH_H
#define __RATGRPH_H

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>

#include <vector>
#include <queue>
#include <string>
#include <map>
#include <set>
#include <algorithm>

extern "C" {
#include "findpath.h"
}

#include "RNAutils.h"
#include "RNAstruc.h"
#include "BHGbuilder.h"

using namespace std;

struct pq_rate {
  int conns;
  int lm;
  int energy;
  double inrate;
  double outrate;

  pq_rate(int c, int l, int e, double inr=0.0, double outr=0.0) {
    conns = c;
    lm = l;
    energy = e;
    inrate = inr;
    outrate = outr;
  }

  bool operator<(const pq_rate &scnd) const {
    if (conns == scnd.conns) {
      return energy < scnd.energy;
    } else return conns > scnd.conns;
  }
};

class pq_comp
{
  char type;  // R - rates, C - connections
public:
  pq_comp(const char& par='C') {
    type = par;
  }

  bool operator() (const pq_rate& lhs, const pq_rate&rhs) const {
    switch (type) {
    case 'R': if (lhs.outrate/lhs.inrate != rhs.outrate/rhs.inrate) {
      return lhs.outrate/lhs.inrate < rhs.outrate/rhs.inrate;
    }
    case 'C':
    default:
      return lhs < rhs;
    }
  }
};


void SwapMinima(int dim_rates, double *rate_matrix, int a, int b);
void MxRShorten(double **shorten, int fulldim, int gdim);

class RateGraph
{
  // rates
  vector < map<int, double> > rates;  // stored in a sparse matrix format due to space.

  // removed:
  set <int> removed;
  map <int, int> lm_to_pos;
  vector<int> pos_to_lm;

  // int
  int edge_count;

  // priority queue for removal
  priority_queue<pq_rate>  to_rem;

  // rate_matrix for shur removal (temporary)
  double *rate_matrix;
  double *prob_matrix;
  unsigned int dim_rates;  // dimension of the matrix rate_matrix

  // local minima
  string seq;
  vector<RNAlocmin> lms;
  set<string> filter;

  // outrates and inrates:
  vector<double> inrates;
  vector<double> outrates;

  // redirection of local minima

  // just for debug and shit:
  vector < map<int, int> > lengths;

public:
  RateGraph(DSU &dsu, double temp, int maxkeep, double minimal_rate, char order);
  ~RateGraph();

  int ReadFilter(char *filename);
  int ConstructQueue(char order, int number_remove, bool leave_trans);

  void PrintDot(FILE *dot);
  void PrintDot(char *filename, bool to_eps);

private:
  int RemoveOne(int rem_lm, double minimal_rate);
public:
  int RemoveX(int x, int stop_fraction = 50, bool reeval = false, int maximal = 8000, double minimal_rate = 0.0);
  int RemoveShur(int x, int step, double minimal_rate);

  void PrintRates(char *filename);
  void PrintRates(FILE *filname);

  void PrintOutput(char *filename);
  void PrintOutput(FILE *output);

  void PrintProb(double t0, char *filename);
  void PrintProb(double t0, FILE *output);
  void CreateProb(double t0);
  void CreateRates();

  int Size() {return lms.size();};

  int RemoveSmall(double *shorten, int dim, double minimal_rate);
};
#endif
