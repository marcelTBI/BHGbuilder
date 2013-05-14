#ifndef __RNAstruc_H
#define __RNAstruc_H

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <vector>
#include <queue>
#include <string>
#include <map>
#include <set>

#include "RNAutils.h"

using namespace std;

enum LMtype { NORMAL, NORM_CF, EE_DSU, EE_COMP }; // normal type, normal which was not in first list, exceeds energy in DSUeval, exceeds energy in connect components
enum SDtype { DIRECT, LDIRECT, NOT_SURE, COMP };   // direct saddle - but not sure if lowest, for sure lowest direct saddle, not sure -- only with outer option, saddle from component join (direct, but maybe not lowest, principially same as DIRECT)

struct RNAstruc {
  int energy;
  short *structure; // in short * format
  char *str_ch;     // in normal char format

  bool operator<(const RNAstruc &second) const {
    if (energy==second.energy) {
      // comparator for structures (in notation as used in energy_of_move)
      int i=1;
      char l=0,r=0;
      while (i<=structure[0]) {
        l = (structure[i]==0?'.':(structure[i]<structure[structure[i]]?')':'('));
        r = (second.structure[i]==0?'.':(second.structure[i]<second.structure[second.structure[i]]?')':'('));
        if (l != r) break;
        i++;
      }
      //fprintf(stderr, "%s %c %s\n", pt_to_str(lhs).c_str(), (i<=lhs[0] && l<r) ? '<':'>', pt_to_str(rhs).c_str());
      if (i>structure[0]) return false;
      return l<r;
    } else return energy<second.energy;
  }

  bool operator==(const RNAstruc &second) const {
    if (energy==second.energy) {
      // comparator for structures (in notation as used in energy_of_move)
      int i=1;
      char l=0,r=0;
      while (i<=structure[0]) {
        l = (structure[i]==0?'.':(structure[i]<structure[structure[i]]?')':'('));
        r = (second.structure[i]==0?'.':(second.structure[i]<second.structure[second.structure[i]]?')':'('));
        if (l != r) break;
        i++;
      }
      //fprintf(stderr, "%s %c %s\n", pt_to_str(lhs).c_str(), (i<=lhs[0] && l<r) ? '<':'>', pt_to_str(rhs).c_str());
      if (i>structure[0]) return true;
      return false;
    } else return false;
  }

  /*RNAstruc(int energy, short *str, bool usable = true) {
    this->energy = energy;
    this->structure = str;
    this->usable = usable;
    if (usable) str_ch = pt_to_char(str);
    else str_ch = NULL;
  }*/

  RNAstruc() {
    structure  = NULL;
    str_ch = NULL;
  }

  void freeMEM();
  void recompute_str();
};

struct RNAlocmin: public RNAstruc {
  LMtype type;

  RNAlocmin():RNAstruc() {
    type = NORMAL;
  }
};

// saddle has 2lm that connects
struct RNAsaddle: public RNAstruc {
  int lm1;
  int lm2;
  SDtype type;

  RNAsaddle(int lm1, int lm2, SDtype type = DIRECT):RNAstruc() {
    this->lm1 = min(lm1, lm2);
    this->lm2 = max(lm1, lm2);
    this->type = type;
  }
};

// just for reverse ordering
struct RNAsaddle_comp {
  bool operator()(const RNAsaddle &first, const RNAsaddle &second) const {
    if (first.lm1 == second.lm1) return first.lm2 < second.lm2;
    else return first.lm1 < second.lm1;
  }
};

// options structure;
struct Opt {
  bool noLP;
  bool shifts;
  bool saddle_conn;

  int flood_num;
  int flood_height;

  // for clustering
  bool debug;
  int maxkeep;
  int num_threshold;
  bool outer;
  float repre_portion;
  bool fbarrier;

  Opt(bool noLP, bool shifts, bool saddle, int num, int height, bool debug, int maxkeep, int num_threshold, bool outer, float repre_portion, bool fbarrier);
};

// just for reverse ordering
struct RNAstruc_rev {
  bool operator()(const RNAstruc &first, const RNAstruc &second) const {
    if (first.energy==second.energy) {
      // comparator for structures (in notation as used in energy_of_move)
      int i=1;
      char l=0,r=0;
      while (i<=first.structure[0]) {
        l = (first.structure[i]==0?'.':(first.structure[i]<first.structure[first.structure[i]]?')':'('));
        r = (second.structure[i]==0?'.':(second.structure[i]<second.structure[second.structure[i]]?')':'('));
        if (l != r) break;
        i++;
      }
      //fprintf(stderr, "%s %c %s\n", pt_to_str(lhs).c_str(), (i<=lhs[0] && l<r) ? '<':'>', pt_to_str(rhs).c_str());
      if (i>first.structure[0]) return false;
      return l>r;
    } else return first.energy>second.energy;
  }

};

// entry for priority queue
struct lm_pair {
  int i, j;
  int hd;

  lm_pair(int i, int j, int hd) {
    this->i = min(i,j);
    this->j = max(i,j);
    this->hd = hd;
  }

  bool operator <(const lm_pair &sec) const {
    if (hd != sec.hd) {
      return hd<sec.hd;
    } else {
      if (i != sec.i) {
        return i<sec.i;
      } else {
        return j<sec.j;
      }
    }
  }
};
// maybe we dont need this for map...
struct pq_setcomp {
  long operator() (const lm_pair &l, const lm_pair &r) const {
    if (l.i == r.i) return l.j<r.j;
    else return l.i < r.i;
  }
};

//edges in graph
struct edge {
  int i;
  int j;

  edge(int ii, int jj) {
    i = ii;
    j = jj;
  }

  bool operator<(const edge &second) const {
    if (i==second.i) {
      return j<second.j;
    } else return i<second.i;
  }

  int goesTo(int src) const { if (i==src) return j; else return i;}

};

struct edgeLL : public edge {
  int saddle;
  int en;

  edgeLL(int ii, int jj, int energy, int saddle):edge(ii,jj) {
    if (i > j) swap(i,j);
    this->saddle = saddle;
    this->en = energy;
  }
};

// energy comparator
struct edgeLL_compen {
  bool operator()(const edgeLL &first, const edgeLL &second) const {
    if (first.en==second.en) {
      if (first.i==second.i) {
        return first.j<second.j;
      } else return first.i<second.i;
    } else return first.en<second.en;
  }
};

struct edgeSS : public edge {
  edgeSS(int i, int j):edge(i,j){
    if (i > j) swap(i,j);
  }
};

struct edgeLS : public edge {
  edgeLS(int lm, int saddle):edge(lm,saddle) {}
};

// edge that has whole refolding path inside (just one refolding path)
class edgeAdv : public edge {
public:
  // highest saddle height
  int max_height;

  // saddle numbers and their energies of the best path
  vector<int> saddles;
  vector<int> energies;

  // rate - should not always correspond to the best path...
  double rate; // now unused

  edgeAdv(int ii, int jj, int energy, int saddle):edge(ii,jj) {
    if (i > j) swap(i,j);
    saddles.push_back(saddle);
    energies.push_back(energy);
    max_height = energy;
  }

  int length() const { return saddles.size(); }
};

// comparator according to lowest energy (and length)
struct edge_comp {
  bool operator ()(const edgeAdv &a, const edgeAdv &b) const {
    if (a.max_height == b.max_height) {
      if (a.length() == b.length()) return a<b;
      else return a.length() < b.length();
    } else return a.max_height < b.max_height;
  }
};

// modes for rates generation
enum mode_rates {VERTEX_CONTR, EDGE_CONTR};

struct Graph {

  // maximal_node
  int max_node;

  // lm for printing:
  int number_lm;

  // mode for joining
  mode_rates mode;

  // lm - just shallow copy for now
  vector<RNAlocmin> &LM;

  // edges -- sparse format adjacency[x]=y only if x>y
  map<std::pair<int, int>, edgeAdv> adjacency;


public:
  // constructor
  Graph(set<edgeLL> &edges, vector<RNAlocmin> &LM, mode_rates mode);
private:
  // helpers
  int Join(const edgeAdv &src, const edgeAdv &dst, int joining_node, edgeAdv &res);
  bool AddEdges(const edgeAdv &found, edgeAdv &res, mode_rates);
public:
  // methods
  int RemoveLastPoint();
  void PrintDot(char *filename, bool dot_prog, bool print, char *file_print);
  void PrintRates(FILE *rates, double temp);
};

// component structure
struct Component {
  // minimal (LM) and maximal (saddle) part
  int min_lm;
  int min_energy;
  int max_saddle;
  int max_energy;
  // local minima and saddles
  vector<int> LMs;
  vector<int> saddles;

  Component() {
    min_lm = -1;
    max_saddle = -1;
    min_energy = INT_MAX;
    max_energy = INT_MIN;
  }

  void AddLM (int num, int energy) {
    if (energy < min_energy) {
      min_lm = num;
      min_energy = energy;
    }
    LMs.push_back(num);
  }

  void AddSadd (int num, int energy){
    if (energy > max_energy) {
      max_saddle = num;
      max_energy = energy;
    }
    saddles.push_back(num);
  }
};

class SimplePath {
public:
  vector<int> points;
  vector<int> energies;

  map<int, int> points_map;

  int max_energy;

  // closed?
  bool closed;

  bool operator<(const SimplePath &left) const {
    if (max_energy == left.max_energy) return points.size() < left.points.size();
    else return max_energy < left.max_energy;
  }

public:
  SimplePath();
  SimplePath(const SimplePath &path);

  void AddPoint(int num, int energy);
  void Close();

  void AddLast(int num, int energy);
  void RemoveLast();

  void Score();

  bool ContainsNode(int num);
  int FindNode(int num);
  void Print(bool whole_path, bool force_print = true, FILE *out = stdout);
};

#endif
