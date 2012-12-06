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
        l = (structure[i]==0?'.':(structure[i]<structure[structure[i]]?'(':')'));
        r = (second.structure[i]==0?'.':(second.structure[i]<second.structure[second.structure[i]]?'(':')'));
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
        l = (structure[i]==0?'.':(structure[i]<structure[structure[i]]?'(':')'));
        r = (second.structure[i]==0?'.':(second.structure[i]<second.structure[second.structure[i]]?'(':')'));
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

// options structure;
struct Opt {
  bool noLP;
  bool shifts;
  bool saddle_conn;

  int flood_num;
  int flood_height;

  Opt(bool noLP, bool shifts, bool saddle, int num, int height);
};

// just for reverse ordering
struct RNAstruc_rev {
  bool operator()(const RNAstruc &first, const RNAstruc &second) const {
    if (first.energy==second.energy) {
      // comparator for structures (in notation as used in energy_of_move)
      int i=1;
      char l=0,r=0;
      while (i<=first.structure[0]) {
        l = (first.structure[i]==0?'.':(first.structure[i]<first.structure[first.structure[i]]?'(':')'));
        r = (second.structure[i]==0?'.':(second.structure[i]<second.structure[second.structure[i]]?'(':')'));
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
class pq_entry {
public:
  int i;       // i is always smaller
  int j;
  int hd;     // hamming dist

public:
  bool operator<(const pq_entry &second) const {
    if (hd==second.hd) {
      if (i==second.i) {
        return j<second.j;
      } else return i<second.i;
    } else return hd>second.hd;
  }

private:
  pq_entry() {};

public:
  pq_entry(int i, int j, int hd);
  ~pq_entry();
};

class saddle {

};


// maybe we dont need this for map...
struct pq_setcomp {
  long operator() (const pq_entry &l, const pq_entry &r) const {
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

struct edgeSS : public edge {
  edgeSS(int i, int j):edge(i,j){
    if (i > j) swap(i,j);
  }
};

struct edgeLS : public edge {
  edgeLS(int lm, int saddle):edge(lm,saddle) {}
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
