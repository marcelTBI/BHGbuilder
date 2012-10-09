#ifndef __DSUEVAL_H
#define __DSUEVAL_H

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

class SimplePath;

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
};

struct RNAstruc2 : public RNAstruc {
  int conn1; // LM it connects
  int conn2;
};

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

// maybe we dont need this for map...
struct pq_setcomp {
  long operator() (const pq_entry &l, const pq_entry &r) const {
    if (l.i == r.i) return l.j<r.j;
    else return l.i < r.i;
  }
};

struct edgeLM {
  int i;
  bool LMi;
  int j;
  bool LMj;

  int en;   // energy and ->
  int sadd; // number of saddle (just if LMi and LMj are true)

  bool component; // if info comes from join components

  edgeLM(int ii, int jj, bool LMii, bool LMjj) {
    i = ii;
    j = jj;
    en = INT_MAX;
    sadd = -1;
    component = false;
    LMi = LMii;
    LMj = LMjj;

    // bad i-j positions
    if (((LMi && LMj) || (!LMi && !LMj)) && (i > j)) {
      swap(i, j);
    }

    // bad position saddle/LM
    if (!LMi && LMj) {
      swap(i,j);
      swap(LMi, LMj);
    }
  }

  void AddSaddle(int energy, int saddle) {
    en = energy;
    sadd = saddle;
  }

  void MarkComp() {
    component = true;
  }

  int goesTo(int src) const { if (i==src) return j; else return i;}

  bool operator<(const edgeLM &second) const {
    if (i==second.i) {
      return j<second.j;
    } else return i<second.i;
  }
};
// energy comparator
struct edgeLM_compen {

  bool operator()(const edgeLM &first, const edgeLM &second) const {
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

class DSU {

private:
  // sequence
  char *seq;
  short *s0;
  short *s1;

  // structures
  vector<RNAstruc> LM;  // contains memory
  int gl_maxen;

  // lists
  priority_queue<pq_entry> TBDlist;
  map<pq_entry, RNAstruc, pq_setcomp> UBlist; // its free in the end unsually

  // output
  vector<std::pair<RNAstruc, pq_entry> > UBoutput; // contains memory

  // -- linkCP
    // vertex sets
    map<RNAstruc, int> vertex_l;  // points to number in LM
    map<RNAstruc, int> vertex_s;  // points to number in saddles

    vector<RNAstruc> saddles; // contains memory (same as UBoutput)

    // edge sets
    set<edgeLM> edges_l; // LM to LM
    set<edgeLM> edges_s; // saddle to saddle
    set<edgeLM> edges_ls; // i is LM, j is saddle

    // edges for graph search
    vector< set<edgeLM> > edgesV_l;

    // components
    vector<Component> comps;
    map<int, int> LM_to_comp;

private:
  DSU() {};

public:
  DSU(FILE *input); // read seq + structs from input (barriers/RNAlocmin output)
  ~DSU();

public:
  // big ones
  int ComputeUB(int maxkeep, int num_threshold, bool outer, bool debug);          // compute all UB (outer - add to UB also outer structures? - we will not have only direct saddles then)
  int LinkCP(bool shifts, bool noLP, bool debug);       // construct vertex and edge set from UBoutput

  // helpers
  int CreateList(int hd_threshold, bool debug);  // create TBDlist
  int FindNum(int energy, short *str);           // find number of structure
  bool InsertUB(int i, int j, int energy, short *saddle, bool debug); // insert into UBlist
  void PrintUBoutput(); // print UBlist to stdout
    // link cp
  void PrintDot(char *filename, bool dot_prog, bool print, char *file_print, bool visual); // print dot file to filename, dot_prog - use dot or neato?; print - print dot output to file_print, visual - use tree for visualisation
  void PrintMatrix(char *filename); // print energy barrier matrix
  int FloodUp(RNAstruc &i, RNAstruc &j, RNAstruc &saddle, bool shifts, bool noLP, bool debug); // flood up from i and j to find direct saddle
  bool FloodSaddle(RNAstruc &saddle_lower, RNAstruc &saddle_higher, bool shifts, bool noLP, bool debug); // flood saddle
  void VisPath(int src, int dest, bool en_barriers, int max_length, bool dot_prog, bool debug);

  vector<SimplePath> ConstructAllPaths(int source, int dest, int max_length, int threshold);
  void ConstructPath(vector<SimplePath> &paths, SimplePath &path, int dest, int max_length, int threshold);
  void FillComps();
  void ConnectComps(int maxkeep, bool debug);
  int AddConnection(int num1, int num2, int energy, short *saddle, UF_set &connected, vector<vector<RNAstruc2*> > &connections);
  void PrintLinkCP(bool full);
  void Color(int lm, int color, Component &cmp, vector<int> &LM_tmp, vector<int> &sadd_tmp);

  // small
  int Size() {return LM.size();}
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
