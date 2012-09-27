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

using namespace std;

class SimplePath;

struct RNAstruc {
  int energy;
  short *structure;
  char *str_ch; // must not be filled

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
  int j;
  int en;

  int sadd; // number of saddle

  edgeLM(int ii, int jj, int enn, int saddl = -1) {
    i = min(ii, jj);
    j = max(ii, jj);
    en = enn;
    sadd = saddl;
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


struct Component {
  int min_lm;
  int max_energy;
  int min_energy;
  vector<int> LMs;

  Component() {
    min_lm = -1;
    max_energy = INT_MIN;
    min_energy = INT_MAX;
  }

  void UpdateEn(int energy) {
    max_energy = max(max_energy, energy);
  }

  void AddLM(int num, int energy) {
    if (energy < min_energy) {
      min_lm = num;
      min_energy = energy;
    }
    max_energy = max(max_energy, energy);
    LMs.push_back(num);
  }
};

class DSU {

private:
  // sequence
  char *seq;
  short *s0;
  short *s1;

  // structures
  vector<RNAstruc> LM;
  int gl_maxen;

  // lists
  priority_queue<pq_entry> TBDlist;
  map<pq_entry, RNAstruc, pq_setcomp> UBlist;

  // output
  vector<std::pair<RNAstruc, pq_entry> > UBoutput;

  // -- linkCP
    // vertex sets (Vl is vector<short*> str)
    map<RNAstruc, int> vertex_l;  // points to number in LM
    map<RNAstruc, int> vertex_s;  // points to number in saddles

    vector<RNAstruc> saddles;

    // edge sets
    set<edgeLM> edges_l;
    set<std::pair<int, int> > edges_s;
    set<std::pair<int, int> > edges_ls; // first int is l, second is s

    // edges for graph search
    vector< set<edgeLM> > edgesV_l;

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
  void VisPath(int src, int dest, bool en_barriers, int max_length, bool dot_prog, bool debug);

  vector<SimplePath> ConstructAllPaths(int source, int dest, int max_length, int threshold);
  void ConstructPath(vector<SimplePath> &paths, SimplePath &path, int dest, int max_length, int threshold);
  void PrintLinkCP();
  void Color(int lm, int color, Component &cmp, vector<int> &LM_tmp);

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
