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

char* my_getline(FILE *fp);             // reads line no matter how long
string pt_to_str(const short *pt);      // convert structure short* to string
bool str_eq(const short *lhs, const short *rhs); // are structures equal?
int en_fltoi(float en); // convert energy from float to int
int HammingDist(const short* struct1, const short* struct2); // ahmming distance beetween 2 structs


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

  bool operator<(const edgeLM &second) const {
    if (i==second.i) {
      return j<second.j;
    } else return i<second.i;
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


private:
  DSU() {};

public:
  DSU(FILE *input); // read seq + structs from input (barriers/RNAlocmin output)
  ~DSU();

public:
  // big ones
  int ComputeUB(int maxkeep, int num_threshold, bool debug);          // compute all UB
  int LinkCP(bool shifts, bool noLP, bool debug);       // construct vertex and edge set from UBoutput

  // helpers
  int CreateList(int hd_threshold, bool debug);  // create TBDlist
  int FindNum(int energy, short *str);           // find number of structure
  bool InsertUB(int i, int j, int energy, short *saddle, bool debug); // insert into UBlist
  void PrintUBoutput();
  void PrintDot(char *filename, bool dot_prog, bool print, char *file_print, bool landmap);
  void PrintMatrix(char *filename);
  int FloodUp(RNAstruc &i, RNAstruc &j, RNAstruc &saddle, bool shifts, bool noLP, bool debug); // flood up to find direct saddle
};

#endif
