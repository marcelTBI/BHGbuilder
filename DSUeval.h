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
#include "RNAstruc.h"

using namespace std;

class DSU {

private:
  // sequence
  char *seq;
  short *s0;
  short *s1;

  // structures
  vector<RNAlocmin> LM;  // contains memory (after split - should be in both programs)
  int gl_maxen;

  // lists
  priority_queue<pq_entry> TBDlist;
  map<pq_entry, RNAsaddle, pq_setcomp> UBlist; // its free in the end

  // output
  vector<std::pair<RNAsaddle, pq_entry> > UBoutput; // contains memory

  // -- linkCP
    // vertex sets
    map<RNAstruc, int> vertex_l;  // points to number in LM
    //map<RNAstruc, int> vertex_s;  // points to number in saddles

    vector<RNAsaddle> saddles; // contains memory (replaces saddles)

    //vector<RNAstruc> saddles; // contains memory (same as UBoutput)

    // edge sets
    set<edgeLL> edges_l; // LM to LM
    set<edgeSS> edges_s; // saddle to saddle
    set<edgeLS> edges_ls; // i is LM, j is saddle

    // edges for graph search
    vector< set<edgeLL> > edgesV_l;

    // components
    vector<Component> comps;
    map<int, int> LM_to_comp;
    map<int, int> saddle_to_comp;

  // -- Height-first Search


private:
  DSU() {};

public:
  DSU(FILE *input); // read seq + structs from input (barriers/RNAlocmin output)
  ~DSU();

public:
  // big ones
  int ComputeUB(int maxkeep, int num_threshold, bool outer, bool noLP, bool shifts, bool debug);          // compute all UB (outer - add to UB also outer structures? - we will not have only direct saddles then)
  int LinkCP(Opt opt, bool debug);       // construct vertex and edge set from UBoutput

  // helpers
  int CreateList(int hd_threshold, bool debug);  // create TBDlist
  int FindNum(int energy, short *str);           // find number of structure
  bool InsertUB(int i, int j, int energy, short *saddle, bool outer, bool debug); // insert into UBlist
  int AddLMtoDSU(short *tmp_str, int tmp_en, int hd_threshold, LMtype type, bool debug);
    // link cp
  void PrintMatrix(char *filename); // print energy barrier matrix
  int FloodUp(RNAlocmin &i, RNAlocmin &j, RNAsaddle &saddle, Opt &opt, bool debug); // flood up from i and j to find direct saddle
  bool FloodSaddle(RNAsaddle &saddle_lower, RNAsaddle &saddle_higher, Opt &opt, bool debug); // flood saddle

  // visualisation
  void VisPath(int src, int dest, bool en_barriers, int max_length, bool dot_prog, bool debug);
  vector<SimplePath> ConstructAllPaths(int source, int dest, int max_length, int threshold);
  void ConstructPath(vector<SimplePath> &paths, SimplePath &path, int dest, int max_length, int threshold);
  void PrintDot(char *filename, bool dot_prog, bool print, char *file_print, bool visual); // print dot file to filename, dot_prog - use dot or neato?; print - print dot output to file_print, visual - use tree for visualisation

  // print text
  void PrintLinkCP(bool fix = true);
  void PrintComps(bool fill = true);
  void PrintUBoutput(); // print UBlist to stdout

    // find components
  void FillComps();
  void Color(int lm, int color, Component &cmp, vector<int> &LM_tmp, vector<int> &sadd_tmp);

  // connect components
  void ConnectComps(int maxkeep, bool debug);
  int AddLMtoComp(short *structure, int energy, bool debug, UF_set &connected, vector<vector<RNAsaddle*> > &connections);
  int AddConnection(int num1, int num2, int energy, short *saddle, UF_set &connected, vector<vector<RNAsaddle*> > &connections);

  // sort results and fix connections, output TRUE if recomputed connections
  bool SortFix();

  // small
  int Size() {return LM.size();}

  // graph techniques
  // -- Height-first Search  -- returns energy barriers to get i-th minima from start minima + distances in graph
  vector<std::pair<int, int> > HeightSearch(int start, vector< set<edgeLL> > &edgesV_l);

  // evaluation
  void EHeights(FILE *heights, bool full);
};
#endif
