#ifndef __DSUEVAL_H
#define __DSUEVAL_H

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

#include "RNAutils.h"
#include "RNAstruc.h"

using namespace std;

enum type1 {INTER_CLUSTER, REPRESENT, CRIT_EDGE, NEW_FOUND};
const char type1_str[][14] = {"INTER_CLUSTER", "REPRESENT", "CRIT_EDGE", "NEW_FOUND"};
const int type1_len = 4;


struct TBDentry {
  int i, j; // two LM indexes
  type1 type_clust;
  bool fiber;
  //bool processed;

  TBDentry(int i, int j, type1 type, bool fiber) {
    this->i = i;
    this->j = j;
    this->type_clust = type;
    this->fiber = fiber;
  }

  bool operator <(const TBDentry &scnd) const {
    // first processed
    //if (processed != scnd.processed) return processed < scnd.processed;
    //else
    if (fiber != scnd.fiber) return fiber < scnd.fiber; //fibers do first
    else if (type_clust != scnd.type_clust) return type_clust > scnd.type_clust;
    else if (i != scnd.i) return i > scnd.i;
    else return j > scnd.j;
  }
};

struct TBD {
  // sets of lm to be processed
private:
  priority_queue<TBDentry> tbd;

  // set that already is processsed
public:
  set<std::pair<int, int> > done;
  int sizes[type1_len];


public:
  TBD();
  bool insert(int i, int j, type1 type, bool fiber);
  TBDentry get_first();
  int size();
  void join(TBD &second);
};

class DSU {

private:
  // sequence
  char *seq;
  short *s0;
  short *s1;

  // time of run
  clock_t time;
  int stop_after;

  // structures
  vector<RNAlocmin> LM;  // contains memory (after split - should be in both programs)
  int number_lm;  // number of lm in the beginning
  int gl_maxen;

  // lists
  priority_queue<lm_pair> TBDlist;


  // output
  vector<std::pair<RNAsaddle, lm_pair> > UBoutput; // contains memory

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
  DSU(FILE *input, bool noLP, bool shifts, int time_max = 0); // read seq + structs from input (barriers/RNAlocmin output)
  ~DSU();

public:
  // big ones
    // obsolete
  int ComputeUB(int maxkeep, int num_threshold, bool outer, bool noLP, bool shifts, bool debug);          // compute all UB (outer - add to UB also outer structures? - we will not have only direct saddles then)
  int LinkCPLM(Opt opt, bool debug);       // construct vertex and edge set from UBoutput (lm to *)
  int LinkCPsaddle(Opt opt, bool debug);       // construct vertex and edge set from UBoutput (saddle to saddle)

  // clustering
  int Cluster(Opt &opt, int kmax, TBD &list);
  //int ComputeTBD(TBD &pqueue, int maxkeep, int num_threshold, bool outer, bool noLP, bool shifts, bool debug);
  vector<std::pair<RNAsaddle, lm_pair> > ComputeTBD2(TBD &pqueue, int maxkeep, int num_threshold, bool outer, bool noLP, bool shifts, bool debug);
  int AddLMtoTBD(short *tmp_str, int tmp_en, LMtype type, bool debug);
  int JoinClusters(Opt &opt, UF_set_child &ufset, set<int> &represents, TBD &output, int i, int j);

  // helpers
  int CreateList(int hd_threshold, bool debug);  // create TBDlist
  int FindNum(int energy, short *str);           // find number of structure
  bool InsertUB(map<lm_pair, RNAsaddle, pq_setcomp> &UBlist, int i, int j, int energy, short *saddle, bool outer, bool debug); // insert into UBlist
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
  void PrintLinkCP(FILE *output = stdout, bool fix = true);
  void PrintLM(FILE *output = stdout, bool fix = true);
  void PrintSaddles(FILE *output = stdout, bool fix = true);
  void PrintComps(FILE *output = stdout, bool fill = true);
  void PrintUBoutput(FILE *output = stdout); // print UBlist to stdout

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
  void ERank(FILE *rank, bool barrier, bool out_conns = false);
};
#endif
