#ifndef __BHGBUILD_H
#define __BHGBUILD_H

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
#include "findpath_pk.h"

#include "RNAutils.h"
#include "RNAstruc.h"
#include "pknots.h"

using namespace std;

enum type1 {INTER_CLUSTER, REPRESENT, CRIT_EDGE, NEW_FOUND, EXPERIM};
const char type1_str[][14] = {"INTER_CLUSTER", "REPRESENT", "CRIT_EDGE", "NEW_FOUND", "EXPERIM"};
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
  vector <vector<bool> > done;
  int sizes[type1_len];


public:
  TBD();
  bool insert(int i, int j, type1 type, bool fiber);
  TBDentry get_first();
  int size();
  void join(TBD &second);
private:
  void ResizeDone(int new_size);
};

struct PKNOT_TYPES {
  bool H;
  bool K;
  bool L;
  bool M;

  PKNOT_TYPES() {
    H = K = L = M = false;
  }

  PKNOT_TYPES(BPAIR_TYPE bpt) {
    H = K = L = M = false;
    Add(bpt);
  }

  void Add(BPAIR_TYPE bpt) {
    switch (bpt) {
    case P_H: H = true; break;
    case P_K: K = true; break;
    case P_L: L = true; break;
    case P_M: M = true; break;
    default: break;
    }
  }

  const char* Print() {
    static string res = "";
    res[0] = H?'H':'_';
    res[1] = K?'K':'_';
    res[2] = L?'L':'_';
    res[3] = M?'M':'_';
    res[4] = '\0';
    return res.c_str();
  }
};

struct HeightData {
  int height; // energies
  int bdist; // base pair distances
  int gdist; // graph distances
  PKNOT_TYPES ptype;


  HeightData(int en, int gdist, int bdist, PKNOT_TYPES ptyp) {
    height = en;
    this->bdist = bdist;
    this->gdist = gdist;
    this->ptype = ptyp;
  }
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

  // temperature
  double _kT;

  // computing with pseudoknots?
  bool pknots;

  // structures
  vector<RNAlocmin> LM;  // contains memory (after split - should be in both programs)
  int number_lm;  // number of lm in the beginning
  int gl_maxen;

  // output
  set<RNAsaddle, RNAsaddle_comp> UBlist; // set of saddles (cannot be two with same connection)

  // -- linkCP
    // vertex sets
    map<RNAstruc, int> vertex_l;  // points to number in LM
    //map<RNAstruc, int> vertex_s;  // points to number in saddles

    vector<RNAsaddle> saddles; // contains memory

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

  // UFset for no-conn in ComputeTBD
  UF_set conectivity;

  // mapping for not full matrices
  vector<int> mapping;
  vector<int> mapping_rev;

private:
  DSU() {};

public:
  DSU(FILE *input, bool noLP, bool shifts, bool pknots, int time_max, int max_lm, bool just_read, bool debug); // read seq + structs from input (barriers/RNAlocmin output)
  ~DSU();

public:
  // big ones
    // obsolete
  int LinkCPLM(Opt opt, bool debug);       // construct vertex and edge set from saddles (lm to *)
  int LinkCPsaddle(Opt opt, bool debug);       // construct vertex and edge set from saddles (saddle to saddle)

  // clustering
  int Cluster(Opt &opt, int kmax);
  void ComputeTBD(TBD &pqueue, int maxkeep, int num_threshold, bool outer, bool noLP, bool shifts, bool debug, vector<RNAsaddle> *output_saddles = NULL, int conn_neighs = 0);
  int AddLMtoTBD(short *tmp_str, int tmp_en, LMtype type, bool debug);
  int JoinClusters(Opt &opt, UF_set_child &ufset, set<int> &represents, TBD &output, int i, int j);
  void GetRepre(TBD &output, set<int> &represents, set<int> &children, Opt &opt);

  // helpers
  int FindNum(int energy, short *str);           // find number of structure
  void FindNumbers(int begin, int end, path_t *path, vector<int> &lm_numbers, bool shifts, bool noLP, bool debug); // find all numbers of LM on path by bisection
  void FindNumbers(int begin, int end, path_pk *path, vector<int> &lm_numbers, bool shifts, bool noLP, bool debug); // and the PK version
  bool InsertUB(RNAsaddle saddle, bool debug); // insert into UBlist
    // link cp
  int FloodUp(RNAlocmin &i, RNAlocmin &j, RNAsaddle &saddle, Opt &opt, bool debug); // flood up from i and j to find direct saddle
  bool FloodSaddle(RNAsaddle &saddle_lower, RNAsaddle &saddle_higher, Opt &opt, bool debug); // flood saddle

  // visualisation
  void VisPath(int src, int dest, bool en_barriers, int max_length, bool dot_prog, bool debug);
  vector<SimplePath> ConstructAllPaths(int source, int dest, int max_length, int threshold);
  void ConstructPath(vector<SimplePath> &paths, SimplePath &path, int dest, int max_length, int threshold);
  void PrintDot(char *filename, bool dot_prog, bool print, char *file_print, bool visual, bool print_energies); // print dot file to filename, dot_prog - use dot or neato?; print - print dot output to file_print, visual - use tree for visualisation

  // print text
  void PrintLinkCP(FILE *output = stdout, bool fix = true);
  void PrintLM(FILE *output = stdout, bool fix = true);
  void PrintSaddles(FILE *output = stdout, bool fix = true);
  void PrintComps(FILE *output = stdout, bool fill = true);
  void PrintBarr(FILE *output = stdout);

  // print files:
  void PrintMatrix(char *filename, bool full, char *filter_file, char type); // print matrices (E - energy, D - distance, G - graph distance, B - base pair distance)
  void PrintRates(char *filename, bool full, double temp, char mode);

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

  void SetkT(double temp) {_kT = 0.00198717*(273.15 + temp);}

  // graph techniques
  // -- Height-first Search  -- returns energy barriers to get i-th minima from start minima + distances in graph
  vector<HeightData> HeightSearch(int start, vector< set<edgeLL> > &edgesV_l);

  // prints optimal path between start and stop to filename.
  void GetPath(int start, int stop,  vector< set<edgeLL> > &edgesV_l, char *filename, int maxkeep, bool rate_path);
  void GetPath(int start, int stop, int maxkeep, bool rate_path);

  // evaluation
  void EHeights(FILE *heights, bool full, bool only_norm);
  void ERank(FILE *rank, bool barrier, bool out_conns = false);

  void Histo(FILE *histo);
};
#endif
