#include <stdio.h>
#include <string.h>
#include <stdlib.h>

extern "C" {
  #include "utils.h"
  #include "fold_vars.h"
  #include "pair_mat.h"

  #include "fold_dsu.h"
  #include "findpath.h"
  #include "move_set.h"
}

#include "DSUeval.h"
#include "RNAutils.h"
#include "hash_util.h"

#include <algorithm>

using namespace std;

// union-find set array
vector<int> parent;
unsigned int num_unions = 0;

// and union-find set functions
int find(int x) {
  if (x != parent[x] && parent[x] != parent[parent[x]])
    parent[x] = find(parent[x]);
  return parent[x];
}

void union_set(int x, int y) {
  int u, v;
  u = find(x);
  v = find(y);
  if (u != v) {
    parent[u] = v;
    num_unions++;
  }
}

bool connected_all() {
  return (num_unions == parent.size()-1);
}

bool joint(int x, int y) {
  return find(x) == find(y);
}

void enlarge_parent() {
  parent.push_back(parent.size());
}

//============================ MAIN functions

pq_entry::pq_entry(int i, int j, int hd)
{
  this->i = min(i, j);
  this->j = max(i, j);
  this->hd = hd;
}

pq_entry::~pq_entry() {
}

DSU::DSU(FILE *input) {

  // NULL::
  seq = NULL;
  s0 = NULL;
  s1 = NULL;
  gl_maxen = INT_MIN;

  if (!input) return ;

  // read seq
  char *line;
  int num = 0;

  line = my_getline(input);
  char *seq2 = strtok(line, " ");
  if (!isSeq(seq2)) {
    free(line);
    return ;
  }
  seq = (char*) malloc((strlen(seq2)+1)*sizeof(char));
  strcpy(seq, seq2);
  free(line);

  //init
  make_pair_matrix();
  s0 = encode_sequence(seq, 0);
  s1 = encode_sequence(seq, 1);

  //read structs
  line = my_getline(input);
  while (line) {
    short *tmp = NULL;
    float energy_tmp;
    bool has_energy = false;

    char *p;
    p = strtok(line, " ");
    while (p !=NULL && (!has_energy || !tmp)) {
      // is struct?
      if (isStruct(p)) {
        if (tmp) free(tmp); // only one struct per line!
        tmp = make_pair_table(p);
        has_energy = false;
      } else {
        // is energy?
        if (sscanf(p, "%f", &energy_tmp)==1) {
          has_energy = true;
        }
      }
      p = strtok(NULL, " ");
    }

    // add info:
    if (tmp) {
      RNAstruc struc;
      struc.structure = tmp;
      char *ch = (char*) malloc((tmp[0]+1)*sizeof(char));
      strcpy(ch, pt_to_str(tmp).c_str());
      struc.str_ch = ch;
      int en = (has_energy?en_fltoi(energy_tmp):energy_of_structure_pt(seq, tmp, s0, s1, 0));
      struc.energy = en;
      LM.push_back(struc);
      gl_maxen = max(gl_maxen, en);
    } else {
      // line without struct???
      fprintf(stderr, "WARNING: line %d without struct: %s\n", num, line);
    }
    free(line);
    line = my_getline(input);
    num++;

    sort(LM.begin(), LM.end());
  }
}

DSU::~DSU() {
  for (unsigned int i=0; i<LM.size(); i++) {
    if (LM[i].structure) free(LM[i].structure);
    if (LM[i].str_ch) free(LM[i].str_ch);
  }
  if (seq) free(seq);
  if (s0) free(s0);
  if (s1) free(s1);
  LM.clear();

  for (map<pq_entry, RNAstruc, pq_setcomp>::iterator it=UBlist.begin(); it!=UBlist.end(); it++) {
    if (it->second.structure) free(it->second.structure);
    if (it->second.str_ch) free(it->second.str_ch);
  }

  for (unsigned int i=0; i<UBoutput.size(); i++) {
    if (UBoutput[i].first.structure) free(UBoutput[i].first.structure);
  }
}

int DSU::CreateList(int hd_threshold, bool debug)
{
  for (unsigned int i=0; i<LM.size(); i++) {
    for (unsigned int j=i+1; j<LM.size(); j++) {
      int hd = HammingDist(LM[i].structure, LM[j].structure);
      if (hd < hd_threshold) {
        TBDlist.push(pq_entry(i, j, hd));
        if (debug) {
          fprintf(stderr, "%s insert que\n%s (%4d,%4d) %3d\n", pt_to_str(LM[i].structure).c_str(), pt_to_str(LM[j].structure).c_str(), i, j, hd);
        }
      }
    }
  }

  return 0;
}

int DSU::FindNum(int en_par, short *str_par)
{
  //fprintf(stderr, "%d %s\n", en_par, pt_to_str(str_par).c_str());
  for (unsigned int i=0; i<LM.size(); i++) {
    if (LM[i].energy == en_par && str_eq(str_par, LM[i].structure)) return i;
  }
  return -1;
}

bool DSU::InsertUB(int i, int j, int energy_par, short *saddle_par, bool debug)
{
  pq_entry pq(i, j, 0);
  map<pq_entry, RNAstruc, pq_setcomp>::iterator it = UBlist.find(pq);
  if (it==UBlist.end()) {
    RNAstruc saddle;
    saddle.energy = energy_par;
    saddle.structure = saddle_par;
    saddle.str_ch = NULL;
    if (debug) fprintf(stderr, "UBins: (%3d, %3d) insert saddle  %6.2f %s\n", i, j, saddle.energy/100.0, pt_to_str(saddle.structure).c_str());
    UBlist.insert(make_pair(pq, saddle));
    return true;
    //fprintf(stderr, "cannot find (UB  ): (%3d, %3d)\n", num1, num2);
  } else {
    // update it
    if (energy_par < it->second.energy) {
      if (debug) fprintf(stderr, "UBupd: (%3d, %3d) from %6.2f to %6.2f %s\n", i, j, it->second.energy/100.0, energy_par/100.0, pt_to_str(it->second.structure).c_str());
      it->second.energy = energy_par;
      if (it->second.structure) free(it->second.structure);
      it->second.structure = saddle_par;
      it->second.str_ch = NULL;
      return true;
    }
  }
  free(saddle_par);
  return false;
}

int DSU::ComputeUB(int maxkeep, int num_threshold, bool outer, bool debug)
{
  int dbg_count = 0;
  int cnt = 0;
  int hd_threshold = INT_MAX;
  while (!TBDlist.empty()) {
    // get pair
    pq_entry pq = TBDlist.top();
    TBDlist.pop();

    // apply threhold
    if (pq.hd > hd_threshold) break;
    if (cnt>num_threshold) {
      if (hd_threshold==INT_MAX) {
        hd_threshold = pq.hd;
        fprintf(stderr, "hd threshold set to %4d\n", hd_threshold);
      }
    } else {
      cnt++;
    }

    // get path
    if (debug) fprintf(stderr, "path between (%3d, %3d) hd=%3d:\n", pq.i, pq.j, pq.hd);
    path_t *path = get_path(seq, LM[pq.i].str_ch, LM[pq.j].str_ch, maxkeep);

    // variables for outer insertion
    double max_energy= -1e8;
    path_t *max_path = path;

    // variables for inner loops and insertions
    path_t *tmp = path;
    path_t *last = NULL;
    short *last_str = NULL;
    int last_en;
    int last_num = -1;

    // loop through whole path
    while (tmp && tmp->s) {
      dbg_count++;
      // debug??
      if (debug) fprintf(stderr, "%s %6.2f\n", tmp->s, tmp->en);

      // update max_energy
      if (max_energy < tmp->en) {
        max_energy = tmp->en;
        max_path = tmp;
      }
      // find adaptive walk
      short *tmp_str = make_pair_table(tmp->s);
      //int tmp_en = move_rand(seq, tmp_str, s0, s1, 0);
      int tmp_en = move_deepest(seq, tmp_str, s0, s1, 0, 0, 0);

      // do the stuff if we have 2 structs and they are not equal
      if (last && !str_eq(last_str, tmp_str)) {
        // not equal LM - we can update something in UBlist
          // find LM num:
        int num1 = (last_num!=-1?last_num:FindNum(last_en, last_str));
        int num2 = FindNum(tmp_en, tmp_str);
        last_num = num2;

        // update UBlist
        if (num1==-1 || num2==-1) {
          if (num2==-1) {
            if (gl_maxen < tmp_en) fprintf(stderr, "exceeds en.: %s %6.2f\n", pt_to_str(tmp_str).c_str(), tmp_en/100.0);
            else {
                                   fprintf(stderr, "cannot find: %s %6.2f\n", pt_to_str(tmp_str).c_str(), tmp_en/100.0);
              // maybe add to list of minima and count with them later... TODO
            }
          }
        } else {
          // store (maybe) better saddle to UB
          int en_tmp = en_fltoi(max(last->en, tmp->en));
          short *saddle = (last->en > tmp->en ? make_pair_table(last->s) : make_pair_table(tmp->s));
          InsertUB(num1, num2, en_tmp, saddle, debug);
        }
      }

      // move one next
      if (last_str) free(last_str);
      last_en = tmp_en;
      last_str = tmp_str;
      last = tmp;
      tmp++;
    } // crawling path

    // insert saddle between outer structures
    if (outer) InsertUB(pq.i, pq.j, en_fltoi(max_energy), make_pair_table(max_path->s), debug);

    // free stuff
    if (last_str) free(last_str);
    free_path(path);
  } // all doing while

  // now just resort UBlist to something sorted according energy
  UBoutput.reserve(UBlist.size());
  for (map<pq_entry, RNAstruc, pq_setcomp>::iterator it=UBlist.begin(); it!=UBlist.end(); it++) {
    UBoutput.push_back(make_pair(it->second, it->first));
  }
  sort(UBoutput.begin(), UBoutput.end());
  UBlist.clear();

  // check if everything has been found:
  if (UBoutput.size() != (LM.size()*(LM.size()-1))/2) {
    fprintf(stderr, "WARNING: All connections have not been found: %d/%d found (%d missing)\n", (int)UBoutput.size(), (int)(LM.size()*(LM.size()-1))/2, (int)((LM.size()*(LM.size()-1))/2 - UBoutput.size()));
  }

  return 0;
}

void DSU::PrintUBoutput()
{
  printf("     %s\n", seq);
  for (unsigned int i=0; i<UBoutput.size(); i++) {
    printf("%4d (%4d,%4d) saddle: %s %6.2f\n", i+1, UBoutput[i].second.i+1, UBoutput[i].second.j+1, pt_to_str(UBoutput[i].first.structure).c_str(), UBoutput[i].first.energy/100.0);
  }
}

int DSU::LinkCP(bool shifts, bool noLP, bool debug)
{
  //edgesV_ls.resize(LM.size());
  edgesV_l.resize(LM.size());

  for (unsigned int i=0; i<UBoutput.size(); i++) {
    RNAstruc stru = UBoutput[i].first;
    pq_entry pq = UBoutput[i].second;

    // flood them up!
    int res = FloodUp(LM[pq.i], LM[pq.j], stru, shifts, noLP, debug);
    if (res == 1) {
      // maybe other color ?
    }// we have found better DS -> true ds ???

    // update vertex/edge sets
    map<RNAstruc, int>::iterator it;
    int saddle_num;
    if ((it = vertex_s.find(stru)) != vertex_s.end()) {
      saddle_num = it->second;
    } else {
      saddle_num = saddles.size();
      // update vertex
      vertex_s.insert(make_pair(stru, saddle_num));
      saddles.push_back(stru);
    }
    // resize
    //edgesV_sl.resize(saddle_num+1);

    // edges ls
    edges_ls.insert(make_pair(pq.i, saddle_num));
    edges_ls.insert(make_pair(pq.j, saddle_num));

    //edgesV_ls[pq.i].insert(saddle_num);
    //edgesV_sl[saddle_num].insert(pq.i);

    //edgesV_ls[pq.j].insert(saddle_num);
    //edgesV_sl[saddle_num].insert(pq.j);

    // edges ll
    edgeLM e(pq.i, pq.j, stru.energy, saddle_num);
    edges_l.insert(e);

    edgesV_l[pq.i].insert(e);
    edgesV_l[pq.j].insert(e);
  }

  return 0;
}

// variables and data strucutres for flooding
bool foundDS;
RNAstruc curr;
int currNum;
int threshold;
bool gl_debug;
unordered_map<RNAstruc, int, hash_fncts, hash_eq> flood_hash;
unordered_map<RNAstruc, int, hash_fncts, hash_eq>::iterator it;
priority_queue<RNAstruc, vector<RNAstruc>, RNAstruc_rev> flood_queue;

int funct(struct_en *moved, struct_en *current) {

  if (moved->energy > curr.energy && moved->energy < threshold) {
    RNAstruc tmp;
    tmp.energy = moved->energy;
    tmp.structure = moved->structure;
    if ((it = flood_hash.find(tmp))!=flood_hash.end()) {
      // found DS!
      if (it->second != currNum) {
        foundDS = true;
        copy_arr(current->structure, moved->structure);
        current->energy = moved->energy;

        if (gl_debug) {
          fprintf(stderr, "FOUND saddle: %s %6.2f\n", pt_to_str(current->structure).c_str(), current->energy/100.0);
        }

        return 1;
      } // else nothing
    } else {
      // insert strucure we havent see yet (and its in range)
      RNAstruc to_ins;
      to_ins.energy = moved->energy;
      to_ins.structure = allocopy(moved->structure);
      flood_hash.insert(make_pair(to_ins, currNum));
      flood_queue.push(to_ins);
    }
  }

  return 0;
}

int DSU::FloodUp(RNAstruc &i, RNAstruc &j, RNAstruc &saddle, bool shifts, bool noLP, bool debug)
{
  // threshold set
  threshold = saddle.energy;
  gl_debug = debug;

  // debug output
  if (debug) {
    fprintf(stderr, "\nbegin1      : %s %6.2f\n", pt_to_str(i.structure).c_str(), i.energy/100.0);
    fprintf(stderr, "begin2      : %s %6.2f\n", pt_to_str(j.structure).c_str(), j.energy/100.0);
    fprintf(stderr, "saddle      : %s %6.2f\n", pt_to_str(saddle.structure).c_str(), saddle.energy/100.0);
  }

  // init
  RNAstruc tmpi, tmpj;
  tmpi.energy = i.energy;
  tmpj.energy = j.energy;
  tmpi.structure = allocopy(i.structure);
  tmpj.structure = allocopy(j.structure);
  flood_hash.insert(make_pair(tmpi, 1));
  flood_hash.insert(make_pair(tmpj, 2));
  flood_queue.push(tmpi);
  flood_queue.push(tmpj);

  // end
  int res = 0;

  // main loop
  curr = flood_queue.top();
  flood_queue.pop();
  currNum = flood_hash[curr];
  while (curr.energy < threshold) {

    // debug output
    if (debug) {
      fprintf(stderr, "browsing    : %s %6.2f\n", pt_to_str(curr.structure).c_str(), curr.energy/100.0);
    }

    // start browsing
    foundDS = false;
    curr.energy = browse_neighs(seq, curr.structure, s0, s1, 0, shifts, noLP, funct);

    // found DS - quitting
    if (foundDS) {
      break;
    }

    // get another struct
    if (flood_queue.empty()) break;
    curr = flood_queue.top();
    flood_queue.pop();
    currNum = flood_hash[curr];
  }

  // we did end sucesfully
  if (foundDS) {
    copy_arr(saddle.structure, curr.structure);
    saddle.energy = curr.energy;
    res = 1;
  }

  // free
  while (!flood_queue.empty()) flood_queue.pop();
  for (it=flood_hash.begin(); it!=flood_hash.end(); it++) {
    if (it->first.structure) free(it->first.structure);
  }
  flood_hash.clear();
  return res;
}

void DSU::PrintDot(char *filename, bool dot_prog, bool print, char *file_print, bool visual)
{
  // landmap not supported yet

  float color = 0.7;
  //open file
  FILE *dot;
  dot = fopen(filename, "w");
  if (dot) {
    fprintf(dot, "Graph G {\n\tnode [width=0.1, height=0.1, shape=circle];\n");
    //nodes LM:
    for (unsigned int i=0; i<LM.size(); i++) {
      fprintf(dot, "\"%d\" [label=\"%d\"]\n", i+1, i+1);
    }
    fprintf(dot, "\n");

    if (visual) {
      for (unsigned int i=0; i<LM.size(); i++) {
        enlarge_parent();
      }
      set<edgeLM, edgeLM_compen> tmp;
      tmp.insert(edges_l.begin(), edges_l.end());
      for (set<edgeLM>::iterator it=tmp.begin(); it!=tmp.end(); it++) {
        if (!joint(it->i, it->j)) {
          fprintf(dot, "\"%d\" -- \"%d\" [label=\"%.2f\"]\n", (it->i)+1, (it->j)+1, it->en/100.0);
          union_set(it->i, it->j);
        }
      }
      fprintf(dot, "\n");

    } else {
      //nodes saddle:
      for (unsigned int i=0; i<saddles.size(); i++) {
        fprintf(dot, "\"S%d\" [label=\"S%d\", color=\"0.0 0.0 %.1f\", fontcolor=\"0.0 0.0 %.1f\"]\n", i+1, i+1, color, color);
      }
      fprintf(dot, "\n");
      // edges l-l
      for (set<edgeLM>::iterator it=edges_l.begin(); it!=edges_l.end(); it++) {
        fprintf(dot, "\"%d\" -- \"%d\" [label=\"%.2f\"]\n", (it->i)+1, (it->j)+1, it->en/100.0);
      }
      fprintf(dot, "\n");
      // edges l-s
      for (set<std::pair<int, int> >::iterator it=edges_ls.begin(); it!=edges_ls.end(); it++) {
        fprintf(dot, "\"%d\" -- \"S%d\" [color=\"0.0 0.0 %.1f\", fontcolor=\"0.0 0.0 %.1f\"]\n", (it->first)+1, (it->second)+1, color, color);
      }

      fprintf(dot, "\n");
      // edges s-s //TODO!
      /*for (set<edgeLM>::iterator it=edges_l.begin(); it!=edges_l.end(); it++) {
        fprintf(dot, "\"%d\" -- \"%d\" [label=\"%.2f\"]\n",it->i,it->j, it->en/100.0);
      }*/
    }
    fprintf(dot, "}\n");
  }

  fclose(dot);

  // start neato/dot:
  if (dot && print && file_print) {
    char syst[200];
    sprintf(syst, "%s -Tps < %s > %s", (dot_prog ? "dot" : "neato"), filename, file_print);
    int res = system(syst);
    printf("%s returned %d\n", syst, res);
  }
}

bool EN_BARRIERS = true;

struct pq_path {
  int lm;
  int dist;
  int en_barr;

  int value() {
    return (EN_BARRIERS ? en_barr : dist);
  }

  bool operator<(const pq_path &second) const {
    if (EN_BARRIERS) {
      if (en_barr != second.en_barr) {
        return en_barr<second.en_barr;
      }
    }
    if (dist == second.dist) {
      return lm<second.lm;
    }
    return dist<second.dist;
  }

  pq_path(int lm, int dist, int en_barr) {
    this->lm = lm;
    this->dist = dist;
    this->en_barr = en_barr;
  }
};

void DSU::VisPath(int src, int dest, bool en_barriers, int max_length, bool dot_prog, bool debug)
{
  EN_BARRIERS = en_barriers;

  char filename[50];
  char file_print[50];
  sprintf(filename, "path%d_%d.dot", src+1, dest+1);
  sprintf(file_print, "path%d_%d.eps", src+1, dest+1);

  float color = 0.5;
  //open file
  FILE *dot;
  dot = fopen(filename, "w");
  if (dot) {
    fprintf(dot, "Graph G {\n\tnode [width=0.1, height=0.1, shape=circle];\n");

    // actual nodes:
      // forward pass:
    vector <int> LM_tmp(LM.size(), INT_MAX);
    priority_queue <pq_path> pq_tmp;
    pq_path to_insert(src, 0, LM[src].energy);
    LM_tmp[src] = to_insert.value();
    pq_tmp.push(to_insert);

    bool found_dst = false;
    int maximum = INT_MAX;
    while (!pq_tmp.empty()) {
      pq_path point = pq_tmp.top();
      pq_tmp.pop();

      if (point.lm == dest) {
        maximum = min(point.value(), maximum);
        found_dst = true;
        if (!EN_BARRIERS) {
          break;
        }
      }

      if (LM_tmp[point.lm] < point.value()) continue; // we have found better and this is obsolete

      for (set<edgeLM>::iterator it=edgesV_l[point.lm].begin(); it!=edgesV_l[point.lm].end(); it++) {
        int en_barr = it->en;
        int goesTo = it->goesTo(point.lm);
        if (found_dst && (maximum < en_barr)) continue; // we dont want higher energy (really dirty programming :/)

        // insert new elements into pq (if they are better)
        pq_path to_insert(goesTo, point.dist+1, max(point.en_barr, en_barr));
        if (LM_tmp[goesTo] <= to_insert.value()) continue;
        pq_tmp.push(to_insert);
        LM_tmp[goesTo] = to_insert.value();
      }
    }

    // output
    set<int> LM_out;
    LM_out.insert(dest);
    LM_out.insert(src);
    set<edgeLM> edge_out;

      // backward pass -- find all paths with same dist/same en_barrier
    if (EN_BARRIERS) {
      // collect all paths (NP-complete)
      vector<SimplePath> paths = ConstructAllPaths(src, dest, max_length, LM_tmp[dest]);

      int en_barr;
      if (paths.size()>0) en_barr = paths[0].max_energy;
      else fprintf(stderr, "WARNING: Cannot reach %d from %d!\n", src+1, dest+1);

      for (unsigned int i=0; i<paths.size(); i++) {
        // stopping condition:
        if (paths[i].max_energy != en_barr) break;

        // print them?
        if (debug) paths[i].Print(true);

        // output them
        for (unsigned int j=0; j<paths[i].points.size(); j++) {
          LM_out.insert(paths[i].points[j]);
          if (j>0) edge_out.insert(edgeLM(paths[i].points[j-1], paths[i].points[j], paths[i].energies[j-1]));
        }
      }

    } else {
      int max = LM_tmp[dest];
      if (max == INT_MAX) {
        fprintf(stderr, "WARNING: Cannot reach %d from %d!\n", src+1, dest+1);
      } else {
        queue<pq_path> que;
        que.push(pq_path(dest, max, INT_MAX));

        while (!que.empty()) {
          pq_path point = que.front();
          que.pop();

          for (set<edgeLM>::iterator it=edgesV_l[point.lm].begin(); it!=edgesV_l[point.lm].end(); it++) {
            int goesTo = it->goesTo(point.lm);

            // if this point is on shortest path:
            if (LM_tmp[goesTo] == point.dist-1) {
              LM_out.insert(goesTo);
              edge_out.insert(*it);

              if (point.dist>1) que.push(pq_path(goesTo, point.dist-1, INT_MAX));
            } // else nothing
          }
        }
      }
    }

    // nodes
    for (set<int>::iterator it=LM_out.begin(); it!=LM_out.end(); it++) {
      if ((*it) == dest || (*it) == src) {
        fprintf(dot, "\"%d\" [label=\"%d\"]\n", (*it)+1, (*it)+1);
      } else fprintf(dot, "\"%d\" [label=\"%d\", color=\"0.0 0.0 %.1f\", fontcolor=\"0.0 0.0 %.1f\"]\n", (*it)+1, (*it)+1, color, color);
    }
    fprintf(dot, "\n");

    // edges:
    for (set<edgeLM>::iterator it=edge_out.begin(); it!=edge_out.end(); it++) {
      fprintf(dot, "\"%d\" -- \"%d\" [label=\"%.2f\", color=\"0.0 0.0 %.1f\", fontcolor=\"0.0 0.0 %.1f\"]\n", (it->i)+1, (it->j)+1, it->en/100.0, color, color);
    }

    fprintf(dot, "}\n");
  }

  fclose(dot);

  // start neato/dot:
  char syst[200];
  sprintf(syst, "%s -Tps < %s > %s", (dot_prog ? "dot" : "neato"), filename, file_print);
  system(syst);
  //printf("%s returned %d", syst, res);
}

void DSU::PrintMatrix(char *filename)
{
  FILE *energies;
  energies = fopen(filename, "w");
  if (energies) {
    // create matrix
    vector<vector<float> > matrix;
    matrix.resize(LM.size());
    for (unsigned int i=0; i<LM.size(); i++) {
      matrix[i].resize(LM.size(), INFINITY);
    }

    // fill it with edges:
    for (set<edgeLM>::iterator it=edges_l.begin(); it!=edges_l.end(); it++) {
      matrix[it->i][it->j] = it->en/100.0;
      matrix[it->j][it->i] = it->en/100.0;
    }

    // resolve i==i
    for (unsigned int i=0; i<matrix.size(); i++) {
      matrix[i][i] = LM[i].energy/100.0;
    }

    // print
    for (unsigned int i=0; i<matrix.size(); i++) {
      for (unsigned int j=0; j<matrix[i].size(); j++) {
        fprintf(energies, "%6.2g ", matrix[i][j]);
      }
      fprintf(energies, "\n");
    }
  }
  fclose(energies);
}

vector<SimplePath> DSU::ConstructAllPaths(int source, int dest, int max_length, int threshold)
{
  vector<SimplePath> paths;

  SimplePath path;
  path.AddLast(source, INT_MIN);

  // construct all path recursively
  ConstructPath(paths, path, dest, max_length-1, threshold);

  // score them
  for (unsigned int i=0; i<paths.size(); i++) {
    paths[i].Score();
  }

  // sort
  sort(paths.begin(), paths.end());

  return paths;
}

void DSU::ConstructPath(vector<SimplePath> &paths, SimplePath &path, int dest, int max_length, int threshold)
{
  int num = path.points[path.points.size()-1];

  if (num == dest) {
    paths.push_back(path);
    //if (paths.size()%100==1) printf("found paths: %d\n", (int)paths.size());
    return ;
  }

  if (max_length == 0) return;

  // all edges
  for (set<edgeLM>::iterator it=edgesV_l[num].begin(); it!=edgesV_l[num].end(); it++) {
    int goesTo = it->goesTo(num);
    if (!path.ContainsNode(goesTo) && it->en <= threshold) {
      path.AddLast(goesTo, it->en);
      ConstructPath(paths, path, dest, max_length-1, threshold);
      path.RemoveLast();
    }
  }
}

void DSU::PrintLinkCP()
{
  vector<Component> comps;
  vector<int> LM_tmp(LM.size(), -1);

  // find components:
  for (unsigned int i=0; i<LM.size(); i++) {
    if (LM_tmp[i]==-1) {
      Component cmp;
      Color(i, comps.size(), cmp, LM_tmp);
      sort(cmp.LMs.begin(), cmp.LMs.end());
      comps.push_back(cmp);
    }
  }

  // print info about them:
  printf("number of components: %d:\n", (int)comps.size());
  for (unsigned int i=0; i<comps.size(); i++) {
    printf("minimal: %4d (%7.2f), max_en: %7.2f :", comps[i].min_lm+1, comps[i].min_energy/100.0, comps[i].max_energy/100.0);
    for (unsigned int j=0; j<comps[i].LMs.size(); j++) {
      printf(" %4d", comps[i].LMs[j]+1);
    }
    printf("\n");
  }
}

void DSU::Color(int lm, int color, Component &cmp, vector<int> &LM_tmp)
{
  // add to component and colour the node
  cmp.AddLM(lm, LM[lm].energy);
  LM_tmp[lm] = color;

  // proceed to non-coloured
  for (set<edgeLM>::iterator it=edgesV_l[lm].begin(); it!=edgesV_l[lm].end(); it++) {
    int goesTo = it->goesTo(lm);
    if (LM_tmp[goesTo] == -1) {
      Color(goesTo, color, cmp, LM_tmp);
    }
  }
}

SimplePath::SimplePath()
{
  closed = false;
  max_energy = INT_MIN;
}

void SimplePath::Close()
{
  closed = true;
}

void SimplePath::Score()
{
  for (unsigned int i=0; i<energies.size(); i++) {
    max_energy = max(max_energy, energies[i]);
  }
}

void SimplePath::AddLast(int num, int energy)
{
  points.push_back(num);
  if (energy > INT_MIN) energies.push_back(energy);
  points_map.insert(make_pair(num, points.size()));
}

void SimplePath::RemoveLast()
{
  points_map.erase(points[points.size()-1]);
  points.pop_back();
  energies.pop_back();
}

SimplePath::SimplePath(const SimplePath &path)
{
  points.assign(path.points.begin(), path.points.end());
  energies.assign(path.energies.begin(), path.energies.end());
  points_map.insert(path.points_map.begin(), path.points_map.end());
  closed = path.closed;
  max_energy = path.max_energy;
}

void SimplePath::AddPoint(int num, int energy)
{
  int h;
  if ((h = FindNode(num)) != -1) {
    points.erase(points.begin()+h+1, points.end());
  } else {
    points.push_back(num);
    points_map.insert(make_pair(num, points.size()));
    max_energy = max(max_energy, energy);
  }
}

int SimplePath::FindNode(int num)
{
  map<int, int>::iterator it;
  if ((it=points_map.find(num))==points_map.end()) return -1;
  else return it->second;
}

bool SimplePath::ContainsNode(int num)
{
  return (bool) points_map.count(num);
}

void SimplePath::Print(bool whole_path, bool force_print, FILE *out)
{
  if (!force_print && !closed) return;
  //fprintf(out, "(%8.4f) ", simple_prob);
  fprintf(out, "(%8.2f) ", max_energy/100.0);
  fprintf(out, "%4d:", (int)points.size());

  if (whole_path) {
    for (unsigned int i=0; i<points.size(); i++) {
      fprintf(out, "%4d ", points[i]+1);
    }
  }

  fprintf(out, "\n");
}
