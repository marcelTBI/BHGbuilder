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
#include "hash_util.h"

#include <algorithm>

using namespace std;

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

      // insert new LM
      vertex_l[struc] = LM.size();
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

  for (unsigned int i=0; i<saddles.size(); i++) {
    if (saddles[i].structure) free(saddles[i].structure);
    if (saddles[i].str_ch) free(saddles[i].str_ch);
  }

  // always should be empty
  for (map<pq_entry, RNAstruc, pq_setcomp>::iterator it=UBlist.begin(); it!=UBlist.end(); it++) {
    if (it->second.structure) free(it->second.structure);
    if (it->second.str_ch) free(it->second.str_ch);
  }

  /*for (unsigned int i=0; i<UBoutput.size(); i++) {
    if (UBoutput[i].first.structure) free(UBoutput[i].first.structure);
    if (UBoutput[i].first.str_ch) free(UBoutput[i].first.str_ch);
  }*/
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
  RNAstruc tmp;
  tmp.structure = str_par;
  tmp.energy = en_par;
  map<RNAstruc, int>::iterator it;
  if ((it = vertex_l.find(tmp))!=vertex_l.end()) {
    return it->second;
  }
  return -1;

  //fprintf(stderr, "%d %s\n", en_par, pt_to_str(str_par).c_str());
  /*for (unsigned int i=0; i<LM.size(); i++) {
    if (LM[i].energy == en_par && str_eq(str_par, LM[i].structure)) return i;
  }
  return -1;  // old version*/
}

bool DSU::InsertUB(int i, int j, int energy_par, short *saddle_par, bool outer, bool debug)
{
  pq_entry pq(i, j, 0);
  map<pq_entry, RNAstruc, pq_setcomp>::iterator it = UBlist.find(pq);
  if (it==UBlist.end()) {
    RNAstruc saddle;
    saddle.energy = energy_par;
    saddle.structure = saddle_par;
    saddle.str_ch = NULL;
    saddle.type = (outer?NOT_SURE:DIRECT);
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
      it->second.type = (outer?NOT_SURE:DIRECT);
      return true;
    }
  }
  free(saddle_par);
  return false;
}

int DSU::ComputeUB(int maxkeep, int num_threshold, bool outer, bool noLP, bool shifts, bool debug)
{
  int dbg_count = 0;
  int cnt = 0;
  int hd_threshold = INT_MAX;

  // go through all pairs in queue
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
      int tmp_en = move_deepest(seq, tmp_str, s0, s1, 0, shifts, noLP);

      // do the stuff if we have 2 structs and they are not equal
      if (last && !str_eq(last_str, tmp_str)) {
        // not equal LM - we can update something in UBlist
          // find LM num:
        int num1 = (last_num!=-1?last_num:FindNum(last_en, last_str));
        int num2 = FindNum(tmp_en, tmp_str);

        // update UBlist
        if (num1==-1 || num2==-1) {
          if (num2==-1) {
            if (gl_maxen < tmp_en) {
              fprintf(stderr, "exceeds en.: %s %6.2f\n", pt_to_str(tmp_str).c_str(), tmp_en/100.0);
              num2 = AddLMtoDSU(tmp_str, tmp_en, hd_threshold, EE_DSU, debug);
            } else {
              fprintf(stderr, "cannot find: %s %6.2f\n", pt_to_str(tmp_str).c_str(), tmp_en/100.0);
              // add to list of minima and count with them later...
              num2 = AddLMtoDSU(tmp_str, tmp_en, hd_threshold, NORMAL, debug);
            }
          }
        }
        // again check if we can add better saddle
        if (num1!=-1 && num2!=-1) {
          // store (maybe) better saddle to UB
          int en_tmp = en_fltoi(max(last->en, tmp->en));
          short *saddle = (last->en > tmp->en ? make_pair_table(last->s) : make_pair_table(tmp->s));
          InsertUB(num1, num2, en_tmp, saddle, false, debug);
        }

        // change last_num
        last_num = num2;
      }

      // move one next
      if (last_str) free(last_str);
      last_en = tmp_en;
      last_str = tmp_str;
      last = tmp;
      tmp++;
    } // crawling path

    // insert saddle between outer structures
    if (outer) InsertUB(pq.i, pq.j, en_fltoi(max_energy), make_pair_table(max_path->s), true, debug);

    // free stuff
    if (last_str) free(last_str);
    free_path(path);
  } // all doing while

  // now just resort UBlist to something sorted according energy
  UBoutput.reserve(UBlist.size());
  for (map<pq_entry, RNAstruc, pq_setcomp>::iterator it=UBlist.begin(); it!=UBlist.end(); it++) {
    if (it->second.str_ch) free(it->second.str_ch);
    it->second.str_ch = pt_to_char(it->second.structure);
    UBoutput.push_back(make_pair(it->second, it->first));
  }
  sort(UBoutput.begin(), UBoutput.end());
  UBlist.clear();

  // check if everything has been found:
  if (UBoutput.size() != (LM.size()*(LM.size()-1))/2) {
    fprintf(stderr, "All connections have not been found: %d/%d found (%d missing)\n", (int)UBoutput.size(), (int)(LM.size()*(LM.size()-1))/2, (int)((LM.size()*(LM.size()-1))/2 - UBoutput.size()));
  }

  return 0;
}

int DSU::AddLMtoDSU(short *tmp_str, int tmp_en, int hd_threshold, int type, bool debug)
{
  RNAstruc rna;
  rna.energy = tmp_en;
  rna.structure = allocopy(tmp_str);
  rna.str_ch = pt_to_char(tmp_str);
  rna.type = type;

  // add pairs to TBDlist
  if (type == NORMAL) {
    for (unsigned int i=0; i<LM.size(); i++) {
      int hd = HammingDist(LM[i].structure, tmp_str);
      if (hd < hd_threshold) {
        TBDlist.push(pq_entry(i, LM.size(), hd));
        if (debug) {
          fprintf(stderr, "%s insert que\n%s (%4d,%4d) %3d\n", pt_to_str(LM[i].structure).c_str(), pt_to_str(rna.structure).c_str(), i, (int)LM.size(), hd);
        }
      }
    }
  }

  // insert LM and return its number
  LM.push_back(rna);
  vertex_l[rna]=LM.size()-1;
  return LM.size()-1;
}

void DSU::PrintUBoutput()
{
  printf("                 %s\n", seq);
  for (unsigned int i=0; i<UBoutput.size(); i++) {
    printf("%4d (%4d,%4d) %s %6.2f\n", i+1, UBoutput[i].second.i+1, UBoutput[i].second.j+1, pt_to_str(UBoutput[i].first.structure).c_str(), UBoutput[i].first.energy/100.0);
  }
}

int DSU::LinkCP(Opt opt, bool debug)
{
  //edgesV_ls.resize(LM.size());
  edgesV_l.resize(LM.size());

  // create lm-saddle and lm-lm edges
  for (unsigned int i=0; i<UBoutput.size(); i++) {
    RNAstruc stru = UBoutput[i].first;
    pq_entry pq = UBoutput[i].second;

    // flood them up!
    int res = FloodUp(LM[pq.i], LM[pq.j], stru, opt, debug);
    if (res == 1 || res == 2) {  // we have found better saddle (1) or reached threshold (2) so better saddle is not possible
      stru.type = LDIRECT; // so our saddle is for sure lowest direct saddle
    }

    // update vertex/edge sets
    map<RNAstruc, int>::iterator it;
    int saddle_num;
    if ((it = vertex_s.find(stru)) != vertex_s.end()) {
      saddle_num = it->second;
      stru.freeMEM();

      // also add missing edges:
      set<int> miss;
      for (set<edgeLM>::iterator it2 = edges_ls.begin(); it2!=edges_ls.end(); it2++) {
        if (it2->j == saddle_num) {
          miss.insert(it2->i);
        }
      }
      miss.insert(pq.i);
      miss.insert(pq.j);

      // edges ll
      for (set<int>::iterator it2 = miss.begin(); it2!=miss.end(); it2++) {
        set<int>::iterator it3 = it2;
        for (it3++; it3!=miss.end(); it3++) {
          edgeLM e(*it2, *it3, true, true);
          e.AddSaddle(stru.energy, saddle_num);
          edges_l.insert(e);

          edgesV_l[*it2].insert(e);
          edgesV_l[*it3].insert(e);
        }
      }

      //fprintf(stderr, "saddle_num = %d\n", saddle_num);
    } else {
      saddle_num = saddles.size();
      // update vertex
      vertex_s.insert(make_pair(stru, saddle_num));
      saddles.push_back(stru);

      // edges ll
      edgeLM e(pq.i, pq.j, true, true);
      e.AddSaddle(stru.energy, saddle_num);
      edges_l.insert(e);

      edgesV_l[pq.i].insert(e);
      edgesV_l[pq.j].insert(e);
    }

    // edges ls
    edges_ls.insert(edgeLM(pq.i, saddle_num, true, false));
    edges_ls.insert(edgeLM(pq.j, saddle_num, true, false));
  }

  // create saddle-saddle edges
    // create saddle set list
  if (opt.saddle_conn) {
    set<std::pair<int, int> > saddle_pairs;
    int last_lm = -1;
    vector<int> tmp;
    for (set<edgeLM>::iterator it = edges_ls.begin(); it!=edges_ls.end(); it++) {

      if (it->i == last_lm) {
        // add every one into saddle_pairs.
        for (unsigned int i=0; i<tmp.size(); i++) {
          saddle_pairs.insert(make_pair(min(tmp[i], it->j), max(tmp[i], it->j)));
        }
      } else {
        // clear
        tmp.clear();
      }
      // add next one
      tmp.push_back(it->j);
      last_lm = it->i;
    }
      // process this list
    for (set<std::pair<int, int> >::iterator it=saddle_pairs.begin(); it!=saddle_pairs.end(); it++) {
      if (debug) {
        fprintf(stderr, "saddles to check: %4d %4d", it->first, it->second);
      }
      RNAstruc &first = saddles[it->first];
      RNAstruc &second = saddles[it->second];
      if (first.energy > second.energy) {
        RNAstruc &tmp = first;
        first = second;
        second = tmp;
      }
      // include them
      if (FloodSaddle(first, second, opt, debug)) {
        edgeLM e(it->first, it->second, false, false);
        edges_s.insert(e);
      }
    }
  }

  return 0;
}

void DSU::PrintDot(char *filename, bool dot_prog, bool print, char *file_print, bool visual)
{
  // landmap not supported yet
  // start find_union stuff
  UF_set connected;

  int color = 180;
  //open file
  FILE *dot;
  dot = fopen(filename, "w");
  if (dot) {
    fprintf(dot, "Graph G {\n\tnode [width=0.1, height=0.1, shape=circle];\n");
    //nodes LM:
    for (unsigned int i=0; i<LM.size(); i++) {
      switch (LM[i].type) {
        case NORMAL: fprintf(dot, "\"%d\" [label=\"%d\"]\n", i+1, i+1); break;
        case EE_DSU: fprintf(dot, "\"%d\" [label=\"%d\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, rgb(0, 0, 255), rgb(0, 0, 255)); break;
        case EE_COMP: fprintf(dot, "\"%d\" [label=\"%d\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, rgb(255, 0, 0), rgb(255, 0, 0)); break;
      }
    }
    fprintf(dot, "\n");

    // visualisation option (not finished -- currently it prints out only without saddles)
    if (visual) {
      connected.enlarge_parent(LM.size());
      set<edgeLM, edgeLM_compen> tmp;
      tmp.insert(edges_l.begin(), edges_l.end());
      for (set<edgeLM>::iterator it=tmp.begin(); it!=tmp.end(); it++) {
        if (!connected.joint(it->i, it->j)) {
          fprintf(dot, "\"%d\" -- \"%d\" [label=\"%.2f\"]\n", (it->i)+1, (it->j)+1, it->en/100.0);
          connected.union_set(it->i, it->j);
        }
      }
      fprintf(dot, "\n");

    } else {

      //nodes saddle:
      for (unsigned int i=0; i<saddles.size(); i++) {
        switch (saddles[i].type) {
          case DIRECT: fprintf(dot, "\"S%d\" [label=\"S%d\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, rgb(color, color, color), rgb(color, color, color)); break;
          case LDIRECT: fprintf(dot, "\"S%d\" [label=\"S%d\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, rgb(color+30, color+30, color+30), rgb(color+30, color+30, color+30)); break;
          case NOT_SURE: fprintf(dot, "\"S%d\" [label=\"S%d\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, rgb(color-30, color-30, color-30), rgb(color-30, color-30, color-30)); break;
          case COMP: fprintf(dot, "\"S%d\" [label=\"S%d\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, rgb(255, color, color), rgb(255, color, color)); break;
        }
      }
      fprintf(dot, "\n");
      // edges l-l
      for (set<edgeLM>::iterator it=edges_l.begin(); it!=edges_l.end(); it++) {
        fprintf(dot, "\"%d\" -- \"%d\" [label=\"%.2f\", color=\"%s\", fontcolor=\"%s\"]\n", (it->i)+1, (it->j)+1, it->en/100.0, (it->component?rgb(255, 0, 0):rgb(0, 0, 0)), (it->component?rgb(255, 0, 0):rgb(0, 0, 0)));
      }
      fprintf(dot, "\n");
      // edges l-s
      for (set<edgeLM>::iterator it=edges_ls.begin(); it!=edges_ls.end(); it++) {
        fprintf(dot, "\"%d\" -- \"S%d\" [color=\"%s\", fontcolor=\"%s\"]\n", (it->i)+1, (it->j)+1, (it->component?rgb(255, color, color):rgb(color, color, color)), (it->component?rgb(255, color, color):rgb(color, color, color)));
      }

      fprintf(dot, "\n");
      // edges s-s
      for (set<edgeLM>::iterator it=edges_s.begin(); it!=edges_s.end(); it++) {
        fprintf(dot, "\"S%d\" -- \"S%d\" [color=\"%s\", fontcolor=\"%s\"]\n",(it->i)+1, (it->j)+1, rgb(color, color, color),rgb(color, color, color));
      }
    }
    fprintf(dot, "\n}\n");
  }

  fclose(dot);

  // start neato/dot:
  if (dot && print && file_print) {
    char syst[200];
    sprintf(syst, "%s -Tps < %s > %s", (dot_prog ? "dot" : "neato"), filename, file_print);
    int res = system(syst);
    //printf("%s returned %d\n", syst, res);
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
    vector<int> LM_tmp(LM.size(), INT_MAX);
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
          if (j>0) {
            edgeLM e(paths[i].points[j-1], paths[i].points[j], true, true);
            e.AddSaddle(paths[i].energies[j-1], -1);
            edge_out.insert(e);
          }
        }
      }

    } else {
      // find the shortest path
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

void DSU::FillComps()
{
  comps.clear();
  LM_to_comp.clear();
  saddle_to_comp.clear();
  vector<int> LM_tmp(LM.size(), -1);
  vector<int> sadd_tmp(saddles.size(), -1);

  // find components:
  for (unsigned int i=0; i<LM.size(); i++) {
    if (LM_tmp[i]==-1) {
      Component cmp;
      Color(i, comps.size(), cmp, LM_tmp, sadd_tmp);
      sort(cmp.LMs.begin(), cmp.LMs.end());
      sort(cmp.saddles.begin(), cmp.saddles.end());
      // assign maps
      for (unsigned int j=0; j<cmp.LMs.size(); j++) {
        LM_to_comp[cmp.LMs[j]] = comps.size();
      }
      for (unsigned int j=0; j<cmp.saddles.size(); j++) {
        saddle_to_comp[cmp.saddles[j]] = comps.size();
      }
      comps.push_back(cmp);
    }
  }
}

int DSU::AddConnection(int num1, int num2, int energy, short *saddle, UF_set &connected, vector<vector<RNAstruc2*> > &connections) {

  int comp1 = LM_to_comp[num1];
  int comp2 = LM_to_comp[num2];

  // we can connect components:
  if (comp1 != comp2) {

    if (connections[comp1][comp2] == NULL || connections[comp1][comp2]->energy > energy) {
      RNAstruc2 *tmpRNA = (RNAstruc2 *) space(sizeof(RNAstruc2));
      tmpRNA->conn1 = num1;
      tmpRNA->conn2 = num2;
      tmpRNA->energy = energy;
      tmpRNA->structure = saddle;
      tmpRNA->str_ch = NULL;

      if (connections[comp1][comp2]) {
        free(connections[comp1][comp2]->structure);
        free(connections[comp1][comp2]);
      }
      // assign new one
      connections[comp1][comp2] = tmpRNA;

      // add to UnionFindSet
      connected.union_set(comp1, comp2);
      return 1;
    }
  }
  free(saddle);
  return 0;
}


// small struct for use in ConnectComps function
struct que_tmp {
  char *str_ch1;
  char *str_ch2;

  int hd;

  // debug info
  int i,j;
  int ci, cj;

  bool operator<(const que_tmp &second) const {
    return hd>second.hd;
  };

  que_tmp(char *ch1, char *ch2, int hd) {
    this->hd = hd;
    str_ch1 = ch1;
    str_ch2 = ch2;
  }

};

void DSU::ConnectComps(int maxkeep, bool debug)
{
  if (comps.size() == 0) FillComps();

  // just debug
  LM_to_comp[-1]=-1;

  // create queue of saddle pairs
  priority_queue<que_tmp> queue;
  for (unsigned int i=0; i<comps.size(); i++) {
    for (unsigned int j=i+1; j<comps.size(); j++) {
      const RNAstruc &first = (comps[i].max_saddle==-1?LM[comps[i].LMs[0]]:saddles[comps[i].max_saddle]);
      const RNAstruc &second = (comps[j].max_saddle==-1?LM[comps[j].LMs[0]]:saddles[comps[j].max_saddle]);

      que_tmp tmp(first.str_ch, second.str_ch, HammingDist(first.structure, second.structure));
      tmp.i = comps[i].LMs[0];
      tmp.j = comps[j].LMs[0];
      tmp.ci = i;
      tmp.cj = j;
      queue.push(tmp);
    }
  }

  int cnt = 0;
  int dbg_count = 0;

  // here we will store output:
  UF_set connected;
  connected.enlarge_parent(comps.size());
  vector<vector<RNAstruc2*> > connections(comps.size(), vector<RNAstruc2*>(comps.size(), NULL));

  // go throught the queue - take pair *
  while (!queue.empty() && !connected.connected_all()) {
    // get pair
    que_tmp top = queue.top();
    queue.pop();

    // apply threshold ??
    cnt++;

    // get path
    if (debug) fprintf(stderr, "path between (LM: %d, %d) (comp: %d, %d) (hd: %d):\n", top.i, top.j, top.ci, top.cj, top.hd);
    path_t *path = get_path(seq, top.str_ch1, top.str_ch2, maxkeep);

    /*// variables for outer insertion
    double max_energy= -1e8;
    path_t *max_path = path;*/

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
      if (debug) fprintf(stderr, "%s %6.2f ", tmp->s, tmp->en);

      /*// update max_energy
      if (max_energy < tmp->en) {
        max_energy = tmp->en;
        max_path = tmp;
      }*/
      // find adaptive walk
      short *tmp_str = make_pair_table(tmp->s);
      //int tmp_en = move_rand(seq, tmp_str, s0, s1, 0);
      int tmp_en = move_deepest(seq, tmp_str, s0, s1, 0, 0, 0);

      if (debug) fprintf(stderr, "%s %6.2f (LM: %4d) (comp: %4d)\n", pt_to_str(tmp_str).c_str(), tmp_en/100.0, FindNum(tmp_en, tmp_str)+1, LM_to_comp[FindNum(tmp_en, tmp_str)]);

      // do the stuff if we have 2 structs and they are not equal
      if (last && !str_eq(last_str, tmp_str)) {
        // not equal LM - we can update something in UBlist
          // find LM num:
        int num1 = (last_num!=-1?last_num:FindNum(last_en, last_str));
        int num2 = FindNum(tmp_en, tmp_str);

        // check if ok
        if (num1==-1 || num2==-1) {
          if (num2==-1) {
            if (gl_maxen < tmp_en) {
              fprintf(stderr, "exceeds en.: %s %6.2f\n", pt_to_str(tmp_str).c_str(), tmp_en/100.0);
              num2 = AddLMtoComp(tmp_str, tmp_en, debug, connected, connections);
            } else {
              fprintf(stderr, "WARNING! cannot find: %s %6.2f\n", pt_to_str(tmp_str).c_str(), tmp_en/100.0);
            }
          }
        }
        // again check
        if (num1!=-1 && num2!=-1) {
          // find component
          if (LM_to_comp.count(num1)==0 || LM_to_comp.count(num2)==0) {
            fprintf(stderr, "Not found component of LM num %d or %d (HUGE ERROR!)\n", num1, num2);
          } else {
            // add connection
            short *saddle = (last->en > tmp->en ? make_pair_table(last->s) : make_pair_table(tmp->s));
            AddConnection(num1, num2, en_fltoi(max(last->en, tmp->en)), saddle, connected, connections);
          }
        }

        // update last_num
        last_num = num2;
      }

      // move one next
      if (last_str) free(last_str);
      last_en = tmp_en;
      last_str = tmp_str;
      last = tmp;
      tmp++;
    } // crawling path

    // insert saddle between outer structures
    //short *saddle = make_pair_table(max_path->s);
    //AddConnection(?, ?, en_fltoi(max_energy), saddle, connected, connections);

    // free stuff
    if (last_str) free(last_str);
    free_path(path);
  }
  // write connections to graph:
  edgesV_l.resize(LM.size());
  for (unsigned int i=0; i<connections.size(); i++) {
    for (unsigned int j=0; j<connections[i].size(); j++) {

      if (!connections[i][j]) continue;
      // saddle
      RNAstruc2* saddle = connections[i][j];
      saddle->type = COMP;
      vertex_s.insert(make_pair(*saddle, saddles.size()));
      saddle->str_ch = pt_to_char(saddle->structure);
      saddles.push_back(*saddle);

      // edge l-l
      edgeLM e(saddle->conn1, saddle->conn2, true, true);
      e.AddSaddle(saddle->energy, saddles.size()-1);
      e.MarkComp();
      edges_l.insert(e);
      edgesV_l[saddle->conn1].insert(e);
      edgesV_l[saddle->conn2].insert(e);
      // edge l-s
      edgeLM e2(saddle->conn1, saddles.size()-1, true, false);
      e2.MarkComp();
      edges_ls.insert(e2);
      e2.i = saddle->conn2;
      edges_ls.insert(e2);

      free(saddle);
    }
  }

  //FillComps();
}

int DSU::AddLMtoComp(short *structure, int energy, bool debug, UF_set &connected, vector<vector<RNAstruc2*> > &connections)
{
  // resize connected and connections (we have new component...)
  int numComp = connected.size();
  connected.enlarge_parent();
  connections.resize(numComp+1);
  for (unsigned int i=0; i< connections.size(); i++) {
    connections[i].resize(numComp+1, NULL);
  }

  // add new LM to LMs
  int numLM = AddLMtoDSU(structure, energy, 0, EE_COMP, debug);


  // add new component to map
  LM_to_comp[numLM] = numComp;

  return numLM;
}

void DSU::PrintLinkCP(bool full)
{
  if (comps.size() == 0 || full) FillComps();

  // print info about comps:
  for (unsigned int i=0; i<comps.size(); i++) {
    printf("%4d %4d (%7.2f)", i, comps[i].min_lm+1, comps[i].min_energy/100.0);
    if (comps[i].max_saddle!=-1)  printf(" %4dS (%7.2f) :", comps[i].max_saddle+1, comps[i].max_energy/100.0);
    else                          printf("                 :");
    for (unsigned int j=0; j<comps[i].LMs.size(); j++) {
      printf(" %4d", comps[i].LMs[j]+1);
    }
    printf("      ");
    for (unsigned int j=0; j<comps[i].saddles.size(); j++) {
      printf(" %3dS", comps[i].saddles[j]+1);
    }
    printf("\n");
  }
  printf("\n");

  // full listing of structures...
  if (full) {
    //printf("Local minima (%4d):\n", (int)LM.size());
    for (unsigned int i=0; i<LM.size(); i++) {
      char type[][10] = { "NORMAL", "EE_DSU", "EE_COMP"};
      printf("%4d  %s (%7.2f) %8s\n", i+1, LM[i].str_ch, LM[i].energy/100.0, type[LM[i].type]);
    }
    //printf("Saddles (%4d):\n", (int)saddles.size());
    printf("\n");
    // collect saddle info
    vector<set<int> > saddle_connLM (saddles.size());
    vector<set<int> > saddle_connSadd (saddles.size());
    for (set<edgeLM>::iterator it=edges_ls.begin(); it!=edges_ls.end(); it++) {
      saddle_connLM[it->j].insert(it->i);
    }
    for (set<edgeLM>::iterator it=edges_s.begin(); it!=edges_s.end(); it++) {
      saddle_connSadd[it->j].insert(it->i);
      saddle_connSadd[it->i].insert(it->j);
    }
    for (unsigned int i=0; i<saddles.size(); i++) {
      char type[][10] = { "DIRECT", "LDIRECT", "NOT_SURE", "COMP" };
      printf("%4dS %s (%7.2f) %8s", i+1, saddles[i].str_ch, saddles[i].energy/100.0, type[saddles[i].type]);
      for (set<int>::iterator it=saddle_connLM[i].begin(); it!=saddle_connLM[i].end(); it++) printf(" %4d", (*it)+1);
      for (set<int>::iterator it=saddle_connSadd[i].begin(); it!=saddle_connSadd[i].end(); it++) printf(" %3dS", (*it)+1);
      printf("\n");
    }
  }

}

void DSU::Color(int lm, int color, Component &cmp, vector<int> &LM_tmp, vector<int> &sadd_tmp)
{
  // add to component and colour the node
  cmp.AddLM(lm, LM[lm].energy);
  LM_tmp[lm] = color;

  // proceed to non-coloured
  for (set<edgeLM>::iterator it=edgesV_l[lm].begin(); it!=edgesV_l[lm].end(); it++) {
    int goesTo = it->goesTo(lm);

    // info about saddle:
    if (it->sadd != -1 && sadd_tmp[it->sadd] == -1) {
      sadd_tmp[it->sadd] = color;
      cmp.AddSadd(it->sadd, it->en);
    }
    // info about LM
    if (LM_tmp[goesTo] == -1) {
      // colour it!
      Color(goesTo, color, cmp, LM_tmp, sadd_tmp);
    }
  }
}
