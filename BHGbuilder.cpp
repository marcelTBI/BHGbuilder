#include <stdio.h>
#include <string.h>
#include <stdlib.h>

extern "C" {
  #include "utils.h"
  #include "fold_vars.h"
  #include "pair_mat.h"

  #include "fold.h"
  //#include "findpath.h"
}

#include "move_set_pk.h"
#include "findpath_pk.h"

#include "BHGbuilder.h"
#include "hash_util.h"

#include <algorithm>

using namespace std;

DSU::DSU(FILE *input, bool noLP, bool shifts, bool pknots, int time_max, int max_lm, bool just_read, bool debug) {

  // NULL::
  seq = NULL;
  s0 = NULL;
  s1 = NULL;
  gl_maxen = INT_MIN;

  // pseudoknots?
  this->pknots = pknots;

  // time:
  if (time_max) {
    time = clock();
  }
  stop_after = time_max;

  if (!input) return ;

  // read seq
  char *line;
  int num = 0;

  line = my_getline(input);
  char *seq2 = strtok(line, " ");
  if (!isSeq(seq2)) {
    free(line);
    line = my_getline(input);
    seq2 = strtok(line, " ");
    if (!isSeq(seq2)) {
      free(line);
      return ;
    }
  }
  seq = (char*) malloc((strlen(seq2)+1)*sizeof(char));
  strcpy(seq, seq2);
  free(line);

  //init
  make_pair_matrix();
  s0 = encode_sequence(seq, 0);
  s1 = encode_sequence(seq, 1);

  // if we read .dsu file
  if (just_read) {
        //read structs
    line = my_getline(input);
    bool saddle_reading = false;
    while (line) {
      bool empty_line = false;
      short *tmp = NULL;
      float energy_tmp;
      int type;

      char *p;
      for (int i=0; i<4; i++) {
        p = strtok(i==0?line:NULL, " ");
        switch (i) {
        case 1:
          if (!p) {empty_line = true; break;}
          if (isStruct(p)) {
            if (tmp) free(tmp); // only one struct per line!
            tmp = pknots?make_pair_table_PK(p):make_pair_table(p);
            break;
          }
        case 2:
          sscanf(p, "%f", &energy_tmp);
          break;
        case 3:
          for (int i=0; i<4; i++) {
            if (!saddle_reading) {
              if (strcmp(p, LMtype_string[i])==0) {
                type = i;
                break;
              }
            } else {
              if (strcmp(p, SDtype_string[i])==0) {
                type = i;
                break;
              }
            }
          }
          break;
        }
        if (empty_line) break;
      }

      if (empty_line) {
        saddle_reading = true;
        free(line);
        line = my_getline(input);

        edgesV_l.resize(LM.size());
        continue;
      }

      vector<int> LM_c;
      vector<int> saddle_c;
      int number;
      if (saddle_reading) {
        while ((p=strtok(NULL, " "))) {
          if (p[strlen(p)-1]=='S') {
            sscanf(p, "%dS", &number);
            saddle_c.push_back(number-1);
          } else {
            sscanf(p, "%d", &number);
            LM_c.push_back(number-1);
          }
        }
      }

      if (saddle_reading && LM_c.size()<2) {
        fprintf(stderr, "File reading error -- too few LM connected by saddle %d\n", (int)saddles.size()+1);
        exit(EXIT_FAILURE);
      }

      // add info:
      if (tmp) {

        if (!saddle_reading) {
          RNAlocmin struc;
          struc.structure = tmp;
          struc.str_ch = pt_to_chars_pk(struc.structure);
          struc.energy = en_fltoi(energy_tmp);
          struc.type = (LMtype)type;

          // insert new LM
          vertex_l[struc] = LM.size();
          LM.push_back(struc);

          gl_maxen = max(gl_maxen, struc.energy);
        } else {
          RNAsaddle struc(LM_c[0], LM_c[1], (SDtype)type);
          struc.structure = tmp;
          struc.str_ch = pt_to_chars_pk(struc.structure);
          struc.energy = en_fltoi(energy_tmp);
          // insert new sadle
          saddles.push_back(struc);

          // edge l-l
          for (int i=0; i<(int)LM_c.size(); i++) {
            for (int j=i+1; j<(int)LM_c.size(); j++) {
              edgeLL e(LM_c[i], LM_c[j], struc.energy, saddles.size()-1);
              edges_l.insert(e);
              edgesV_l[LM_c[i]].insert(e);
              edgesV_l[LM_c[j]].insert(e);
            }
          }
          // edge l-s
          for (int i=0; i<(int)LM_c.size(); i++) {
            edgeLS e2(LM_c[i], saddles.size()-1);
            edges_ls.insert(e2);
          }

          // edges s-s
          for (int i=0; i<(int)saddle_c.size(); i++) {
            edgeSS e(saddle_c[i], saddles.size()-1);
            edges_s.insert(e);
          }

        }
      } else {
        free(tmp);
      }

      free(line);
      line = my_getline(input);
      num++;
    }

  } else {
    //read structs
    line = my_getline(input);

    while (line) {
      short *tmp = NULL;
      float energy_tmp;
      bool has_energy = false;

      //if (debug) fprintf(stderr, "%s\n", line);


      char *p;
      p = strtok(line, " ");
      while (p !=NULL && (!has_energy || !tmp)) {
        // is struct?
        if (isStruct(p)) {
          if (tmp) free(tmp); // only one struct per line!
          tmp = pknots?make_pair_table_PK(p):make_pair_table(p);
          //if (debug) fprintf(stderr, "%s\n", pt_to_str_pk(tmp).c_str());

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
        // check if its true local minima
        int en;
        if (pknots) {
          Structure str(seq, tmp, s0, s1);
          en = move_gradient_pk(seq, &str, s0, s1, shifts, 0);
          copy_arr(tmp, str.str);
        } else {
          en = move_gradient(seq, tmp, s0, s1, 0, shifts, noLP);
        }

        // if we don't have it yet, add it:
        if (FindNum(en, tmp)==-1) {

          RNAlocmin struc;
          struc.structure = tmp;
          struc.str_ch = pt_to_chars_pk(struc.structure);
          struc.energy = en;
          //int en = (has_energy?en_fltoi(energy_tmp):energy_of_structure_pt(seq, tmp, s0, s1, 0));

          // insert new LM
          vertex_l[struc] = LM.size();
          LM.push_back(struc);

          gl_maxen = max(gl_maxen, en);
          if (max_lm>0 && max_lm==(int)LM.size()) {
            free(line);
            break;
          }
        } else {
          free(tmp);
        }
      } else {
        // line without struct???
        fprintf(stderr, "WARNING: line %d without struct: %s\n", num, line);
      }
      free(line);
      line = my_getline(input);
      num++;

    }
  }

  if (debug) PrintLM(stderr, false);

  // sort them!
  SortFix();
	number_lm = (int)LM.size();
	//printf("------------------------------------------------------------%d\n", number_lm);
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

  if (pknots) freeP();

  for (unsigned int i=0; i<saddles.size(); i++) {
    if (saddles[i].structure) free(saddles[i].structure);
    if (saddles[i].str_ch) free(saddles[i].str_ch);
  }
}

bool DSU::InsertUB(RNAsaddle saddle, bool debug)
{
  set<RNAsaddle, RNAsaddle_comp>::iterator it = UBlist.find(saddle);
  if (it==UBlist.end()) {
    if (debug) fprintf(stderr, "UBins: (%3d, %3d) insert saddle  %6.2f %s\n", saddle.lm1, saddle.lm2, saddle.energy/100.0, pt_to_str(saddle.structure).c_str());
    UBlist.insert(saddle);
    return true;
    //fprintf(stderr, "cannot find (UB  ): (%3d, %3d)\n", num1, num2);
  } else {
    // update it
    if (saddle.energy < it->energy) {
      if (debug) fprintf(stderr, "UBupd: (%3d, %3d) from %6.2f to %6.2f %s\n", saddle.lm1, saddle.lm2, it->energy/100.0, saddle.energy/100.0, pt_to_str(it->structure).c_str());

      if (it->structure) free(it->structure);
      UBlist.erase(it);
      UBlist.insert(saddle);
      return true;
    }
  }
  free(saddle.structure);
  return false;
}

int DSU::LinkCPLM(Opt opt, bool debug)
{
  //edgesV_ls.resize(LM.size());
  edgesV_l.resize(LM.size());

  fprintf(stderr, "Computing lm-* edges.\n");

  int trueds = 0;
  // create lm-saddle and lm-lm edges
  for (unsigned int i=0; i<saddles.size(); i++) {
    RNAsaddle stru = saddles[i];

    // flood them up!
    if (opt.flood_num>0) {
      int res = FloodUp(LM[stru.lm1], LM[stru.lm2], stru, opt, debug);
      if (res == 1 || res == 2) {  // we have found better saddle (1) or reached threshold (2) so better saddle is not possible
        stru.type = LDIRECT; // so our saddle is for sure lowest direct saddle
        trueds ++;
      }
    }

    // update vertex/edge sets
    /*map<RNAstruc, int>::iterator it;
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
    } else {*/

    // update vertex
    //vertex_s.insert(make_pair(stru, saddle_num));
    saddles[i] = stru;

    // edges ll
    edgeLL e(stru.lm1, stru.lm2, stru.energy, i);
    edges_l.insert(e);

    edgesV_l[stru.lm1].insert(e);
    edgesV_l[stru.lm2].insert(e);

    // edges ls
    edges_ls.insert(edgeLS(stru.lm1, i));
    edges_ls.insert(edgeLS(stru.lm2, i));
  }

  fprintf(stderr, "Recomputed %d(%d) edges (%d are true direct saddles)\n", (int)edges_l.size(), (int)saddles.size(), trueds);
  return 0;
}

int DSU::LinkCPsaddle(Opt opt, bool debug) {       // construct vertex and edge set from saddles (saddle to saddle)

  fprintf(stderr, "Computing saddle-saddle edges.\n");
  set<std::pair<int, int> > saddle_pairs;
  int last_lm = -1;
  vector<int> tmp;
  for (set<edgeLS>::iterator it = edges_ls.begin(); it!=edges_ls.end(); it++) {

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
    RNAsaddle &first = ((saddles[it->first].energy <= saddles[it->second].energy) ? saddles[it->first] : saddles[it->second]);
    RNAsaddle &second = ((saddles[it->first].energy > saddles[it->second].energy) ? saddles[it->first] : saddles[it->second]);
    // include them
    if (FloodSaddle(first, second, opt, debug)) {
      edgeSS e(it->first, it->second);
      edges_s.insert(e);
    }
  }
  fprintf(stderr, "Recomputed %d saddle-saddle edges (%d computations done)\n", (int)edges_s.size(), (int)saddle_pairs.size());

  return 0;
}

void DSU::PrintDot(char *filename, bool dot_prog, bool print, char *file_print, bool visual, bool print_energies)
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
      char energy[20] = "";
      if (print_energies) sprintf(energy, "\\n%.2f", LM[i].energy/100.0);
      switch (LM[i].type) {
        case NORMAL:
        case NORM_CF: fprintf(dot, "\"%d\" [label=\"%d%s\"]\n", i+1, i+1, energy); break;
        case EE_DSU: fprintf(dot, "\"%d\" [label=\"%d%s\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, energy, rgb(0, 0, 255), rgb(0, 0, 255)); break;
        case EE_COMP: fprintf(dot, "\"%d\" [label=\"%d%s\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, energy, rgb(255, 0, 0), rgb(255, 0, 0)); break;
      }
    }
    fprintf(dot, "\n");

    // visualisation option (not finished -- currently it prints out only without saddles)
    if (visual) {
      connected.enlarge_parent(LM.size());
      set<edgeLL, edgeLL_compen> tmp;
      tmp.insert(edges_l.begin(), edges_l.end());
      for (set<edgeLL>::iterator it=tmp.begin(); it!=tmp.end(); it++) {
        if (!connected.joint(it->i, it->j)) {
          fprintf(dot, "\"%d\" -- \"%d\" [label=\"%.2f\"]\n", (it->i)+1, (it->j)+1, it->en/100.0);
          connected.union_set(it->i, it->j);
        }
      }
      fprintf(dot, "\n");

    } else {

      //nodes saddle:
      for (unsigned int i=0; i<saddles.size(); i++) {
        char energy[20] = "";
        if (print_energies) sprintf(energy, "\\n%.2f", saddles[i].energy/100.0);
        switch (saddles[i].type) {
          /*case DIRECT: fprintf(dot, "\"S%d\" [label=\"S%d (%d %d)\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, saddles[i].lm1+1, saddles[i].lm2+1, rgb(color, color, color), rgb(color, color, color)); break;
          case LDIRECT: fprintf(dot, "\"S%d\" [label=\"S%d (%d %d)\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, saddles[i].lm1+1, saddles[i].lm2+1, rgb(color+30, color+30, color+30), rgb(color+30, color+30, color+30)); break;
          case NOT_SURE: fprintf(dot, "\"S%d\" [label=\"S%d (%d %d)\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, saddles[i].lm1+1, saddles[i].lm2+1, rgb(color-30, color-30, color-30), rgb(color-30, color-30, color-30)); break;
          case COMP: fprintf(dot, "\"S%d\" [label=\"S%d (%d %d)\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, saddles[i].lm1+1, saddles[i].lm2+1, rgb(255, color, color), rgb(255, color, color)); break;*/

          case DIRECT: fprintf(dot, "\"S%d\" [label=\"S%d%s\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, energy, rgb(color, color, color), rgb(color, color, color)); break;
          case LDIRECT: fprintf(dot, "\"S%d\" [label=\"S%d%s\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, energy, rgb(color+30, color+30, color+30), rgb(color+30, color+30, color+30)); break;
          case NOT_SURE: fprintf(dot, "\"S%d\" [label=\"S%d%s\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, energy, rgb(color-30, color-30, color-30), rgb(color-30, color-30, color-30)); break;
          case COMP: fprintf(dot, "\"S%d\" [label=\"S%d%s\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, energy, rgb(255, color, color), rgb(255, color, color)); break;
        }
      }
      fprintf(dot, "\n");
      // edges l-l
      for (set<edgeLL>::iterator it=edges_l.begin(); it!=edges_l.end(); it++) {
        bool component = (saddles[it->saddle].type == COMP);
        fprintf(dot, "\"%d\" -- \"%d\" [label=\"%.2f\", color=\"%s\", fontcolor=\"%s\"]\n", (it->i)+1, (it->j)+1, it->en/100.0, (component?rgb(255, 0, 0):rgb(0, 0, 0)), (component?rgb(255, 0, 0):rgb(0, 0, 0)));
      }
      fprintf(dot, "\n");
      // edges l-s
      for (set<edgeLS>::iterator it=edges_ls.begin(); it!=edges_ls.end(); it++) {
        bool component = (saddles[it->j].type == COMP);
        fprintf(dot, "\"%d\" -- \"S%d\" [color=\"%s\", fontcolor=\"%s\"]\n", (it->i)+1, (it->j)+1, (component?rgb(255, color, color):rgb(color, color, color)), (component?rgb(255, color, color):rgb(color, color, color)));
      }

      fprintf(dot, "\n");
      // edges s-s
      for (set<edgeSS>::iterator it=edges_s.begin(); it!=edges_s.end(); it++) {
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
    system(syst);
    //printf("%s returned %d\n", syst, res);
  }
}

static int EN_BARRIERS = true;

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

      for (set<edgeLL>::iterator it=edgesV_l[point.lm].begin(); it!=edgesV_l[point.lm].end(); it++) {
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
    set<edgeLL> edge_out;

      // backward pass -- find all paths with same dist/same en_barrier
    if (EN_BARRIERS) {
      // collect all paths (NP-complete)
      vector<SimplePath> paths = ConstructAllPaths(src, dest, max_length, LM_tmp[dest]);

      int en_barr;
      if (paths.size()>0) {
        en_barr = paths[0].max_energy;

        for (unsigned int i=0; i<paths.size(); i++) {
          // stopping condition:
          if (paths[i].max_energy != en_barr) break;

          // print them?
          if (debug) paths[i].Print(true);

          // output them
          for (unsigned int j=0; j<paths[i].points.size(); j++) {
            LM_out.insert(paths[i].points[j]);
            if (j>0) {
              edgeLL e(paths[i].points[j-1], paths[i].points[j], paths[i].energies[j-1], -1);
              edge_out.insert(e);
            }
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

          for (set<edgeLL>::iterator it=edgesV_l[point.lm].begin(); it!=edgesV_l[point.lm].end(); it++) {
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
    for (set<edgeLL>::iterator it=edge_out.begin(); it!=edge_out.end(); it++) {
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

vector<vector<std::pair<int, int> > > matrix;
bool generated = false;

void DSU::PrintMatrix(char *filename, bool full, char* filter_file, char type)
{
  // firest read the filter file
  if (filter_file) {

    FILE *filter = fopen(filter_file, "r");
    if (filter) {
      mapping.clear();

      char *line = my_getline(filter);
      while (line) {
        //fprintf(stderr, "working: %s\n",line);
        bool empty_line = false;
        short *tmp = NULL;
        int num;
        int energy_tmp;

        char *p;
        for (int i=0; i<3; i++) {
          p = strtok(i==0?line:NULL, " ");
          switch (i) {
          case 0:
            if (!p) {empty_line = true; break;}
            sscanf(p, "%d", &num);
            break;
          case 1:
            if (!p) {empty_line = true; break;}
            if (isStruct(p)) {
              if (tmp) free(tmp); // only one struct per line!
              tmp = pknots?make_pair_table_PK(p):make_pair_table(p);
              break;
            }
          case 2:
            sscanf(p, "%d", &energy_tmp);
            break;
          }
          if (empty_line) {
            free(line);
            break;
          }
        }

        // write down the data to mapping
        if (tmp) {
          RNAstruc tmpstruc;
          tmpstruc.structure = tmp;
          tmpstruc.energy = energy_tmp;

          if (vertex_l.count(tmpstruc)==0) {
            fprintf(stderr, "WARNING: structure %5d ""%s"" %6.2f  not found in DSU\n", num, pt_to_str_pk(tmp).c_str(), tmpstruc.energy/100.0);
          } else {
            int num_real = vertex_l[tmpstruc];
            mapping.push_back(num_real);
          }
        }

        if (tmp) free(tmp); tmp = NULL;
        free(line);
        line = my_getline(filter);
      }
      fclose(filter);

      mapping_rev.clear();
      mapping_rev.resize(LM.size());
      for (int i=0; i<(int)mapping.size(); i++) {
        mapping_rev[mapping[i]] = i;
      }

    } else {
      fprintf(stderr, "WARNING: cannot open filter file ""%s""!", filter_file);
    }
  }


  if (!full && mapping.size() == 0) {
    mapping_rev.resize(LM.size(), 0);
    fprintf(stderr, "Mapping: \n");
    for (int i=0; i<(int)LM.size(); i++) {
      if (LM[i].type==NORMAL) {
        mapping.push_back(i);
        mapping_rev[i] = mapping.size()-1;
        fprintf(stderr, "%4d[%4d] ", i+1, (int)mapping.size());
      }
    }
    //fprintf(stderr, "(%d mapped to %d) \n", (int)LM.size(), (int)mapping.size());
  }

  fprintf(stderr, "Generating the matrix %c\n", type);
  int size = (full?LM.size():mapping.size());
  FILE *energies;
  energies = fopen(filename, "w");
  if (energies) {
    // create matrix
    if (type != 'D' && !generated) {
      matrix.clear();
      for (int i=0; i<size; i++) {
        int to_search = full?i:mapping[i];
        matrix.push_back(HeightSearch(to_search, edgesV_l));
        fprintf(stderr, "done %d/%d\r", i, size);
        // discard those we dont want
        if (!full) {
          for (int j=0; j<(int)mapping.size(); j++) {
            matrix[i][j] = matrix[i][mapping[j]];
          }
        }
        if ((int)matrix[i].size()!=size) {
          matrix[i].resize(size);
        }
      }
      // resolve i==i
      for (unsigned int i=0; i<matrix.size(); i++) {
        matrix[i][i] = make_pair(LM[i].energy, 0);
      }
      generated = true;
    }

    // print
    //fprintf(stderr, "Printing the matrix type %c\n", type);
    for (unsigned int i=0; i<matrix.size(); i++) {
      fprintf(energies, "%6d %s ", mapping[i]+1, LM[mapping[i]].str_ch);
      for (unsigned int j=0; j<matrix[i].size(); j++) {
        switch (type) {
        case 'E': fprintf(energies, "%6.2f ", matrix[i][j].first/100.0); break;
        case 'D': fprintf(energies, "%5d ", HammingDist(LM[i].structure, LM[j].structure)); break;
        case 'G': fprintf(energies, "%5d ", matrix[i][j].second); break;
        }
      }
      fprintf(energies, "\n");
    }
  }
  fclose(energies);
}

void DSU::PrintRates(char *filename, bool full, double temp, char mode)
{
  int size = (full?LM.size():number_lm);
  FILE *rates;
  rates = fopen(filename, "w");
  if (rates) {
    mode_rates mod;
    switch (mode) {
      case 'V': mod = VERTEX_CONTR; break;
      case 'E': mod = EDGE_CONTR_MIN; break;
      case 'M': mod = EDGE_CONTR_MAX; break;
      case 'F': mod = NO_CONTR; break;
      case 'S': mod = VERTEX_CONTR_SUM; break;
      default: mod = VERTEX_CONTR;
    }
    Graph graph(edges_l, LM, mod);
    fprintf(stderr, "graph created...%d LM\n", (int)LM.size());
    if (mod!=NO_CONTR) {
      for(int i=LM.size()-1; i>=size; i--) {
        /*char filename[20];
        char filename_eps[20];
        sprintf(filename, "smth%d.dot", i);
        sprintf(filename_eps, "smth%d.eps", i);
        graph.PrintDot(filename, true, true, filename_eps);*/

        // call appropriate contraction
        if (mod == VERTEX_CONTR || mod == VERTEX_CONTR_SUM) {
          graph.RemoveLastPoint();
        } else if (mod == EDGE_CONTR_MAX || mod == EDGE_CONTR_MIN) {
          while (!graph.RemoveLowestEdge()) ;
        }
        if (i%100 ==0) fprintf(stderr, "removed point %d\n", i);
      }
      char filename[20];
		  char filename_eps[20];
		  sprintf(filename, "reduced%c.dot", mode);
		  sprintf(filename_eps, "reduced%c.eps", mode);
		  graph.PrintDot(filename, true, true, filename_eps);
    }
    graph.PrintRates(rates, temp);
  }
  fclose(rates);
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
  for (set<edgeLL>::iterator it=edgesV_l[num].begin(); it!=edgesV_l[num].end(); it++) {
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
  fprintf(stderr, "Filling components.");
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
  fprintf(stderr, "..done, found %4d components\n", (int)comps.size());
}

void DSU::Color(int lm, int color, Component &cmp, vector<int> &LM_tmp, vector<int> &sadd_tmp) {
  // add to component and colour the node
  cmp.AddLM(lm, LM[lm].energy);
  LM_tmp[lm] = color;

  // proceed to non-coloured
  for (set<edgeLL>::iterator it=edgesV_l[lm].begin(); it!=edgesV_l[lm].end(); it++) {
    int goesTo = it->goesTo(lm);
    // info about saddle:
    if (it->saddle != -1 && sadd_tmp[it->saddle] == -1) {
      sadd_tmp[it->saddle] = color;
      cmp.AddSadd(it->saddle, it->en);
    }
    // info about LM
    if (LM_tmp[goesTo] == -1) {
      // colour it!
      Color(goesTo, color, cmp, LM_tmp, sadd_tmp);
    }
  }
}

int DSU::AddConnection(int num1, int num2, int energy, short *saddle, UF_set &connected, vector<vector<RNAsaddle*> > &connections) {

  int comp1 = LM_to_comp[num1];
  int comp2 = LM_to_comp[num2];

  // we can connect components:
  if (comp1 != comp2) {

    if (connections[comp1][comp2] == NULL || connections[comp1][comp2]->energy > energy) {
      RNAsaddle *tmpRNA = new RNAsaddle(num1, num2, COMP);
      tmpRNA->energy = energy;
      tmpRNA->structure = saddle;
      tmpRNA->str_ch = NULL;

      if (connections[comp1][comp2]) {
        free(connections[comp1][comp2]->structure);
        delete connections[comp1][comp2];
      }
      // assign new one
      connections[comp1][comp2] = tmpRNA;

      // add to UnionFindSet
      connected.union_set(comp1, comp2);
      return 1;
    }
  }
  // if it didn't work then release memory
  free(saddle);
  return 0;
}


// small struct for use in ConnectComps function
struct que_tmp {
  char *str_ch1;
  char *str_ch2;

  short *str1;
  short *str2;

  int hd;

  // debug info
  int i,j;
  int ci, cj;

  bool operator<(const que_tmp &second) const {
    return hd>second.hd;
  };

  que_tmp(short *st1, short *st2, char *ch1, char *ch2, int hd) {
    this->hd = hd;
    str_ch1 = ch1;
    str_ch2 = ch2;
    str1 = st1;
    str2 = st2;
  }

};

void DSU::ConnectComps(int maxkeep, bool debug)
{
  if (comps.size() == 0) FillComps();

  // just debug
  LM_to_comp[-1]=-1;


  if (comps.size()==1) return ;
  fprintf(stderr, "Connecting components (%d components).\n", (int)comps.size());
  int lastLM = LM.size();

  // create queue of saddle pairs
  priority_queue<que_tmp> queue;
  for (unsigned int i=0; i<comps.size(); i++) {
    for (unsigned int j=i+1; j<comps.size(); j++) {
      RNAstruc first = LM[comps[i].LMs[0]];
      RNAstruc second = LM[comps[j].LMs[0]];
      if (comps[i].max_saddle!=-1) first = saddles[comps[i].max_saddle];
      if (comps[j].max_saddle!=-1) second = saddles[comps[j].max_saddle];

      que_tmp tmp(first.structure, second.structure, first.str_ch, second.str_ch, HammingDist(first.structure, second.structure));
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
  vector<vector<RNAsaddle*> > connections(comps.size(), vector<RNAsaddle*>(comps.size(), NULL));

  // go throught the queue - take pair *
  while (!queue.empty() && !connected.connected_all()) {
    // get pair
    que_tmp top = queue.top();
    queue.pop();

    // apply threshold ??
    cnt++;

    // get path
    if (debug) fprintf(stderr, "path between (LM: %d, %d) (comp: %d, %d) (hd: %d):\n", top.i, top.j, top.ci, top.cj, top.hd);
    if (pknots) {
      path_pk *path = get_path_light_pk(seq, top.str1, top.str2, maxkeep);

      /*// variables for outer insertion
      double max_energy= -1e8;
      path_t *max_path = path;*/

      // variables for inner loops and insertions
      path_pk *tmp = path;
      path_pk *last = NULL;
      short *last_str = NULL;
      int last_en = tmp->en;
      int last_num = -1;

      // loop through whole path
      while (tmp && tmp->structure) {
        dbg_count++;
        // debug??
        if (debug) fprintf(stderr, "%s %6.2f ", pt_to_str_pk(tmp->structure).c_str(), tmp->en);

        /*// update max_energy
        if (max_energy < tmp->en) {
          max_energy = tmp->en;
          max_path = tmp;
        }*/
        // find adaptive walk
        short *tmp_str = allocopy(tmp->structure);
        //int tmp_en = move_rand(seq, tmp_str, s0, s1, 0);

        int tmp_en;
        Structure str(seq, tmp_str, s0, s1);
        tmp_en = move_gradient_pk(seq, &str, s0, s1, 0, 0);
        copy_arr(tmp_str, str.str);

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
              short *saddle = (last->en > tmp->en ? allocopy(last->structure) : allocopy(tmp->structure));
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
      free_path_pk(path);
    } else {
      path_t *path = get_path(seq, top.str_ch1, top.str_ch2, maxkeep);

      /*// variables for outer insertion
      double max_energy= -1e8;
      path_t *max_path = path;*/

      // variables for inner loops and insertions
      path_t *tmp = path;
      path_t *last = NULL;
      short *last_str = NULL;
      int last_en = tmp->en;
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
        int tmp_en = move_gradient(seq, tmp_str, s0, s1, 0, 0, 0);

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
  }
  // write connections to graph:
  edgesV_l.resize(LM.size());
  for (unsigned int i=0; i<connections.size(); i++) {
    for (unsigned int j=0; j<connections[i].size(); j++) {

      if (!connections[i][j]) continue;
      // saddle
      RNAsaddle* saddle = connections[i][j];
      saddle->type = COMP;
      //vertex_s.insert(make_pair(*saddle, saddles.size()));
      saddle->str_ch = pt_to_chars_pk(saddle->structure);
      saddles.push_back(*saddle);

      // edge l-l
      edgeLL e(saddle->lm1, saddle->lm2, saddle->energy, saddles.size()-1);
      edges_l.insert(e);
      edgesV_l[saddle->lm1].insert(e);
      edgesV_l[saddle->lm2].insert(e);
      // edge l-s
      edgeLS e2(saddle->lm1, saddles.size()-1);
      edges_ls.insert(e2);
      e2.i = saddle->lm2;
      edges_ls.insert(e2);

      // free the pointer
      free(saddle);
    }
  }

  fprintf(stderr, "LM resized from %d to %d\n", lastLM, (int)LM.size());

  //FillComps();
}

int DSU::AddLMtoComp(short *structure, int energy, bool debug, UF_set &connected, vector<vector<RNAsaddle*> > &connections)
{
  // resize connected and connections (we have new component...)
  int numComp = connected.size();
  connected.enlarge_parent();
  connections.resize(numComp+1);
  for (unsigned int i=0; i< connections.size(); i++) {
    connections[i].resize(numComp+1, NULL);
  }

  // add new LM to LMs
  int numLM = AddLMtoTBD(structure, energy, EE_COMP, debug);


  // add new component to map
  LM_to_comp[numLM] = numComp;

  return numLM;
}

void DSU::PrintComps(FILE *output, bool fill)
{
  // fix the sorting disorder
  if (comps.size()==0 || fill) FillComps();

  // print info about comps:
  for (unsigned int i=0; i<comps.size(); i++) {
    fprintf(output, "%4d %4d (%7.2f)", i, comps[i].min_lm+1, comps[i].min_energy/100.0);
    if (comps[i].max_saddle!=-1)  fprintf(output, " %4dS (%7.2f) :", comps[i].max_saddle+1, comps[i].max_energy/100.0);
    else                          fprintf(output, "                 :");
    for (unsigned int j=0; j<comps[i].LMs.size(); j++) {
      fprintf(output, " %4d", comps[i].LMs[j]+1);
    }
    fprintf(output, "      ");
    for (unsigned int j=0; j<comps[i].saddles.size(); j++) {
      fprintf(output, " %3dS", comps[i].saddles[j]+1);
    }
    fprintf(output, "\n");
  }
  fprintf(output, "\n");
}

void DSU::PrintLinkCP(FILE *output, bool fix)
{
  PrintLM(output, true);
  PrintSaddles(output, false);
}

void DSU::PrintLM(FILE *output, bool fix)
{
  if (fix) SortFix();
  fprintf(output, "      %s\n", seq);

  //printf("Local minima (%4d):\n", (int)LM.size());
  for (unsigned int i=0; i<LM.size(); i++) {
    fprintf(output, "%4d  %s %7.2f %8s\n", i+1, LM[i].str_ch, LM[i].energy/100.0, LMtype_string[LM[i].type]);
  }
}

void DSU::PrintBarr(FILE *output)
{
  fprintf(output, "     %s\n", seq);

  //get heights
  vector<vector<std::pair<int, int> > > res(number_lm);
  for (int i=0; i<number_lm; i++) {
    res[i] = HeightSearch(i, edgesV_l);
  }

  // get minimal height:
  vector<int> saddle_en(number_lm);
  vector<int> father(number_lm);
  for (int i=number_lm-1; i>=1; i--) {
    int minim = INT_MAX;
    int index;
    for (int j=i-1; j>=0; j--) {
      if (res[i][j].first <= minim) {
        minim = res[i][j].first;
        index = j;
      }
    }
    saddle_en[i] = minim;
    father[i] = index+1;
  }

  saddle_en[0] = LM[0].energy+0.1;
  father[0] = 0;


  for (int i=0; i<number_lm; i++) {
    fprintf(output, "%4d %s %6.2f %4d %6.2f\n", i+1, LM[i].str_ch, LM[i].energy/100.0, father[i], (saddle_en[i]-LM[i].energy)/100.0);
  }
}


void DSU::PrintSaddles(FILE *output, bool fix)
{
  if (fix) SortFix();
  //printf("Saddles (%4d):\n", (int)saddles.size());
  fprintf(output, "\n");
  // collect saddle info
  vector<set<int> > saddle_connLM (saddles.size());
  vector<set<int> > saddle_connSadd (saddles.size());
  for (set<edgeLS>::iterator it=edges_ls.begin(); it!=edges_ls.end(); it++) {
    saddle_connLM[it->j].insert(it->i);
  }
  for (set<edgeSS>::iterator it=edges_s.begin(); it!=edges_s.end(); it++) {
    saddle_connSadd[it->j].insert(it->i);
    saddle_connSadd[it->i].insert(it->j);
  }
  for (unsigned int i=0; i<saddles.size(); i++) {
    fprintf(output, "%4dS %s %7.2f %8s", i+1, saddles[i].str_ch, saddles[i].energy/100.0, SDtype_string[saddles[i].type]);
    for (set<int>::iterator it=saddle_connLM[i].begin(); it!=saddle_connLM[i].end(); it++) fprintf(output, " %4d", (*it)+1);
    for (set<int>::iterator it=saddle_connSadd[i].begin(); it!=saddle_connSadd[i].end(); it++) fprintf(output, " %3dS", (*it)+1);
    fprintf(output, "\n");
  }
}

bool DSU::SortFix()
{
  //return ;
  vector<int> mappingLM(LM.size());
  vector<int> mappingSD(saddles.size());
  {
    // create mapping (and sort LMs)
    vector<std::pair<RNAlocmin, int> > temp(LM.size());
    for (unsigned int i=0; i<LM.size(); i++) {
      temp[i] = make_pair(LM[i], i);
    }
    sort(temp.begin(), temp.end());
    for (unsigned int i=0; i<temp.size(); i++) {
      LM[i]=temp[i].first;
      mappingLM[temp[i].second]=i;
      vertex_l[temp[i].first]=i;
    }
  }
  {
    // create mapping (and sort saddles)
    vector<std::pair<RNAsaddle, int> > temp;
    for (unsigned int i=0; i<saddles.size(); i++) {
      temp.push_back(make_pair(saddles[i], i));
    }
    sort(temp.begin(), temp.end());
    for (unsigned int i=0; i<temp.size(); i++) {
      saddles[i]=temp[i].first;
      saddles[i].lm1 = mappingLM[saddles[i].lm1];
      saddles[i].lm2 = mappingLM[saddles[i].lm2];
      if (saddles[i].lm1 > saddles[i].lm2) swap(saddles[i].lm1, saddles[i].lm2);
      mappingSD[temp[i].second]=i;
    }
  }

  // now process all edges:
  set<edgeLL> edges_l_t;
  for (set<edgeLL>::iterator it=edges_l.begin(); it!=edges_l.end(); it++) {
    edges_l_t.insert(edgeLL(mappingLM[it->i], mappingLM[it->j], saddles[mappingSD[it->saddle]].energy, mappingSD[it->saddle]));
  }
  edges_l = edges_l_t;

  set<edgeLS> edges_ls_t;
  for (set<edgeLS>::iterator it=edges_ls.begin(); it!=edges_ls.end(); it++) {
    edges_ls_t.insert(edgeLS(mappingLM[it->i], mappingSD[it->j]));
  }
  edges_ls = edges_ls_t;

  set<edgeSS> edges_s_t;
  for (set<edgeSS>::iterator it=edges_s.begin(); it!=edges_s.end(); it++) {
    edges_s_t.insert(edgeSS(mappingSD[it->i], mappingSD[it->j]));
  }

  // process edgesV_l:
  vector<set<edgeLL> > edges_temp(edgesV_l.size());
  for (unsigned int i=0; i<edgesV_l.size(); i++) {
    for (set<edgeLL>::iterator it=edgesV_l[i].begin(); it!=edgesV_l[i].end(); it++) {
      edges_temp[i].insert(edgeLL(mappingLM[it->i], mappingLM[it->j], saddles[mappingSD[it->saddle]].energy, mappingSD[it->saddle]));
    }
  }
  for (unsigned int i=0; i<edgesV_l.size(); i++) {
    edgesV_l[mappingLM[i]]=edges_temp[i];
  }

  // recompute components:
  if (comps.size()>0) {
    FillComps();
    return true;
  }
  return false;
}
