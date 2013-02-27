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

struct cluster_pair {
  int i, j;
  int d;

  cluster_pair(int i, int j, int d) {
    this->i = i;
    this->j = j;
    this->d = d;
  }

  bool operator <(const cluster_pair &sec) const {
    if (d != sec.d) {
      return d<sec.d;
    } else {
      if (i != sec.i) {
        return i<sec.i;
      } else {
        return j<sec.j;
      }
    }
  }
};

TBD::TBD()
{
  for (int i=0; i<type1_len; i++) {
    sizes[i] = 0;
  }
}

bool TBD::insert(int i, int j, type1 type, bool fiber)
{
  sizes[type]++;

  if (done.count(make_pair(i,j))>0) return false;
  else {
    //fprintf(stderr, "inserting %d %d %s\n", i, j, type1_str[type]);
    tbd.push(TBDentry(i, j, type, fiber));
    return true;
  }
}

int TBD::size() {
  return tbd.size();
}

TBDentry TBD::get_first()
{
  if (size()==0) return TBDentry(-1,-1,NEW_FOUND,-1);
  TBDentry tbde = tbd.top();
  tbd.pop();
  done.insert(make_pair(tbde.i, tbde.j));
  return tbde;
}


int DSU::Cluster(int kmax, TBD &output)
{
  // create data structures
  vector<cluster_pair> to_cluster;
  to_cluster.reserve(LM.size());
  UF_set_child ufset;
  ufset.enlarge_parent(LM.size());

  // fill it
  for (unsigned int i=0; i<LM.size(); i++) {
    for (unsigned int j=i+1; j<LM.size(); j++) {
      to_cluster.push_back(cluster_pair(i,j,HammingDist(LM[i].structure, LM[j].structure)));
    }
  }
  sort(to_cluster.begin(), to_cluster.end());

  // process:
  int last_hd = to_cluster[0].d;
  for (unsigned int i=0; i<to_cluster.size(); i++) {
    cluster_pair &cp = to_cluster[i];

    if (cp.d!=last_hd) {
      // do something, cause we are on higher level...
    }

    // see if we are not joint yet:
    if (!ufset.joint(cp.i, cp.j)) {

      //fprintf(stderr, "clustering %d %d (%d)\n", cp.i, cp.j, cp.d);
      // try to connect

      if (ufset.count(ufset.find(cp.i)) + ufset.count(ufset.find(cp.j)) > kmax) { // cannot connect them, need to insert all edges into the TBD
        // insert crit edge:
        output.insert(cp.i, cp.j, CRIT_EDGE, false);

        set<int> first = ufset.get_children(ufset.find(cp.i));
        set<int> second = ufset.get_children(ufset.find(cp.j));

        // insert rpresentative edge
        output.insert(*first.begin(), *second.begin(), REPRESENT, false);

        // insert all inter edges:
        for (set<int>::iterator it=first.begin(); it!=first.end(); it++) {
          set<int>::iterator it2 = it; it2++;
          for (;it2!=first.end(); it2++) {
            output.insert(*it, *it2, INTER_CLUSTER, false);
          }
        }

        for (set<int>::iterator it=second.begin(); it!=second.end(); it++) {
          set<int>::iterator it2 = it; it2++;
          for (;it2!=second.end(); it2++) {
            output.insert(*it, *it2, INTER_CLUSTER, false);
          }
        }

        // now make from this group only one vertex (maybe wrong)
        ufset.union_set(cp.i, cp.j);
        ufset.make_single(cp.i);


      } else {
        // connect them
        ufset.union_set(cp.i, cp.j);

      }
    }
    last_hd = cp.d;
  }

  // now we have just one cluster, we have to add all intercluster connections that are left:
  int father = ufset.find(0);
  set<int> first = ufset.get_children(father);
  // insert all inter edges:
  for (set<int>::iterator it=first.begin(); it!=first.end(); it++) {
    set<int>::iterator it2 = it; it2++;
    for (;it2!=first.end(); it2++) {
      output.insert(*it, *it2, INTER_CLUSTER, false);
    }
  }


  fprintf(stderr, "output size = %d = (%d, %d, %d)\n", output.size(), output.sizes[0], output.sizes[1], output.sizes[2]);

  return 0;
}

int DSU::ComputeTBD(TBD &pqueue, int maxkeep, int num_threshold, bool outer, bool noLP, bool shifts, bool debug)
{
  int dbg_count = 0;
  int cnt = 0;
  int hd_threshold = INT_MAX;
  int norm_cf = 0;

  // go through all pairs in queue
  while (pqueue.size()>0) {
    // check time:
    if (stop_after && (time  - clock())/(double)CLOCKS_PER_SEC > stop_after) {
      break;
    }

    // get next
    TBDentry tbd = pqueue.get_first();

    // apply threshold
    if (cnt>num_threshold) {
      break;
    } else {
      cnt++;
    }

    // get path
    if (debug) fprintf(stderr, "path between (%3d, %3d) type=%s fiber=%d:\n", tbd.i, tbd.j, type1_str[tbd.type_clust], tbd.fiber);
    path_t *path = get_path(seq, LM[tbd.i].str_ch, LM[tbd.j].str_ch, maxkeep);

    // variables for outer insertion
    double max_energy= -1e8;
    path_t *max_path = path;

    // variables for inner loops and insertions
    path_t *tmp = path;
    path_t *last = NULL;
    short *last_str = NULL;
    int last_en;
    int last_num = -1;

    // how long this path is:
    int path_length = 1;

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
        path_length++;
        // not equal LM - we can update something in UBlist
          // find LM num:
        int num1 = (last_num!=-1?last_num:FindNum(last_en, last_str));
        int num2 = FindNum(tmp_en, tmp_str);

        // update UBlist
        if (num1==-1 || num2==-1) {
          if (num2==-1) {
            if (gl_maxen <= tmp_en) {
              //fprintf(stderr, "exceeds en.: %s %6.2f\n", pt_to_str(tmp_str).c_str(), tmp_en/100.0);
              num2 = AddLMtoTBD(tmp_str, tmp_en, EE_DSU, debug);

            } else {
              fprintf(stderr, "cannot find: %s %6.2f\n", pt_to_str(tmp_str).c_str(), tmp_en/100.0);
              // add to list of minima and count with them later...
              num2 = AddLMtoTBD(tmp_str, tmp_en, NORM_CF, debug);
              norm_cf++;
            }
          }
        }
        // again check if we can add better saddle
        if (num1!=-1 && num2!=-1) {
          // store (maybe) better saddle to UB
          int en_tmp = en_fltoi(max(last->en, tmp->en));
          short *saddle = (last->en > tmp->en ? make_pair_table(last->s) : make_pair_table(tmp->s));
          InsertUB(num1, num2, en_tmp, saddle, false, debug);


          // try to insert new things into TBD:
          bool do_insert = path_length>2;
          if (path_length==2) {
            path_t *check = tmp;
            check++;
            if (check != NULL) {
              do_insert = true;
            }
          }
          if (do_insert) {
            pqueue.insert(num1, num2, NEW_FOUND, true);
          }
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
    if (outer) InsertUB(tbd.i, tbd.j, en_fltoi(max_energy), make_pair_table(max_path->s), true, debug);

    // free stuff
    if (last_str) free(last_str);
    free_path(path);
  } // all doing while

  // now just resort UBlist to something sorted according energy
  UBoutput.reserve(UBlist.size());
  for (map<pq_entry, RNAsaddle, pq_setcomp>::iterator it=UBlist.begin(); it!=UBlist.end(); it++) {
    if (it->second.str_ch) free(it->second.str_ch);
    it->second.str_ch = pt_to_char(it->second.structure);
    UBoutput.push_back(make_pair(it->second, it->first));
  }
  sort(UBoutput.begin(), UBoutput.end());
  UBlist.clear();

  // check if everything has been found:
  fprintf(stderr, "Found: %d connections\nLM resized from: %d to %d (%d missing, %d above the energy threshold)\n", (int)UBoutput.size(), number_lm, (int)LM.size(), norm_cf, (int)LM.size()-number_lm-norm_cf);

  return 0;
}

int DSU::AddLMtoTBD(short *tmp_str, int tmp_en, LMtype type, bool debug)
{
  RNAlocmin rna;
  rna.energy = tmp_en;
  rna.structure = allocopy(tmp_str);
  rna.str_ch = pt_to_char(tmp_str);
  rna.type = type;

  // insert LM and return its number
  LM.push_back(rna);
  vertex_l[rna]=LM.size()-1;
  return LM.size()-1;
}

