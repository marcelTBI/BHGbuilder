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

TBD::TBD()
{
  for (int i=0; i<type1_len; i++) {
    sizes[i] = 0;
  }
}

bool TBD::insert(int i, int j, type1 type, bool fiber)
{
  sizes[type]++;
  ResizeDone(max(i,j)+1);

  if (done[i][j]) return false;
  else {
    //fprintf(stderr, "inserting %d %d %s\n", i, j, type1_str[type]);
    tbd.push(TBDentry(i, j, type, fiber));
    return true;
  }
}

int TBD::size() {
  return tbd.size();
}

void TBD::ResizeDone(int new_size) {
  if (new_size>(int)done.size()) {
    // resize existing:
    for (unsigned int i=0; i<done.size(); i++) {
      done[i].resize(new_size, false);
    }

    done.resize(new_size, vector<bool> (new_size, false));
  }
}

TBDentry TBD::get_first()
{
  if (size()==0) return TBDentry(-1,-1,NEW_FOUND,-1);
  TBDentry tbde = tbd.top();
  tbd.pop();
  ResizeDone(max(tbde.i, tbde.j)+1);
  while ((size() > 0) && (done[tbde.i][tbde.j])) {  // maybe we dont need this
    tbde = tbd.top();
    ResizeDone(max(tbde.i, tbde.j)+1);
    tbd.pop();
  }
  if (size() == 0) return TBDentry(-1,-1,NEW_FOUND,-1);
  ResizeDone(max(tbde.i, tbde.j)+1);
  done[tbde.i][tbde.j] = true;
  return tbde;
}

void TBD::join(TBD &second) {   // can be efficient
  //assert(second.size()==0);
  ResizeDone(second.done.size()+1);
  for (unsigned int i=0; i<second.done.size(); i++) {
    for (unsigned int j=0; j<second.done[i].size(); j++) {
      done[i][j] = max(second.done[i][j], done[i][j]);
    }
  }
}

int DSU::Cluster(Opt &opt, int kmax, TBD &output)
{
  // create data structures
  vector<lm_pair> to_cluster;
  to_cluster.reserve(LM.size());
  UF_set_child ufset;
  ufset.enlarge_parent(LM.size());

  // representative nodes
  set<int> represents;

  // fill it
  for (unsigned int i=0; i<LM.size(); i++) {
    for (unsigned int j=i+1; j<LM.size(); j++) {
      to_cluster.push_back(lm_pair(i,j,HammingDist(LM[i].structure, LM[j].structure)));
    }
  }
  sort(to_cluster.begin(), to_cluster.end());

  // process:
  int last_hd = to_cluster[0].hd;
  for (unsigned int i=0; i<to_cluster.size(); i++) {
    lm_pair &cp = to_cluster[i];

    if (cp.hd!=last_hd) {
      // do something, cause we are on higher level...
    }

    // see if we are not joint yet:
    if (!ufset.joint(cp.i, cp.j)) {

      //fprintf(stderr, "clustering %d %d (%d)\n", cp.i, cp.j, cp.d);
      // try to connect
      int father1 = ufset.find(cp.i);
      int father2 = ufset.find(cp.j);
      if (ufset.count(father1) + ufset.count(father2) > kmax) { // cannot connect them, need to insert all edges into the TBD

        // join clusters
        JoinClusters(opt, ufset, represents, output, cp.i, cp.j);

      } else {
        // connect them
        ufset.union_set(cp.i, cp.j);

      }
    }
    last_hd = cp.hd;
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
  // and its represent node
  represents.insert(father);

  // and finally add represen edges:
  for (set<int>::iterator it=represents.begin(); it!=represents.end(); it++) {
    set<int>::iterator it2 = it; it2++;
    for (;it2!=represents.end(); it2++) {
      output.insert(*it, *it2, REPRESENT, false);
    }
  }


  fprintf(stderr, "output size = %d (%d, %d, %d)\n", output.size(), output.sizes[0], output.sizes[1], output.sizes[2]);

  // now finish:
  vector<RNAsaddle> UBout;
  UBout = ComputeTBD2(output, opt.maxkeep, opt.num_threshold, opt.outer, opt.noLP, opt.shifts, opt.debug);
  UBoutput.insert(UBoutput.end(), UBout.begin(), UBout.end());
  sort(UBoutput.begin(), UBoutput.end());

  return 0;
}

void DSU::GetRepre(vector<RNAsaddle> &UBoutput, TBD &output, set<int> &represents, set<int> &children, Opt &opt) {

  // get intercluster connections
  TBD cluster;
  for (set<int>::iterator it=children.begin(); it!=children.end(); it++) {
    set<int>::iterator it2 = it; it2++;
    for (;it2!=children.end(); it2++) {
      cluster.insert(*it, *it2, INTER_CLUSTER, false);
    }
  }

  // insert representative minima:
  represents.insert(*children.begin());

  // do we want exactly number of represents?
  int more = true; // now hardcoded, and I dont think it would be different ;-)
  int rsize = represents.size();

  // collect saddles for intercluster connections
  vector<RNAsaddle> input;
  input = ComputeTBD2(cluster, opt.maxkeep, opt.num_threshold, opt.outer, opt.noLP, opt.shifts, opt.debug);

  // insert new representatives
  if (opt.fbarrier) {
    // insert according to highest barrier
    vector<std::pair<int, int> > sort_by_barr;
    for (unsigned int i=0; i<input.size(); i++) {
      int barrier = input[i].energy - max(LM[input[i].lm1].energy, LM[input[i].lm2].energy);
      sort_by_barr.push_back(make_pair(barrier, i));
    }
    sort(sort_by_barr.begin(), sort_by_barr.end());


    int many = max(1, (int)(opt.repre_portion*children.size()));
    if (more) {
      // insert as loong as we need them
      int pos = sort_by_barr.size()-1;
      while (pos>=0 && (int)represents.size() - rsize < many*2) {
        represents.insert(input[sort_by_barr[pos].second].lm1);
        represents.insert(input[sort_by_barr[pos].second].lm2);
        pos--;
      }
    } else {
      // get repre
      for (int i=0; i<many; i++) {
        int pos = sort_by_barr.size()-1-i;
        if (pos<0) break;
        represents.insert(input[sort_by_barr[pos].second].lm1);
        represents.insert(input[sort_by_barr[pos].second].lm2);
      }
    }

  } else {
    // insert according to highest saddl
    int many = max(1, (int)(opt.repre_portion/2.0*children.size()));

    if (more) {
      // insert as long as we need them
      int pos = input.size()-1;
      while (pos>=0 && (int)represents.size() - rsize < many*2) {
        represents.insert(input[pos].lm1);
        represents.insert(input[pos].lm2);
        pos--;
      }
    } else {
      // insert just approx.
      for (int i=0; i<many; i++) {
        int pos = input.size()-1-i;
        if (pos<0) break;
        represents.insert(input[pos].lm1);
        represents.insert(input[pos].lm2);
      }
    }
    if (opt.debug) fprintf(stderr, "cluster size: %5d, acquiring %3d represents, repre size: %4d\n", (int)children.size(), many*2, (int)represents.size());
  }

  // add them to global output.
  UBoutput.insert(UBoutput.end(), input.begin(), input.end());
  // join queues
  output.join(cluster);
}

int DSU::JoinClusters(Opt &opt, UF_set_child &ufset, set<int> &represents, TBD &output, int i, int j) {

  // insert crit edge:
  output.insert(i, j, CRIT_EDGE, false);

  set<int> childreni = ufset.get_children(ufset.find(i));
  set<int> childrenj = ufset.get_children(ufset.find(j));

  // do computation + represent LM generation for each of 2 clusters:
  GetRepre(UBoutput, output, represents, childreni, opt);
  GetRepre(UBoutput, output, represents, childrenj, opt);

  // now make from this group only one vertex (maybe wrong)
  ufset.union_set(i, j);
  ufset.make_single(i);

  if (opt.debug) fprintf(stderr, "repre size = %d\n", represents.size());

  return 0;
}

vector<RNAsaddle> DSU::ComputeTBD2(TBD &pqueue, int maxkeep, int num_threshold, bool outer, bool noLP, bool shifts, bool debug)
{
  int dbg_count = 0;
  int cnt = 0;
  int norm_cf = 0;

  // partial results
  map<lm_pair, RNAsaddle, pq_setcomp> UBlist;

  // go through all pairs in queue
  while (pqueue.size()>0) {
    // check time:
    double time_secs = ((clock()  - time)/(double)CLOCKS_PER_SEC);
    if (stop_after && (time_secs > stop_after)) {
      fprintf(stderr, "Time threshold reached (%d secs.), processed %d/%d\n", stop_after, cnt, pqueue.size()+cnt);
      break;
    }

    // apply threshold
    if (cnt>num_threshold) {
      fprintf(stderr, "Number threshold reached, processed %d/%d\n", cnt, pqueue.size()+cnt);
      break;
    } else {
      cnt++;
    }

    // get next
    TBDentry tbd = pqueue.get_first();
    if (tbd.i==-1) break;


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
          InsertUB(UBlist, num1, num2, en_tmp, saddle, false, debug);


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
    if (outer) InsertUB(UBlist, tbd.i, tbd.j, en_fltoi(max_energy), make_pair_table(max_path->s), true, debug);

    // free stuff
    if (last_str) free(last_str);
    free_path(path);
  } // all doing while

  // now just resort UBlist to something sorted according energy
  vector<RNAsaddle> UBoutput;
  UBoutput.reserve(UBlist.size());
  for (map<lm_pair, RNAsaddle, pq_setcomp>::iterator it=UBlist.begin(); it!=UBlist.end(); it++) {
    if (it->second.str_ch) free(it->second.str_ch);
    it->second.str_ch = pt_to_char(it->second.structure);
    UBoutput.push_back(it->second);
  }
  sort(UBoutput.begin(), UBoutput.end());
  UBlist.clear();

  return UBoutput;
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

