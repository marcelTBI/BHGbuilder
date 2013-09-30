#include <stdio.h>
#include <string.h>
#include <stdlib.h>

extern "C" {
  #include "utils.h"
  #include "fold_vars.h"
  #include "pair_mat.h"

  #include "fold.h"
  #include "findpath.h"
  #include "move_set.h"
}

#include "BHGbuilder.h"

TBD::TBD()
{
  for (int i=0; i<type1_len; i++) {
    sizes[i] = 0;
  }
}

bool TBD::insert(int i, int j, type1 type, bool fiber)
{
  ResizeDone(max(i,j)+1);

  if (i>j) swap(i, j);

  if (done[i][j]) return false;
  else {
    sizes[type]++;
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
  TBDentry tbde = tbd.top(); tbd.pop();
  ResizeDone(max(tbde.i, tbde.j)+1);
  while (done[tbde.i][tbde.j]) {  // maybe we dont need this
    if (size()==0) return TBDentry(-1,-1,NEW_FOUND,-1);
    tbde = tbd.top(); tbd.pop();
    ResizeDone(max(tbde.i, tbde.j)+1);
  }
  done[tbde.i][tbde.j] = true;
  return tbde;
}

void TBD::join(TBD &second) {   // can be more efficient
  //assert(second.size()==0);
  ResizeDone(second.done.size()+1);
  for (unsigned int i=0; i<second.done.size(); i++) {
    for (unsigned int j=0; j<second.done[i].size(); j++) {
      done[i][j] = max(second.done[i][j], done[i][j]);
    }
  }
}

//debug:
int debug_c = 0;
int debug_c2 = 0;


int DSU::FindNum(int en_par, short *str_par)
{
  debug_c2++;
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

vector<vector<int> > histo(100, vector<int>(100, 0));

int DSU::Cluster(Opt &opt, int kmax)
{
  // pqueue for pairs of LM
  TBD output;

  // if no-conn flag:
  if (opt.no_conn) conectivity.enlarge_parent(LM.size());

  if (kmax>0) {
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
        // do something?, cause we are on higher level...
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

    // and finally add represent edges:
    for (set<int>::iterator it=represents.begin(); it!=represents.end(); it++) {
      set<int>::iterator it2 = it; it2++;
      for (;it2!=represents.end(); it2++) {
        output.insert(*it, *it2, REPRESENT, false);
      }
    }
  } else {
    // now we don't do clustering, we have to add all intercluster connections
    for (unsigned int i=0; i<LM.size(); i++) {
      for (unsigned int j=i+1; j<LM.size(); j++) {
        output.insert(i, j, INTER_CLUSTER, false);
      }
    }
  }

  fprintf(stderr, "output size = %d (%d, %d, %d)\n", output.size(), output.sizes[0], output.sizes[1], output.sizes[2]);

  // now finish:
  ComputeTBD(output, opt.maxkeep, opt.num_threshold, opt.outer, opt.noLP, opt.shifts, opt.debug);

  // now just resort UBlist to something sorted according energy
  saddles.reserve(UBlist.size());
  for (set<RNAsaddle, RNAsaddle_comp>::iterator it=UBlist.begin(); it!=UBlist.end(); it++) {
    RNAsaddle saddle = *it;
    if (it->str_ch) free(it->str_ch);
    saddle.str_ch = pt_to_char(it->structure);
    saddles.push_back(saddle);
  }
  sort(saddles.begin(), saddles.end());
  UBlist.clear();
/*
  for (int i=0; i<saddles.size(); i++) {
    fprintf(stderr, "%d %d %.2f\n", saddles[i].lm1, saddles[i].lm2, saddles[i].energy/100.0);
  }*/

  // debug
  if (opt.debug) {
    fprintf(stderr, "found %d, not found %d\n", debug_c, debug_c2);
    for (int i=0; i<(int)histo.size(); i++) {
      if (histo[i][0]) {
        fprintf(stderr, "%5d(%5d) |", i, histo[i][0]);
        for (int j=1; j<min(50, (int)histo[i].size()); j++) {
          fprintf(stderr, "%5d", histo[i][j]);
        }
        fprintf(stderr, "\n");
      }
    }
  }


  return 0;
}

int DSU::JoinClusters(Opt &opt, UF_set_child &ufset, set<int> &represents, TBD &output, int i, int j) {

  // insert crit edge:
  output.insert(i, j, CRIT_EDGE, false);

  set<int> childreni = ufset.get_children(ufset.find(i));
  set<int> childrenj = ufset.get_children(ufset.find(j));

  // do computation + represent LM generation for each of 2 clusters:
  GetRepre(output, represents, childreni, opt);
  GetRepre(output, represents, childrenj, opt);

  // now make from this group only one vertex (maybe wrong)
  ufset.union_set(i, j);
  ufset.make_single(i);

  // message:
  fprintf(stderr, "Joining clusters, reduced to dimension %6d/%6d \n", ufset.dimension(), ufset.size());

  //fprintf(stderr, "repre size = %d\n", (int)represents.size());

  return 0;
}

void DSU::GetRepre(TBD &output, set<int> &represents, set<int> &children, Opt &opt) {

  // get intercluster connections
  if (opt.debug) fprintf(stderr, "Children: ");
  TBD cluster;
  for (set<int>::iterator it=children.begin(); it!=children.end(); it++) {
    if (opt.debug) fprintf(stderr, "%d ", *it);
    set<int>::iterator it2 = it; it2++;
    for (;it2!=children.end(); it2++) {
      cluster.insert(*it, *it2, INTER_CLUSTER, false);
    }
  }
  if (opt.debug) fprintf(stderr, "\n");

  // insert representative minima:
  represents.insert(*children.begin());

  // do we want exactly number of represents?
  int more = true; // now hardcoded, and I dont think it would be different ;-)
  int rsize = represents.size();

  // collect saddles for intercluster connections
  vector<RNAsaddle> input;
  ComputeTBD(cluster, opt.maxkeep, opt.num_threshold, opt.outer, opt.noLP, opt.shifts, opt.debug, &input);

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
    // insert according to highest saddle
    int many = max(1, (int)(opt.repre_portion/2.0*children.size()));
    sort(input.begin(), input.end());

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

  // join queues
  output.join(cluster);
}

void DSU::FindNumbers(int begin, int end, path_t *path, vector<int> &lm_numbers, bool shifts, bool noLP, bool debug)
{
  // first resolve small case:
  if (end-begin<4) {
    bool begins = true;
    for (int i=begin+1; i<end; i++) {

      // get the minimum
      short *tmp_str = make_pair_table(path[i].s);
      int tmp_en = move_deepest(seq, tmp_str, s0, s1, 0, shifts, noLP);
      // speedup
      if (begins && tmp_en == LM[lm_numbers[begin]].energy && str_eq(LM[lm_numbers[begin]].structure, tmp_str)) {
        lm_numbers[i] = lm_numbers[begin];
      } else {
        begins = false;
        if (tmp_en == LM[lm_numbers[end]].energy && str_eq(LM[lm_numbers[end]].structure, tmp_str)) {
          lm_numbers[i] = lm_numbers[end];
        }
      }

      if (lm_numbers[i]==-1) {
        lm_numbers[i] = FindNum(tmp_en, tmp_str);

        // update UBlist
        if (lm_numbers[i]==-1) {
          if (gl_maxen < tmp_en) {
            //fprintf(stderr, "exceeds en.: %s %6.2f\n", pt_to_str(tmp_str).c_str(), tmp_en/100.0);
            lm_numbers[i] = AddLMtoTBD(tmp_str, tmp_en, EE_DSU, debug);

          } else {
            if (debug) fprintf(stderr, "cannot find: %s %6.2f\n", pt_to_str(tmp_str).c_str(), tmp_en/100.0);
            // add to list of minima and count with them later...
            lm_numbers[i] = AddLMtoTBD(tmp_str, tmp_en, NORM_CF, debug);
          }
        }
      } else debug_c++;

      free(tmp_str);
    }
    return ;
  }

  // da middle one
  int pivot = (end+begin)/2;

  short *tmp_str = make_pair_table(path[pivot].s);
  //fprintf(stderr, "%s\n", pt_to_str(tmp_str).c_str());
  int tmp_en = move_deepest(seq, tmp_str, s0, s1, 0, shifts, noLP);

  //fprintf(stderr, "%s\n", pt_to_str(tmp_str).c_str());

  // speed up:
  if (tmp_en == LM[lm_numbers[begin]].energy && str_eq(LM[lm_numbers[begin]].structure, tmp_str)) {
    lm_numbers[pivot] = lm_numbers[begin];
  } else {
    if (tmp_en == LM[lm_numbers[end]].energy && str_eq(LM[lm_numbers[end]].structure, tmp_str)) {
      lm_numbers[pivot] = lm_numbers[end];
    }
  }

  // normal behaviour
  if (lm_numbers[pivot]==-1) {
    lm_numbers[pivot] = FindNum(tmp_en, tmp_str);

    // update UBlist
    if (lm_numbers[pivot]==-1) {
      if (gl_maxen < tmp_en) {
        //fprintf(stderr, "exceeds en.: %s %6.2f\n", pt_to_str(tmp_str).c_str(), tmp_en/100.0);
        lm_numbers[pivot] = AddLMtoTBD(tmp_str, tmp_en, EE_DSU, debug);

      } else {
        if (debug) fprintf(stderr, "cannot find: %s %6.2f\n", pt_to_str(tmp_str).c_str(), tmp_en/100.0);
        // add to list of minima and count with them later...
        lm_numbers[pivot] = AddLMtoTBD(tmp_str, tmp_en, NORM_CF, debug);
      }
    }
  } else debug_c++;

  free(tmp_str);

  // continue recursion:
  if (lm_numbers[pivot]!=lm_numbers[begin] && pivot-begin>1) FindNumbers(begin, pivot, path, lm_numbers, shifts, noLP, debug);
  if (lm_numbers[pivot]!=lm_numbers[end] && end-pivot>1) FindNumbers(pivot, end, path, lm_numbers, shifts, noLP, debug);

  // return maximal energy
  return ;
}

void DSU::ComputeTBD(TBD &pqueue, int maxkeep, int num_threshold, bool outer, bool noLP, bool shifts, bool debug, vector<RNAsaddle> *output_saddles)
{
  int cnt = 0;

  // go through all pairs in queue
  while (pqueue.size()>0) {
    // check time:
    double time_secs = ((clock()  - time)/(double)CLOCKS_PER_SEC);
    if (stop_after && (time_secs > stop_after)) {
      fprintf(stderr, "Time threshold reached (%d secs.), processed %d/%d\n", stop_after, cnt, pqueue.size()+cnt);
      break;
    }

    // just visualisation
    if (!output_saddles && cnt%1000==0) {
      fprintf(stderr, "Finding path: %7d/%7d\n", cnt, pqueue.size()+cnt);
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

    // check no-conn
    if (conectivity.size() > 0 && !tbd.fiber && conectivity.joint(tbd.i, tbd.j)) continue;


    // get path
    if (debug) fprintf(stderr, "path between (%3d, %3d) type=%s fiber=%d:\n", tbd.i, tbd.j, type1_str[tbd.type_clust], tbd.fiber);
    //2fprintf(stderr, "depth: %d\n%s\n%s\n%s\n", maxkeep, seq, LM[tbd.i].str_ch, LM[tbd.j].str_ch);
    path_t *path = get_path(seq, LM[tbd.i].str_ch, LM[tbd.j].str_ch, maxkeep);

    // variables for outer insertion
    double max_energy= -1e8;
    path_t *max_path = path;

    // variables for inner loops and insertions

    // get the length of path for speed up
    int length = 0;
    for (path_t *tmp = path; tmp && tmp->s; tmp++) {
      length ++;
    }

    // create vector of known LM numbers on path (where 0 and length-1 are known)
    vector<int> lm_numbers(length, -1);
    lm_numbers[0] = tbd.i;
    lm_numbers[length-1] = tbd.j;

    // bisect the path and find new LMs:
    FindNumbers(0, length-1, path, lm_numbers, shifts, noLP, debug);


    // debug
    if (debug) {
      int diff = 1;
      int last_num = lm_numbers[0];
      for (int i=0; i<length; i++) {
        fprintf(stderr, "path[%3d]= %4d (%s %6.2f)\n", i, lm_numbers[i], path[i].s, path[i].en);
        if (lm_numbers[i]!=last_num && lm_numbers[i]!=-1) {
          diff++;
          last_num=lm_numbers[i];
        }
      }
      histo[length][diff]++;
      histo[length][0]++;
    }

    // now process the array of found numbers:
    int last_num = lm_numbers[0];
    for (int i=1; i<length; i++) {
      if (lm_numbers[i]!=-1 && lm_numbers[i]!=last_num) {

        // save saddle
        RNAsaddle saddle(last_num, lm_numbers[i], DIRECT);
        saddle.energy = en_fltoi(max(path[i-1].en, path[i].en));
        saddle.str_ch = NULL;
        saddle.structure = (path[i-1].en > path[i].en ? make_pair_table(path[i-1].s) : make_pair_table(path[i].s));
        bool inserted = InsertUB(saddle, debug);

        // ???
        if (output_saddles && inserted) {
          output_saddles->push_back(saddle);
        }

        // try to insert new things into TBD:
        if (lm_numbers[i]!=lm_numbers[length-1] || lm_numbers[i-1]!=lm_numbers[0]) {
          // check no-conn
          if (conectivity.size() > 0) conectivity.union_set(tbd.i, tbd.j);
          pqueue.insert(lm_numbers[i-1], lm_numbers[i], NEW_FOUND, true);
        }
        last_num = lm_numbers[i];
      }
    }

   /* // loop through whole path
    while (tmp && tmp->s) {
      dbg_count++;
      // debug??
      if (debug) fprintf(stderr, "%s %6.2f", tmp->s, tmp->en);

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

        if (debug) fprintf(stderr, " %d\n", num2);



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
            }
          }
        }
        // again check if we can add better saddle
        if (num1!=-1 && num2!=-1) {
          // store (maybe) better saddle to UB
          RNAsaddle saddle(num1, num2, DIRECT);
          saddle.energy = en_fltoi(max(last->en, tmp->en));
          saddle.str_ch = NULL;
          saddle.structure = (last->en > tmp->en ? make_pair_table(last->s) : make_pair_table(tmp->s));

          bool inserted = InsertUB(saddle, debug);

          if (output_saddles && inserted) {
            output_saddles->push_back(saddle);
          }

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
      } else if (debug) fprintf(stderr, "\n");

      // move one next
      if (last_str) free(last_str);
      last_en = tmp_en;
      last_str = tmp_str;
      last = tmp;
      tmp++;
    } // crawling path*/

    // insert saddle between outer structures
    if (outer) {
      RNAsaddle tmp(tbd.i, tbd.j, NOT_SURE);
      tmp.energy = en_fltoi(max_energy);
      tmp.str_ch = NULL;
      tmp.structure = make_pair_table(max_path->s);

      bool inserted = InsertUB(tmp, debug);

      if (output_saddles && inserted) {
        output_saddles->push_back(tmp);
      }
    }

    // free stuff
    //if (last_str) free(last_str);
    free_path(path);
  } // all doing while
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

