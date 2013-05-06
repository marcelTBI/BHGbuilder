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

#include "DSUeval.h"
#include "hash_util.h"

#include <algorithm>

using namespace std;

Opt::Opt(bool noLP, bool shifts, bool saddle, int num, int height, bool debug, int maxkeep, int num_threshold, bool outer, float repre_portion, bool fbarrier)
{
  this->noLP = noLP;
  this->shifts = shifts;
  this->saddle_conn = saddle;
  this->flood_num = num;
  this->flood_height = height;
  this->debug = debug;
  this->maxkeep = maxkeep;
  this->num_threshold = num_threshold;
  this->outer = outer;
  this->repre_portion = repre_portion;
  this->fbarrier = fbarrier;
}

void RNAstruc::freeMEM() {
  if (structure) free(structure);
  if (str_ch) free(str_ch);
}

void RNAstruc::recompute_str()
{
  if (str_ch) free(str_ch);
  if (!structure) str_ch = NULL;
  else str_ch = pt_to_char(structure);
}

Graph::Graph(set<edgeLL> &edges, vector<RNAlocmin> &LM)
{

  this->max_node = this->number_lm = LM.size();
  adjacency.resize(max_node);
  for (int i=0; i<max_node; i++) {
    adjacency[i].resize(max_node);
  }

  // adjacency creation
  for (set<edgeLL>::iterator it=edges.begin(); it!=edges.end(); it++) {
    int i=min(it->i, it->j);
    int j=max(it->i, it->j);
    adjacency[i][j].insert(edgeAdv(i, j, it->en, it->saddle));
  }

  // just shallow copy
  this->LM = LM;
  /*for (unsigned int i=0; i<LM.size(); i++) {
    this->LM[i] = LM[i];
  }*/
}

int Graph::Join(edgeAdv &src, edgeAdv &dst, int joining_node, edgeAdv &res)
{
  res = src;

  // in and out node
  int src_node, dst_node;
  if (src.i == joining_node) {
    src_node = src.j;
  } else {
    src_node = src.i;
    if (src.j != joining_node) {fprintf(stderr, "ERROR: joining node not found(src)\n"); return -1;}
  }
  if (dst.i == joining_node) {
    dst_node = dst.j;
  } else {
    dst_node = dst.i;
    if (dst.j != joining_node) {fprintf(stderr, "ERROR: joining node not found(dst)\n"); return -1;}
  }

  // edges that are parrallel
  if (dst_node == src_node) return 1;

  // start and end point
  res.i = min(dst_node, src_node);
  res.j = max(dst_node, src_node);

  // energies & saddles
  res.saddles.insert(res.saddles.end(), dst.saddles.begin(), dst.saddles.end());
  res.energies.insert(res.energies.end(), dst.energies.begin(), dst.energies.end());

  // max_height
  for(unsigned int i=0; i<dst.energies.size(); i++) {
    res.max_height = max(res.max_height, dst.energies[i]);
  }

  return 0;
}

int Graph::RemovePoint(int point, int keep) {

  // sets the number of active nodes
  if (number_lm > point) number_lm = point;

  int count = 0;
  // find candidates + delete old ones
  vector<edgeAdv> candidates;
  for (int i=0; i<max_node; i++) {
    for (multiset<edgeAdv, edge_comp>::iterator it = adjacency[min(i,point)][max(i,point)].begin(); it!=adjacency[min(i,point)][max(i,point)].end(); it++) {
      edgeAdv res = *it;
      candidates.push_back(res);
    }
    count-= adjacency[min(i,point)][max(i,point)].size();
    adjacency[min(i,point)][max(i,point)].clear();
  }

  // get there "keep" new ones:
  for (unsigned int i=0; i<candidates.size(); i++) {
    for (unsigned int j=i+1; j<candidates.size(); j++) {
      edgeAdv res(candidates[i]);
      int status = Join(candidates[i], candidates[j], point, res);
      if (status==0) {
        adjacency[res.i][res.j].insert(res);
        count++;
        if ((int)adjacency[res.i][res.j].size()>keep) {
          multiset<edgeAdv, edge_comp>::iterator it = adjacency[res.i][res.j].end(); it--;
          adjacency[res.i][res.j].erase(it);
          count--;
        }
      }
    }
  }

  return count; //  returns change in edge count
}

void Graph::PrintDot(char *filename, bool dot_prog, bool print, char *file_print)
{
  //open file
  FILE *dot;
  dot = fopen(filename, "w");
  if (dot) {
    fprintf(dot, "Graph G {\n\tnode [width=0.1, height=0.1, shape=circle];\n");
    //nodes LM:
    for (int i=0; i<number_lm; i++) {
      switch (LM[i].type) {
        case NORMAL:
        case NORM_CF: fprintf(dot, "\"%d\" [label=\"%d\"]\n", i+1, i+1); break;
        case EE_DSU: fprintf(dot, "\"%d\" [label=\"%d\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, rgb(0, 0, 255), rgb(0, 0, 255)); break;
        case EE_COMP: fprintf(dot, "\"%d\" [label=\"%d\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, rgb(255, 0, 0), rgb(255, 0, 0)); break;
      }
    }
    fprintf(dot, "\n");

    // edges l-l
    for (int i=0; i<max_node; i++) {
      for (int j=0; j<max_node; j++) {
        for (multiset<edgeAdv>::iterator it=adjacency[i][j].begin(); it!=adjacency[i][j].end(); it++) {
          char length[10]="";
          bool component = (it->length()>1);
          if (component) sprintf(length, "(%d)", it->length());
          fprintf(dot, "\"%d\" -- \"%d\" [label=\"%.2f%s\", color=\"%s\", fontcolor=\"%s\"]\n", (it->i)+1, (it->j)+1, it->max_height/100.0, length, (component?rgb(255, 0, 0):rgb(0, 0, 0)), (component?rgb(255, 0, 0):rgb(0, 0, 0)));
        }
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


void Graph::PrintRates(FILE *rates, double temp, mode_rates mode)
{
  double _kT = 0.00198717*(273.15 + temp);

  // create matrix
  vector<vector<double> > mat_rates(number_lm);
  for (int i=0; i<number_lm; i++) {
    mat_rates[i].resize(number_lm, 0.0);
  }

  // fill rates matrix
  for (int i=0; i<number_lm; i++) {
    for (int j=i+1; j<number_lm; j++) {
      for (multiset<edgeAdv, edge_comp>::iterator it=adjacency[i][j].begin(); it!=adjacency[i][j].end(); it++) {
        //fprintf(stderr, "%d %d %f(%d)\n", i, j, it->max_height/100.0, it->length());
        switch (mode) {
          case JUST_BEST:
            mat_rates[i][j] = 1.0*exp(-(it->max_height-LM[i].energy)/100.0/_kT);
            mat_rates[j][i] = 1.0*exp(-(it->max_height-LM[j].energy)/100.0/_kT);
            break;
          case ADDITIVE:
            mat_rates[i][j] += 1.0*exp(-(it->max_height-LM[i].energy)/100.0/_kT);
            mat_rates[j][i] += 1.0*exp(-(it->max_height-LM[j].energy)/100.0/_kT);
            break;
        }
        if (mode==JUST_BEST) break;
        //fprintf(stderr, "%d %d %f(%d)\n", i, j, it->max_height/100.0, it->length());
      }
    }
  }

  // print rate matrix
  for (int i=0; i<number_lm; i++) {
    for (int j=0; j<number_lm; j++) {
      fprintf(rates, "%11.5g ", mat_rates[i][j]);
    }
    fprintf(rates, "\n");
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
