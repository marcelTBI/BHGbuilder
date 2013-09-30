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
  #include "read_epars.h"
}

#include "BHGbuilder.h"
#include "hash_util.h"

#include <algorithm>

using namespace std;

Opt::Opt(gengetopt_args_info &args_info)
{
  // read parameter file
  if (args_info.paramFile_given) {
    read_parameter_file(args_info.paramFile_arg);
  }

  // dangle setup
  if (args_info.dangles_given) {
    dangles = args_info.dangles_arg;
  }

  this->noLP = args_info.noLP_flag;
  this->shifts = args_info.shift_flag;
  this->saddle_conn = !args_info.noSaddle_flag;
  this->flood_num = args_info.floodMax_arg;
  this->flood_height = args_info.floodHeight_arg;
  this->debug = args_info.debug_flag;
  this->maxkeep = args_info.depth_arg;
  this->num_threshold = args_info.num_threshold_arg;
  this->outer = args_info.outer_flag;
  this->repre_portion = args_info.cluster_repre_arg;
  this->fbarrier = !args_info.cluster_fsaddle_flag;
  this->no_conn = args_info.no_conn_flag;

  //update_fold_params();
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

inline double rate(int en_from, int en_to, double _kT) {
  return 1.0*exp(-(en_to-en_from)/100.0/_kT);
}

void edgeAdv::FillRate(mode_rates mode, double _kT, vector<RNAlocmin> &LM)
{
  switch (mode) {
  case VERTEX_CONTR_SUM: {
    double cum_rate_toi = 0.0;
    double cum_rate_toj = 0.0;
    //  cum_rate = 1/rat + 1/rat + ...
    for (int k=0; k<(int)saddles.size(); k++) {
      int lms_from = (k==0?i:lms[k-1]);
      int lms_to = (k==(int)saddles.size()-1?j:lms[k-1]);
      cum_rate_toj += 1/rate(LM[lms_from].energy, energies[k], _kT);
      cum_rate_toi += 1/rate(LM[lms_to].energy, energies[k], _kT);
    }
    rate_toi = 1/cum_rate_toi;
    rate_toj = 1/cum_rate_toj;
    break;
    }
  case VERTEX_CONTR:
  case NO_CONTR:
  default:
    rate_toj = rate(LM[i].energy, max_height, _kT);
    rate_toi = rate(LM[j].energy, max_height, _kT);
    break;
  }
}

Graph::Graph(set<edgeLL> &edges, vector<RNAlocmin> &LM, mode_rates mode)
  :LM(LM), lowest(edge_comp_adv(LM, mode==EDGE_CONTR_MAX)) //shallow copy
{
  this->mode = mode;
  this->max_node = this->number_lm = LM.size();

  // adjacency creation
  for (set<edgeLL>::iterator it=edges.begin(); it!=edges.end(); it++) {
    int i=min(it->i, it->j);
    int j=max(it->i, it->j);
    //edgeAdv ea(i, j, it->en, it->saddle);
    edgeAdv &edge = (adjacency.insert(make_pair(make_pair(i,j), edgeAdv(i, j, it->en, it->saddle)))).first->second;
    if (mode == EDGE_CONTR_MAX || mode == EDGE_CONTR_MIN) lowest.push(edge);
  }

  // edge contraction works with nodes:
  if (mode == EDGE_CONTR_MAX || mode == EDGE_CONTR_MIN) {
    ufset.enlarge_parent(number_lm);
  }
}

// JOIN TWO EDGES TO MAKE ONE LONGER
int Graph::Join(const edgeAdv &src, const edgeAdv &dst, int joining_node, edgeAdv &res)
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
  res.lms.push_back(joining_node);
  res.lms.insert(end(res.lms), begin(dst.lms), end(dst.lms));

  // max_height
  res.max_height = max(res.max_height, dst.max_height);
  /*
  for(unsigned int i=0; i<dst.energies.size(); i++) {
    res.max_height = max(res.max_height, dst.energies[i]);
  }*/

  return 0;
}

static edge_comp ec;
// two parallel edges come together
bool Graph::AddEdges(const edgeAdv &found, edgeAdv &res, mode_rates mode)
{
  bool changed = false;
  if (!ec(found, res)) {
    changed = true;
  }

  return changed;
}

int Graph::RemoveLowestEdge() {

  int count = 0;

  // get new lowest edge:
  edgeAdv edge = lowest.top(); lowest.pop();

  // check if we will reduce nodes:
  if (ufset.joint(edge.i, edge.j)) {
    // skip
    return 0;
  }

  // fix:
  int i = ufset.find(edge.i);
  int j = ufset.find(edge.j);
  edge.i = min(i,j);
  edge.j = max(i,j);

  ufset.union_set(edge.i, edge.j);

  // remove that edge
  adjacency.erase(make_pair(edge.i, edge.j));

  // add edges from j to i
  for (int i=0; i<max_node; i++) {
    int point = edge.j;
    map<std::pair<int, int>, edgeAdv>::iterator it;
    if ((it = adjacency.find(make_pair(min(point, i), max(point,i))))!=adjacency.end()) {
      edgeAdv ev = it->second;
      adjacency.erase(it);
      int goes = ev.goesTo(point);
      ev.i = min(goes, edge.i);
      ev.j = max(goes, edge.i);

      if ((it = adjacency.find(make_pair(ev.i, ev.j))) != adjacency.end()) {
        bool changed = AddEdges(it->second, ev, mode);
        if (changed) adjacency.insert(make_pair(make_pair(ev.i, ev.j), ev)); // insert it back modified
      } else {
        adjacency.insert(make_pair(make_pair(ev.i, ev.j), ev));
      }
      count++;

    }
  }
  number_lm--;
  if (count==0) return -1;
  return count;
}


int Graph::RemoveLastPoint() {

  // remove last point;
  int point = number_lm-1;

  int count = 0;

  // candidates:
  vector<edgeAdv> candidates;
  for (int i=0; i<number_lm; i++) {
    map<std::pair<int, int>, edgeAdv>::iterator it;
    if ((it = adjacency.find(make_pair(min(point, i), max(point,i))))!=adjacency.end()) {
      candidates.push_back(it->second);
      adjacency.erase(it);
    }
  }

  for (int i=0; i<(int)candidates.size(); i++) {
    for (int j=i+1; j<(int)candidates.size(); j++) {
      edgeAdv res = candidates[i];
      int status = Join(candidates[i], candidates[j], point, res);
      if (status==0) {
        // if we have not found already some edge:
        map<std::pair<int, int>, edgeAdv>::iterator it;
        if  ((it = adjacency.find(make_pair(res.i, res.j)))==adjacency.end()) {
          adjacency.insert(make_pair(make_pair(res.i, res.j), res));
        } else {
          bool changed = AddEdges(it->second, res, mode);
          if (changed) {
            adjacency.insert(make_pair(make_pair(res.i, res.j), res));
          }
        }

        count++;
      }
    }
  }

  // last was erased
  number_lm--;

/*
  // try to join every possible pair of edges:
  for (set<edgeAdv>::iterator it=adjacency[point].begin(); it!=adjacency[point].end(); it++) {
    set<edgeAdv>::iterator it2 = it; it2++;
    set<edgeAdv>::iterator found;
    for (; it2!=adjacency[point].end(); it2++) {
      edgeAdv res(*it);
      int status = Join(*it, *it2, point, res);
      if (status==0) {
        // if we have not found already some edge:
        if ((found = adjacency[res.j].find(res))==adjacency[res.j].end()) {
          adjacency[res.j].insert(res);
        } else {
          bool changed = AddEdges(*found, res, mode);
          if (changed) {
            adjacency[res.j].erase(found);
            adjacency[res.j].insert(res);
          }
        }

        count++;
      }
    }
  }*/

  return count; //  returns change in edge count
}

void Graph::PrintDot(char *filename, bool dot_prog, bool print, char *file_print)
{
  fprintf(stderr, "prinitng graph...\n");
  //open file
  FILE *dot;
  dot = fopen(filename, "w");
  if (dot) {
    fprintf(dot, "Graph G {\n\tnode [width=0.1, height=0.1, shape=circle];\n");
    if (mode == VERTEX_CONTR) {
      //nodes LM:
      for (int i=0; i<number_lm; i++) {
        switch (LM[i].type) {
          case NORMAL:
          case NORM_CF: fprintf(dot, "\"%d\" [label=\"%d\"]\n", i+1, i+1); break;
          case EE_DSU: fprintf(dot, "\"%d\" [label=\"%d\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, rgb(0, 0, 255), rgb(0, 0, 255)); break;
          case EE_COMP: fprintf(dot, "\"%d\" [label=\"%d\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, rgb(255, 0, 0), rgb(255, 0, 0)); break;
        }
      }
    } else if (mode == EDGE_CONTR_MAX || mode == EDGE_CONTR_MIN) {
      // nodes LM:
      for (int i=0; i<max_node; i++) {
        if (ufset.count(i)>0) {
          set<int> childs = ufset.get_children(i);

          fprintf(dot, "\"%d\" [label=\"%d", i+1, i+1);

          set<int>::iterator it = childs.begin(); it++;
          for (;it!=childs.end(); it++) {
            fprintf(dot, " %d", (*it)+1);
          }

          switch (LM[i].type) {
            case NORMAL:
            case NORM_CF: fprintf(dot, "\"]\n"); break;
            case EE_DSU: fprintf(dot, "\", color=\"%s\", fontcolor=\"%s\"]\n", rgb(0, 0, 255), rgb(0, 0, 255)); break;
            case EE_COMP: fprintf(dot, "\", color=\"%s\", fontcolor=\"%s\"]\n", rgb(255, 0, 0), rgb(255, 0, 0)); break;
          }
        }
      }
    }

    fprintf(dot, "\n");

    // edges
    for (map<std::pair<int, int>, edgeAdv>::iterator it = adjacency.begin(); it!=adjacency.end(); it++) {
      char length[10]="";
     bool component = (it->second.length()>1);
     if (component) sprintf(length, "(%d)", it->second.length());
     fprintf(dot, "\"%d\" -- \"%d\" [label=\"%.2f%s\", color=\"%s\", fontcolor=\"%s\"]\n", (it->second.i)+1, (it->second.j)+1, it->second.max_height/100.0, length, (component?rgb(255, 0, 0):rgb(0, 0, 0)), (component?rgb(255, 0, 0):rgb(0, 0, 0)));
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


void Graph::PrintRates(FILE *rates, double temp)
{
  double _kT = 0.00198717*(273.15 + temp);
  fprintf(stderr, "prinitng rates...\n");

  // create matrix
  vector<vector<double> > mat_rates(number_lm);
  for (int i=0; i<number_lm; i++) {
    mat_rates[i].resize(number_lm, 0.0);
  }

  // fill rates matrix
  /*for (int i=0; i<number_lm; i++) {
    for (set<edgeAdv>::iterator it=adjacency[i].begin(); it!=adjacency[i].end(); it++) {
      //fprintf(stderr, "%d %d %f(%d)\n", i, j, it->max_height/100.0, it->length());
      mat_rates[it->i][it->j] = 1.0*exp(-(it->max_height-LM[it->i].energy)/100.0/_kT);
      mat_rates[it->j][it->i] = 1.0*exp(-(it->max_height-LM[it->j].energy)/100.0/_kT);
    }
  }*/

  if (mode == VERTEX_CONTR || mode == NO_CONTR || mode == VERTEX_CONTR_SUM) {
    for (map<std::pair<int, int>, edgeAdv>::iterator it = adjacency.begin(); it!=adjacency.end(); it++) {
      it->second.FillRate(mode, _kT, LM);
      mat_rates[it->second.i][it->second.j] = it->second.rate_toj;
      mat_rates[it->second.j][it->second.i] = it->second.rate_toi;
    }
  } else if (mode == EDGE_CONTR_MAX || mode == EDGE_CONTR_MIN) {
    map<int, int> inverted = ufset.get_invert();
    for (map<std::pair<int, int>, edgeAdv>::iterator it = adjacency.begin(); it!=adjacency.end(); it++) {
      int i = inverted[it->second.i];
      int j = inverted[it->second.j];
      mat_rates[i][j] = 1.0*exp(-(it->second.max_height-LM[it->second.i].energy)/100.0/_kT);
      mat_rates[j][i] = 1.0*exp(-(it->second.max_height-LM[it->second.j].energy)/100.0/_kT);
    }
  }

  // print rate matrix
  for (int i=0; i<number_lm; i++) {
    for (int j=0; j<number_lm; j++) {
      fprintf(rates, "%11.5g ", mat_rates[i][j]);
    }
    fprintf(rates, "\n");
  }

  // print map of contracted nodes if EDGE_CONTR was specified:
  if (mode == EDGE_CONTR_MAX || mode == EDGE_CONTR_MIN) {
    FILE *contr = fopen(mode == EDGE_CONTR_MAX?"edge_contr.mapM":"edge_contr.mapE", "w");
    if (contr) {
      int j=1;
      for (int i=0; i<ufset.size(); i++) {
        if (ufset.count(i)>0) {
          set<int> edges = ufset.get_children(i);
          fprintf(contr, "%4d   ", j);
          j++;
          for (set<int>::iterator it = edges.begin(); it!=edges.end(); it++) {
            fprintf(contr, "%d ", (*it)+1);
          }
          fprintf(contr, "\n");
        }
      }
      fclose(contr);
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
