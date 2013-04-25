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
