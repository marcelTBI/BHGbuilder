#include <stdio.h>
#include <string.h>
#include <stdlib.h>

extern "C" {
  #include "fold.h"
  #include "findpath.h"
  #include "move_set.h"
}

#include "BHGbuilder.h"
#include "RNAutils.h"
#include "hash_util.h"

#include <algorithm>

using namespace std;

// variables and data strucutres for flooding (LM to LM)
bool foundDS;
RNAstruc curr;
int currNum;
int threshold;
bool gl_debug;
unordered_map<RNAstruc, int, hash_fncts, hash_eq> flood_hash;
unordered_map<RNAstruc, int, hash_fncts, hash_eq>::iterator it;
priority_queue<RNAstruc, vector<RNAstruc>, RNAstruc_rev> flood_queue;

int funct(struct_en *moved, struct_en *current) {

  if (moved->energy >= curr.energy && moved->energy < threshold) {
    RNAstruc tmp;
    tmp.energy = moved->energy;
    tmp.structure = moved->structure;
    if ((it = flood_hash.find(tmp))!=flood_hash.end()) {
      // found DS!
      if (it->second != currNum) {
        if (!foundDS || current->energy > moved->energy) {
          foundDS = true;
          copy_arr(current->structure, moved->structure);
          current->energy = moved->energy;
        }

        if (gl_debug) {
          fprintf(stderr, "FOUND saddle: %s %6.2f\n", pt_to_str(current->structure).c_str(), current->energy/100.0);
        }

        return 0;
      } // else nothing
    } else {
      // insert strucure we havent see yet (and it is in range)
      RNAstruc to_ins;
      to_ins.energy = moved->energy;
      to_ins.structure = allocopy(moved->structure);
      flood_hash.insert(make_pair(to_ins, currNum));
      flood_queue.push(to_ins);
    }
  }

  return 0;
}

int DSU::FloodUp(RNAlocmin &i, RNAlocmin &j, RNAsaddle &saddle, Opt &opt, bool debug)
{
  // threshold set
  if (opt.flood_height == 0) {
    if (saddle.type == NOT_SURE) threshold = INT_MAX; // we cannot specify the threshold
    else threshold = saddle.energy;
  } else {
    threshold = min(saddle.energy, min(i.energy, j.energy) + opt.flood_height);
  }

  gl_debug = debug;

  // debug output
  if (gl_debug) {
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
  foundDS = false;

  // main loop
  curr = flood_queue.top();
  flood_queue.pop();
  currNum = flood_hash[curr];
  while (curr.energy < threshold && (int)flood_hash.size() < opt.flood_num) {

    // debug output
    if (debug) {
      fprintf(stderr, "browsing    : %s %6.2f\n", pt_to_str(curr.structure).c_str(), curr.energy/100.0);
    }

    // start browsing
    curr.energy = browse_neighs(seq, curr.structure, s0, s1, 0, opt.shifts, opt.noLP, funct);

    // found DS - quitting
    if (foundDS) {
      break;
    }

    // get another struct
    if (flood_queue.empty()) {
      break;
    }
    curr = flood_queue.top();
    flood_queue.pop();
    currNum = flood_hash[curr];
  }

  // if we have searched up to known direct saddle threshold
  if (flood_queue.empty()) {
    if (saddle.type != NOT_SURE) {
      res = 2;
    }
  }

  //int size = flood_hash.size();

  // we did end sucesfully
  if (foundDS) {
    copy_arr(saddle.structure, curr.structure);
    saddle.energy = curr.energy;
    saddle.recompute_str();
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

// variables and data strucutres for flooding (saddle to saddle)
bool found_saddle;
RNAstruc to_find;
//int threshold;
//bool gl_debug;
unordered_set<RNAstruc, hash_fncts, hash_eq> flood_set;
unordered_set<RNAstruc, hash_fncts, hash_eq>::iterator it_set;
//priority_queue<RNAstruc, vector<RNAstruc>, RNAstruc_rev> flood_queue;

int funct_saddle(struct_en *moved, struct_en *current) {

  if (moved->energy >= curr.energy) {
    RNAstruc tmp;
    tmp.energy = moved->energy;
    tmp.structure = moved->structure;
    if (to_find == tmp) {
      // found it!!!
      found_saddle = true;
      return 1;
    }
    if (moved->energy < threshold) {

      if ((it_set = flood_set.find(tmp))!=flood_set.end()) {
        // nothing...
        return 0;
      } else {
        // add it
        RNAstruc to_ins;
        to_ins.energy = moved->energy;
        to_ins.structure = allocopy(moved->structure);
        flood_set.insert(to_ins);
        flood_queue.push(to_ins);
        return 0;
      }
    }
  }

  return 0;
}

bool DSU::FloodSaddle(RNAsaddle &saddle_lower, RNAsaddle &saddle_higher, Opt &opt, bool debug)
{
  // threshold set
  threshold = saddle_higher.energy; // not very good to use floodHeight restriction here :-)
  gl_debug = debug;

  // debug output
  if (debug) {
    fprintf(stderr, "\nsadle1      : %s %6.2f\n", pt_to_str(saddle_lower.structure).c_str(), saddle_lower.energy/100.0);
    fprintf(stderr, "sadle2      : %s %6.2f\n", pt_to_str(saddle_higher.structure).c_str(), saddle_higher.energy/100.0);
  }

  // init
  to_find = saddle_higher;
  RNAstruc tmp;
  tmp.energy = saddle_lower.energy;
  tmp.structure = allocopy(saddle_lower.structure);
  flood_set.insert(tmp);
  while (!flood_queue.empty()) flood_queue.pop();
  flood_queue.push(tmp);

  // end
  bool res = false;
  found_saddle = false;

  // main loop
  curr = flood_queue.top();
  flood_queue.pop();
  while (curr.energy < threshold && (int)flood_set.size() < opt.flood_num) {

    // debug output
    if (debug) {
      fprintf(stderr, "browsing2   : %s %6.2f\n", pt_to_str(curr.structure).c_str(), curr.energy/100.0);
    }

    // start browsing
    curr.energy = browse_neighs(seq, curr.structure, s0, s1, 0, opt.shifts, opt.noLP, funct_saddle);

    // found DS - quitting
    if (found_saddle) {
      break;
    }

    // get another struct
    if (flood_queue.empty()) break;
    curr = flood_queue.top();
    flood_queue.pop();
  }

  // we did end sucesfully
  if (found_saddle) {
    res = true;
  }

  // free
  while (!flood_queue.empty()) flood_queue.pop();
  for (it_set=flood_set.begin(); it_set!=flood_set.end(); it_set++) {
    if (it_set->structure) free(it_set->structure);
  }
  flood_set.clear();
  return res;
}
