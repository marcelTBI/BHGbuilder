#include "DSUeval.h"

struct pq_height {
  int height;
  int distance;
  int number;
  int from;
  int saddle_num;

  bool operator<(const pq_height &right) const {
    if (height == right.height) {
      if (distance == right.distance) return number > right.number;
      else return distance > right.distance;
    } else return height > right.height;
  }

  pq_height(int h, int d, int n, int f, int s) {
    height = h;
    distance = d;
    number = n;
    from = f;
    saddle_num = s;
  }
};

vector<std::pair<int, int> > DSU::HeightSearch(int start, vector< set<edgeLL> > &edgesV_l)
{
  // define + init
  vector<int> heights(LM.size(), INT_MAX);
  vector<int> distance(LM.size(), INT_MAX);
  vector<bool> done(LM.size(), false);

  priority_queue<pq_height> que;

  // starting point -- all points from start
  done[start] = true;
  distance[start] = 0;
  heights[start] = LM[start].energy;
  for (set<edgeLL>::iterator it=edgesV_l[start].begin(); it!=edgesV_l[start].end(); it++) {
    que.push(pq_height(it->en, 1, it->goesTo(start), start, it->saddle));
  }

  // main loop -- dijkstra-like (take one with lowest energy, proceed it)
  while (!que.empty()) {
    // get next one to do
    pq_height pq = que.top(); que.pop();
    if (done[pq.number]) continue;

    //fprintf(stderr, "from %d to %d, en %d dist %d\n", pq.from+1, pq.number+1, pq.height, pq.distance);

    // write him:
    done[pq.number] = true;
    distance[pq.number] = pq.distance;
    heights[pq.number] = pq.height;

    // push next ones:
    for (set<edgeLL>::iterator it=edgesV_l[pq.number].begin(); it!=edgesV_l[pq.number].end(); it++) {
      int to = it->goesTo(pq.number);
      // if we already had him
      if (done[to]) {
        if (max(it->en, pq.height) < heights[to]) fprintf(stderr, "WRONG from: %d to: %d enLM %d enHeight %d\n", pq.number+1, to+1, pq.height, heights[to]);
        continue;
      }
      // push next ones
      //fprintf(stderr, "adding(%d): from %d to %d, en %d dist %d\n", (int)!done[to], pq.number+1, to+1, max(it->en, pq.height), pq.distance+1);
      if (!done[to]) {
        que.push(pq_height(max(it->en, pq.height), pq.distance+1, to, pq.number, it->saddle));
      }
    }
  }



  vector<std::pair<int, int> > res(LM.size());
  for (unsigned int i=0; i<heights.size(); i++) {
    //fprintf(stderr, "%d en: %d dist: %d prev: %d\n", i+1, heights[i], distance[i], previous[i]+1);
    res[i] = make_pair(heights[i], distance[i]);
  }

  return res;
}

void DSU::GetPath(int start, int stop, int maxkeep)
{
  char filename[100];
  sprintf(filename, "path%d_%d.path", start+1, stop+1);
  GetPath(start, stop, edgesV_l, filename, maxkeep);
}

void DSU::GetPath(int start, int stop,  vector< set<edgeLL> > &edgesV_l, char *filename, int maxkeep)
{
  // define + init
  vector<int> heights(LM.size(), INT_MAX);
  vector<int> distance(LM.size(), INT_MAX);
  vector<int> previous(LM.size(), -1);
  vector<int> saddle_num(LM.size(), -1);
  vector<bool> done(LM.size(), false);

  priority_queue<pq_height> que;

  // starting point -- all points from start
  done[start] = true;
  distance[start] = 0;
  heights[start] = LM[start].energy;
  for (set<edgeLL>::iterator it=edgesV_l[start].begin(); it!=edgesV_l[start].end(); it++) {
    que.push(pq_height(it->en, 1, it->goesTo(start), start, it->saddle));
  }

  // main loop -- dijkstra-like (take one with lowest energy, proceed it)
  while (!que.empty()) {
    // get next one to do
    pq_height pq = que.top(); que.pop();
    if (done[pq.number]) continue;

    //fprintf(stderr, "from %d to %d, en %d dist %d\n", pq.from+1, pq.number+1, pq.height, pq.distance);

    // write him:
    done[pq.number] = true;
    distance[pq.number] = pq.distance;
    previous[pq.number] = pq.from;
    heights[pq.number] = pq.height;
    saddle_num[pq.number] = pq.saddle_num;

    // end when we have reached destination
     if (pq.number==stop) break;

    // push next ones:
    for (set<edgeLL>::iterator it=edgesV_l[pq.number].begin(); it!=edgesV_l[pq.number].end(); it++) {
      int to = it->goesTo(pq.number);
      // if we already had him
      if (done[to]) {
        if (max(it->en, pq.height) < heights[to]) fprintf(stderr, "WRONG from: %d to: %d enLM %d enHeight %d\n", pq.number+1, to+1, pq.height, heights[to]);
        continue;
      }
      // push next ones
      //fprintf(stderr, "adding(%d): from %d to %d, en %d dist %d\n", (int)!done[to], pq.number+1, to+1, max(it->en, pq.height), pq.distance+1);
      if (!done[to]) {
        que.push(pq_height(max(it->en, pq.height), pq.distance+1, to, pq.number, it->saddle));
      }
    }
  }

  // retreive path:
  vector<int> lms;
  vector<int> sdd;
  int number = stop;
  while (number!=start) {
    lms.push_back(number);
    sdd.push_back(saddle_num[number]);
    number = previous[number];
  }
  lms.push_back(number);


  // write it down:
  FILE *fil;
  fil = fopen(filename, "w");
  int count = 0;
  if (fil) {
    fprintf(fil,"        %s\n", seq);
    for (int i=(int)lms.size()-1; i>0; i--) {
      fprintf(fil, "%6d  %s %7.2f %3d\n", lms[i]+1, LM[lms[i]].str_ch, LM[lms[i]].energy/100.0, HammingDist(LM[lms[i]].structure, saddles[sdd[i-1]].structure));
      if (maxkeep) {
        path_t *tmp = get_path(seq, LM[lms[i]].str_ch, saddles[sdd[i-1]].str_ch, maxkeep);
        path_t *path = tmp+1;
        while (path && path->s && (path+1) && (path+1)->s) {
          fprintf(fil, "        %s %7.2f\n", path->s, path->en);
          path++;
          count++;
        }
        free_path(tmp);
      }
      fprintf(fil, "%6dS %s %7.2f %3d\n", sdd[i-1]+1, saddles[sdd[i-1]].str_ch, saddles[sdd[i-1]].energy/100.0, HammingDist(LM[lms[i-1]].structure, saddles[sdd[i-1]].structure));
      if (maxkeep) {
        path_t *tmp = get_path(seq, saddles[sdd[i-1]].str_ch, LM[lms[i-1]].str_ch, maxkeep);
        path_t *path = tmp+1;
        while (path && path->s && (path+1) && (path+1)->s) {
          fprintf(fil, "        %s %7.2f\n", path->s, path->en);
          path++;
          count++;
        }
        free_path(tmp);
      }
    }
    fprintf(fil, "%6d  %s %7.2f\n", lms[0]+1, LM[lms[0]].str_ch, LM[lms[0]].energy/100.0);
    if (maxkeep)  fprintf(fil, "Path from %6d to %6d: %4d local minima, %d structures, %6.2f kcal/mol highest point.\n", start+1, stop+1, (int)lms.size(), count + (int)lms.size() + (int)sdd.size(), heights[stop]/100.0);
    else          fprintf(fil, "Path from %6d to %6d: %4d local minima, %6.2f kcal/mol highest point.\n", start+1, stop+1, (int)lms.size(), heights[stop]/100.0);
    fclose(fil);
  } else {
    fprintf(stderr, "Unable to open file %s!\n", filename);
  }
}

void DSU::EHeights(FILE *heights, bool full)
{
  if (!full) {
    for (unsigned int i=0; i< saddles.size(); i++) {
      //fprintf(heights, "%s %.2f %s %.2f %s %.2f\n", LM[saddles[i].lm1].str_ch, LM[saddles[i].lm1].energy/100.0, LM[saddles[i].lm2].str_ch, LM[saddles[i].lm2].energy/100.0, saddles[i].str_ch, saddles[i].energy/100.0);
      fprintf(heights, "%4d %4d %.2f\n", saddles[i].lm1+1, saddles[i].lm2+1, saddles[i].energy/100.0);
    }
  } else {
    vector<vector<std::pair<int, int> > > res(LM.size());
    for (unsigned int i=0; i<LM.size(); i++) {
      res[i] = HeightSearch(i, edgesV_l);
    }
    for (unsigned int i=0; i<res.size(); i++) {
      for (unsigned int j=i+1; j<res[i].size(); j++) {
        // print lines: energy heights: from_node, to_node, energy height, distance
        fprintf(heights, "%4d %4d %.2f %3d\n", i+1, j+1, res[i][j].first/100.0, res[i][j].second);
      }
    }
  }
}


void DSU::ERank(FILE *rank, bool barr, bool out_conns)
{
  // saddle point energies + distances
  vector<vector<std::pair<int, int> > > res(number_lm);
  for (int i=0; i<number_lm; i++) {
    res[i] = HeightSearch(i, edgesV_l);
  }

  struct energy_pair {
    int i,j;
    int barrier;
    bool operator<(const energy_pair &scnd) const {
      return  barrier>scnd.barrier;
    }
  };

  int n = number_lm;

  priority_queue<energy_pair> saddles;
  for (int i=0; i<n; i++) {
    //if (nodes[i].father!=-1) continue;
    for (int j=i+1; j<n; j++) {
      //if (nodes[j].father!=-1) continue;
      if (res[i][j].first<1e8) {
        energy_pair ep;
        ep.barrier = res[i][j].first;
        if (barr) {
          ep.barrier -= max(LM[i].energy, LM[j].energy);
        }
        ep.i=i;
        ep.j=j;
        saddles.push(ep);
      }
    }
  }

  // variables
  int count = 1;
  int total_en = 0.0;

  // build a tree.
  UF_set ufset;
  ufset.enlarge_parent(n);
  while (!saddles.empty()) {
    energy_pair ep = saddles.top();
    saddles.pop();

    //fprintf(stderr, "%4d %4d %7.2f\n", ep.i, ep.j, ep.barrier);

    if (!ufset.joint(ep.i, ep.j)) {

      int i=ufset.find(ep.i);
      int j=ufset.find(ep.j);

      int father = min(i, j);
      int child = max(i, j);

      fprintf(rank, "%4d %8.2f", count++, ep.barrier/100.0);
      total_en += ep.barrier;
      if (out_conns) {
        fprintf(rank, " %4d %4d", ep.i+1, ep.j+1);
      }
      fprintf(rank, "\n");

      // finally join them
      ufset.union_set(father, child);
    }
  }
  fprintf(stderr, "average energy = %6.2f/%4d = %11.8f\n", total_en/100.0, (count-1), total_en/100.0/(double)(count-1));
}

