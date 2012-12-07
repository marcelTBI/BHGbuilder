#include "DSUeval.h"

struct pq_height {
  int height;
  int distance;
  int number;
  int from;

  bool operator<(const pq_height &right) const {
    if (height == right.height) {
      if (distance == right.distance) return number > right.number;
      else return distance > right.distance;
    } else return height > right.height;
  }

  pq_height(int h, int d, int n, int f) {
    height = h;
    distance = d;
    number = n;
    from = f;
  }
};

vector<std::pair<int, int> > DSU::HeightSearch(int start, vector< set<edgeLL> > &edgesV_l)
{
  // define + init
  vector<int> heights(LM.size(), INT_MAX);
  vector<int> distance(LM.size(), INT_MAX);
  vector<int> previous(LM.size(), -1);
  vector<bool> done(LM.size(), false);

  priority_queue<pq_height> que;

  // starting point
  done[start] = true;
  distance[start] = 0;
  heights[start] = LM[start].energy;
  for (set<edgeLL>::iterator it=edgesV_l[start].begin(); it!=edgesV_l[start].end(); it++) {
    que.push(pq_height(it->en, 1, it->goesTo(start), start));
  }

  // main loop
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
        que.push(pq_height(max(it->en, pq.height), pq.distance+1, to, pq.number));
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
        fprintf(heights, "%4d %4d %.2f %3d\n", i+1, j+1, res[i][j].first/100.0, res[i][j].second);
      }
    }
  }
}
