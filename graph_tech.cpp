#include "BHGbuilder.h"

static char optimality = 'B';

class pq_height {
public:
  int height;
  int distance;
  int number;
  int from;
  int saddle_num;
  double rate;

  // probability/time dependent likelihood (the higher the better)
  double criterion;
  double time_min;

  int bdist;
  PKNOT_TYPES pt;

public:
  bool operator<(const pq_height &right) const {

    switch (optimality) {
      case 'B': if (height != right.height) return height > right.height;
      case 'L': if (bdist != right.bdist) return bdist > right.bdist;
      case 'M': if (distance != right.distance) return distance > right.distance;
      case 'R': if (rate != right.rate) return rate < right.rate;
      case 'T':
      case 'P': if (criterion != right.criterion) return criterion < right.criterion;
    }

    // default ordering:
    if (height == right.height) {
      if (bdist == right.bdist) {
        if (distance == right.distance) {
          return rate < right.rate;
        } else return distance > right.distance;
      } else return bdist > right.bdist;
    } else return height > right.height;
  }

  pq_height(int h, int d, int n, int f, int s, double rat = 0.0, int bist = 0, double other = 0.0, double time = 0.0, PKNOT_TYPES pt = PKNOT_TYPES()) {
    height = h;
    distance = d;
    number = n;
    from = f;
    saddle_num = s;
    rate = rat;
    bdist = bist;
    criterion = other;
    time_min = time;
    this->pt = pt;
  }
};


vector<HeightData> DSU::HeightSearch(int start, vector< set<edgeLL> > &edgesV_l)
{
  optimality = 'B';

  // define + init
  vector<int> heights(LM.size(), INT_MAX);
  vector<int> distance(LM.size(), INT_MAX);
  vector<int> bdist(LM.size(), INT_MAX);
  vector<bool> done(LM.size(), false);
  vector<PKNOT_TYPES> pknot_types(LM.size());

  priority_queue<pq_height> que;

  // starting point -- all points from start
  done[start] = true;
  distance[start] = 0;
  heights[start] = LM[start].energy;
  bdist[start] = 0;
  for (set<edgeLL>::iterator it=edgesV_l[start].begin(); it!=edgesV_l[start].end(); it++) {
    int to = it->goesTo(start);
    int sadd = it->saddle;
    int hd = HammingDist(LM[to].structure, saddles[sadd].structure) + HammingDist(saddles[sadd].structure, LM[start].structure);
    PKNOT_TYPES pt(Identify_AllPK(LM[start].structure));
    pt.Add(Identify_AllPK(LM[to].structure));
    pt.Add(Identify_AllPK(saddles[sadd].structure));
    que.push(pq_height(it->en, 1, it->goesTo(start), start, it->saddle, 0.0, hd, 0.0, 0.0, pt));
  }

  // main loop -- dijkstra-like (take one with lowest energy, proceed it)
  while (!que.empty()) {
    // get next one to do
    pq_height pq = que.top(); que.pop();
    if (done[pq.number]) continue;

    //fprintf(stderr, "from %6d to %6d, en %6.2f dist %6d\n", pq.from+1, pq.number+1, pq.height/100.0, pq.distance);

    // write him:
    done[pq.number] = true;
    distance[pq.number] = pq.distance;
    heights[pq.number] = pq.height;
    bdist[pq.number] = pq.bdist;
    pknot_types[pq.number] = pq.pt;

    // push next ones:
    for (set<edgeLL>::iterator it=edgesV_l[pq.number].begin(); it!=edgesV_l[pq.number].end(); it++) {
      int to = it->goesTo(pq.number);
      int sadd = it->saddle;
      // if we already had him
      if (done[to]) {
        if (max(it->en, pq.height) < heights[to]) fprintf(stderr, "WRONG from: %d to: %d enLM %d enHeight %d\n", pq.number+1, to+1, pq.height, heights[to]);
        continue;
      }
      // push next ones
      //fprintf(stderr, "adding(%d): from %d to %d, en %d dist %d\n", (int)!done[to], pq.number+1, to+1, max(it->en, pq.height), pq.distance+1);
      if (!done[to]) {
        int hd = HammingDist(LM[to].structure, saddles[sadd].structure) + HammingDist(saddles[sadd].structure, LM[pq.number].structure);
        //fprintf(stderr, "%4d %4d %4d \n", pq.number, to, hd);
        PKNOT_TYPES pt = pq.pt;
        pt.Add(Identify_AllPK(LM[to].structure));
        pt.Add(Identify_AllPK(saddles[sadd].structure));
        que.push( pq_height(max(it->en, pq.height), pq.distance+1, to, pq.number, it->saddle, 0.0, pq.bdist+hd, 0.0, 0.0, pt) );
      }
    }
  }

  // construct the result: vector of heights and distances from the start point to all other points.
  vector<HeightData> res;
  for (unsigned int i=0; i<heights.size(); i++) {
    //fprintf(stderr, "%d en: %d dist: %d prev: %d\n", i+1, heights[i], distance[i], previous[i]+1);
    //res[i] = make_pair(heights[i], distance[i]);
    res.push_back(HeightData(heights[i], distance[i], bdist[i], pknot_types[i]));
  }

  return res;
}

void DSU::GetPath(int start, int stop, int maxkeep, char optimality, int max_saddle, double time_fold, int bulk)
{
  char filename[100];
  if (bulk) {
    sprintf(filename, "pathB%d.%spath%c%s", bulk, optimality=='T'?std::to_string((int)time_fold).c_str():"", optimality, max_saddle<1000?"_":"");
  } else {
    sprintf(filename, "path%d_%d.%spath%c%s", start+1, stop+1, optimality=='T'?std::to_string((time_fold<=1e9?(int)time_fold:time_fold)).c_str():"", optimality, max_saddle<1000?"_":"");
  }
  GetPath(start, stop, edgesV_l, filename, maxkeep, optimality, max_saddle, time_fold);
}

void DSU::GetPath(int start, int stop,  vector< set<edgeLL> > &edgesV_l, char *filename, int maxkeep, char optim, int max_saddle, double time_fold)
{
  /*bool swapped = false;
  if (LM[start].energy > LM[stop].energy && !rate_path) {
    swapped = true;
    swap(start, stop);
  }*/

  optimality = optim;

  switch (optimality) {
    case 'B':
    case 'R':
    case 'L':
    case 'M':
    case 'T':
    case 'P':  break;
    default: fprintf(stderr, "ERROR: weird optimality criterion (%c).\n", optimality); return;
  }

  // calculate outgoing rates
  vector<double> outgoing_rates(LM.size());
  if (optimality == 'P' || optimality == 'T') {
    for (unsigned int i=0; i<outgoing_rates.size(); i++) {
      for (set<edgeLL>::iterator it=edgesV_l[i].begin(); it!=edgesV_l[i].end(); it++) {
        outgoing_rates[i] += 1.0*exp(-(it->en-LM[i].energy)/100.0/_kT);
      }
    }
  }

  // define + init
  vector<int> heights(LM.size(), INT_MAX);
  vector<int> distance(LM.size(), INT_MAX);
  vector<int> previous(LM.size(), -1);
  vector<int> saddle_num(LM.size(), -1);
  vector<bool> done(LM.size(), false);
  vector<double> rates(LM.size(), 0.0);

  priority_queue<pq_height> que;

  // starting point -- all points from start
  done[start] = true;
  distance[start] = 0;
  rates[start] = 0.0;
  heights[start] = LM[start].energy;
  for (set<edgeLL>::iterator it=edgesV_l[start].begin(); it!=edgesV_l[start].end(); it++) {
    //  add it if it does not contradict the max_saddle
    if (it->en > max_saddle) continue;

    // calculate other info
    double rate = 1.0*exp(-(it->en-heights[start])/100.0/_kT);
    int to = it->goesTo(start);
    int sadd = it->saddle;
    int hd = HammingDist(LM[to].structure, saddles[sadd].structure) + HammingDist(saddles[sadd].structure, LM[start].structure);

    double crit = 0.0;
    double time = 0.0;
    if (optimality == 'P') {
      crit = log(rate/outgoing_rates[start]);
    }
    if (optimality == 'T') {
      crit = log(rate);
      time = min(outgoing_rates[to], outgoing_rates[start]);
      crit += -(time)*time_fold;
    }
    que.push(pq_height(it->en, 1, it->goesTo(start), start, it->saddle, rate, hd, crit, time));
  }

  double end_rate = 0.0;
  double end_crit = 0.0;

  // main loop -- dijkstra-like (take one with lowest energy, proceed it)
  while (!que.empty()) {
    // get next one to do
    pq_height pq = que.top(); que.pop();
    if (done[pq.number]) continue;

    //fprintf(stderr, "from %6d to %6d, en %6.2f dist %6d\n", pq.from+1, pq.number+1, pq.height/100.0, pq.distance);

    // write him:
    done[pq.number] = true;
    distance[pq.number] = pq.distance;
    previous[pq.number] = pq.from;
    heights[pq.number] = pq.height;
    saddle_num[pq.number] = pq.saddle_num;
    rates[pq.number] = pq.rate;

    // end when we have reached destination
     if (pq.number==stop) {
      end_rate = pq.rate;
      end_crit = pq.criterion;
      break;
    }

    // push next ones:
    for (set<edgeLL>::iterator it=edgesV_l[pq.number].begin(); it!=edgesV_l[pq.number].end(); it++) {
      //  add it if it does not contradict the max_saddle
      if (it->en > max_saddle) continue;

      int to = it->goesTo(pq.number);
      // if we already had him
      if (done[to] && optimality == 'B') {
        if (max(it->en, pq.height) < heights[to]) fprintf(stderr, "WRONG from: %6d to: %6d enLM %6.2f enHeight %6d\n", pq.number+1, to+1, pq.height/100.0, heights[to]);
        continue;
      }
      // push next ones
      //fprintf(stderr, "adding(%d): from %d to %d, en %d dist %d\n", (int)!done[to], pq.number+1, to+1, max(it->en, pq.height), pq.distance+1);
      if (!done[to]) {
        double rate = 0.0;
        double r23 = 1.0*exp(-(it->en-LM[pq.number].energy)/100.0/_kT);
        double r21 = 1.0*exp(-(saddles[pq.saddle_num].energy-LM[pq.number].energy)/100.0/_kT);
        //fprintf(stderr, "%10.8g %10.8g %10.8g => %10.8g\n", pq.rate, r23, r21, pq.rate*r23/(r23+r21));
        rate = pq.rate*r23/(r23+r21); /// combine

        int sadd = it->saddle;
        int hd = HammingDist(LM[to].structure, saddles[sadd].structure) + HammingDist(saddles[sadd].structure, LM[pq.number].structure);

        double crit = 0.0;
        double time = 0.0;
        if (optimality == 'P') {
          crit = pq.criterion + log(r23/outgoing_rates[pq.number]);
        }
        if (optimality == 'T') {
          crit = pq.criterion + log(r23);
          time = outgoing_rates[to];
          if (pq.time_min > time) {
            crit += -(pq.time_min - time)*time_fold;
          } else time = pq.time_min;
        }

        que.push( pq_height(max(it->en, pq.height), pq.distance+1, to, pq.number, it->saddle, rate, pq.bdist+hd, crit, time) );
      }
    }
  }

  // retreive path:
  vector<int> lms;
  vector<int> sdd;
  vector<double> rats;
  int number = stop;
  while (number!=start) {
    lms.push_back(number);
    sdd.push_back(saddle_num[number]);
    rats.push_back(rates[number]);
    //fprintf(stderr, "%14.8g %4d\n", rates[number], number);
    number = previous[number];
  }
  lms.push_back(number);
  rats.push_back(0.0);

  // swap back:
  /*if (swapped) {
    swap(start, stop);
    reverse(lms.begin(), lms.end());
    reverse(sdd.begin(), sdd.end());
    // recompute rates:
    rats.clear();

    rats.push_back(0.0);

    double rate = 1.0*exp(-(saddles[sdd[sdd.size()-1]].energy-LM[lms[lms.size()-1]].energy)/100.0/_kT);
    rats.push_back(rate);
    for (int i=sdd.size()-2; i>=0; i--) {
      double r23 = 1.0*exp(-(saddles[sdd[i+1]].energy-LM[lms[i]].energy)/100.0/_kT);
      double r21 = 1.0*exp(-(saddles[sdd[i]].energy-LM[lms[i]].energy)/100.0/_kT);
      rate = rate*r23/(r23+r21);
      rats.push_back(rate);
      end_rate = rate;
    }
    reverse(rats.begin(),rats.end());
  }*/

  int cumulative_height = 0;

  // write it down:
  FILE *fil;
  fil = fopen(filename, "w");
  int count = 0;
  if (fil) {
    fprintf(fil,"        %s\n", seq);
    for (int i=(int)lms.size()-1; i>0; i--) {
      fprintf(fil, "%6d  %s %7.2f %3d %3d %14.8g\n", lms[i]+1, LM[lms[i]].str_ch, LM[lms[i]].energy/100.0, HammingDist(LM[lms[i]].structure, LM[lms[0]].structure), HammingDist(LM[lms[i]].structure, saddles[sdd[i-1]].structure), rats[i]);
      cumulative_height += saddles[sdd[i-1]].energy - LM[lms[i]].energy;
      if (maxkeep) {
        if (pknots) {
          path_pk *tmp = get_path_pk(seq, LM[lms[i]].str_ch, saddles[sdd[i-1]].str_ch, maxkeep);
          path_pk *path = tmp+1;
          while (path && path->s && (path+1) && (path+1)->s) {
            fprintf(fil, " inter  %s %7.2f %3d\n", path->s, path->en/100.0, HammingDist(path->s, LM[lms[0]].structure));
            path++;
            count++;
          }
          free_path_pk(tmp);

        } else {
          path_t *tmp = get_path(seq, LM[lms[i]].str_ch, saddles[sdd[i-1]].str_ch, maxkeep);
          path_t *path = tmp+1;
          while (path && path->s && (path+1) && (path+1)->s) {
            fprintf(fil, " inter  %s %7.2f %3d\n", path->s, path->en, HammingDist(path->s, LM[lms[0]].structure));
            path++;
            count++;
          }
          free_path(tmp);
        }
      }
      fprintf(fil, "%6dS %s %7.2f %3d %3d\n", sdd[i-1]+1, saddles[sdd[i-1]].str_ch, saddles[sdd[i-1]].energy/100.0, HammingDist(saddles[sdd[i-1]].structure, LM[lms[0]].structure), HammingDist(LM[lms[i-1]].structure, saddles[sdd[i-1]].structure));
      if (maxkeep) {
        if (pknots) {
          path_pk *tmp = get_path_pk(seq, saddles[sdd[i-1]].str_ch, LM[lms[i-1]].str_ch, maxkeep);
          path_pk *path = tmp+1;
          while (path && path->s && (path+1) && (path+1)->s) {
            fprintf(fil, " inter  %s %7.2f %3d\n", path->s, path->en/100.0, HammingDist(path->s, LM[lms[0]].structure));
            path++;
            count++;
          }
          free_path_pk(tmp);
        } else {
          path_t *tmp = get_path(seq, saddles[sdd[i-1]].str_ch, LM[lms[i-1]].str_ch, maxkeep);
          path_t *path = tmp+1;
          while (path && path->s && (path+1) && (path+1)->s) {
            fprintf(fil, " inter  %s %7.2f %3d\n", path->s, path->en, HammingDist(path->s, LM[lms[0]].structure));
            path++;
            count++;
          }
          free_path(tmp);
        }
      }
    }
    fprintf(fil, "%6d  %s %7.2f   0   0 %14.8g\n", lms[0]+1, LM[lms[0]].str_ch, LM[lms[0]].energy/100.0, rats[0]);
    if (maxkeep)  fprintf(fil, "Path from %6d to %6d: %4d local minima, %d structures, %6.2f kcal/mol highest point.\n", start+1, stop+1, (int)lms.size(), count + (int)lms.size() + (int)sdd.size(), heights[stop]/100.0);
    else          fprintf(fil, "Path from %6d to %6d: %4d local minima, %6.2f kcal/mol highest point.\n", start+1, stop+1, (int)lms.size(), heights[stop]/100.0);
    fprintf(fil, "Path has rate of %14.8g, cumulative energy barrier %6.2f, criterion (log-prob or log (time dependent likelihood)): %14.8g", end_rate, cumulative_height/100.0, end_crit);
    fclose(fil);
  } else {
    fprintf(stderr, "Unable to open file %s!\n", filename);
  }
}

void DSU::EHeights(FILE *heights, bool full, bool only_norm)
{
  if (!full) {
    for (unsigned int i=0; i< saddles.size(); i++) {
      //fprintf(heights, "%s %.2f %s %.2f %s %.2f\n", LM[saddles[i].lm1].str_ch, LM[saddles[i].lm1].energy/100.0, LM[saddles[i].lm2].str_ch, LM[saddles[i].lm2].energy/100.0, saddles[i].str_ch, saddles[i].energy/100.0);
      fprintf(heights, "%4d %4d %.2f\n", saddles[i].lm1+1, saddles[i].lm2+1, saddles[i].energy/100.0);
    }
  } else {
    vector<vector<HeightData> > res(LM.size());
    for (unsigned int i=0; i<LM.size(); i++) {
      if (i % 1000 == 0) fprintf(stderr, "searching... %6d/%7d\n", i, (int)res.size());
      if (only_norm && LM[i].type != NORMAL) continue;

      // search:
      res[i] = HeightSearch(i, edgesV_l);
    }
    for (unsigned int i=0; i<res.size(); i++) {
      for (unsigned int j=i+1; j<res[i].size(); j++) {
        if (only_norm && LM[i].type == NORMAL && LM[j].type == NORMAL) {
          // print lines: energy heights: from_node, to_node, energy height, distance
          fprintf(heights, "%4d %4d %.2f %3d %3d\n", i+1, j+1, res[i][j].height/100.0, res[i][j].gdist, res[i][j].bdist);
        }
      }
    }
  }
}


void DSU::ERank(FILE *rank, bool barr, bool out_conns)
{
  // saddle point energies + distances
  vector<vector<HeightData> > res(number_lm);
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
      if (res[i][j].height<1e8) {
        energy_pair ep;
        ep.barrier = res[i][j].height;
        if (barr) {
          ep.barrier -= (LM[i].energy + LM[j].energy)/2;
        }
        ep.i=i;
        ep.j=j;
        saddles.push(ep);
      }
    }
  }

  // variables
  int count = 1;
  int total_en = 0;

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

    if (saddles.size()<10) {
      fprintf(stderr, "%4d %8.2f %4d %4d\n", count++, ep.barrier/100.0, ep.i+1, ep.j+1);
    }
  }
  fprintf(stderr, "average energy = %6.2f/%4d = %11.8f\n", total_en/100.0, (count-1), total_en/100.0/(double)(count-1));
}

void DSU::Histo(FILE *histo)
{
  vector<int> histogram(LM.size(), 0);
  vector<int> histogramf(LM.size(), 0);

  int filter = -5920;

  // build
  for (unsigned int i=0; i<saddles.size(); i++) {
    if (LM[saddles[i].lm1].energy<=filter) histogramf[saddles[i].lm1]++;
    if (LM[saddles[i].lm2].energy<=filter) histogramf[saddles[i].lm2]++;
    histogram[saddles[i].lm1]++;
    histogram[saddles[i].lm2]++;
  }

  // and count them
  map<int, int> hst_map;
  map<int, int> hst_mapf;
  for (unsigned int i=0; i<histogram.size(); i++) {
    hst_map[histogram[i]]++;
  }
  for (unsigned int i=0; i<histogramf.size(); i++) {
    hst_mapf[histogramf[i]]++;
  }

  // print:
  int count = 0;
  int countf = 0;
  for (map<int, int>::iterator it=hst_map.begin(); it!=hst_map.end(); it++) {
    count += it->second;
    countf += hst_mapf[it->first];
    fprintf(histo, "%6d %6d %6d %6d %6d %6d %6d %6d\n", it->first, it->second, hst_mapf[it->first], count, countf, (int)LM.size()-count, (int)LM.size()-countf, count-countf);
  }

}

