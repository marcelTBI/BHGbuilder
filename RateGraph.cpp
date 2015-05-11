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
  #include "mxccm.h"
}

#include "RateGraph.h"
#include "hash_util.h"

#include <algorithm>

using namespace std;

double rate(int en_from, int en_to, double _kT) {
  double res = 1.0*exp(-(en_to-en_from)/100.0/_kT);
  if (res > 1.0) return 1.0;
  return res;
}

void print_path(path_t *path) {
  fprintf(stderr, "path:\n");
  while (path && path->s) {
    fprintf(stderr, " inter  %s %7.2f\n", path->s, path->en);
    path++;
  }
}

path_t *get_whole_path(DSU &dsu, int lm1, int lm2, int saddle, int maxkeep, int &count1, int &count2)
{
  path_t *tmp1 = get_path(dsu.seq, dsu.LM[lm1].str_ch, dsu.saddles[saddle].str_ch, maxkeep);
  path_t *tmp2 = get_path(dsu.seq, dsu.saddles[saddle].str_ch, dsu.LM[lm2].str_ch, maxkeep);

  //print_path(tmp1);
  //print_path(tmp2);

  path_t *path = tmp2;
  count2 = 0;
  count1 = 0;
  while (path && path->s) {
    count2++;
    path++;
  }

  path = tmp1;
  while (path && path->s) {
    count1++;
    path++;
  }

  //fprintf(stderr, "%4d %4d %d %d\n", lm1+1, lm2+1, count1, count2);


  // realloc
  tmp1 = (path_t*) realloc(tmp1, sizeof(path_t)*(count1+count2));

  // and relocate:
  path_t *path2 = tmp2+1;
  int i= count1;
  while (path2 && path2->s) {
    tmp1[i].s = path2->s;
    tmp1[i].en = path2->en;
    path2++;
    i++;
  }
  tmp1[i].s = NULL;
  free(tmp2->s);
  free(tmp2);

  //print_path(tmp1);

  count1--;
  count2--;

  return tmp1;
}

std::pair<double, double> rate_whole(path_t *path, double _kT)
{
  // create array of rates:
  vector<double> rats_up;
  vector<double> rats_dw;

  int last_en = en_fltoi(path->en);
  //path_t *tmp = path;
  path++;
  while (path && path->s) {
    int en = en_fltoi(path->en);
    rats_up.push_back(rate(last_en, en, _kT));
    rats_dw.push_back(rate(en, last_en, _kT));
    last_en = en;
    path++;
  }

  // now get a global rate out of small ones:
  double rate_up = rats_up[0];
  double rate_dw = rats_dw[0];
  for (unsigned int i=1; i<rats_up.size(); i++) {
    rate_up = rate_up*rats_up[i]/(rate_dw+rats_up[i]);
    rate_dw = rate_dw*rats_dw[i]/(rate_dw+rats_up[i]);
  }

  // debug:
  //double drat1 = rate(en_fltoi(tmp[0].en), en_fltoi(tmp[1].en), _kT);
  //double drat2 = rate(en_fltoi(tmp[2].en), en_fltoi(tmp[1].en), _kT);

  return make_pair(rate_up, rate_dw);
}

RateGraph::RateGraph(DSU &dsu, double temp, int maxkeep, double minimal_rate, char order)
{
  double _kT = 0.00198717*(273.15 + temp);

  rates.resize(dsu.Size());
  if (maxkeep) lengths.resize(dsu.Size());

  inrates.resize(dsu.Size(), 0.0);
  outrates.resize(dsu.Size(), 0.0);

  // construct rates data structr
  for (unsigned int i=0; i<dsu.edgesV_l.size(); i++) {
    for (auto j=dsu.edgesV_l[i].begin(); j!=dsu.edgesV_l[i].end(); j++) {

      double rate1, rate2;

      // extract info
      int en = j->en;
      int lm1 = j->i;
      int lm2 = j->j;
      int saddle = j->saddle;

      if (!maxkeep) {

        rate1 = rate(dsu.LM[lm1].energy, en, _kT);
        rate2 = rate(dsu.LM[lm2].energy, en, _kT);

        // construct Rate-s
        //Rate rt1 (rate1, j->saddle);
        //Rate rt2 (rate2, j->saddle);
      } else {

        // first get whole path
        int count1, count2;
        auto whole_path = get_whole_path(dsu, lm1, lm2, saddle, maxkeep, count1, count2);

        // then convert it to both rates:
        std::pair<double, double> rats = rate_whole(whole_path, _kT);
        free_path(whole_path);
        rate1 = rats.first;
        rate2 = rats.second;

        if (rate1 > minimal_rate && rate2 > minimal_rate) {
          lengths[lm1][lm2] = count1;
          lengths[lm2][lm1] = count2;
        }
      }

      // add 'em
      if (rate1 > minimal_rate && rate2 > minimal_rate) {
        rates[lm1][lm2] = rate1;
        rates[lm2][lm1] = rate2;

        inrates[lm1] += rate2;
        inrates[lm2] += rate1;
        outrates[lm1] += rate1;
        outrates[lm2] += rate2;

        edge_count++;
      }

    }
  }

  // copy local minima:
  lms.resize(dsu.LM.size());
  for (unsigned int i=0; i<dsu.LM.size(); i++) {
    lms[i] = dsu.LM[i];
    lms[i].structure = allocopy(dsu.LM[i].structure);
    lms[i].str_ch = NULL;
    lms[i].recompute_str();
  }

  // and structure:
  seq = dsu.seq;

  // maps:
  pos_to_lm.resize(lms.size());
  for (unsigned int i=0; i<lms.size(); i++) {
    pos_to_lm[i] = i;
    lm_to_pos[i] = i;
  }

  // null it.
  rate_matrix = NULL;
}

int RateGraph::ReadFilter(char *filename)
{
  FILE *file_filt;
  file_filt = fopen(filename, "r");
  if (file_filt) {
    char *line;
    line = my_getline(file_filt);
    while (line) {
      char *p = strtok(line, " ");
      while (p) {
        if (isStruct(p)) {
          //fprintf(stderr, "%s \n",p);
          filter.insert(p);
          break;
        }
        p = strtok(NULL, " ");
      }
      free(line);
      line = my_getline(file_filt);
    }
    free(line);

    fclose(file_filt);
  } else {
    fprintf(stderr, "WARNING: cannot open filter file ""%s""!\n", filename);
  }
  int size = filter.size();
  return size;
}

RateGraph::~RateGraph()
{
  if (rate_matrix) free(rate_matrix);
  for (unsigned int i=0; i<lms.size(); i++) {
    free(lms[i].structure);
    if (lms[i].str_ch) free(lms[i].str_ch);
  }
}

int RateGraph::ConstructQueue(char order, int number_remove, bool leave_trans)
{
  // tmp pq for ordering
  priority_queue<pq_rate, vector<pq_rate>, pq_comp> pq_tmp(order);
  int added = 0;
  while (!to_rem.empty()) to_rem.pop();

  switch (order) {
  case 'R':
  case 'C':
    // make a priority queue with just top "number_remove"
    for (unsigned int i=0; i<rates.size(); i++) {
      if (filter.size() == 0 || filter.count(lms[i].str_ch)==0) {
        pq_tmp.push(pq_rate(rates[i].size(), i, lms[i].energy, inrates[i], outrates[i]));
      }
    }
    // add only first number_remove to the pq_tmp:
    while ((int)to_rem.size() != number_remove && !pq_tmp.empty()) {
      pq_rate pq = pq_tmp.top();
      //fprintf(stderr, "%5d %5d %10g\n", pq.lm, pq.conns, pq.outrate/pq.inrate);
      to_rem.push(pq);
      pq_tmp.pop();
    }

    // now include also connections <= 2 if not flagged out:
    if (!leave_trans) {
      while (!pq_tmp.empty()) {
        pq_rate pq = pq_tmp.top();
        if (pq.conns<=2) to_rem.push(pq);
        pq_tmp.pop();
      }
    }

    break;
  case 'E':
    // the same, but insert there only high energy ones. assume sorted order.
    for (unsigned int i=rates.size()-1; i>=rates.size(); i--) {
      if (filter.size() == 0 || filter.count(lms[i].str_ch)==0) {
        if (added < number_remove) {
          to_rem.push(pq_rate(rates[i].size(), i, lms[i].energy, inrates[i], outrates[i]));
          added++;
        }
        // add conns <= 2
        if (rates[i].size()<=2 && !leave_trans) {
          to_rem.push(pq_rate(rates[i].size(), i, lms[i].energy, inrates[i], outrates[i]));
        }
      }
    }
    break;
  default: return 1;
  }
  //int res = to_rem.size();
  return to_rem.size();
}

int RateGraph::RemoveOne(int rem_lm, double minimal_rate)
{
  int count = 0;

  // calculate sum
  double sum = 0.0;
  for (auto it=rates[rem_lm].begin(); it!=rates[rem_lm].end(); it++) {
    sum += it->second;
  }

  // go through doubles
  for (auto it=rates[rem_lm].begin(); it!=rates[rem_lm].end(); it++) {

    double rR_1 = it->second;
    int lm1 = it->first;
    auto rat_it = rates[lm1].find(rem_lm);
    double r1_R = rat_it->second;
    // this can happen due to minimal rate:
    if (rat_it == rates[lm1].end()) {

    } else {
      rates[lm1].erase(rat_it);
      count--;
    }

    auto it2 = it; it2++;
    for (; it2!=rates[rem_lm].end(); it2++) {
      // we have rates rem_lm -> lm1; rem_lm> lm2
      // get others
      double rR_2 = it2->second;
      int lm2 = it2->first;
      auto rat_it = rates[lm2].find(rem_lm);
      double r2_R = rat_it->second;
      //rates[lm2].erase(rat_it);

      //join 'em
      double r1_2 = r1_R*rR_2/sum;
      double r2_1 = r2_R*rR_1/sum;

      // add them to rates:
      auto f12 = rates[lm1].find(lm2);
      auto f21 = rates[lm2].find(lm1);

      bool not_enough = false;
      if (f12!=rates[lm1].end()) {
        f12->second += r1_2;
      } else {
        if (r1_2 > minimal_rate) {
          count++;
          rates[lm1][lm2] = r1_2;
        } else {
          not_enough = true;
        }
      }
      if (f21!=rates[lm2].end()) {
        f21->second += r2_1;
      } else {
        if (r2_1 > minimal_rate && !not_enough) {
          count++;
          rates[lm2][lm1] = r2_1;
        } else{
          rates[lm2].erase(lm1);
          count--;
        }
      }
    }
  }

  // erase edges:
  count -= (int)rates[rem_lm].size();
  rates[rem_lm].clear();

  // return how many edges inserted
  return count;
}

// remove nodes one by one.
int RateGraph::RemoveX(int x, int stop_fraction, bool reeval, int maximal, double minimal_rate)
{
  clock_t time = clock();
  while ((int)removed.size()!=x) {

    // get lowest
    pq_rate pqr = to_rem.top(); to_rem.pop();

    // check it:
    if (reeval && pqr.conns != (int)rates[pqr.lm].size()) {
      pqr.conns = (int)rates[pqr.lm].size();
      to_rem.push(pqr);
      continue;
    }

    // stopping criterion:
    int conns = (int)rates[pqr.lm].size();
    if (conns * stop_fraction > (int)(rates.size()-removed.size()) && (int)(rates.size()-removed.size()) <= maximal) {
      to_rem.push(pqr);
      break;
    }

    // debug:
    if (removed.size()%10000 == 0 || (conns>100 && removed.size()%1000==0) || (conns>200 && removed.size()%100==0) || (conns>500)) {
      fprintf(stderr, "removing node %8d (%8d remaining) (conns=%4d, energy=%6.2f) (%6.2f secs.)\n", (int)removed.size(), (int)(rates.size()-removed.size()), conns, pqr.energy/100.0, (clock()-time)/(double)CLOCKS_PER_SEC);
      time = clock();
    }

    /*// debug
    char name[100];
    sprintf(name, "test%d.dot", (int)rates.size()-removed.size());
    PrintDot(name, true);
*/
    // remove it:
    int new_edges = RemoveOne(pqr.lm, minimal_rate);
    removed.insert(pqr.lm);
    edge_count += new_edges;
  }

  /*// debug
  char name[100];
  sprintf(name, "test%d.dot", (int)rates.size()-count);
  PrintDot(name, true);
*/
  // reduce size of rates vector:
  // update maps:
  if (removed.size() > 0) {
    lm_to_pos.clear();
    pos_to_lm.clear();
    auto it_rem = removed.begin();
    for(unsigned int i=0, j=0; i<rates.size(); i++) {
      if (it_rem==removed.end() || *it_rem != (int)i) {
        lm_to_pos[i] = j;
        j++;
      } else {
        if (it_rem!=removed.end()) it_rem++;
      }
    }
    pos_to_lm.resize(lm_to_pos.size());
    for (auto it=lm_to_pos.begin(); it!=lm_to_pos.end(); it++) {
      pos_to_lm[it->second] = it->first;
    }

    // rewrite rates vector:
    for (unsigned int i=0; i<rates.size(); i++) {
      map<int, double> tmp;
      for (auto it=rates[i].begin(); it!=rates[i].end(); it++) {
        tmp.insert(make_pair(lm_to_pos[it->first], it->second));
      }
      rates[i] = tmp;
    }
    // and remove empty lines
    int new_i = 0;
    for(unsigned int i=0; i<rates.size()-removed.size(); i++, new_i++) {
      while(removed.count(new_i)) new_i++;
      rates[i] = rates[new_i];
    }
    rates.resize(rates.size()-removed.size());
  }

  /*// print the conversion
  for (unsigned int i=0; i<pos_to_lm.size(); i++) {
    fprintf(stderr, "%8d -> %8d\n", pos_to_lm[i]+1, i+1);
  }*/

  // return how many removed
  return removed.size();
}

void SwapMinima(int dim_rates, double *rate_matrix, int a, int b)
{
  for (int i=0; i<dim_rates; i++) {
    //fprintf(stderr, "swapping %5d %5d %d (%d)\n", a, b, dim_rates, i );
    // change lines:
    swap(rate_matrix[a*dim_rates+i], rate_matrix[b*dim_rates+i]);
  }
  for (int i=0; i<dim_rates; i++) {
    // change columns
    swap(rate_matrix[i*dim_rates+a], rate_matrix[i*dim_rates+b]);
  }
}

void MxRShorten(double **shorten, int fulldim, int gdim)
{
  //does: shortened = GG - GB*BB^(-1)*BG, where matrix tmp_rates is split as:
  //tmp_rates = (GG | GB)
  //            (BG | BB)
  // GG has dimension gdim*gdim; tmp_rates fulldim*fulldim
  // create matrices:

  int bdim = fulldim - gdim;
  int i,j;

  //double *gg = (double *)calloc(gdim*gdim,sizeof(double));
  double *bg = (double *)calloc(bdim*gdim,sizeof(double));
  double *bb = (double *)calloc(bdim*bdim,sizeof(double));
  double *gb = (double *)calloc(gdim*bdim,sizeof(double));

  // first we need to fix the diagonal entries tmp_rates[i][i] = sum_j tmp_rates[i][j]
  double *tmp_rates = *shorten;
  for (i = 0; i < fulldim; i++) tmp_rates[fulldim*i+i] = 0.0;
  for (i = 0; i < fulldim; i++) {
    double tmp = 0.00;
    // calculate row sum
    for(j = 0; j < fulldim; j++)  tmp += tmp_rates[fulldim*i+j];
    tmp_rates[fulldim*i+i] = -tmp;
  }

  // fill the matrices: (row = i; column = j)
  for (i=0; i<bdim; i++) {
    for (j=0; j<gdim; j++) {
      bg[gdim*i+j] = tmp_rates[fulldim*(i+gdim)+j];
    }
  }

  for (i=0; i<gdim; i++) {
    for (j=0; j<bdim; j++) {
      gb[bdim*i+j] = tmp_rates[fulldim*i+j+gdim];
    }
  }

  for (i=0; i<bdim; i++) {
    for (j=0; j<bdim; j++) {
      bb[bdim*i+j] = tmp_rates[fulldim*(i+gdim)+j+gdim];
    }
  }

  /*MxFPrintD(tmp_rates, "Q", my_dim, my_dim, stderr);
  MxFPrintD(gg, "GG", dim, dim, stderr);
  MxFPrintD(bg, "BG", bdim, dim, stderr);
  MxFPrintD(gb, "GB", dim, bdim, stderr);
  MxFPrintD(bb, "BB", bdim, bdim, stderr);
*/
  // result2 = gb*bb^(-1)*bg
  minv(bb, bdim);
  //MxFPrintD(bb, "BBinv", bdim, bdim, stderr);
  double *result = (double *)calloc(gdim*bdim,sizeof(double));
  mmul_singular(result, gb, bb, gdim, bdim, bdim, 0);
  //MxFPrintD(result, "gb*bb-1", dim, bdim, stderr);
  double *result2 = (double *)calloc(gdim*gdim,sizeof(double));
  mmul_singular(result2, result, bg, gdim, bdim, gdim, 1);

  //if (opt.want_verbose) MxFPrintD(result2, "gb*bb-1*bg", gdim, gdim, stderr);

  // result2 = gg - result2
  for (i=0; i<gdim; i++) {
    for (j=0; j<gdim; j++) {
      result2[gdim*i+j] = tmp_rates[fulldim*i+j] - result2[gdim*i+j];
    }
  }

  //MxFPrintD(result2, "matrix after shortening", dim ,dim, stderr);
  free(result);
  free(*shorten);
  free(gb);
  free(bg);
  free(bb);
  *shorten = result2;
}

int RateGraph::RemoveSmall(double *shorten, int dim, double minimal_rate)
{
  int count = 0;

  if (minimal_rate == 0.0) return 0;

  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      if (shorten[i*dim + j]<minimal_rate) {
        count++;
        shorten[i*dim + j] = 0.0;
        shorten[j*dim + i] = 0.0;
      }
    }
  }

  return count;
}

struct sort_struc {
  int number;
  bool to_remove;

  sort_struc(int number, bool to_remove) {
    this->number = number;
    this->to_remove = to_remove;
  }

  const bool operator<(const sort_struc &second) const {
    if (to_remove == second.to_remove) {
      return number < second.number;
    } else {
      return to_remove < second.to_remove;
    }
  }
};

int RateGraph::RemoveShur(int x, int step, double minimal_rate)
{
  //fprintf(stderr, "calling Shur %5d %5d\n", x, step );
  // get set to remove:
  /*set<int> to_remove;
  while ((int)to_remove.size() != x) {
    pq_rate tmp = to_rem.top(); to_rem.pop();
    to_remove.insert(lm_to_pos[tmp.lm]);
  }*/

  //build rate matrix:
  dim_rates = rates.size();
  if (rate_matrix) free(rate_matrix);
  rate_matrix = (double *) malloc(sizeof(double)*dim_rates*dim_rates);
  for (unsigned int i=0; i<dim_rates*dim_rates; i++) rate_matrix[i] = 0.0;

  for (unsigned int i=0; i<dim_rates; i++) {
    double sum = 0.0;
    for (auto a=rates[i].begin(); a!=rates[i].end(); a++) {
      rate_matrix[i*dim_rates + a->first] = a->second;
      sum += a->second;
    }
    rate_matrix[i*dim_rates + i] = -sum;
  }
  rates.clear();


  PrintRates("before_reorder.rat");

  // now resort matrix to have the minima to remove on bottom.
  vector<sort_struc> tmp_sort;
  for (unsigned int i=0; i<dim_rates; i++) {
    tmp_sort.push_back(sort_struc(i, false));
  }
  int to_remove = 0;
  while (to_remove != x) {
    pq_rate tmp = to_rem.top(); to_rem.pop();
    tmp_sort[lm_to_pos[tmp.lm]].to_remove = true;
    to_remove++;
  }
  sort(tmp_sort.begin(), tmp_sort.end());
  // make inverse:
  vector<int> tmp_inv(tmp_sort.size());
  for (unsigned int i=0; i<tmp_sort.size(); i++) {
    tmp_inv[tmp_sort[i].number] = i;
  }

  // sort it by swapping (in place ordering):
  for (unsigned int i=0; i<dim_rates; i++) {
    while ((int)i!=tmp_inv[i]) {
      SwapMinima(dim_rates, rate_matrix, i, tmp_inv[i]);

      swap(pos_to_lm[i], pos_to_lm[tmp_inv[i]]);
      swap(lm_to_pos[pos_to_lm[tmp_inv[i]]], lm_to_pos[pos_to_lm[i]]);

      swap(tmp_inv[i], tmp_inv[tmp_inv[i]]);
    }
  }

/*
  // sort it temporary:
  {
    double *tmp_rates = (double*) malloc(dim_rates*dim_rates*sizeof(double));

    for (int i=0; i<dim_rates; i++) {
      for (int j=0; j<dim_rates; j++) {
        tmp_rates[i*dim_rates+j] = rate_matrix[tmp_sort[i].number*dim_rates+tmp_sort[j].number];
      }
    }
    swap(rate matrix, tmp_rates);
    free(tmp_rates);
  }

  // adjust mapping:
  for (int i=0; i<dim_rates; i++) {
    swap(pos_to_lm[i], pos_to_lm[tmp_sort[i].number]);
    swap(lm_to_pos[pos_to_lm[*it]], lm_to_pos[pos_to_lm[i]]);
  }*/

  //PrintRates("after_reorder.rat");

  //fprintf(stderr, "calling Shur %5d %5d\n", x, step );

  clock_t time = clock();

  // finally run removal on that matrix
  int res_dim = dim_rates-x;
  while ((int)dim_rates > res_dim) {
    // reduce of:
    int reduce = min(step, (int)dim_rates-res_dim);
    // now Shur it!
    MxRShorten(&rate_matrix, dim_rates, dim_rates-reduce);
    dim_rates -= reduce;

    int removed = RemoveSmall(rate_matrix, dim_rates, minimal_rate);

    // report:
    fprintf(stderr, "reduced dimension from %8d to %8d (%6.2f secs.) (%d removed)\n", dim_rates+reduce, dim_rates, (clock()-time)/(double)CLOCKS_PER_SEC, removed);
    time = clock();
  }

  // add those, that we have removed
  for (unsigned int i=res_dim; i<dim_rates+x; i++) {
    removed.insert(pos_to_lm[i]);
    lm_to_pos.erase(pos_to_lm[i]);
  }
  pos_to_lm.resize(res_dim);

  // lastly, rewrite the "rates" data structure:
  rates.clear();
  rates.resize(res_dim);
  for (int i=0; i<res_dim; i++) {
    for (int j=0; j<res_dim; j++) {
      if (rate_matrix[i*res_dim+j]>0 && i!=j) {
        rates[i][j] = rate_matrix[i*res_dim+j];
      }
    }
  }

  return x;
}

void RateGraph::PrintRates(char *filename)
{
  FILE *rate_file = fopen(filename, "w");
  if (rate_file) {
    PrintRates(rate_file);
    fclose(rate_file);
  }
}

void RateGraph::PrintRates(FILE *filname)
{
  if (!rate_matrix) {
      //build rate matrix:
    dim_rates = rates.size();
    rate_matrix = (double *) malloc(sizeof(double)*dim_rates*dim_rates);
    for (unsigned int i=0; i<dim_rates*dim_rates; i++) rate_matrix[i] = 0.0;

    for (unsigned int i=0; i<dim_rates; i++) {
      double sum = 0.0;
      for (auto a=rates[i].begin(); a!=rates[i].end(); a++) {
        rate_matrix[i*dim_rates + a->first] = a->second;
        sum += a->second;
      }
      rate_matrix[i*dim_rates + i] = -sum;
    }
  }

  // just print the rate matrix:
  for (unsigned int i=0; i<dim_rates; i++) {
    for (unsigned int j=0; j<dim_rates; j++) {
        fprintf(filname, "%11.5g ", rate_matrix[dim_rates*i+j]);
    }
    fprintf(filname, "\n");
  }

  // print lengths:
  if (lengths.size()>0) {
    FILE *lengths_file = fopen("lenghts.rat", "w");
    if (lengths_file) {
      vector<vector<int> > len_mat(dim_rates);
      for (unsigned int i=0; i<dim_rates; i++) {
        len_mat[i].resize(dim_rates);
        for (auto a=lengths[i].begin(); a!=lengths[i].end(); a++) {
          len_mat[i][a->first] = a->second;
        }
      }

      for (unsigned int i=0; i<dim_rates; i++) {
        for (unsigned int j=0; j<dim_rates; j++) {
            fprintf(lengths_file, "%11d ", len_mat[i][j]);
        }
        fprintf(lengths_file, "\n");
      }
      fclose(lengths_file);
    }
  }
}

void RateGraph::PrintDot(char *filename, bool to_eps)
{
  FILE *dot = fopen(filename, "w");
  if (dot) {
    PrintDot(dot);
    fclose(dot);
  }

  if (to_eps) {
    // start neato/dot:
    char syst[200];
    char filename2[200];
    strcpy(filename2, filename);
    filename2[strlen(filename2)-4]='\0';
    sprintf(syst, "%s -Tps < %s > %s.eps", "dot", filename, filename2);
    system(syst);
  }
}

void RateGraph::PrintDot(FILE *dot)
{
  fprintf(dot, "Digraph G {\n\tnode [width=0.1, height=0.1, shape=circle];\n");
  //nodes LM:
  for (unsigned int i=0; i<rates.size(); i++) {
    fprintf(dot, "\"%d\" [label=\"%d\"]\n", i+1, pos_to_lm[i]+1/*, lms[pos_to_lm[i]].energy*/);
  }
  fprintf(dot, "\n");

  for (unsigned int i=0; i<rates.size(); i++) {
    for (auto it=rates[i].begin(); it!=rates[i].end(); it++) {
      // edges l-l
      char length[20] = "";
      if (lengths.size()) sprintf(length, " (%d)", lengths[i][it->first]);
      fprintf(dot, "\"%d\" -> \"%d\" [label=\"%6.4g%s\"]\n", i+1, (it->first)+1, it->second, length);
    }
  }
  fprintf(dot, "\n}\n");
}

void RateGraph::PrintOutput(char *filename)
{
  FILE *outpu = fopen(filename, "w");
  if (outpu) {
    PrintOutput(outpu);
    fclose(outpu);
  }
}

void RateGraph::PrintOutput(FILE *output)
{
  fprintf(output, "      %s\n", seq.c_str());
  for (unsigned int i=0; i<pos_to_lm.size(); i++) {
    fprintf(output, "%5d %s %6.2f %5d\n", pos_to_lm[i]+1, lms[pos_to_lm[i]].str_ch, lms[pos_to_lm[i]].energy/100.0, i+1);
  }
}
