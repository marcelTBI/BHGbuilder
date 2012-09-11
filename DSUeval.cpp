#include <stdio.h>
#include <string.h>
#include <stdlib.h>

extern "C" {
  #include "pair_mat.h"
  #include "findpath.h"
  #include "move_set.h"
  #include "fold.h"
}

#include "DSUeval.h"
#include "hash_util.h"

#include <algorithm>

using namespace std;


/* reads a line no matter how long*/
char* my_getline(FILE *fp)
{
  char s[512], *line, *cp;
  line = NULL;
  do {
    if(fgets(s, 512, fp) == NULL) break;
    cp = strchr(s, '\n');
    if(cp != NULL) *cp = '\0';
    if(line == NULL) line = (char *) calloc(strlen(s) + 1, sizeof(char));
    else line = (char *) realloc(line, strlen(s) + strlen(line) + 1);
    strcat (line, s);
  } while (cp == NULL);
  return (line);
}

// pt to str
string pt_to_str(const short *pt)
{
  string str;
  str.resize(pt[0]);
  for (int i=1; i<=pt[0]; i++) {
    if (pt[i]==0) str[i-1]='.';
    else if (pt[i]<i) str[i-1]=')';
    else str[i-1]='(';
  }
  return str;
}

// structure equality
bool str_eq(const short *lhs, const short *rhs) {
  int i=1;
  while (i<=lhs[0] && lhs[i]==rhs[i]) {
    i++;
  }
  if (i>lhs[0]) return true;
  else return false;
}

int en_fltoi(float en)
{
  if (en < 0.0) return (int)(en*100 - 0.5);
  else return (int)(en*100 + 0.5);
}

inline bool isStruct(char *p)
{
  // check first two chars - should be enough
  if (strlen(p)<2) return false;
  if ((p[0]=='.' || p[0]=='(' || p[0]==')') && (p[1]=='.' || p[1]=='(' || p[1]==')')) return true;
  else return false;
}

int HammingDist(const short* struct1, const short* struct2)
{
  int match = 0;
  int str1_par = 0;
  int str2_par = 0;
  for (int i=1; i<=struct1[0]; i++) {
    if (struct1[i]!=0 && struct1[i]>i) {  //'('
      str1_par++;
      // count '(' that does match
      if (struct1[i]==struct2[i]) {
        match++;
      }
    }
    if (struct2[i]!=0 && struct2[i]>i) str2_par++;
  }

  // return all pairs minus those that matches
  return str1_par+str2_par-(2*match);
}

inline bool isSeq(char *p)
{
  if (strlen(p)<2) return false;
  // check first two chars - should be enough
  switch (p[0]){
    case 'A':
    case 'C':
    case 'G':
    case 'T':
    case 'U':
    case 'a':
    case 'c':
    case 'g':
    case 't':
    case 'u': switch (p[1]){
      case 'A':
      case 'C':
      case 'G':
      case 'T':
      case 'U':
      case 'a':
      case 'c':
      case 'g':
      case 't':
      case 'u': return true;
    }
    default : return false;
  }
}

pq_entry::pq_entry(int i, int j, int hd)
{
  this->i = min(i, j);
  this->j = max(i, j);
  this->hd = hd;
}

pq_entry::~pq_entry() {
}

DSU::DSU(FILE *input) {

  // NULL::
  seq = NULL;
  s0 = NULL;
  s1 = NULL;
  gl_maxen = INT_MIN;

  if (!input) return ;

  // read seq
  char *line;

  line = my_getline(input);
  char *seq2 = strtok(line, " ");
  if (!isSeq(seq2)) {
    free(line);
    return ;
  }
  seq = (char*) malloc((strlen(seq2)+1)*sizeof(char));
  strcpy(seq, seq2);
  free(line);

  //init
  make_pair_matrix();
  s0 = encode_sequence(seq, 0);
  s1 = encode_sequence(seq, 1);

  //read structs
  line = my_getline(input);
  while (line) {
    short *tmp = NULL;
    float energy_tmp;
    bool has_energy = false;

    char *p;
    p = strtok(line, " ");
    while (p !=NULL && (!has_energy || !tmp)) {
      // is struct?
      if (isStruct(p)) {
        if (tmp) free(tmp); // only one struct per line!
        tmp = make_pair_table(p);
        has_energy = false;
      } else {
        // is energy?
        if (sscanf(p, "%f", &energy_tmp)==1) {
          has_energy = true;
        }
      }
      p = strtok(NULL, " ");
    }
    // add info:
    RNAstruc struc;
    struc.structure = tmp;
    char *ch = (char*) malloc((tmp[0]+1)*sizeof(char));
    strcpy(ch, pt_to_str(tmp).c_str());
    struc.str_ch = ch;
    int en = (has_energy?en_fltoi(energy_tmp):energy_of_structure_pt(seq, tmp, s0, s1, 0));
    struc.energy = en;
    LM.push_back(struc);
    gl_maxen = max(gl_maxen, en);
    free(line);
    line = my_getline(input);
  }
}

DSU::~DSU() {
  for (unsigned int i=0; i<LM.size(); i++) {
    if (LM[i].structure) free(LM[i].structure);
    if (LM[i].str_ch) free(LM[i].str_ch);
  }
  if (seq) free(seq);
  if (s0) free(s0);
  if (s1) free(s1);
  LM.clear();

  for (map<pq_entry, RNAstruc, pq_setcomp>::iterator it=UBlist.begin(); it!=UBlist.end(); it++) {
    if (it->second.structure) free(it->second.structure);
    if (it->second.str_ch) free(it->second.str_ch);
  }

  for (unsigned int i=0; i<UBoutput.size(); i++) {
    if (UBoutput[i].first.structure) free(UBoutput[i].first.structure);
  }
}

int DSU::CreateList(int hd_threshold, bool debug)
{
  for (unsigned int i=0; i<LM.size(); i++) {
    for (unsigned int j=i+1; j<LM.size(); j++) {
      int hd = HammingDist(LM[i].structure, LM[j].structure);
      if (hd < hd_threshold) {
        TBDlist.push(pq_entry(i, j, hd));
        if (debug) {
          fprintf(stderr, "%s insert que\n%s (%4d,%4d) %3d\n", pt_to_str(LM[i].structure).c_str(), pt_to_str(LM[j].structure).c_str(), i, j, hd);
        }
      }
    }
  }

  return 0;
}

int DSU::FindNum(int en_par, short *str_par)
{
  //fprintf(stderr, "%d %s\n", en_par, pt_to_str(str_par).c_str());
  for (unsigned int i=0; i<LM.size(); i++) {
    if (LM[i].energy == en_par && str_eq(str_par, LM[i].structure)) return i;
  }
  return -1;
}

bool DSU::InsertUB(int i, int j, int energy_par, short *saddle_par, bool debug)
{
  pq_entry pq(i, j, 0);
  map<pq_entry, RNAstruc, pq_setcomp>::iterator it = UBlist.find(pq);
  if (it==UBlist.end()) {
    RNAstruc saddle;
    saddle.energy = energy_par;
    saddle.structure = saddle_par;
    saddle.str_ch = NULL;
    if (debug) fprintf(stderr, "UBins: (%3d, %3d) insert saddle  %6.2f %s\n", i, j, saddle.energy/100.0, pt_to_str(saddle.structure).c_str());
    UBlist.insert(make_pair(pq, saddle));
    return true;
    //fprintf(stderr, "cannot find (UB  ): (%3d, %3d)\n", num1, num2);
  } else {
    // update it
    if (energy_par < it->second.energy) {
      if (debug) fprintf(stderr, "UBupd: (%3d, %3d) from %6.2f to %6.2f %s\n", i, j, it->second.energy/100.0, energy_par/100.0, pt_to_str(it->second.structure).c_str());
      it->second.energy = energy_par;
      if (it->second.structure) free(it->second.structure);
      it->second.structure = saddle_par;
      it->second.str_ch = NULL;
      return true;
    }
  }
  free(saddle_par);
  return false;
}

int DSU::ComputeUB(int maxkeep, int num_threshold, bool debug)
{
  int dbg_count = 0;
  int cnt = 0;
  int hd_threshold = INT_MAX;
  while (!TBDlist.empty()) {
    // get pair
    pq_entry pq = TBDlist.top();
    TBDlist.pop();

    // apply threhold
    if (pq.hd > hd_threshold) break;
    if (cnt>num_threshold) {
      if (hd_threshold==INT_MAX) {
        hd_threshold = pq.hd;
        fprintf(stderr, "hd threshold set to %4d\n", hd_threshold);
      }
    } else {
      cnt++;
    }

    // get path
    if (debug) fprintf(stderr, "path between (%3d, %3d) hd=%3d:\n", pq.i, pq.j, pq.hd);
    path_t *path = get_path(seq, LM[pq.i].str_ch, LM[pq.j].str_ch, maxkeep);

    // variables for outer insertion
    double max_energy= -1e8;
    path_t *max_path = path;

    // variables for inner loops and insertions
    path_t *tmp = path;
    path_t *last = NULL;
    short *last_str = NULL;
    int last_en;
    int last_num = -1;

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
      int tmp_en = move_rand(seq, tmp_str, s0, s1, 0);

      // do the stuff if we have 2 structs and they are not equal
      if (last && !str_eq(last_str, tmp_str)) {
        // not equal LM - we can update something in UBlist
          // find LM num:
        int num1 = (last_num!=-1?last_num:FindNum(last_en, last_str));
        int num2 = FindNum(tmp_en, tmp_str);
        last_num = num2;

        // update UBlist
        if (num1==-1 || num2==-1) {
          /*if (num1==-1) {   // double outputs :/
            if (gl_maxen < last_en) fprintf(stderr, "exceeds en. (last): %s %6.2f\n", pt_to_str(last_str).c_str(), last_en/100.0);
            else                    fprintf(stderr, "cannot find (last): %s %6.2f\n", pt_to_str(last_str).c_str(), last_en/100.0);
          }*/
          if (num2==-1) {
            if (gl_maxen < tmp_en) fprintf(stderr, "exceeds en.: %s %6.2f\n", pt_to_str(tmp_str).c_str(), tmp_en/100.0);
            else {
                                   fprintf(stderr, "cannot find: %s %6.2f\n", pt_to_str(tmp_str).c_str(), tmp_en/100.0);
              // maybe add to list of minima and count with them later... TODO
            }
          }
        } else {
          // store (maybe) better saddle to UB
          int en_tmp = en_fltoi(max(last->en, tmp->en));
          short *saddle = (last->en > tmp->en ? make_pair_table(last->s) : make_pair_table(tmp->s));
          InsertUB(num1, num2, en_tmp, saddle, debug);
        }
      }

      // move one next
      if (last_str) free(last_str);
      last_en = tmp_en;
      last_str = tmp_str;
      last = tmp;
      tmp++;
    } // crawling path

    // insert saddle between outer structures
    InsertUB(pq.i, pq.j, en_fltoi(max_energy), make_pair_table(max_path->s), debug);

    // free stuff
    if (last_str) free(last_str);
    free_path(path);
  } // all doing while

  // now just resort UBlist to something sorted according energy
  UBoutput.reserve(UBlist.size());
  for (map<pq_entry, RNAstruc, pq_setcomp>::iterator it=UBlist.begin(); it!=UBlist.end(); it++) {
    UBoutput.push_back(make_pair(it->second, it->first));
  }
  sort(UBoutput.begin(), UBoutput.end());
  UBlist.clear();

  // check if everything has been found:
  if (UBoutput.size() != (LM.size()*(LM.size()-1))/2) {
    fprintf(stderr, "WARNING: All connections have not been found: %d/%d found (%d missing)\n", (int)UBoutput.size(), (int)(LM.size()*(LM.size()-1))/2, (int)((LM.size()*(LM.size()-1))/2 - UBoutput.size()));
  }

  return 0;
}

void DSU::PrintUBoutput()
{
  printf("     %s\n", seq);
  for (unsigned int i=0; i<UBoutput.size(); i++) {
    printf("%4d (%4d,%4d) saddle: %s %6.2f\n", i+1, UBoutput[i].second.i, UBoutput[i].second.j, pt_to_str(UBoutput[i].first.structure).c_str(), UBoutput[i].first.energy/100.0);
  }
}

int DSU::LinkCP(bool shifts, bool noLP, bool debug)
{
  for (unsigned int i=0; i<UBoutput.size(); i++) {
    RNAstruc stru = UBoutput[i].first;
    pq_entry pq = UBoutput[i].second;

    // flood them up!
    int res = FloodUp(LM[pq.i], LM[pq.j], stru, shifts, noLP, debug);
    if (res == 1) {
      // maybe other color ?
    }// we have found better DS -> true ds ???

    // update vertex/edge sets
    map<RNAstruc, int>::iterator it;
    int saddle_num;
    if ((it = vertex_s.find(stru)) != vertex_s.end()) {
      saddle_num = it->second;
    } else {
      saddle_num = saddles.size();
      // update vertex
      vertex_s.insert(make_pair(stru, saddle_num));
      saddles.push_back(stru);
    }
    // edges ls
    edges_ls.insert(make_pair(pq.i, saddle_num));
    edges_ls.insert(make_pair(pq.j, saddle_num));
    // edges ll
    edgeLM e;
    e.i = pq.i;
    e.j = pq.j;
    e.en = stru.energy;
    edges_l.insert(e);
  }

  return 0;
}

// variables and data strucutres for flooding
bool foundDS;
RNAstruc curr;
int currNum;
int threshold;
bool gl_debug;
unordered_map<RNAstruc, int, hash_fncts, hash_eq> flood_hash;
unordered_map<RNAstruc, int, hash_fncts, hash_eq>::iterator it;
priority_queue<RNAstruc, vector<RNAstruc>, RNAstruc_rev> flood_queue;

int funct(struct_en *moved, struct_en *current) {

  if (moved->energy > curr.energy && moved->energy < threshold) {
    RNAstruc tmp;
    tmp.energy = moved->energy;
    tmp.structure = moved->structure;
    if ((it = flood_hash.find(tmp))!=flood_hash.end()) {
      // found DS!
      if (it->second != currNum) {
        foundDS = true;
        copy_arr(current->structure, moved->structure);
        current->energy = moved->energy;

        if (gl_debug) {
          fprintf(stderr, "FOUND saddle: %s %6.2f\n", pt_to_str(current->structure).c_str(), current->energy/100.0);
        }

        return 1;
      } // else nothing
    } else {
      // insert strucure we havent see yet (and its in range)
      RNAstruc to_ins;
      to_ins.energy = moved->energy;
      to_ins.structure = allocopy(moved->structure);
      flood_hash.insert(make_pair(to_ins, currNum));
      flood_queue.push(to_ins);
    }
  }

  return 0;
}

int DSU::FloodUp(RNAstruc &i, RNAstruc &j, RNAstruc &saddle, bool shifts, bool noLP, bool debug)
{
  // threshold set
  threshold = saddle.energy;
  gl_debug = debug;

  // debug output
  if (debug) {
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

  // main loop
  curr = flood_queue.top();
  flood_queue.pop();
  currNum = flood_hash[curr];
  while (curr.energy < threshold) {

    // debug output
    if (debug) {
      fprintf(stderr, "browsing    : %s %6.2f\n", pt_to_str(curr.structure).c_str(), curr.energy/100.0);
    }

    // start browsing
    foundDS = false;
    curr.energy = browse_neighs(seq, curr.structure, s0, s1, 0, shifts, noLP, funct);

    // found DS - quitting
    if (foundDS) {
      break;
    }

    // get another struct
    if (flood_queue.empty()) break;
    curr = flood_queue.top();
    flood_queue.pop();
    currNum = flood_hash[curr];
  }

  // we did end sucesfully
  if (foundDS) {
    copy_arr(saddle.structure, curr.structure);
    saddle.energy = curr.energy;
    res = 1;
  }

  // free
  for (it=flood_hash.begin(); it!=flood_hash.end(); it++) {
    if (it->first.structure) free(it->first.structure);
  }
  flood_hash.clear();
  while (!flood_queue.empty()) flood_queue.pop();
  return res;
}

void DSU::PrintDot(char *filename, bool dot_prog, bool print, char *file_print, bool landmap)
{
  // landmap not supported yet

  float color = 0.7;
  //open file
  FILE *dot;
  dot = fopen(filename, "w");
  if (dot) {
    fprintf(dot, "Graph G {\n\tnode [width=0.1, height=0.1, shape=circle];\n");
    //nodes LM:
    for (unsigned int i=0; i<LM.size(); i++) {
      fprintf(dot, "\"%d\" [label=\"%d\"]\n", i+1, i+1);
    }
    fprintf(dot, "\n");
    //nodes saddle:
    for (unsigned int i=0; i<saddles.size(); i++) {
      fprintf(dot, "\"S%d\" [label=\"S%d\", color=\"0.0 0.0 %.1f\", fontcolor=\"0.0 0.0 %.1f\"]\n", i+1, i+1, color, color);
    }
    fprintf(dot, "\n");
    // edges l-l
    for (set<edgeLM>::iterator it=edges_l.begin(); it!=edges_l.end(); it++) {
      fprintf(dot, "\"%d\" -- \"%d\" [label=\"%.2f\"]\n", (it->i)+1, (it->j)+1, it->en/100.0);
    }
    fprintf(dot, "\n");
    // edges l-s
    for (set<std::pair<int, int> >::iterator it=edges_ls.begin(); it!=edges_ls.end(); it++) {
      fprintf(dot, "\"%d\" -- \"S%d\" [color=\"0.0 0.0 %.1f\", fontcolor=\"0.0 0.0 %.1f\"]\n", (it->first)+1, (it->second)+1, color, color);
    }

    fprintf(dot, "\n");
    // edges s-s //TODO!
    /*for (set<edgeLM>::iterator it=edges_l.begin(); it!=edges_l.end(); it++) {
      fprintf(dot, "\"%d\" -- \"%d\" [label=\"%.2f\"]\n",it->i,it->j, it->en/100.0);
    }*/
    fprintf(dot, "}\n");
  }

  fclose(dot);

  // start neato/dot:
  if (dot && print && file_print) {
    char syst[200];
    sprintf(syst, "%s -Tps < %s > %s", (dot_prog ? "dot" : "neato"), filename, file_print);
    int res = system(syst);
    printf("%s returned %d", syst, res);
  }
}

void DSU::PrintMatrix(char *filename)
{
  FILE *energies;
  energies = fopen(filename, "w");
  if (energies) {
    // create matrix
    vector<vector<float> > matrix;
    matrix.resize(LM.size());
    for (unsigned int i=0; i<LM.size(); i++) {
      matrix[i].resize(LM.size(), INFINITY);
    }

    // fill it with edges:
    for (set<edgeLM>::iterator it=edges_l.begin(); it!=edges_l.end(); it++) {
      matrix[it->i][it->j] = it->en/100.0;
      matrix[it->j][it->i] = it->en/100.0;
    }

    // resolve i==i
    for (unsigned int i=0; i<matrix.size(); i++) {
      matrix[i][i] = LM[i].energy/100.0;
    }

    // print
    for (unsigned int i=0; i<matrix.size(); i++) {
      for (unsigned int j=0; j<matrix[i].size(); j++) {
        fprintf(energies, "%6.2g ", matrix[i][j]);
      }
      fprintf(energies, "\n");
    }
  }
  fclose(energies);
}
