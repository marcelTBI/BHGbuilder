#include <stdio.h>
#include <string.h>

extern "C" {
  #include "pair_mat.h"
}

#include "RNAutils.h"

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

char *pt_to_char(const short *pt)
{
  char *str = (char*) malloc(pt[0]*(sizeof(char)+1));
  for (int i=1; i<=pt[0]; i++) {
    if (pt[i]==0) str[i-1]='.';
    else if (pt[i]<i) str[i-1]=')';
    else str[i-1]='(';
  }
  str[pt[0]]='\0';
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

bool isStruct(const char *p)
{
  // check first two chars - should be enough
  if (strlen(p)<2) return false;
  if ((p[0]=='.' || p[0]=='(' || p[0]==')') && (p[1]=='.' || p[1]=='(' || p[1]==')')) return true;
  else return false;
}

int HammingDist(char* struct1, const short* struct2) {
  short* s1 = make_pair_table(struct1);

  int bpdist = HammingDist(s1, struct2);

  free(s1);

  return bpdist;
}

int HammingDist(char* struct1, char* struct2) {
  short* s1 = make_pair_table(struct1);
  short* s2 = make_pair_table(struct2);

  int bpdist = HammingDist(s1, s2);

  free(s1);
  free(s2);

  return bpdist;
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

bool isSeq(const char *p)
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

inline char to16(unsigned short num) {
  if (num<10) return '0'+num;
  else return 'a'+num-10;
}

//convert rgb to #rgb
char *rgb(unsigned short red, unsigned short green, unsigned short blue)
{
  static char rgbchar[8];
  rgbchar[0]='#';

  red %= 256;
  blue %= 256;
  green %= 256;

  rgbchar[1]=to16(red/16);
  rgbchar[2]=to16(red%16);
  rgbchar[3]=to16(green/16);
  rgbchar[4]=to16(green%16);
  rgbchar[5]=to16(blue/16);
  rgbchar[6]=to16(blue%16);
  rgbchar[7]='\0';

  return rgbchar;
}

char *rgb_d(double red, double green, double blue) {
  return rgb((unsigned short) red*255, (unsigned short) green*255, (unsigned short) blue*255);
}

// UNION FIND set functions
UF_set::UF_set() {
  parent.clear();
  num_unions = 0;
}

int UF_set::find(int x, bool fix) {
  if (x >= (int) parent.size()) {
    if (!fix) return -1;
    else enlarge_parent(x+1);
  }
  if (x != parent[x] && parent[x] != parent[parent[x]])
    parent[x] = find(parent[x]);
  return parent[x];
}

void UF_set::union_set(int x, int y) {
  int u, v;
  u = find(x, true);
  v = find(y, true);
  if (u != v) {
    parent[u] = v;
    num_unions++;
  }
}

bool UF_set::connected_all() {
  return (num_unions == parent.size()-1);
}

bool UF_set::joint(int x, int y) {
  return find(x) == find(y);
}

void UF_set::enlarge_parent() {
  parent.push_back(parent.size());
}

void UF_set::enlarge_parent(int cnt) {
  while (size()<cnt) enlarge_parent();
}

int UF_set::size() {
  return parent.size();
}

void UF_set::clear() {
  parent.clear();
  num_unions = 0;
}

vector<int> UF_set::get_parents() {
  vector<int> res;
  for (int i=0; i<(int)parent.size(); i++) {
    if (find(i)==i) res.push_back(i);
  }

  return res;
}

map<int, int> UF_set::get_invert() {
  map<int, int> res;
  int count = 0;
  for (int i=0; i<(int)parent.size(); i++) {
    if (find(i)==i) res[i] = count++;
  }

  return res;
}

// UF_set_child

set<int> UF_set_child::get_children(int which)
{
  if (which < (int)children.size() && which>=0) return children[which];
  else return set<int>();
}

int UF_set_child::count(int which)
{
  if (which < (int)children.size() && which>=0) return children[which].size();
  else return -1;
}

void UF_set_child::make_single(int which)
{
  // just try to simulate that we have only one vertex, so count = 1 and children = this
  if (which < (int)children.size() && which>=0) {
    which = find(which);
    reduced += children[which].size()-1;
    children[which].clear();
    children[which].insert(which);
  }
}

UF_set_child::UF_set_child():
  UF_set()
{
  clear();
}

void UF_set_child::enlarge_parent() {
  UF_set::enlarge_parent();
  set<int> s;
  s.insert(size()-1);
  children.push_back(s);
}

void UF_set_child::enlarge_parent(int cnt) {
  for (int i=0; i<cnt; i++) enlarge_parent();
}

void UF_set_child::clear() {
  UF_set::clear();
  children.clear();
  reduced = 0;
}

int UF_set_child::dimension() {
  return size()-reduced;
}

void UF_set_child::union_set(int x, int y) {
  int u, v;
  u = find(x);
  v = find(y);
  if (v>u) swap(u,v);
  if (u != v && u>=0 && v>=0) {
    UF_set::union_set(u,v);
    children[v].insert(children[u].begin(), children[u].end());
    children[u].clear();
  }
}
