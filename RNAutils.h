#ifndef __RNAUTILS_H
#define __RNAUTILS_H

#include <string>
#include <vector>

using namespace std;

// file with simple, but widely used utilities in RNA struct programs.

// reads line no matter how long - returns solid pointer - must be freed
char* my_getline(FILE *fp);
// convert structure short* to string
string pt_to_str(const short *pt);
char *pt_to_char(const short *pt); // must be freed!!!
// are structures equal?
bool str_eq(const short *lhs, const short *rhs);
// convert energy from float to int
int en_fltoi(float en);
// hamming distance beetween 2 structs
int HammingDist(const short* struct1, const short* struct2);
// is string an RNA structure?
bool isStruct(const char *p);
// is string an RNA sequence?
bool isSeq(const char *p);

//convert rgb to #rgb
char *rgb(unsigned short red, unsigned short green, unsigned short blue);
char *rgb_d(double red, double green, double blue);

// class for encapsulating Union-Find set operations
class UF_set {
private:
  vector<int> parent;
  unsigned int num_unions;

public:
  // and union-find set functions
  UF_set();
  int find(int x);
  void union_set(int x, int y);
  bool connected_all() ;
  bool joint(int x, int y);
  void enlarge_parent();
  void enlarge_parent(int count);
  int size();
  void clear();
};

#endif
