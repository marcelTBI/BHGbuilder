#ifndef __RNAUTILS_H
#define __RNAUTILS_H

#include <string>

using namespace std;

// file with simple, but widely used utilities in RNA struct programs.

// reads line no matter how long - returns solid pointer - must be freed
char* my_getline(FILE *fp);
// convert structure short* to string
string pt_to_str(const short *pt);
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

#endif
