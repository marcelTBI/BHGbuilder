#include <stdio.h>
#include <string.h>

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
