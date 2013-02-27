#include "DSUeval.h"

struct cluster_pair {
  int i, j;
  int d;

  cluster_pair(int i, int j, int d) {
    this->i = i;
    this->j = j;
    this->d = d;
  }

  bool operator <(const cluster_pair &sec) const {
    if (d != sec.d) {
      return d<sec.d;
    } else {
      if (i != sec.i) {
        return i<sec.i;
      } else {
        return j<sec.j;
      }
    }
  }
};

TBD::TBD()
{
  for (int i=0; i<type1_len; i++) {
    sizes[i] = 0;
  }
}

bool TBD::insert(int i, int j, type1 type, bool fiber)
{
  sizes[type]++;

  if (done.count(make_pair(i,j))>0) return false;
  else {
    fprintf(stderr, "inserting %d %d %s\n", i, j, type1_str[type]);
    tbd.push(TBDentry(i, j, type, fiber));
    return true;
  }
}

int TBD::size() {
  return tbd.size();
}

TBDentry TBD::get_first()
{
  if (size()==0) return TBDentry(-1,-1,NEW_FOUND,-1);
  TBDentry tbde = tbd.top();
  tbd.pop();
  done.insert(make_pair(tbde.i, tbde.j));
  return tbde;
}


int DSU::Cluster(int kmax, TBD &output)
{
  // create data structures
  vector<cluster_pair> to_cluster;
  to_cluster.reserve(LM.size());
  UF_set_child ufset;
  ufset.enlarge_parent(LM.size());

  // fill it
  for (unsigned int i=0; i<LM.size(); i++) {
    for (unsigned int j=i+1; j<LM.size(); j++) {
      to_cluster.push_back(cluster_pair(i,j,HammingDist(LM[i].structure, LM[j].structure)));
    }
  }
  sort(to_cluster.begin(), to_cluster.end());

  // process:
  int last_hd = to_cluster[0].d;
  for (unsigned int i=0; i<to_cluster.size(); i++) {
    cluster_pair &cp = to_cluster[i];

    //fprintf(stderr, "clustering %d %d (%d)\n", cp.i, cp.j, cp.d);

    if (cp.d!=last_hd) {
      // do something, cause we are on higher level...
    }

    // see if we are not joint yet:
    if (!ufset.joint(cp.i, cp.j)) {
      // try to connect
      if (ufset.count(ufset.find(cp.i)) + ufset.count(ufset.find(cp.j)) > kmax) { // cannot connect them, need to insert all edges into the TBD
        // insert crit edge:
        output.insert(cp.i, cp.j, CRIT_EDGE, false);

        set<int> first = ufset.get_children(ufset.find(cp.i));
        set<int> second = ufset.get_children(ufset.find(cp.j));

        // insert rpresentative edge
        output.insert(*first.begin(), *second.begin(), REPRESENT, false);

        // insert all inter edges:
        for (set<int>::iterator it=first.begin(); it!=first.end(); it++) {
          set<int>::iterator it2 = it; it2++;
          for (;it2!=first.end(); it2++) {
            output.insert(*it, *it2, INTER_CLUSTER, false);
          }
        }

        for (set<int>::iterator it=second.begin(); it!=second.end(); it++) {
          set<int>::iterator it2 = it; it2++;
          for (;it2!=second.end(); it2++) {
            output.insert(*it, *it2, INTER_CLUSTER, false);
          }
        }

        // now make from this group only one vertex (maybe wrong)
        ufset.union_set(cp.i, cp.j);
        ufset.make_single(cp.i);


      } else {
        // connect them
        ufset.union_set(cp.i, cp.j);

      }
    }
    last_hd = cp.d;
  }

  fprintf(stderr, "output size = %d = (%d, %d, %d)\n", output.size(), output.sizes[0], output.sizes[1], output.sizes[2]);

  return 0;
}
