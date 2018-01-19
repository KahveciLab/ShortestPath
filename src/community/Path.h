#ifndef PATH_H_
#define PATH_H_
#include "Poly.h"
#include <vector>
#include <iostream>
#include <boost/dynamic_bitset.hpp>
using namespace std;
using namespace boost;

class Path{
public:
   int nextLenIndex;
   dynamic_bitset<> edges;
   vector<Poly> polyTerms;

public:
	Path(): nextLenIndex(-1){};
	void setEdgeSize(int edgeSize){
		edges.resize(edgeSize);
	}
};

#endif /* PATH_H_ */