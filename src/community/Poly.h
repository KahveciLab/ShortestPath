#ifndef POLY_H_
#define POLY_H_
#include <vector>
#include <iostream>
#include <boost/dynamic_bitset.hpp>
using namespace std;
using namespace boost;

class Poly{
public:
    double prob;
    dynamic_bitset<> usedEdges;
    dynamic_bitset<> possPaths;
    vector<int> finishedPaths;

public:
	Poly(int numEdges, int numPaths){
		usedEdges.resize(numEdges);
		possPaths.resize(numPaths);
	}
};

#endif /* POLY_H_ */