#ifndef GRAPH_H_
#define GRAPH_H_
#include <map>
#include <vector>
#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include <string>
#include <armadillo>
using namespace arma;
using namespace std;
using namespace boost;

class Graph{
public:
	int nodeSize;
	int edgeSize;
	vector<dynamic_bitset<> > adjacencyList;
	vector<vector<int> > adjacencyMatrix;
	vector<double> edgeProbabilities;
	map<int, string> nodeNameMap;
	vec nonzeroExpectPath;
	typedef Mat<unsigned char> bmat;
	typedef Row<unsigned char> brow;
	vec nonzeroDeterPath;
	map<int, pair<int,int> > edgeToNodeMap;
	uvec pq;

public:
	Graph(int nodeNum, char* fileName): nodeSize(nodeNum), adjacencyList(nodeNum, dynamic_bitset<>(nodeNum)){
		adjacencyMatrix.resize(nodeNum);
		for(int i=0;i<nodeNum;i++){
			adjacencyMatrix[i].resize(nodeNum);
			std::fill(adjacencyMatrix[i].begin(),adjacencyMatrix[i].end(),-1);
		}
		readNetworkFile(fileName);
	}

	void countShortestPaths();
	double countBinaryShortestPaths();
	double countThresholdShortestPaths(double threshold);
	double countSampleShortestPaths();

private:
	void findAllSimplePath(vector<dynamic_bitset<> >& paths, int source, int end);
	void readNetworkFile(char* fileName);
	void dfs(vector<dynamic_bitset<> >& paths, dynamic_bitset<>& path, dynamic_bitset<>& visited, int curr, int end);
};

#endif /* GRAPH_H_ */
