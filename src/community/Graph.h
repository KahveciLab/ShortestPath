#ifndef GRAPH_H_
#define GRAPH_H_
#include <map>
#include <vector>
#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include <string>
#include "Path.h"
#include <armadillo>
using namespace arma;
using namespace std;
using namespace boost;

class Graph{
public:
	int nodeSize;
	int edgeSize;
	vector<dynamic_bitset<> > adjacencyList; //directed adjacency list
	vector<dynamic_bitset<> > originalAdjacencyList; //original undirected adjacency list
	vector<dynamic_bitset<> > undirectedAL; //undirected adjacency list
	vector<vector<int> > adjacencyMatrix;
	vector<double> edgeProbabilities;

	map<int, pair<int,int> > edgeToNodeMap;
	map<int, string> nodeNameMap;

	vector<double> edgeBetweenness;
	vector<vector<int> > communityV;
	double maxV; //modularity value
	vector<int> member; //indicate membership

public:
	Graph(int nodeNum, char* fileName): nodeSize(nodeNum), adjacencyList(nodeNum, dynamic_bitset<>(nodeNum)), undirectedAL(nodeNum, dynamic_bitset<>(nodeNum)), originalAdjacencyList(nodeNum, dynamic_bitset<>(nodeNum)){
		adjacencyMatrix.resize(nodeNum);
		for(int i=0;i<nodeNum;i++){
			adjacencyMatrix[i].resize(nodeNum);
			std::fill(adjacencyMatrix[i].begin(),adjacencyMatrix[i].end(),-1);
		}
		maxV=0;

		readNetworkFile(fileName);
		edgeBetweenness.resize(edgeSize);
		std::fill(edgeBetweenness.begin(),edgeBetweenness.end(),0);

		member.resize(nodeSize);
		std::fill(member.begin(),member.end(),0);
	}

	void printGraph(){
		dynamic_bitset<>::size_type it;
		for(int i=0;i<nodeSize;i++){
			it=adjacencyList[i].find_first();
			while(it!=dynamic_bitset<>::npos){
				int index=adjacencyMatrix[i][it];
				cout<<"edge between "<<nodeNameMap[i]<<" and "<<nodeNameMap[it]<<" pro: "<<edgeProbabilities[index]<<endl;
				it=adjacencyList[i].find_next(it);
			}
		}
		cout<<"Node and edge: "<<nodeSize<<"\t"<<edgeSize<<endl;
	}

	void calculateBetweenness(vector<int>& nodes, dynamic_bitset<>& edgeSameComponent);
	double findCommunity();
	double calModularity(vector<vector<int> >& components);

private:
	void findAllSimplePath(vector<vector<int> >& res, int source, int end);
	void readNetworkFile(char* fileName);
	void dfs(vector<vector<int> >& res, vector<int>& path, dynamic_bitset<>& visited, int curr, int end);
	void dfs(vector<vector<int> >& components, dynamic_bitset<>& visited, int curr, int index);
	void collapse(double val, vector<vector<int> >& compressEdges, Path paths[], vector<int>& finishedPaths, dynamic_bitset<>& edgeSameComponent, vector<double>& tempBetweenness);
	double multiplyPoly(dynamic_bitset<>& within_edges, dynamic_bitset<>& connect_edges, dynamic_bitset<>& other_edges);
	void cal_coeff(dynamic_bitset<>& edges, rowvec& coeff);
};

#endif /* GRAPH_H_ */
