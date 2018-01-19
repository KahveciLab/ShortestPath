#include "Graph.h"
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <unordered_map>
#include <iterator>

void Graph::readNetworkFile(char* fileName){
	int sourceNode, destNode, mappingCounter = 0;
	ifstream in_stream;
	string sourceName, destName, edgeProb,str1;
	int edgeIndex=0;
	in_stream.open(fileName);
	vector<string> strs;
	map<string,int> nodeNameToNoMap;

	while(!in_stream.eof()) {
		getline(in_stream, str1);
		if(str1.size()>0){
			strs.clear();
			boost::split(strs,str1,boost::is_any_of("\t"));
			sourceName=strs[0];
			destName=strs[1];
			edgeProb=strs[2];
			if(edgeProb.compare("-")==0)
			continue;

			if(nodeNameToNoMap.find(sourceName) == nodeNameToNoMap.end()) {
				nodeNameToNoMap.insert(pair<string,int>(sourceName, mappingCounter));
				nodeNameMap.insert(pair<int,string>(mappingCounter,sourceName));
				++mappingCounter;

			}

			if(nodeNameToNoMap.find(destName) == nodeNameToNoMap.end()) {
				nodeNameToNoMap.insert(pair<string,int>(destName, mappingCounter));
				nodeNameMap.insert(pair<int,string>(mappingCounter,destName));
				++mappingCounter;
			}

			sourceNode = nodeNameToNoMap[sourceName];
			destNode = nodeNameToNoMap[destName];

			if(adjacencyList[sourceNode][destNode]==0){
				adjacencyMatrix[sourceNode][destNode] = edgeIndex;
				adjacencyMatrix[destNode][sourceNode] = edgeIndex;

				adjacencyList[sourceNode][destNode] = 1;
				undirectedAL[sourceNode][destNode]=1;
				undirectedAL[destNode][sourceNode]=1;
				originalAdjacencyList[sourceNode][destNode] = 1;
				originalAdjacencyList[destNode][sourceNode] = 1;

				edgeProbabilities.push_back(strtof(edgeProb.c_str(), 0));
				edgeToNodeMap.insert(pair<int, pair<int,int> >(edgeIndex, pair<int,int>(sourceNode,destNode)));
				edgeIndex++;
			}
		}
	}

	in_stream.close();
	edgeSize=edgeIndex;
	strs.clear();
	nodeNameToNoMap.clear();
}

void Graph::cal_coeff(dynamic_bitset<>& edges, rowvec& coeff){
	coeff(0)=1;
	dynamic_bitset<>::size_type it;
	it=edges.find_first();
	int size=1;
	while(it!=dynamic_bitset<>::npos){
		double p=edgeProbabilities[it];
		rowvec temp(coeff.cols(0,size-1));
		coeff *= (1-p);
		temp *=p;
		coeff.cols(1,size) +=temp;
		size++;
		it=edges.find_next(it);
	}
}

double Graph::multiplyPoly(dynamic_bitset<>& within_edges, dynamic_bitset<>& connect_edges, dynamic_bitset<>& other_edges){
	if(within_edges.none() && connect_edges.none())
	return 0;
	else{
		//multiply the polynomials for edges within the community
		int num_within=within_edges.count()+1;
		rowvec res1(num_within,fill::zeros);
		cal_coeff(within_edges,res1);

		int num_con=connect_edges.count()+1;
		rowvec res2(num_con,fill::zeros);
		cal_coeff(connect_edges,res2);

		int num_others=other_edges.count()+1;
		rowvec res3(num_others,fill::zeros);
		cal_coeff(other_edges,res3);

		double res=0;
		for(int i=0;i<num_within;i++)
		for(int j=0;j<num_con;j++)
		{
			if(i==0 && j==0)
			continue;
			for(int k=0;k<num_others;k++){
				double coeff=res1(i)*res2(j)*res3(k);
				int lc=i;
				int m=i+j+k;
				int dc=2*i+j;
				res+=coeff*((i*1.0)/m - (dc*dc*1.0)/(4*m*m));
			}
		}
		return res;
	}
}

double Graph::calModularity(vector<vector<int> >& components){
	double res=0;

	for(int i=0;i<components.size();i++){
		dynamic_bitset<> within_edges(edgeSize);
		dynamic_bitset<> connect_edges(edgeSize);
		dynamic_bitset<> other_edges(edgeSize);

		for(int j=0;j<components[i].size();j++)
		{
			int node=components[i][j];
			dynamic_bitset<>::size_type it;
			it=originalAdjacencyList[node].find_first();
			while(it!=dynamic_bitset<>::npos){
				if(member[it]==member[node])
				within_edges[adjacencyMatrix[node][it]]=1;

				else connect_edges[adjacencyMatrix[node][it]]=1;
				it=originalAdjacencyList[node].find_next(it);
			}
		}
		other_edges= within_edges;
		other_edges |= connect_edges;
		other_edges = other_edges.flip();
		res+=multiplyPoly(within_edges,connect_edges,other_edges);
	}

	return res;
}

double Graph::findCommunity(){
	vector<vector<int> > components;
	dynamic_bitset<> visited(nodeSize);
	dynamic_bitset<> edgeSameComponent(edgeSize);

	int index=0;
	for(int j=0;j<nodeSize;j++){
		if(visited[j]==0){
			vector<int> temp;
			components.push_back(temp);
			dfs(components,visited,j,index);
			index++;
		}
	}

	int initialComNum=index;
	vector<int> nodes;
	for(int i=0;i<nodeSize;i++)
	  nodes.push_back(i);

	for(int i=0;i<edgeSize;i++){
		edgeSameComponent.reset();
		calculateBetweenness(nodes, edgeSameComponent);
		nodes.clear();

		for(int j=0;j<edgeSize;j++){
			if(edgeSameComponent[j]==1)
			   edgeBetweenness[j] /= edgeProbabilities[j];
		}

		int pos=std::distance(edgeBetweenness.begin(), max_element(edgeBetweenness.begin(), edgeBetweenness.end()));
		int sourceNode=edgeToNodeMap[pos].first;
		int destNode=edgeToNodeMap[pos].second;
		adjacencyList[sourceNode][destNode]=0;
		undirectedAL[sourceNode][destNode]=0;
		undirectedAL[destNode][sourceNode]=0;

		edgeBetweenness[adjacencyMatrix[sourceNode][destNode]]=0;

		for(int j=0;j<nodeSize;j++){
			if(member[j]==member[sourceNode]){
				nodes.push_back(j);
			}
		}
		//recalculate community value, it is based on the original graph
		// The degree of each node is fixed, use the full network to calculate the modularity, even some edges have been removed
		index=0;
		components.clear();
		visited.reset();
		//first dfs

		for(int j=0;j<nodeSize;j++){
			if(visited[j]==0){
				vector<int> temp;
				components.push_back(temp);
				dfs(components,visited,j,index);
				index++;
			}
		}

		if(index>initialComNum){
			double res=0;
			res=calModularity(components);

			if(res<0)
			res=0;
			//            cout<<"ModularityValue: "<<res<<" and community number: "<<index<<endl;

			if(res>maxV){
				maxV=res;
				communityV=components;
			}

			initialComNum=index;
		}
	}

	//print community
	cout<<"---------------print community----------------"<<endl;
	for(int i=0;i<communityV.size();i++){
		cout<<"community "<<i<<": ";
		for(int j=0;j<communityV[i].size();j++){
			member[communityV[i][j]]=i;
			cout<<nodeNameMap[communityV[i][j]]<<"\t";
		}
		cout<<endl;
  }
	cout<<"The expected modularity value is "<<maxV<<endl;
	return maxV;

}

void Graph::dfs(vector<vector<int> >& components, dynamic_bitset<>& visited, int curr, int index){
	components[index].push_back(curr);
	member[curr]=index;
	visited[curr]=1;
	dynamic_bitset<>::size_type it;
	it=undirectedAL[curr].find_first();
	while(it!=dynamic_bitset<>::npos){
		if(visited[it]==0){
			dfs(components,visited,it,index);
		}
		it=undirectedAL[curr].find_next(it);
	}
}

void Graph::calculateBetweenness(vector<int>& nodes, dynamic_bitset<>& edgeSameComponent){
	vector<vector<int> > res;
	for(int ni=0;ni<nodes.size();ni++){
		int i=nodes[ni];
		if(adjacencyList[i].any()){
			for(int nj=0;nj<nodes.size();nj++){
				int j=nodes[nj];
				if(i==j || member[i]!=member[j])
				continue;

				findAllSimplePath(res, i, j);
				if(res.size()>0){
					if(res.size()==1){

						for(int k=0;k<res[0].size();k++){
							int edgeIndex=res[0][k];
							if(edgeSameComponent[edgeIndex]==0){
								edgeBetweenness[edgeIndex]=1;
								edgeSameComponent[edgeIndex]=1;
							}
							else edgeBetweenness[edgeIndex]+=1;
						}

					}
					else{
						vector<double> tempBetweenness;
						tempBetweenness.resize(edgeSize);
						std::fill(tempBetweenness.begin(),tempBetweenness.end(),0);
						double proSum=0;

						//sort path based on the length
						sort(res.begin(),res.end(),[](const vector<int>& a, const vector<int>& b){return a.size() < b.size();});
						int pathNum=res.size();

						//compress path
						//1. build edge matrix ( edge index, path bit set )
						unordered_map<int, dynamic_bitset<> > edgeMatrix;
						for(int k=0;k<res.size();k++){
							for(int l=0;l<res[k].size();l++){
								int index=res[k][l];
								if(edgeMatrix.find(index)==edgeMatrix.end()){
									dynamic_bitset<> temp(pathNum);
									temp[k]=1;
									edgeMatrix.insert({index, temp});
								}
								else edgeMatrix[index][k]=1;
							}
						}

						//2. compress, reNumber the edge Index
						unordered_map<int, int> edgeIndexMap; //old edge index and renumbered edge index
						unordered_map<string, int> bitEdgeMap; //path list and renumbered edge index
						vector<double> edgeProb; // renumbered edge probabilities
						int edgeCount=0;
						vector<vector<int> >compressEdges;
						vector<dynamic_bitset<> > newEdgeMatrix;
						for(std::pair<int, dynamic_bitset<> > element : edgeMatrix){
							string s;
							to_string(element.second, s);
							if(bitEdgeMap.find(s)==bitEdgeMap.end()){
								bitEdgeMap.insert({s, edgeCount});
								edgeIndexMap.insert({element.first, edgeCount});
								vector<int> temp;
								temp.push_back(element.first);
								compressEdges.push_back(temp);
								edgeProb.push_back(edgeProbabilities[element.first]);
								newEdgeMatrix.push_back(element.second);
								edgeCount++;
							}
							else{
								int eIndex=bitEdgeMap[s];
								edgeIndexMap.insert({element.first,eIndex});
								compressEdges[eIndex].push_back(element.first);
								edgeProb[eIndex]*=edgeProbabilities[element.first];
							}
						}
						//                       cout<<"Compress edges from "<<edgeMatrix.size()<<" to "<<compressEdges.size()<<endl;
						edgeMatrix.clear();
						bitEdgeMap.clear();

						//generate Path
						Path paths[pathNum];
						int preIndex=pathNum;
						int preLen=res[pathNum-1].size();
						for(int k=pathNum-1;k>=0;k--){
							paths[k].setEdgeSize(edgeCount);
							if(res[k].size()!=preLen){
								preIndex=k+1;
								preLen=res[k].size();
							}
							paths[k].nextLenIndex=preIndex;
							for(int l=0;l<res[k].size();l++){
								paths[k].edges[edgeIndexMap[res[k][l]]]=1;
							}
						}

						//deal with polynomials
						//                        cout<<"polynomials start"<<endl;
						Poly a(edgeCount, pathNum);
						a.prob=1;
						a.possPaths.set();
						paths[0].polyTerms.push_back(a);
						dynamic_bitset<> tbitset;
						for(int pIndex=0;pIndex<pathNum;pIndex++){
							for(int k=0; k< paths[pIndex].polyTerms.size(); k++){

								dynamic_bitset<>::size_type it;
								it=paths[pIndex].edges.find_first();
								int numPossPath=paths[pIndex].polyTerms[k].possPaths.count();

								while(it!=dynamic_bitset<>::npos){
									if(paths[pIndex].polyTerms[k].usedEdges[it]==0){
										if(numPossPath!=1){
											paths[pIndex].polyTerms[k].usedEdges[it]=1;
											tbitset = paths[pIndex].polyTerms[k].possPaths;
											tbitset ^= newEdgeMatrix[it];
											tbitset &= paths[pIndex].polyTerms[k].possPaths;

											if(tbitset.none()){
												if(paths[pIndex].polyTerms[k].finishedPaths.size()>0){
													proSum+=paths[pIndex].polyTerms[k].prob * (1-edgeProb[it]);//newly written
													collapse(paths[pIndex].polyTerms[k].prob * (1-edgeProb[it]),compressEdges, paths, paths[pIndex].polyTerms[k].finishedPaths, edgeSameComponent, tempBetweenness);
												}
											}
											else if(tbitset.find_next(pIndex)>=paths[pIndex].nextLenIndex && paths[pIndex].polyTerms[k].finishedPaths.size()>0){
												proSum+=paths[pIndex].polyTerms[k].prob * (1-edgeProb[it]);//newly written
												collapse(paths[pIndex].polyTerms[k].prob * (1-edgeProb[it]),compressEdges, paths, paths[pIndex].polyTerms[k].finishedPaths, edgeSameComponent, tempBetweenness);
											}
											else{
												Poly temp(paths[pIndex].polyTerms[k]);
												temp.prob *= (1-edgeProb[it]);
												temp.possPaths = tbitset;
												paths[tbitset.find_next(pIndex)].polyTerms.push_back(temp);
											}
										}

										paths[pIndex].polyTerms[k].prob *= edgeProb[it];
									}
									it=paths[pIndex].edges.find_next(it);
								}

								paths[pIndex].polyTerms[k].possPaths[pIndex]=0;
								paths[pIndex].polyTerms[k].finishedPaths.push_back(pIndex);
								if(paths[pIndex].polyTerms[k].possPaths.none() || paths[pIndex].polyTerms[k].possPaths.find_next(pIndex)>=paths[pIndex].nextLenIndex)
								{
									proSum+=paths[pIndex].polyTerms[k].prob;//newly written
									collapse(paths[pIndex].polyTerms[k].prob,compressEdges, paths, paths[pIndex].polyTerms[k].finishedPaths, edgeSameComponent, tempBetweenness);
								}
								else{
									paths[paths[pIndex].polyTerms[k].possPaths.find_next(pIndex)].polyTerms.push_back(paths[pIndex].polyTerms[k]);
								}

							}

							paths[pIndex].polyTerms.clear();
						}
						for(int k=0;k<edgeSize;k++){
							if(tempBetweenness[k]>0)
							edgeBetweenness[k]+=tempBetweenness[k]/proSum;
						}
					}


				}
				res.clear();
			}
		}
	}

}

void Graph::collapse(double val, vector<vector<int> >& compressEdges, Path paths[], vector<int>& finishedPaths, dynamic_bitset<>& edgeSameComponent, vector<double>& tempBetweenness){
	int num=finishedPaths.size();
	double factor=1.0/num;
	for(int i=0;i<finishedPaths.size();i++){
		dynamic_bitset<>::size_type it;
		int pIndex=finishedPaths[i];
		it=paths[pIndex].edges.find_first();
		while(it!=dynamic_bitset<>::npos){
			for(int j=0;j<compressEdges[it].size();j++){
				int edgeIndex=compressEdges[it][j];
				if(edgeSameComponent[edgeIndex]==0){
					edgeBetweenness[edgeIndex]=0; //newly written
					edgeSameComponent[edgeIndex]=1;
				}
				tempBetweenness[edgeIndex]+=val*factor;//newly written
			}

			it=paths[pIndex].edges.find_next(it);
		}
	}
}

void Graph::findAllSimplePath(vector<vector<int> >& res, int source, int end){
	vector<int> path;
	dynamic_bitset<> visited(nodeSize);
	dfs(res, path, visited, source, end);
}

void Graph::dfs(vector<vector<int> >& res, vector<int>& path, dynamic_bitset<>& visited, int curr, int end){
	visited[curr]=1;
	dynamic_bitset<>::size_type it;
	it=adjacencyList[curr].find_first();
	while(it!=dynamic_bitset<>::npos){
		path.push_back(adjacencyMatrix[curr][it]);
		if(it==end){
			res.push_back(path);
		}
		else if(visited[it]==0){
			dfs(res,path,visited, it, end);
		}
		path.pop_back();
		it=adjacencyList[curr].find_next(it);
	}
	visited[curr]=0;
}
