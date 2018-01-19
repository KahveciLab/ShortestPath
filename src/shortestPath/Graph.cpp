#include "Graph.h"
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <iterator>
#include <random>
#include <cmath>

struct TreeNode{
  int id;
  int distance;
  int weight;
  TreeNode(int d){
    id=d;
    distance=weight=0;
  }
};

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
        adjacencyList[sourceNode][destNode] = 1;

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

void Graph::countShortestPaths(){
  mat expectPath;
  expectPath.zeros(nodeSize,nodeSize);
  dynamic_bitset<>::size_type it;

  for(int i=0;i<nodeSize;i++)
    for(int j=0;j<nodeSize;j++){
      if(i==j){
        continue;
      }
//      cout<<"print "<<i<<"\t"<<j<<endl;
      vector<dynamic_bitset<> > paths;
      findAllSimplePath(paths, i, j);
//      cout<<paths.size()<<endl;
      if(paths.size()>0){
        if(paths.size()==1){ // expected number of shortest paths is the probability of its existence
          double expNum=1;
          it=paths[0].find_first();
          while(it!=dynamic_bitset<>::npos){
            expNum*=edgeProbabilities[it];
            it=paths[0].find_next(it);
          }
          expectPath(i,j)=expNum;
        }
        else{
          //sort path based on the length
          sort(paths.begin(),paths.end(),[](const dynamic_bitset<>& a, const dynamic_bitset<>& b){return a.count() < b.count();});
          int pathNum=paths.size();

          vector<int> group; //the last edge's index of each group
          int initalLen=paths[0].count();
          for(int k=1;k<pathNum;k++){
            if(paths[k].count()>initalLen){
              group.push_back(k);
              initalLen=paths[k].count();
            }
          }
          group.push_back(pathNum);

            //create edge matrix
            vector<dynamic_bitset<> > edgeVector;
            unordered_map<int,int> edgeIndexMap; //original edge index to the matrix edge index
            int edgeCount=0;
            vector<int> edgeGroup;
            int lastIndex=0;
            vector<double> tempEdgeProb;

            for(int m=0;m<group.size();m++){
              for(int k=lastIndex;k<group[m];k++){
                it=paths[k].find_first();
                while(it!=dynamic_bitset<>::npos){
                  if(edgeIndexMap.find(it)==edgeIndexMap.end()){
                    edgeIndexMap.insert({it,edgeCount});
                    dynamic_bitset<> eachEdge(pathNum);
                    eachEdge.set();
                    eachEdge[k]=0;
                    edgeVector.push_back(eachEdge);
                    tempEdgeProb.push_back(edgeProbabilities[it]);
                    edgeCount++;
                  }
                  else{
                    int eIndex=edgeIndexMap[it];
                    edgeVector[eIndex][k]=0;
                  }
                  it=paths[k].find_next(it);
                }
              }
              edgeGroup.push_back(edgeCount);
              lastIndex=group[m];
            }
            paths.clear();

            //print edgeMatrix
            // for(int k=0;k<edgeVector.size();k++)
            // cout<<edgeVector[k]<<endl;

//           cout<<"edge vecor: "<<edgeVector.size()<<"\t group: "<<edgeGroup.size()<<endl;
            //print edge vector
            // for(int k=0;k<edgeGroup.size();k++)
            // cout<<edgeGroup[k]<<"\t";
            // cout<<endl;

            edgeCount =0;
            lastIndex=0;
            vector<dynamic_bitset<> > compressEdgeVector;
            vector<int> compressEdgeGroup;
            vector<double> edgeProb; // renumbered edge probabilities
            for(int m=0;m<edgeGroup.size();m++){
                for(int k=lastIndex;k<edgeGroup[m];k++){
                    if(k==lastIndex){
                      compressEdgeVector.push_back(edgeVector[k]);
                      edgeProb.push_back(tempEdgeProb[k]);
                      edgeCount++;
                    }
                    else{
                      int indexStart;
                      if(compressEdgeGroup.size()==0)
                         indexStart = 0;
                      else indexStart = compressEdgeGroup[m-1];
                      dynamic_bitset<> tdb(pathNum);
                      bool found=false;
                      for(int l=indexStart;l<compressEdgeVector.size();l++){
                         tdb = edgeVector[k];
                         tdb ^= compressEdgeVector[l];
                         if(tdb.none()){
                           edgeProb[l] *= tempEdgeProb[k];
                           found = true;
                           break;
                         }
                      }
                      if(!found){
                        compressEdgeVector.push_back(edgeVector[k]);
                        edgeProb.push_back(tempEdgeProb[k]);
                        edgeCount++;
                      }
                    }
                }
                lastIndex=edgeGroup[m];
                compressEdgeGroup.push_back(edgeCount);
            }
            edgeVector.clear();
            edgeGroup.clear();
            tempEdgeProb.clear();

//            cout<<"after compress: "<<compressEdgeVector.size()<<endl;
//            cout<<"edge vecor: "<<compressEdgeVector.size()<<"\t group: "<<compressEdgeGroup.size()<<"\t prob: "<<edgeProb.size()<<endl;
            //print edgeMatrix
            // for(int k=0;k<compressEdgeVector.size();k++)
            // cout<<compressEdgeVector[k]<<endl;
            //
            //
            // //print edge vector
            // for(int k=0;k<compressEdgeGroup.size();k++)
            // cout<<compressEdgeGroup[k]<<"\t";
            // cout<<endl;

            uvec q0;
            umat edgeMatrix;
            string str;
            vector<int> numVector;
            dynamic_bitset<> init(pathNum);
            init.set();
            compressEdgeVector.push_back(init);
            for(int k=0;k<compressEdgeVector.size();k++){
               lastIndex=0;
               umat eachCol;
               for(int m=0;m<group.size();m++){
                 int numSize = group[m]-lastIndex;

                 int blockNum = ceil(numSize/64.0);
                 umat eachBlock(blockNum,1);
                 if(k==0){
                     numVector.push_back(blockNum);
                 }
                 for(int l=0;l<blockNum;l++){
                   int subsize = (64 < numSize)? 64 : numSize;
                   dynamic_bitset<> tdb = compressEdgeVector[k];
                   tdb.resize(subsize);
//                   cout<<"tdb: "<<tdb<<endl;
                   numSize -= subsize;
                   eachBlock(l) = tdb.to_ulong();
                   compressEdgeVector[k] >>= subsize;
//                   cout<<"com: "<<compressEdgeVector[k]<<endl;
                 }
                 eachCol.insert_rows(eachCol.n_rows, eachBlock);
                 lastIndex = group[m];
               }
               edgeMatrix.insert_cols(edgeMatrix.n_cols,eachCol);
            }
//            edgeMatrix.print("edgeMatrix");
//            cout<<"numVector"<<endl;
            // for(int m=0;m<numVector.size();m++)
            //   cout<<numVector[m]<<"\t";
            //   cout<<endl;
            compressEdgeVector.clear();

            //Then for each group, do edge Multiplication
            umat resM(edgeMatrix.col(edgeMatrix.n_cols-1));
            edgeMatrix.shed_cols(edgeMatrix.n_cols-1,edgeMatrix.n_cols-1);

            mat coeff(1,1,arma::fill::ones);
            umat nresM;
            mat ncoeff;
            urowvec q1;
            double expNum=0;

            lastIndex = 0;
            for(int m=0;m<compressEdgeGroup.size() && resM.n_cols>0;m++){
              for(int k=lastIndex;k<compressEdgeGroup[m];k++){
                ncoeff=coeff;
                ncoeff*=(1-edgeProb[k]);
                coeff*=edgeProb[k];

                nresM.resize(resM.n_rows,resM.n_cols);
                for(int ii=0;ii<resM.n_rows;ii++)
                   for(int jj=0;jj<resM.n_cols;jj++)
                      nresM(ii,jj) = resM(ii,jj) & edgeMatrix(ii,k);

               q1 = any(nresM, 0);
               q0 = find(q1==1);
               nresM = nresM.cols(q0);
               ncoeff = ncoeff.cols(q0);

                resM.insert_cols(resM.n_cols,nresM);
                coeff.insert_cols(coeff.n_cols,ncoeff);

//                cout<<"multiply "<<k<<endl;
//                resM.print("resM");
              }

              if(resM.n_cols==0)
                 continue;

              umat subMa=resM.head_rows(numVector[m]);
              edgeMatrix.shed_rows(0,numVector[m]-1);
              resM.shed_rows(0,numVector[m]-1);
              q1= any(subMa, 0);
              q0=find(q1==0);
              uvec q2=find(q1==1);

              int numBits;
              if(m==0)
                 numBits = group[0];
              else numBits = group[m]-group[m-1];


              brow b(1,q2.n_elem);
              b.zeros();
              for(int jj=0;jj<q2.n_elem;jj++)
                 for(int ii=0;ii<subMa.n_rows;ii++){
                   dynamic_bitset<> tdb(numBits, subMa(ii,q2(jj)));
//                   cout<<"tdb: "<<tdb<<endl;
                   b(jj) += tdb.count();
                 }
            double avg=sum(coeff.cols(q2) % b);
             expNum += avg;
  //           cout<<"avg: "<<avg<<endl;

             resM=resM.cols(q0);
             coeff=coeff.cols(q0);

              lastIndex=compressEdgeGroup[m];
            }
          expectPath(i,j)=expNum;


        }

      }

      cout<<nodeNameMap[i]<<"\t"<<nodeNameMap[j]<<": "<<expectPath(i,j)<<endl;

    }

    pq=arma::find(expectPath>0);
    nonzeroExpectPath = expectPath.elem(pq);

  }

  void Graph::findAllSimplePath(vector<dynamic_bitset<> >& paths, int source, int end){
    dynamic_bitset<> path(edgeSize);
    dynamic_bitset<> visited(nodeSize);
    dfs(paths, path, visited, source, end);
  }

  void Graph::dfs(vector<dynamic_bitset<> >& paths, dynamic_bitset<>& path, dynamic_bitset<>& visited, int curr, int end){
    visited[curr]=1;
    dynamic_bitset<>::size_type it;
    it=adjacencyList[curr].find_first();
    while(it!=dynamic_bitset<>::npos){
      int edgeIndex = adjacencyMatrix[curr][it];
      if(it==end){
        path[edgeIndex]=1;
        paths.push_back(path);
        path[edgeIndex]=0;
      }
      else if(visited[it]==0){
        path[edgeIndex]=1;
        dfs(paths,path,visited,it,end);
        path[edgeIndex]=0;
      }
      it=adjacencyList[curr].find_next(it);
    }
    visited[curr]=0;
  }

  double Graph::countBinaryShortestPaths(){
    dynamic_bitset<>::size_type it;
    mat deterPath;
    deterPath.zeros(nodeSize,nodeSize);
    //for every pair, find all shortest path
    for(int i=0;i<nodeSize;i++){
      vector<TreeNode*> allNodes(nodeSize);

      vector<vector<TreeNode*> > tree;
      int index=0; //the level of the tree
      //create the root node
      TreeNode* root=new TreeNode(i);
      root->weight=1;
      vector<TreeNode*> temp;
      temp.push_back(root);
      tree.push_back(temp);
      allNodes[i]=root;

      while(index<tree.size()){
        vector<TreeNode*> a; //new generated level
        for(int k=0;k<tree[index].size();k++){
          TreeNode* n=tree[index][k];
          int node=n->id;
          it=adjacencyList[node].find_first();
          while(it!=dynamic_bitset<>::npos){
            if(allNodes[it]==NULL){
              TreeNode* newNode=new TreeNode(it);
              newNode->distance=n->distance+1;
              newNode->weight=n->weight;
              a.push_back(newNode);
              allNodes[it]=newNode;
            }
            else if(allNodes[it]->distance==n->distance+1){
              allNodes[it]->weight += n->weight;
            }
            it=adjacencyList[node].find_next(it);
          }
        }
        if(a.size()>0)
        tree.push_back(a);
        index++;
      }

      for(int k=0;k<nodeSize;k++){
        if(allNodes[k]!=NULL)
        {
          if(i!=k)
          deterPath(i,k)=allNodes[k]->weight;
          delete allNodes[k];
        }
      }
      allNodes.clear();
      tree.clear();

    }

    nonzeroDeterPath = deterPath.elem(pq);
    double averError = accu((abs(nonzeroDeterPath- nonzeroExpectPath)) / nonzeroExpectPath);
    averError /= pq.n_elem;
    return averError;
  }

  double Graph::countThresholdShortestPaths(double threshold){
    mat deterPath;
    deterPath.zeros(nodeSize,nodeSize);
    vector<dynamic_bitset<> > al(nodeSize, dynamic_bitset<>(nodeSize));
    for(int i=0; i< edgeProbabilities.size();i++){
      if(edgeProbabilities[i]>threshold){
        int source = edgeToNodeMap[i].first;
        int dest = edgeToNodeMap[i].second;
        if(adjacencyMatrix[source][dest]==i)
        al[source][dest]=1;
        if(adjacencyMatrix[dest][source]==i)
        al[dest][source]=1;
      }
    }

    dynamic_bitset<>::size_type it;
    deterPath.zeros();
    //for every pair, find all shortest path
    for(int i=0;i<nodeSize;i++){
      vector<TreeNode*> allNodes(nodeSize);

      vector<vector<TreeNode*> > tree;
      int index=0; //the level of the tree
      //create the root node
      TreeNode* root=new TreeNode(i);
      root->weight=1;
      vector<TreeNode*> temp;
      temp.push_back(root);
      tree.push_back(temp);
      allNodes[i]=root;

      while(index<tree.size()){
        vector<TreeNode*> a; //new generated level
        for(int k=0;k<tree[index].size();k++){
          TreeNode* n=tree[index][k];
          int node=n->id;
          it=al[node].find_first();
          while(it!=dynamic_bitset<>::npos){
            if(allNodes[it]==NULL){
              TreeNode* newNode=new TreeNode(it);
              newNode->distance=n->distance+1;
              newNode->weight=n->weight;
              a.push_back(newNode);
              allNodes[it]=newNode;
            }
            else if(allNodes[it]->distance==n->distance+1){
              allNodes[it]->weight += n->weight;
            }
            it=al[node].find_next(it);
          }
        }
        if(a.size()>0)
        tree.push_back(a);
        index++;
      }

      for(int k=0;k<nodeSize;k++){
        if(allNodes[k]!=NULL)
        {

          if(i!=k)
          deterPath(i,k)=allNodes[k]->weight;
          delete allNodes[k];
        }
      }
      allNodes.clear();
      tree.clear();

    }
    nonzeroDeterPath = deterPath.elem(pq);
    double averError = accu((abs(nonzeroDeterPath- nonzeroExpectPath)) / nonzeroExpectPath);
    averError /= pq.n_elem;
    return averError;
  }

  double Graph::countSampleShortestPaths(){
    mat deterPath;
    deterPath.zeros(nodeSize,nodeSize);
    std::random_device r;
    std::default_random_engine generator(r());
    vector<dynamic_bitset<> > al(nodeSize, dynamic_bitset<>(nodeSize));
    for(int i=0; i< edgeProbabilities.size();i++){
      std::bernoulli_distribution distribution(edgeProbabilities[i]);
      bool isEx=distribution(generator);
      if(isEx){
        int source = edgeToNodeMap[i].first;
        int dest = edgeToNodeMap[i].second;
        if(adjacencyMatrix[source][dest]==i)
        al[source][dest]=1;
        if(adjacencyMatrix[dest][source]==i)
        al[dest][source]=1;
      }
    }

    dynamic_bitset<>::size_type it;
    deterPath.zeros();
    //for every pair, find all shortest path
    for(int i=0;i<nodeSize;i++){
      vector<TreeNode*> allNodes(nodeSize);

      vector<vector<TreeNode*> > tree;
      int index=0; //the level of the tree
      //create the root node
      TreeNode* root=new TreeNode(i);
      root->weight=1;
      vector<TreeNode*> temp;
      temp.push_back(root);
      tree.push_back(temp);
      allNodes[i]=root;

      while(index<tree.size()){
        vector<TreeNode*> a; //new generated level
        for(int k=0;k<tree[index].size();k++){
          TreeNode* n=tree[index][k];
          int node=n->id;
          it=al[node].find_first();
          while(it!=dynamic_bitset<>::npos){
            if(allNodes[it]==NULL){
              TreeNode* newNode=new TreeNode(it);
              newNode->distance=n->distance+1;
              newNode->weight=n->weight;
              a.push_back(newNode);
              allNodes[it]=newNode;
            }
            else if(allNodes[it]->distance==n->distance+1){
              allNodes[it]->weight += n->weight;
            }
            it=al[node].find_next(it);
          }
        }
        if(a.size()>0)
        tree.push_back(a);
        index++;
      }

      for(int k=0;k<nodeSize;k++){
        if(allNodes[k]!=NULL)
        {
          if(i!=k)
              deterPath(i,k)=allNodes[k]->weight;
          delete allNodes[k];
        }
      }
      allNodes.clear();
      tree.clear();

    }

    nonzeroDeterPath = deterPath.elem(pq);
    double averError = accu((abs(nonzeroDeterPath- nonzeroExpectPath)) / nonzeroExpectPath);
    averError /= pq.n_elem;
    return averError;
  }
