#include "Graph.h"
#include <cstdlib>
#include <cstring>
int main(int argc, char *argv[]){
    int numberNodes = atoi(argv[1]);
	  char fileName[100];

	  strcpy(fileName,argv[2]);
	  Graph g(numberNodes,fileName);
    clock_t begin_time=clock();
    cout<<"----------Print the expected number of paths between pairs of nodes--------"<<endl;
    g.countShortestPaths();
    //The following compute the error rate of binary, threshold, sampling methods
    /*
    cout <<float( clock () - begin_time ) /  CLOCKS_PER_SEC<<"\t";
    double error=0;
    begin_time=clock();
    error = g.countBinaryShortestPaths();
    cout<<error<<"\t"<<float( clock () - begin_time ) /  CLOCKS_PER_SEC<<"\t";
    begin_time=clock();
    error = g.countThresholdShortestPaths(0.6);
    cout<<error<<"\t"<<float( clock () - begin_time ) /  CLOCKS_PER_SEC<<"\t";

    begin_time=clock();
    double sampleError = 0;
    vec paths;
    paths.zeros(g.pq.n_elem);
    for(int ii=1;ii<=1000;ii++){
      g.countSampleShortestPaths();
      paths = paths + g.nonzeroDeterPath;
    }
    paths=paths/1000;
    sampleError = accu((abs(paths- g.nonzeroExpectPath)) / g.nonzeroExpectPath);
    sampleError /= g.pq.n_elem;
    cout<<sampleError<<"\t"<<float( clock () - begin_time ) /  CLOCKS_PER_SEC<<endl;
    */
	return 0;

}
