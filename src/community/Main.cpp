#include "Graph.h"
#include <cstdlib>
#include <cstring>
int main(int argc, char *argv[]){
  int numberNodes = atoi(argv[1]);
	char fileName[100];
	strcpy(fileName,argv[2]);
	Graph g(numberNodes,fileName);
  g.findCommunity();  
}
