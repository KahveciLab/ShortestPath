CC=g++
ST=--std=c++11
LOP=-o

MAIN=*
TAG=shortestPath
LIB=-larmadillo

all : clean code1

code1 :
	$(CC) $(ST) $(LOP) $(TAG) $(MAIN).cpp $(MAIN).h $(LIB)

clean :
	rm -f $(TAG)
