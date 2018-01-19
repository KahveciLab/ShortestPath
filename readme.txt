1. Before running our method, you need to install boost C++ library and armadillo software.
2. To complie our method: make. 
   If you have not installed the separately-compiled Boost libraries, you need include the boost directory in make file.
   ex: g++ --std=c++11 -o a.out *.cpp *.h -larmadillo -I  path/to/boost_1_65_1
3. To run our method:
   input parameters: number of nodes, file path
   ex: ./path 100 ~/dataset



