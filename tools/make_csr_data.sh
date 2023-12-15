#!/bin/bash

g++ -O3 -o make_csr.exe  make_csr.cpp

g++ -O3 -I../includes -fopenmp -o make_reachable_src_dest.exe ./make_reachable_src_dest.cpp

# use PaRMAT to generate edge list based graph
cd PaRMAT/Release/
make clean
make
./PaRMAT -nVertices 100 -nEdges 5000  -a 0.45 -b 0.22 -c 0.22 -sorted -noEdgeToSelf -noDuplicateEdges -output edge_list.graph

#change edge list into csr
../../make_csr.exe edge_list.graph begin.bin adj.bin value.bin

#generate random reachable pairs of src and dest
../../make_reachable_src_dest.exe begin.bin adj.bin value.bin > reachable_src_dest.txt
