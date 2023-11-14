scale_status: peek_adaptive_with_status_array_kbound_thread_32.exe peek_adaptive_with_status_array_kbound_thread_16.exe peek_adaptive_with_status_array_kbound_thread_8.exe peek_adaptive_with_status_array_kbound_thread_4.exe peek_adaptive_with_status_array_kbound_thread_2.exe peek_adaptive_with_status_array_kbound_thread_1.exe

scale_swap: peek_adaptive_with_edge_swap_kbound_thread_32.exe peek_adaptive_with_edge_swap_kbound_thread_16.exe peek_adaptive_with_edge_swap_kbound_thread_8.exe peek_adaptive_with_edge_swap_kbound_thread_4.exe peek_adaptive_with_edge_swap_kbound_thread_2.exe peek_adaptive_with_edge_swap_kbound_thread_1.exe

runSTPruning.exe: ./experiments/runSTPruning.cpp
	g++ -O3 -I./includes -fopenmp -o runSTPruning.exe ./experiments/runSTPruning.cpp

runSTPruningWithKBound.exe: ./experiments/runSTPruning.cpp
	g++ -O3 -I./includes -DTECHNIQUE_BOUND -fopenmp -o runSTPruningWithKBound.exe ./experiments/runSTPruning.cpp

peek_adaptive_with_status_array_kbound_thread_32.exe: ./experiments/PeeKAdaptiveWithStatusArrayPerf.cpp
	g++ -O3 -I./includes -DTECHNIQUE_BOUND -DKSP_NUM_THREADS=32 -DKSP_PARALLEL_DEVIATIONS_L2 -DCOMPACT_KSP_NUM_THREADS=1 -DADAPTIVE_GRAPH_COMPACT_VERTEX -DADAPTIVE_GRAPH_COMPACT_EDGE -fopenmp -o peek_adaptive_with_status_array_kbound_thread_32.exe ./experiments/PeeKAdaptiveWithStatusArrayPerf.cpp

peek_adaptive_with_status_array_kbound_thread_16.exe: ./experiments/PeeKAdaptiveWithStatusArrayPerf.cpp
	g++ -O3 -I./includes -DTECHNIQUE_BOUND -DKSP_NUM_THREADS=16 -DKSP_PARALLEL_DEVIATIONS_L2 -DCOMPACT_KSP_NUM_THREADS=1 -DADAPTIVE_GRAPH_COMPACT_VERTEX -DADAPTIVE_GRAPH_COMPACT_EDGE -fopenmp -o peek_adaptive_with_status_array_kbound_thread_16.exe ./experiments/PeeKAdaptiveWithStatusArrayPerf.cpp

peek_adaptive_with_status_array_kbound_thread_8.exe: ./experiments/PeeKAdaptiveWithStatusArrayPerf.cpp
	g++ -O3 -I./includes -DTECHNIQUE_BOUND -DKSP_NUM_THREADS=8 -DKSP_PARALLEL_DEVIATIONS_L2 -DCOMPACT_KSP_NUM_THREADS=1 -DADAPTIVE_GRAPH_COMPACT_VERTEX -DADAPTIVE_GRAPH_COMPACT_EDGE -fopenmp -o peek_adaptive_with_status_array_kbound_thread_8.exe ./experiments/PeeKAdaptiveWithStatusArrayPerf.cpp

peek_adaptive_with_status_array_kbound_thread_4.exe: ./experiments/PeeKAdaptiveWithStatusArrayPerf.cpp
	g++ -O3 -I./includes -DTECHNIQUE_BOUND -DKSP_NUM_THREADS=4 -DKSP_PARALLEL_DEVIATIONS_L2 -DCOMPACT_KSP_NUM_THREADS=1 -DADAPTIVE_GRAPH_COMPACT_VERTEX -DADAPTIVE_GRAPH_COMPACT_EDGE -fopenmp -o peek_adaptive_with_status_array_kbound_thread_4.exe ./experiments/PeeKAdaptiveWithStatusArrayPerf.cpp

peek_adaptive_with_status_array_kbound_thread_2.exe: ./experiments/PeeKAdaptiveWithStatusArrayPerf.cpp
	g++ -O3 -I./includes -DTECHNIQUE_BOUND -DKSP_NUM_THREADS=2 -DKSP_PARALLEL_DEVIATIONS_L2 -DCOMPACT_KSP_NUM_THREADS=1 -DADAPTIVE_GRAPH_COMPACT_VERTEX -DADAPTIVE_GRAPH_COMPACT_EDGE -fopenmp -o peek_adaptive_with_status_array_kbound_thread_2.exe ./experiments/PeeKAdaptiveWithStatusArrayPerf.cpp

peek_adaptive_with_status_array_kbound_thread_1.exe: ./experiments/PeeKAdaptiveWithStatusArrayPerf.cpp
	g++ -O3 -I./includes -DTECHNIQUE_BOUND -DKSP_NUM_THREADS=1 -DKSP_PARALLEL_DEVIATIONS_L2 -DCOMPACT_KSP_NUM_THREADS=1 -DADAPTIVE_GRAPH_COMPACT_VERTEX -DADAPTIVE_GRAPH_COMPACT_EDGE -fopenmp -o peek_adaptive_with_status_array_kbound_thread_1.exe ./experiments/PeeKAdaptiveWithStatusArrayPerf.cpp

peek_adaptive_with_edge_swap_kbound_thread_32.exe: ./experiments/PeeKAdaptiveWithEdgeSwapPerf.cpp
	g++ -O3 -I./includes -DTECHNIQUE_BOUND -DKSP_NUM_THREADS=32 -DKSP_PARALLEL_DEVIATIONS_L2 -DCOMPACT_KSP_NUM_THREADS=1 -fopenmp -o peek_adaptive_with_edge_swap_kbound_thread_32.exe ./experiments/PeeKAdaptiveWithEdgeSwapPerf.cpp

peek_adaptive_with_edge_swap_kbound_thread_16.exe: ./experiments/PeeKAdaptiveWithEdgeSwapPerf.cpp
	g++ -O3 -I./includes -DTECHNIQUE_BOUND -DKSP_NUM_THREADS=16 -DKSP_PARALLEL_DEVIATIONS_L2 -DCOMPACT_KSP_NUM_THREADS=1 -fopenmp -o peek_adaptive_with_edge_swap_kbound_thread_16.exe ./experiments/PeeKAdaptiveWithEdgeSwapPerf.cpp

peek_adaptive_with_edge_swap_kbound_thread_8.exe: ./experiments/PeeKAdaptiveWithEdgeSwapPerf.cpp
	g++ -O3 -I./includes -DTECHNIQUE_BOUND -DKSP_NUM_THREADS=8 -DKSP_PARALLEL_DEVIATIONS_L2 -DCOMPACT_KSP_NUM_THREADS=1 -fopenmp -o peek_adaptive_with_edge_swap_kbound_thread_8.exe ./experiments/PeeKAdaptiveWithEdgeSwapPerf.cpp

peek_adaptive_with_edge_swap_kbound_thread_4.exe: ./experiments/PeeKAdaptiveWithEdgeSwapPerf.cpp
	g++ -O3 -I./includes -DTECHNIQUE_BOUND -DKSP_NUM_THREADS=4 -DKSP_PARALLEL_DEVIATIONS_L2 -DCOMPACT_KSP_NUM_THREADS=1 -fopenmp -o peek_adaptive_with_edge_swap_kbound_thread_4.exe ./experiments/PeeKAdaptiveWithEdgeSwapPerf.cpp

peek_adaptive_with_edge_swap_kbound_thread_2.exe: ./experiments/PeeKAdaptiveWithEdgeSwapPerf.cpp
	g++ -O3 -I./includes -DTECHNIQUE_BOUND -DKSP_NUM_THREADS=2 -DKSP_PARALLEL_DEVIATIONS_L2 -DCOMPACT_KSP_NUM_THREADS=1 -fopenmp -o peek_adaptive_with_edge_swap_kbound_thread_2.exe ./experiments/PeeKAdaptiveWithEdgeSwapPerf.cpp

peek_adaptive_with_edge_swap_kbound_thread_1.exe: ./experiments/PeeKAdaptiveWithEdgeSwapPerf.cpp
	g++ -O3 -I./includes -DTECHNIQUE_BOUND -DKSP_NUM_THREADS=1 -DKSP_PARALLEL_DEVIATIONS_L2 -DCOMPACT_KSP_NUM_THREADS=1 -fopenmp -o peek_adaptive_with_edge_swap_kbound_thread_1.exe ./experiments/PeeKAdaptiveWithEdgeSwapPerf.cpp

runYen.exe: ./experiments/runYen.cpp
	g++ -O3 -I./includes -fopenmp -o runYen.exe ./experiments/runYen.cpp

runYenWithKBound.exe: ./experiments/runYen.cpp
	g++ -O3 -I./includes -DTECHNIQUE_BOUND -fopenmp -o runYenWithKBound.exe ./experiments/runYen.cpp

stpruningStatistics.exe: ./experiments/stpruningStatistics.cpp
	g++ -O3 -I./includes -fopenmp -o stpruningStatistics.exe ./experiments/stpruningStatistics.cpp

adaptive_graph_compact_vertex_edge.exe: ./experiments/adaptiveGraphCompact.cpp
	g++ -O3 -I./includes -DTECHNIQUE_BOUND -DKSP_NUM_THREADS=4 -DKSP_PARALLEL_DEVIATIONS_L2 -DCOMPACT_KSP_NUM_THREADS=4 -DADAPTIVE_GRAPH_COMPACT_VERTEX -DADAPTIVE_GRAPH_COMPACT_EDGE -fopenmp -o adaptive_graph_compact_vertex_edge.exe ./experiments/adaptiveGraphCompact.cpp

adaptive_graph_compact_vertex.exe: ./experiments/adaptiveGraphCompact.cpp
	g++ -O3 -I./includes -DTECHNIQUE_BOUND -DKSP_NUM_THREADS=4 -DKSP_PARALLEL_DEVIATIONS_L2 -DCOMPACT_KSP_NUM_THREADS=4 -DADAPTIVE_GRAPH_COMPACT_VERTEX -fopenmp -o adaptive_graph_compact_vertex.exe ./experiments/adaptiveGraphCompact.cpp

adaptive_graph_compact_edge.exe: ./experiments/adaptiveGraphCompact.cpp
	g++ -O3 -I./includes -DTECHNIQUE_BOUND -DKSP_NUM_THREADS=4 -DKSP_PARALLEL_DEVIATIONS_L2 -DCOMPACT_KSP_NUM_THREADS=4 -DADAPTIVE_GRAPH_COMPACT_EDGE -fopenmp -o adaptive_graph_compact_edge.exe ./experiments/adaptiveGraphCompact.cpp

clean:
	rm -rf *.o
	rm -rf *.exe


.PHONY: clean all