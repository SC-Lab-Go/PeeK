
---
Software requirement
-----
g++ (8.5.0 tested)

OpenMP (4.5 tested)

Boost C++ Libraries (1.81.0 tested it's only for `PeeKAdaptiveWithEdgeSwapParallel.h`)

---

Installation and usage
------
```bash
To generate CSR format graph data and its reachable pairs of source and target

cd tools/
./make_csr_data.sh

The data is located in `PaRMAT/Release/begin.bin`, `PaRMAT/Release/adj.bin` and `PaRMAT/Release/value.bin`
The reachable pairs of source and target are located in `PaRMAT/Release/reachable_src_dest.txt`


To run PeeK
cd ..
make peek_adaptive_with_edge_swap_kbound_thread_32.exe
./peek_adaptive_with_edge_swap_kbound_thread_32.exe tools/PaRMAT/Release/begin.bin tools/PaRMAT/Release/adj.bin tools/PaRMAT/Release/value.bin tools/PaRMAT/Release/reachable_src_dest.txt


```

----
Reference
-------

If you use PeeK in your project, please cite the following paper.

```python
@inproceedings{feng2023peek,
  title={PeeK: A Prune-Centric Approach for K Shortest Path Computation},
  author={Feng, Wang and Chen, Shiyang and Liu, Hang and Ji, Yuede},
  booktitle={Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis},
  pages={1--14},
  year={2023}
}
```

----
Compared work
-------
This is the open version of OptYen.
https://github.com/AbcAeffchen/k-shortest-path