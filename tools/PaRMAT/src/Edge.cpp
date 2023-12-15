#include "Edge.hpp"
#include <random>
#include <chrono>

Edge::Edge(EdgeIndexType rec_src, EdgeIndexType rec_dst):
		src(rec_src), dst(rec_dst)
{}

Edge::Edge(const Edge &cSource) {
	src = cSource.src;
	dst = cSource.dst;
}

Edge& Edge::operator= (const Edge &cSource) {
	src = cSource.src;
	dst = cSource.dst;
	return *this;
}

bool Edge::selfEdge(){
	return ( src == dst );
}

bool operator< (const Edge& cR1, const Edge& cR2) {
	if( cR1.src < cR2.src )
		return true;
	else if( cR1.src > cR2.src )
		return false;
	else if( cR1.dst < cR2.dst )
		return true;
	else
		return false;
}

bool operator== (const Edge& cR1, const Edge& cR2) {
	return ( cR1.dst == cR2.dst && cR1.src == cR2.src);
}

std::ostream& operator<< (std::ostream &out, Edge &cEdge) {
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine gen(seed);
	std::normal_distribution<double> n_dis(0.5,0.2);
	double randomValue;
	do {
		randomValue = n_dis(gen);
	} while (randomValue < 0.0 || randomValue > 1.0);
	out << cEdge.src << "\t" << cEdge.dst << "\t" << randomValue << "\n";
	return out;
}
