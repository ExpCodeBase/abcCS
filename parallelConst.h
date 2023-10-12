#pragma once
#ifndef MULTITHREAD
#define MULTITHREAD
#include "bigraph.h"
#include "indConst.h"


extern void parallel_build_skyline_index_lv1(BiGraph& g, 
	std::unordered_map<std::vector<int>, std::vector<std::string>, VectorHasher>& skyNode, 
	std::unordered_map<std::vector<int>, skyline_index_hub*, VectorHasher>& skyline_index_lv0,
	std::unordered_map<std::string, std::vector<skyline_index_hub*>>& skyline_index_v2h, 
    std::unordered_map<std::string, std::vector<skyline_index_ccblock*>>& skyline_index_v2ccb, 
	std::unordered_map<std::vector<int>, skyline_index_ccblock*, VectorHasher>& skyline_index_lv1,
    int threads);


#endif // !MULTITHREAD
