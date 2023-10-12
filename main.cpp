#include <ctime>
#include "indConst.h"
#include "bigraph.h"
#include "parallelConst.h"
using namespace std;



int main(int argc, char **argv) {
	if (argc == 1) {
		cout << "error in number of arguments" << endl;
	}
	string exec_type = argv[1];
	if (exec_type == "-skyIndDelta") {
		cout << "start skyIndDelta for " << argv[2] << endl;
		unordered_map<string, vector<vector<int>>> skyLable;
		int delta = 0;
		BiGraph g(argv[2]);
		auto start = chrono::system_clock::now();
		skyIndexDelta(g, skyLable, delta);
		auto end = chrono::system_clock::now();
		chrono::duration<double> time = end - start;
		cout << "skyIndDelta running time: " << time.count() << endl;
	}
    else if (exec_type == "-skyLv1Query") {
		cout << "start skyLv1Query for " << argv[2] << endl;
		
		BiGraph g(argv[2]);
		unordered_map<string, vector<vector<int>>> skyLable;
		unordered_map<vector<int>, vector<string>, VectorHasher> skyNode;
		unordered_map<vector<int>, skyline_index_hub*, VectorHasher> skyline_index_lv0;
		unordered_map<vector<int>, skyline_index_ccblock*, VectorHasher> skyline_index_lv1;
		unordered_map<string, vector<skyline_index_ccblock*>> skyline_index_v2ccb;
		unordered_map<string, vector<skyline_index_hub*>> skyline_index_v2h;
		unordered_map<string, skyline_index_ccblock*> skyline_index_v2par;
		int delta = 0;
		string v_q = argv[5];
		vector<bool> left; vector<bool> right;
		// all the vertices in query result are set as true
		left.resize(g.num_v1, false); right.resize(g.num_v2, false);

		skyIndexDelta(g, skyLable, delta);
		cout << "delta: " << delta << endl;
		build_skyline_index_lv0(g, skyLable, skyNode, skyline_index_lv0, skyline_index_v2h, delta);
		build_skyline_index_lv1(g, skyNode, skyline_index_lv0, skyline_index_v2h, skyline_index_v2ccb, skyline_index_v2par, skyline_index_lv1, delta);
		retrieve_skyline_index_lv1(g, skyline_index_v2ccb, skyline_index_lv1, left, right, stoi(argv[3]), stoi(argv[4]), v_q);
		output_skyline_index(argv[2], skyline_index_v2ccb, skyline_index_lv1, skyline_index_lv0, skyline_index_v2h);
		
	}
	else if (exec_type == "-skyLv1Par") {
		cout << "start skyLv1Parallel for " << argv[2] << endl;
		
		BiGraph g(argv[2]);
		unordered_map<string, vector<vector<int>>> skyLable;
		unordered_map<vector<int>, vector<string>, VectorHasher> skyNode;
		unordered_map<vector<int>, skyline_index_hub*, VectorHasher> skyline_index_lv0;
		unordered_map<vector<int>, skyline_index_ccblock*, VectorHasher> skyline_index_lv1;
		unordered_map<string, vector<skyline_index_ccblock*>> skyline_index_v2ccb;
		unordered_map<string, vector<skyline_index_hub*>> skyline_index_v2h;
		
		int delta = 0;
		int threads = stoi(argv[3]);
		vector<bool> left; vector<bool> right;
		// all the vertices in query result are set as true
		left.resize(g.num_v1, false); right.resize(g.num_v2, false);

		skyIndexDelta(g, skyLable, delta);
		cout << "delta: " << delta << endl;
		build_skyline_index_lv0(g, skyLable, skyNode, skyline_index_lv0, skyline_index_v2h, delta);
		parallel_build_skyline_index_lv1(g, skyNode, skyline_index_lv0, skyline_index_v2h, skyline_index_v2ccb, skyline_index_lv1, threads);

	}	

	else {
		cout << "illegal arguments" << endl;
	}
	return 0;
}