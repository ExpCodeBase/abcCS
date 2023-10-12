#include "parallelConst.h"

using namespace std;

vector<vector<atomic_int>> atomic_index;

void multithread() {

}

void build_path_left(vector<int> build_start, 
	BiGraph& g,
	unordered_map<vector<int>, skyline_index_hub*, VectorHasher>& skyline_index_lv0,
	unordered_map<string, vector<skyline_index_hub*>>& skyline_index_v2h,
    unordered_map<string, vector<skyline_index_ccblock*>>& skyline_index_v2ccb, 
	unordered_map<vector<int>, skyline_index_ccblock*, VectorHasher>& skyline_index_lv1){

    unordered_map<string, skyline_index_ccblock*> skyline_index_v2par;

	vector<bool> left_pushed; vector<bool> right_pushed;
	vector<bool> left_pushed_glb; vector<bool> right_pushed_glb;

	left_pushed.resize(g.num_v1, false); 
	right_pushed.resize(g.num_v2, false);
	left_pushed_glb.resize(g.num_v1, false); 
	right_pushed_glb.resize(g.num_v2, false);

	vector<int> tar;
	skyline_index_hub* cur_bd = skyline_index_lv0[build_start];
	skyline_index_hub* cur_lk = skyline_index_lv0[build_start];

	vector<vector<int>> vt_tbre;
	vector<vector<int>> vt_ngp;

	// build based on hub

	// global state labeling
	skyline_index_hub* cur_st = cur_bd;
	while(cur_st){
		for(auto vt: cur_st->nodeset){
			if(vt[0]=='u'){
				vt.erase(vt.begin());
				int vtid = stoi(vt)-1;
				if(!left_pushed_glb[vtid]) left_pushed_glb[vtid] = true;
			} else {
				vt.erase(vt.begin());
				int vtid = stoi(vt)-1;
				if(!right_pushed_glb[vtid]) right_pushed_glb[vtid] = true;
			}
		}
		cur_st = cur_st->right_c;	
	}
	
	// building
	while(cur_bd){

		fill_n(left_pushed.begin(), left_pushed.size(), false);
		fill_n(right_pushed.begin(), right_pushed.size(), false);

		// delete edges of blocks above
		
		if(cur_bd->left_p){
			for(auto vt: cur_bd->left_p->nodeset){
				if(vt[0]=='u'){
					vt.erase(vt.begin());
					int vtid = stoi(vt)-1;
					for(auto e: g.getV1Neighbors(vtid)){
						vector<int> tmp;
						if(g.isEdge(vtid, e)){
							g.deleteEdge(vtid,e);
							tmp = {vtid, e};
							vt_tbre.push_back(tmp);
						}
					}
				} else {
					vt.erase(vt.begin());
					int vtid = stoi(vt)-1;
					for(auto e: g.getV2Neighbors(vtid)){
						vector<int> tmp;
						if(g.isEdge(e, vtid)){
							g.deleteEdge(e,vtid);
							tmp = {e, vtid};
							vt_tbre.push_back(tmp);
						}
					}
				}
			}
		}
		
		tar = cur_bd->hubid;

		if(tar.size()<3) {
			tar.push_back(0);
		} else {
			tar[2] = 0;
		}


		// state labeling of vertices in nodeset
		for(auto vt: cur_bd->nodeset){
			if(vt[0]=='u'){
				vt.erase(vt.begin());
				int vtid = stoi(vt)-1;
				if(!left_pushed[vtid]) left_pushed[vtid] = true;
			} else {
				vt.erase(vt.begin());
				int vtid = stoi(vt)-1;
				if(!right_pushed[vtid]) right_pushed[vtid] = true;
			}
		}
		

		queue<string> ccQue;
		string now;
		
		// ccQue for BFS to build each tar with cc_id
		for(auto vt: cur_bd->nodeset){
			
			if(vt[0]=='u'){
				vt.erase(vt.begin());
				int vtid = stoi(vt)-1;
				if(left_pushed[vtid]){
					vt = "u" + vt;
					ccQue.push(vt);
				}
			} else {
				vt.erase(vt.begin());
				int vtid = stoi(vt)-1;
				if(right_pushed[vtid]){
					vt = "v" + vt;
					ccQue.push(vt);
				}
			}

			// stopFlag for ccid iteration
			bool stopFlag = true;
			
			
			while(!ccQue.empty()){

				now = ccQue.front();
				ccQue.pop();

				if(!skyline_index_lv1[tar]) skyline_index_lv1[tar] = new skyline_index_ccblock;

				
				skyline_index_lv1[tar]->ccid = tar[2];			
				skyline_index_lv1[tar]->hubid = cur_bd->hubid;
				skyline_index_lv0[cur_bd->hubid]->blockset.push_back(skyline_index_lv1[tar]);
		

				int vtid = 0;
				
				if(now[0]=='u'){
					now.erase(now.begin());
					vtid = stoi(now)-1;
					string nowPu = "u" + now;
					if(left_pushed[vtid]){
						left_pushed[vtid] = false;
						left_pushed_glb[vtid] = false;
						skyline_index_v2ccb[nowPu].push_back(skyline_index_lv1[tar]);
						skyline_index_lv1[tar]->nodeset.push_back(nowPu);
					}
					for(auto e: g.getV1Neighbors(vtid)){
						if(right_pushed[e]){
							string nowNei = "v" + to_string(e+1);
							skyline_index_lv1[tar]->nodeset.push_back(nowNei);
							skyline_index_v2ccb[nowNei].push_back(skyline_index_lv1[tar]);
							right_pushed[e] = false;
							right_pushed_glb[e] = false;
							ccQue.push(nowNei);
							stopFlag = false;
						}
						else if(right_pushed_glb[e]){
							right_pushed_glb[e] = false;
							string nowNei = "v" + to_string(e+1);
							ccQue.push(nowNei);
							for(auto e: skyline_index_v2h[nowNei]){
								if(e->hubid == cur_bd->right_c->hubid){
									skyline_index_v2par[nowNei] = skyline_index_lv1[tar];
									break;
								}
							}
							
							stopFlag = false;
						}
					}
				} else {
					now.erase(now.begin());
					vtid = stoi(now)-1;
					string nowPu = "v" + now;
					if(right_pushed[vtid]){
						right_pushed[vtid] = false;
						right_pushed_glb[vtid] = false;
						skyline_index_v2ccb[nowPu].push_back(skyline_index_lv1[tar]);
						skyline_index_lv1[tar]->nodeset.push_back(nowPu);
					}
					for(auto e: g.getV2Neighbors(vtid)){
						if(left_pushed[e]){
							string nowNei = "u" + to_string(e+1);
							skyline_index_lv1[tar]->nodeset.push_back(nowNei);
							skyline_index_v2ccb[nowNei].push_back(skyline_index_lv1[tar]);
							left_pushed[e] = false;
							left_pushed_glb[e] = false;
							ccQue.push(nowNei);
							stopFlag = false;
						}
						else if(left_pushed_glb[e]){
							left_pushed_glb[e] = false;
							string nowNei = "u" + to_string(e+1);
							ccQue.push(nowNei);

							for(auto e: skyline_index_v2h[nowNei]){
								if(e->hubid == cur_bd->right_c->hubid){
									skyline_index_v2par[nowNei] = skyline_index_lv1[tar];
									break;
								}
							}
							stopFlag = false;
						}
					}
					
				}
			}

			if(!stopFlag){
				tar[2]++;
			}


		}


		// global state labeling
		skyline_index_hub* cur_st1 = cur_bd->right_c;
		while(cur_st1){
			for(auto vt: cur_st1->nodeset){
				if(vt[0]=='u'){
					vt.erase(vt.begin());
					int vtid = stoi(vt)-1;
					if(!left_pushed_glb[vtid]) left_pushed_glb[vtid] = true;
				} else {
					vt.erase(vt.begin());
					int vtid = stoi(vt)-1;
					if(!right_pushed_glb[vtid]) right_pushed_glb[vtid] = true;
				}
			}
			cur_st1 = cur_st1->right_c;	
		}


		// jumping to next level
		cur_bd = cur_bd->right_c;

	}

	
	// link based on v2par

	while(cur_lk){
		for(auto vt: cur_lk->nodeset){
			if(skyline_index_v2par[vt]){
				for(auto e: skyline_index_v2ccb[vt]){
					if(e->hubid == cur_lk->hubid
					&& !e->left_p){
						e->left_p = skyline_index_v2par[vt];
						skyline_index_v2par[vt]->right_c.push_back(e);
					
					}

				}
			}
				
		}

		cur_lk = cur_lk->right_c;
	}


	// restore edges for next path 

	for(auto e: vt_tbre){
		g.addEdge(e[0],e[1]);
	}
	for(auto e: vt_ngp){
		g.addEdge(e[0],e[1]);
	}
	
	vt_tbre.clear();
	vt_ngp.clear();
	tar.clear();
	skyline_index_v2par.clear();

}

void link_path_right(vector<int> link_start,
	BiGraph& g,
	unordered_map<vector<int>, skyline_index_hub*, VectorHasher>& skyline_index_lv0,
	unordered_map<string, vector<skyline_index_hub*>>& skyline_index_v2h,
    unordered_map<string, vector<skyline_index_ccblock*>>& skyline_index_v2ccb, 
	unordered_map<vector<int>, skyline_index_ccblock*, VectorHasher>& skyline_index_lv1){

    unordered_map<string, skyline_index_ccblock*> skyline_index_v2par;

	vector<bool> left_pushed; vector<bool> right_pushed;
	vector<bool> left_pushed_glb; vector<bool> right_pushed_glb;

	left_pushed.resize(g.num_v1, false); 
	right_pushed.resize(g.num_v2, false);
	left_pushed_glb.resize(g.num_v1, false); 
	right_pushed_glb.resize(g.num_v2, false);

	vector<int> tar;
	skyline_index_hub* cur_bd = skyline_index_lv0[link_start];
	skyline_index_hub* cur_lk = skyline_index_lv0[link_start];

	vector<vector<int>> vt_tbre;
	vector<vector<int>> vt_ngp;

	// global state labeling
	skyline_index_hub* cur_st = cur_bd;
	while(cur_st){
		for(auto vt: cur_st->nodeset){
			if(vt[0]=='u'){
				vt.erase(vt.begin());
				int vtid = stoi(vt)-1;
				if(!left_pushed_glb[vtid]) left_pushed_glb[vtid] = true;
			} else {
				vt.erase(vt.begin());
				int vtid = stoi(vt)-1;
				if(!right_pushed_glb[vtid]) right_pushed_glb[vtid] = true;
			}
		}
		cur_st = cur_st->left_c;	
	}


	// building
	while(cur_bd){

		fill_n(left_pushed.begin(), left_pushed.size(), false);
		fill_n(right_pushed.begin(), right_pushed.size(), false);

		// delete edges of blocks above
		
		if(cur_bd->right_p){
			for(auto vt: cur_bd->right_p->nodeset){
				if(vt[0]=='u'){
					vt.erase(vt.begin());
					int vtid = stoi(vt)-1;
					for(auto e: g.getV1Neighbors(vtid)){
						vector<int> tmp;
						if(g.isEdge(vtid, e)){
							g.deleteEdge(vtid,e);
							tmp = {vtid, e};
							vt_tbre.push_back(tmp);
						}
					}
				} else {
					vt.erase(vt.begin());
					int vtid = stoi(vt)-1;
					for(auto e: g.getV2Neighbors(vtid)){
						vector<int> tmp;
						if(g.isEdge(e, vtid)){
							g.deleteEdge(e,vtid);
							tmp = {e, vtid};
							vt_tbre.push_back(tmp);
						}
					}
				}
			}
		}
		
		tar = cur_bd->hubid;

		

		// state labeling of vertices in nodeset
		for(auto vt: cur_bd->nodeset){
			if(vt[0]=='u'){
				vt.erase(vt.begin());
				int vtid = stoi(vt)-1;
				if(!left_pushed[vtid]) left_pushed[vtid] = true;
			} else {
				vt.erase(vt.begin());
				int vtid = stoi(vt)-1;
				if(!right_pushed[vtid]) right_pushed[vtid] = true;
			}
		}
		

		queue<string> ccQue;
		string now;
		
		// ccQue for BFS to build each tar with cc_id
		for(auto vt: cur_bd->nodeset){

			if(vt[0]=='u'){
				vt.erase(vt.begin());
				int vtid = stoi(vt)-1;
				if(left_pushed[vtid]){
					vt = "u" + vt;
					ccQue.push(vt);
					for(auto p: skyline_index_v2ccb[vt]){
						if(p->hubid == tar){
							tar.push_back(p->ccid);
							break;
						}
					}
				}
			} else {
				vt.erase(vt.begin());
				int vtid = stoi(vt)-1;
				if(right_pushed[vtid]){
					vt = "v" + vt;
					ccQue.push(vt);
					for(auto p: skyline_index_v2ccb[vt]){
						if(p->hubid == tar){
							tar.push_back(p->ccid);
							break;
						}
					}
				}
			}
			
			while(!ccQue.empty()){

				now = ccQue.front();
				ccQue.pop();
				int vtid = 0;

				
				
				if(now[0]=='u'){

					now.erase(now.begin());
					vtid = stoi(now)-1;
					string nowPu = "u" + now;

					for(auto e: g.getV1Neighbors(vtid)){
						if(right_pushed[e]){
							string nowNei = "v" + to_string(e+1);
							right_pushed[e] = false;
							right_pushed_glb[e] = false;
							ccQue.push(nowNei);
						}
						else if(right_pushed_glb[e]){
							right_pushed_glb[e] = false;
							string nowNei = "v" + to_string(e+1);
							ccQue.push(nowNei);

							for(auto e: skyline_index_v2h[nowNei]){
								if(e->hubid == cur_bd->left_c->hubid
								&& !skyline_index_v2par[nowNei]){
									skyline_index_v2par[nowNei] = skyline_index_lv1[tar];
									break;
								}
							}
						}
					}
				} else {
					now.erase(now.begin());
					vtid = stoi(now)-1;
					string nowPu = "v" + now;

					for(auto e: g.getV2Neighbors(vtid)){
						if(left_pushed[e]){
							string nowNei = "u" + to_string(e+1);
							left_pushed[e] = false;
							left_pushed_glb[e] = false;
							ccQue.push(nowNei);
						}
						else if(left_pushed_glb[e]){
							left_pushed_glb[e] = false;
							string nowNei = "u" + to_string(e+1);
							ccQue.push(nowNei);
							for(auto e: skyline_index_v2h[nowNei]){
								if(e->hubid == cur_bd->left_c->hubid
								&& !skyline_index_v2par[nowNei]){
									skyline_index_v2par[nowNei] = skyline_index_lv1[tar];
									break;
								}
							}

						}
					}
					
				}
				

			}
			
			if(tar.size()==3) tar.pop_back();
			
		}


		// global state labeling
		skyline_index_hub* cur_st2 = cur_bd->left_c;
		while(cur_st2){
			for(auto vt: cur_st2->nodeset){
				if(vt[0]=='u'){
					vt.erase(vt.begin());
					int vtid = stoi(vt)-1;
					if(!left_pushed_glb[vtid]) left_pushed_glb[vtid] = true;
				} else {
					vt.erase(vt.begin());
					int vtid = stoi(vt)-1;
					if(!right_pushed_glb[vtid]) right_pushed_glb[vtid] = true;
				}
			}
			cur_st2 = cur_st2->left_c;	
		}


		// jumping to next level
		cur_bd = cur_bd->left_c;

	}

	
	// link based on v2par

	while(cur_lk){
		for(auto vt: cur_lk->nodeset){
			if(skyline_index_v2par[vt]){
				for(auto e: skyline_index_v2ccb[vt]){
					if(e->hubid == cur_lk->hubid
					&& !e->right_p){
						e->right_p = skyline_index_v2par[vt];
						skyline_index_v2par[vt]->left_c.push_back(e);
					
					}

				}
			}
				
		}

		cur_lk = cur_lk->left_c;
	}


	// restore edges for next path 

	for(auto e: vt_tbre){
		g.addEdge(e[0],e[1]);
	}
	for(auto e: vt_ngp){
		g.addEdge(e[0],e[1]);
	}
	
	vt_tbre.clear();
	vt_ngp.clear();
	tar.clear();
	skyline_index_v2par.clear();

}

void par_build_path_left(int start, int end,
    vector<vector<int>> build_start_set,
    BiGraph& g,
	unordered_map<vector<int>, skyline_index_hub*, VectorHasher>& skyline_index_lv0,
	unordered_map<string, vector<skyline_index_hub*>>& skyline_index_v2h,
    unordered_map<string, vector<skyline_index_ccblock*>>& skyline_index_v2ccb, 
	unordered_map<vector<int>, skyline_index_ccblock*, VectorHasher>& skyline_index_lv1){
        for(int i = start; i < end; i++){
            vector<int> build_start = build_start_set[i];
            build_path_left(build_start, g, skyline_index_lv0, skyline_index_v2h, skyline_index_v2ccb, 
	                        skyline_index_lv1);
        }

}

void par_link_path_right(int start, int end,
    vector<vector<int>> link_remain_set,
    BiGraph& g,
	unordered_map<vector<int>, skyline_index_hub*, VectorHasher>& skyline_index_lv0,
	unordered_map<string, vector<skyline_index_hub*>>& skyline_index_v2h,
    unordered_map<string, vector<skyline_index_ccblock*>>& skyline_index_v2ccb, 
	unordered_map<vector<int>, skyline_index_ccblock*, VectorHasher>& skyline_index_lv1){
        for(int i = start; i < end; i++){
            vector<int> link_start = link_remain_set[i];
            link_path_right(link_start, g, skyline_index_lv0, skyline_index_v2h, skyline_index_v2ccb, 
	                        skyline_index_lv1);
        }

}


void merge_skyline_index_lv1(
    unordered_map<string, vector<skyline_index_ccblock*>>& skyline_index_v2ccb, 
	unordered_map<vector<int>, skyline_index_ccblock*, VectorHasher>& skyline_index_lv1,
    vector<unordered_map<string, vector<skyline_index_ccblock*>>>& skyline_index_v2ccbsA,
    vector<unordered_map<string, vector<skyline_index_ccblock*>>>& skyline_index_v2ccbsB,
    vector<unordered_map<vector<int>, skyline_index_ccblock*, VectorHasher>>& skyline_index_lv1sA,
    vector<unordered_map<vector<int>, skyline_index_ccblock*, VectorHasher>>& skyline_index_lv1sB){

    for(int i=0; i < skyline_index_lv1sA.size(); i++){
        skyline_index_lv1.merge(skyline_index_lv1sA[i]);
    }
    for(int i=0; i < skyline_index_lv1sB.size(); i++){
        skyline_index_lv1.merge(skyline_index_lv1sB[i]);
    }
    for(int i=0; i < skyline_index_v2ccbsA.size(); i++){
        skyline_index_v2ccb.merge(skyline_index_v2ccbsA[i]);
    }
    for(int i=0; i < skyline_index_v2ccbsB.size(); i++){
        skyline_index_v2ccb.merge(skyline_index_v2ccbsB[i]);
    }

}

void parallel_build_skyline_index_lv1(BiGraph& g, 
	unordered_map<vector<int>, vector<string>, VectorHasher>& skyNode, 
	unordered_map<vector<int>, skyline_index_hub*, VectorHasher>& skyline_index_lv0,
	unordered_map<string, vector<skyline_index_hub*>>& skyline_index_v2h,
    unordered_map<string, vector<skyline_index_ccblock*>>& skyline_index_v2ccb, 
	unordered_map<vector<int>, skyline_index_ccblock*, VectorHasher>& skyline_index_lv1, 
    int threads){

    cout << "multithread thread_num: " << threads << endl;
    vector<BiGraph> bigraphsA;
	vector<BiGraph> bigraphsB;
    vector<vector<vector<int>>> build_start_setA;
    vector<vector<vector<int>>> link_remain_setB;
    vector<unordered_map<vector<int>, skyline_index_hub*, VectorHasher>> skyline_index_lv0sA;
    vector<unordered_map<string, vector<skyline_index_hub*>>> skyline_index_v2hsA;
    vector<unordered_map<string, vector<skyline_index_ccblock*>>> skyline_index_v2ccbsA;
    vector<unordered_map<string, vector<skyline_index_ccblock*>>> skyline_index_v2ccbsB;
    vector<unordered_map<vector<int>, skyline_index_ccblock*, VectorHasher>> skyline_index_lv1sA;
    vector<unordered_map<vector<int>, skyline_index_ccblock*, VectorHasher>> skyline_index_lv1sB;
	vector<thread> threadgroupsA;
	vector<thread> threadgroupsB;


    vector<bool> left_pushed; vector<bool> right_pushed;
	vector<bool> left_pushed_glb; vector<bool> right_pushed_glb;

	left_pushed.resize(g.num_v1, false); 
	right_pushed.resize(g.num_v2, false);
	left_pushed_glb.resize(g.num_v1, false); 
	right_pushed_glb.resize(g.num_v2, false);

	vector<vector<int>> build_start_set;
	vector<vector<int>> link_remain_set;

	// build from left side
	for(auto e: skyNode){
		if(skyline_index_lv0[e.first]->left_p == nullptr){
			build_start_set.push_back(e.first);
		}
	}

    // link remaining from right side
    for(auto e: skyNode){
		if(skyline_index_lv0[e.first]->right_p == nullptr){
			link_remain_set.push_back(e.first);
		}
	}

    int lsize = build_start_set.size();
    int rsize = link_remain_set.size();
    int gap_1 = ceil(lsize * 1.0 / threads);
    int gap_2 = ceil(rsize * 1.0 / threads);

    for (int i = 0; i < threads; i++) {
		bigraphsA.push_back(g);
		bigraphsB.push_back(g);
        build_start_setA.push_back(build_start_set);
        link_remain_setB.push_back(link_remain_set);
        skyline_index_lv0sA.push_back(skyline_index_lv0);
        skyline_index_v2hsA.push_back(skyline_index_v2h);
        skyline_index_v2ccbsA.push_back(skyline_index_v2ccb);
        skyline_index_v2ccbsB.push_back(skyline_index_v2ccb);
        skyline_index_lv1sA.push_back(skyline_index_lv1);
        skyline_index_lv1sB.push_back(skyline_index_lv1);
	}

    for (int i = 0; i < threads; i++) {
        threadgroupsA.push_back(thread(par_build_path_left,
                        i * gap_1,
                        MIN(lsize, i * gap_1 + gap_1),
                        ref(build_start_setA[i]),
                        ref(bigraphsA[i]),
                        ref(skyline_index_lv0sA[i]), 
                        ref(skyline_index_v2hsA[i]), 
                        ref(skyline_index_v2ccbsA[i]), 
                        ref(skyline_index_lv1sA[i])));
	}
    for (int i = 0; i < threadgroupsA.size(); i++) {
		threadgroupsA[i].join();
	}

    for (int i = 0; i < threads; i++) {
        threadgroupsB.push_back(thread(par_link_path_right,
                        i * gap_2,
                        MIN(rsize, i * gap_2 + gap_2),
                        ref(link_remain_setB[i]),
                        ref(bigraphsB[i]),
                        ref(skyline_index_lv0sA[i]), 
                        ref(skyline_index_v2hsA[i]), 
                        ref(skyline_index_v2ccbsB[i]), 
                        ref(skyline_index_lv1sB[i])));

	}
	for (int i = 0; i < threadgroupsB.size(); i++) {
		threadgroupsB[i].join();		
	}

    merge_skyline_index_lv1(skyline_index_v2ccb, skyline_index_lv1, skyline_index_v2ccbsA, 
                            skyline_index_v2ccbsB, skyline_index_lv1sA, skyline_index_lv1sB);


}