#ifndef BP_HPP_
#define BP_HPP_

#include "Instance_Stadium.h"
#include "uti.h"
#include <iostream>     // allows you to read and write to terminal
#include <vector>       // allows you to use vectors
#include <fstream>      // allows you to read and write to files
//#include  <queue>		// priority queue of nodes
#include <limits>
#include <stdlib.h>
#include <ctime>
#include <tuple>
#include <cstdlib>
#include <random>
#include <sstream>
#include <string>
#include <list>
#include <cmath>
#include "math.h"
#include <vector>
#include <algorithm>
#include <iomanip>
#include <ilcplex/ilocplex.h>
#include <set>
#include "gurobi_c++.h"

using namespace std;

struct BDDNode {

	int							index;
	vector<int>					parent_indexes;
	int							segment_capacity;		//number of seat in the segemnt that node belongs to
	int							party;					//party index

	int							capacity;				//number of seats left in the segment
	vector<bool>				party_taken;			//whether the party is selected or not (considered "verctor" is because of merging nodes)
	bool						if_party_from_same_group_selected;	// It is equal to one if at least one party from the current group is selected; 0 otherwise
	vector<bool>				value;					//it is equal to one if "selected" is greater than zero and it is the first party from its group that is selected up to that node 

	vector<int> input_arcs;								//indexes of the arcs that enter the node
	vector<int> output_arcs;							//indexes of the arcs that go out of the node


	BDDNode() {

		index = -1;
		party = -1;
		capacity = 0;
		if_party_from_same_group_selected = false;

	}

};

struct BDDArc {

	int							index;
	int							i;						//origin node
	int							j;						//destination node
	int							party;					//party index corresponding to the arc
	int							party_taken;			//whether the party is selected or not
	bool						value;					//value is equal to one if it is the first party in the corresponding group which is selected

	BDDArc() {
		value = false;
	}

};

struct BDDSolver {

	// objects
	vector<vector<BDDNode* >>				BDD_nodes;
	vector<vector<BDDNode* >>				BDD_nodes_left_to_branch;
	vector<vector<int>>						T_conn_nodes_index;			//including the set of nodes connected to the terminal node in each BDD
	vector<vector<BDDArc* >>				BDD_arcs;
	Stadium_Seating* 						inst;

	// member functions
	void									create_first_BDD_node(int segment_capacity);
	void 									create_BDDs(double & time_limit, ofstream& file1);
	void									merge_nodes(bool & merged, int seg_capacity, vector<vector<BDDNode* >> & BDD_nodes, vector<vector<BDDNode* >> & BDD_nodes_left_to_branch_on, BDDNode* new_BDD_node);

	// print functions
	void									print_node(BDDNode* node);
	void									print_arc(BDDArc* node);

	// solving functions
	void									solve_net_model(double & time_limit, ofstream& file1);

	BDDSolver(Stadium_Seating* _inst) {
		inst = _inst;
	}

};



#endif
