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

struct MDDNode {

	int							index;
	vector<int>					parent_indexes;
	int							segment_capacity;		//number of seat in the segemnt that node belongs to
	vector<int>					size_group;				//it has two elements: party-size and group corresponding to the arc

	int							capacity;				//number of seats left in the segment
	vector<int>					n_party_selected;		//nubmer of parties from a specific size and group seleced
	bool						if_party_from_same_group_selected;	// It is equal to one if at least one party from the current group is selected; 0 otherwise
	vector<bool>				value;					//it is equal to one if "selected" is greater than zero and it is the first party from its group that is selected up to that node 

	vector<int> input_arcs;								//indexes of the arcs that enter the node
	vector<int> output_arcs;							//indexes of the arcs that go out of the node

	MDDNode() {
		index = -1;
		size_group.resize(2, -1);
		capacity = 0;
		if_party_from_same_group_selected = false;
	}

};

struct MDDArc {

	int							index;
	int							i;						//origin node
	int							j;						//destination node
	vector<int>					size_group;				//it has two elements: party-size and group corresponding to the arc
	int							n_party_selected;		//nubmer of parties from a specific size and group seleced
	bool						value;					//value is equal to one if it is the first party in the corresponding group which is selected

	MDDArc() {
		value = false;
		size_group.resize(2);
	}

};

struct MDDSolver {

	// objects
	vector<vector<MDDNode* >>				MDD_nodes;
	vector<vector<MDDNode* >>				MDD_nodes_left_to_branch;
	vector<vector<int>>						T_conn_nodes_index;			//including the set of nodes connected to the terminal node in each MDD
	vector<vector<MDDArc* >>				MDD_arcs;
	Stadium_Seating* 						inst;

	// member functions
	void									create_first_MDD_node(int segment_capacity);
	void 									create_MDDs(double & time_limit, ofstream& file1);
	void									merge_nodes(bool & merged, int seg_capacity, vector<vector<MDDNode* >> & MDD_nodes, vector<vector<MDDNode* >> & MDD_nodes_left_to_branch_on, MDDNode* new_MDD_node);

	// print functions
	void									print_node(MDDNode* node);
	void									print_arc(MDDArc* node);

	// solving functions
	/*void									solve_MDDs(double & time_limit, ofstream& file1);*/
	void									solve_net_model(double & time_limit, ofstream& file1);

	MDDSolver(Stadium_Seating* _inst) {
		inst = _inst;
	}

};

#endif