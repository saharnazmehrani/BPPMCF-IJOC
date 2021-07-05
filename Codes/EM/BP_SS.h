#ifndef BP_HPP_
#define BP_HPP_

#include "Instance_Stadium.h"
#include "MDD2.h"
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
#include <set>
#include "gurobi_c++.h"

using namespace std;

struct Variable {

	int 						index;
	vector<vector<int>>			assigned_group_size;	// equal to 1 if party "p" is assigned to the row
	int							capacity;				// capacity of the column
	double						value;

	Variable() {
		index = -1;
		value = 0;
		capacity = 0;
	}

};

struct SearchNode {

	int 							index;
	double							parent_lower_bound;     // lower bound for parent
	double							LP_lower_bound;         // lower bound for node

	vector<Variable* > 				variables_set_to_one;
	vector<Variable* > 				variables_set_to_zero;

	vector<int>						UB_seg_no_for_group;    //Lower bound for the number of segments must be assigned to each group
	vector<int>						LB_seg_no_for_group;	//Upper bound for the number of segments must be assigned to each group

	SearchNode() {
		index = -1;
		parent_lower_bound = numeric_limits<double>::min();
		LP_lower_bound = numeric_limits<double>::min();
	}

};

bool compare_nodes(SearchNode* node1, SearchNode* node2);

struct BPSolver {

	// objects

	Stadium_Seating* 						inst;
	MDDSolver*                              MDD;
	vector<Variable* > 						variables;
	//priority_queue<SearchNode* > 			search_nodes;
	vector<SearchNode* >					search_nodes; //NOTE: It should be a priority queue - much faster
	double									global_lower_bound;
	double									global_upper_bound;
	vector<Variable* >						best_solution;

	// member function
	void 									create_first_set_of_variables();
	void 									create_first_set_of_variables2(ofstream& file1);
	void									create_first_search_node();

	// print functions
	void 									print_variable(Variable* var_to_print);
	void 									print_search_node(SearchNode* node);
	void 									print_selected_node(SearchNode* node, ofstream& file1);
	void									print_best_solution();

	// solving functions
	void									pick_branch_variable(SearchNode* node, double & time_limit, ofstream& file1, int & branch_var_index, vector<int> & group_branch_index_and_value);
	void									node_feasibility_check(SearchNode* node, double & time_limit, bool & node_is_feasible, int n_node_explored);
	void 									solve_LP(SearchNode* node, double & time_left, double time_limit, clock_t startTime, ofstream& file1, int& n_CG_iter, int& n_price_iter, double& total_price_time);
	void 									solve_IP(SearchNode* node, double & time_limit, ofstream& file1);
	void 									solve(double & time_limit, clock_t startTime, ofstream& file1);
	void									solve_RMP_LP(SearchNode* node, double & time_limit, ofstream& file1, vector<vector<double>> & party_duals, vector<double> & capacity_duals, vector<double> & seg_UB_duals, vector<double> & new_LB_duals);
	void									solve_pricing_problem(SearchNode* node, double & time_limit, ofstream& file1, vector<vector<double>> & party_duals, vector<double> & capacity_duals, vector<double> & seg_UB_duals, vector<double> & new_LB_duals, int capacity, bool & need_another_rmp_lp_solution, int& n_price_iter, double& total_price_time);


	// new functions added in revision 1 from IJOC
	void                                    solve_MDD_net_model(double & time_limit, ofstream& file1);
	void                                    solve_MDD_pricing_problem(SearchNode* node, double & time_limit, ofstream& file1, vector<vector<double>> & party_duals, vector<double> & capacity_duals, vector<double> & seg_LB_duals, vector<double> & seg_UB_duals, int capacity, bool & need_another_rmp_lp_solution, int& n_price_iter, double& total_price_time);
	void                                    solve_LP_dual_stabilization(SearchNode* node, double & time_left, double time_limit, clock_t startTime, ofstream& file1, int& n_CG_iter, int& n_price_iter, double& total_price_time);
	void                                    solve_LP_with_pricing_and_DUALstabil(SearchNode* node, double & time_left, double time_limit, clock_t startTime, ofstream& file1, int& n_CG_iter, int& n_price_iter, double& total_price_time);
	void                                    solve_LP_with_pricing_and_DUALstabil2(SearchNode* node, double & time_left, double time_limit, clock_t startTime, ofstream& file1, int& n_CG_iter, int& n_price_iter, double& total_price_time);
	void									Farkas_node_feasibility_check(SearchNode* node, double & time_limit, clock_t startTime, bool & node_is_feasible, int n_node_explored);
	void                                    create_first_set_of_variables_IP(ofstream& file1);
	void                                    create_first_set_of_variables_bin_completion(ofstream& file1);


	BPSolver(Stadium_Seating* _inst, MDDSolver* _MDD) {
		inst = _inst;
		MDD = _MDD;
		global_lower_bound = 0;
		global_upper_bound = inst->Tn_party;
	}

};

#endif