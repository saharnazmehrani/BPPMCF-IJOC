#include "BDD.h"
#include <ctime>

void BDDSolver::create_BDDs(double & time_limit, ofstream& file1) {

	cout << "creating BDD" << endl;

	clock_t startTime = clock();
	double time_left;

	//creating the nodes in all BDDs (We have BDD per segemnt capacity
	BDD_nodes.resize(inst->n_SeatRow);
	BDD_nodes_left_to_branch.resize(inst->n_SeatRow);
	T_conn_nodes_index.resize(inst->n_SeatRow);		//will include all the nodes which are connected to the terminal node

	for (int cap = 0; cap < inst->n_SeatRow; ++cap) {
		if (inst->no_seg_of_each_size[cap] > 0 && double(clock() - startTime) / (double)CLOCKS_PER_SEC < time_limit) {

			create_first_BDD_node(cap + 1);

			while (BDD_nodes_left_to_branch[cap].size() > 0) {
				if (double(clock() - startTime) / (double)CLOCKS_PER_SEC >= time_limit) {
					cout << "time out and exit!!!!" << endl;
					break;
				}
				else {
					time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
				}

				BDDNode* current_node;
				current_node = BDD_nodes_left_to_branch[cap].front();
				BDD_nodes_left_to_branch[cap].erase(BDD_nodes_left_to_branch[cap].begin());

				BDDNode* new_BDD_node;

				//Branching the current node
				bool did_branching = false;
				if (current_node->party < inst->Tn_party - 1) {
					if (current_node->capacity > 1) {

						//// party is not selected
						new_BDD_node = new BDDNode();
						new_BDD_node->index = BDD_nodes[cap].size();
						new_BDD_node->parent_indexes.push_back(current_node->index);
						new_BDD_node->segment_capacity = cap + 1;
						new_BDD_node->party = current_node->party + 1;
						new_BDD_node->capacity = current_node->capacity;
						new_BDD_node->party_taken.push_back(false);
						if (current_node->if_party_from_same_group_selected && inst->subG_party[current_node->party] == inst->subG_party[current_node->party + 1]) {
							new_BDD_node->if_party_from_same_group_selected = true;
						}
						else {
							new_BDD_node->if_party_from_same_group_selected = false;
						}
						new_BDD_node->value.push_back(false);

						//merging this new node with any current one
						bool merged = false;
						merge_nodes(merged, cap + 1, BDD_nodes, BDD_nodes_left_to_branch, new_BDD_node);
						if (merged == false) {
							BDD_nodes[cap].push_back(new_BDD_node);
							BDD_nodes_left_to_branch[cap].push_back(new_BDD_node);
						}
						did_branching = true;

						//// party is selected
						if (inst->n[current_node->party + 1] <= current_node->capacity) {
							new_BDD_node = new BDDNode();
							new_BDD_node->index = BDD_nodes[cap].size();
							new_BDD_node->parent_indexes.push_back(current_node->index);
							new_BDD_node->segment_capacity = cap + 1;
							new_BDD_node->party = current_node->party + 1;
							new_BDD_node->capacity = current_node->capacity - inst->n[current_node->party + 1];
							new_BDD_node->party_taken.push_back(true);
							new_BDD_node->if_party_from_same_group_selected = true;
							if (current_node->if_party_from_same_group_selected && inst->subG_party[current_node->party] == inst->subG_party[current_node->party + 1]) {
								new_BDD_node->value.push_back(false);
							}
							else {
								new_BDD_node->value.push_back(true);
							}

							//merging this new node with any current one
							bool merged = false;
							merge_nodes(merged, cap + 1, BDD_nodes, BDD_nodes_left_to_branch, new_BDD_node);
							if (merged == false) {
								BDD_nodes[cap].push_back(new_BDD_node);
								BDD_nodes_left_to_branch[cap].push_back(new_BDD_node);
							}

							did_branching = true;
						}

					}
				}

				if (did_branching == false) {
					T_conn_nodes_index[cap].push_back(current_node->index);
				}
				//end of Branching the current node

			}

		}
	}

	/////////////// Creating the Arcs in all BDDs /////////////////////
	///////////////////////////////////////////////////////////////////

	if (double(clock() - startTime) / (double)CLOCKS_PER_SEC < time_limit) {
		time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);

		BDD_arcs.resize(inst->n_SeatRow);
		BDDArc* new_arc;
		for (int cap = 0; cap < inst->n_SeatRow; ++cap) {
			if (inst->no_seg_of_each_size[cap] > 0) {

				for (int jj = 1; jj < (int)BDD_nodes[cap].size(); jj++) {
					for (int aa = 0; aa < BDD_nodes[cap][jj]->parent_indexes.size(); ++aa) {
						new_arc = new BDDArc();
						new_arc->j = jj;
						new_arc->index = BDD_arcs[cap].size();
						new_arc->i = BDD_nodes[cap][jj]->parent_indexes[aa];
						new_arc->party = BDD_nodes[cap][jj]->party;
						new_arc->party_taken = BDD_nodes[cap][jj]->party_taken[aa];
						new_arc->value = BDD_nodes[cap][jj]->value[aa];
						BDD_nodes[cap][new_arc->j]->input_arcs.push_back(new_arc->index);
						BDD_nodes[cap][new_arc->i]->output_arcs.push_back(new_arc->index);

						BDD_arcs[cap].push_back(new_arc);
					}
				}

				//adding Arcs (i,terminal)
				for (int ii = 0; ii < (int)T_conn_nodes_index[cap].size(); ii++) {
					new_arc = new BDDArc();
					new_arc->index = BDD_arcs[cap].size();
					new_arc->j = BDD_nodes[cap].size();
					new_arc->i = T_conn_nodes_index[cap][ii];
					new_arc->party = -1;
					new_arc->party_taken = false;
					new_arc->value = false;

					BDD_nodes[cap][new_arc->i]->output_arcs.push_back(new_arc->index);

					BDD_arcs[cap].push_back(new_arc);
				}

			}
		}

	}

	if (double(clock() - startTime) / (double)CLOCKS_PER_SEC < time_limit) {
		time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
	}
	else {
		time_left = 0;
		time_limit = time_left;
	}
	//time_limit = time_left;

}


void BDDSolver::solve_net_model(double & time_limit, ofstream& file1) {

	//creating BDD nodes and arcs
	create_BDDs(time_limit, file1);
	cout << "BDD is created" << endl;

	int n_arcs = 0;
	int n_nodes = 0;

	//int n_BDDarcs = BDD_arcs[inst->n_SeatRow - 1].size();
	n_arcs = BDD_arcs[inst->n_SeatRow - 1].size();
	n_nodes = BDD_nodes[inst->n_SeatRow - 1].size();

	file1 << n_nodes << ",";								//number of MDD nodes
	file1 << n_arcs << ",";								//number of MDD arcs


	if (time_limit > PRECISION) {
		cout << "=============================================" << endl << "Solving Arc Based Network Flow Model:" << endl << "=============================================";
		GRBEnv env;
		try {

			GRBModel model(env);

			//Defining and adding arc variables
			int count = 0;
			int n_const = 0;

			GRBVar** x;																		//x[s][a]=1 : if arc "a" in the BDD of segment "s" is selected
			x = new GRBVar*[inst->n_seg];
			for (int s = 0; s < inst->n_seg; ++s) {
				int n_BDDarcs = BDD_arcs[inst->no_seat_in_seg[s] - 1].size();
				x[s] = new GRBVar[n_BDDarcs];
				for (int a = 0; a < n_BDDarcs; ++a) {
					x[s][a] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					count++;
				}
			}

			model.update();


			//adding constraints

			// in-flow = out-flow
			for (int s = 0; s < inst->n_seg; ++s) {
				for (int k = 1; k < n_nodes; ++k) {
					GRBLinExpr equ = 0.0;
					for (int in = 0; in < BDD_nodes[inst->no_seat_in_seg[s] - 1][k]->input_arcs.size(); ++in) {
						equ += x[s][BDD_nodes[inst->no_seat_in_seg[s] - 1][k]->input_arcs[in]];
					}
					for (int out = 0; out < BDD_nodes[inst->no_seat_in_seg[s] - 1][k]->output_arcs.size(); ++out) {
						equ -= x[s][BDD_nodes[inst->no_seat_in_seg[s] - 1][k]->output_arcs[out]];
					}
					model.addConstr(equ == 0.0);
					n_const++;
				}
			}


			// root-node single flow
			for (int s = 0; s < inst->n_seg; s++) {
				GRBLinExpr equ = 0.0;
				for (int out = 0; out < BDD_nodes[inst->no_seat_in_seg[s] - 1][0]->output_arcs.size(); ++out) {
					equ += x[s][BDD_nodes[inst->no_seat_in_seg[s] - 1][0]->output_arcs[out]];
				}
				model.addConstr(equ == 1);
				n_const++;
			}


			//Each party must be selected once
			for (int p = 0; p < inst->Tn_party; ++p) {
				GRBLinExpr party_equ = 0.0;
				for (int s = 0; s < inst->n_seg; ++s) {
					int n_BDDarcs = BDD_arcs[inst->no_seat_in_seg[s] - 1].size();
					for (int a = 0; a < n_BDDarcs; ++a) {
						if (BDD_arcs[inst->no_seat_in_seg[s] - 1][a]->party == p) {
							if (BDD_arcs[inst->no_seat_in_seg[s] - 1][a]->party_taken) {
								party_equ += x[s][a];
							}
						}
					}
				}
				model.addConstr(party_equ == 1.0);
				n_const++;
			}


			//adding the new lower bounds from Silvano Martello and Paolo Toth per group
			for (int g = 0; g < inst->n_subGroup; ++g) {
				GRBLinExpr LB_group_equ = 0.0;
				for (int s = 0; s < inst->n_seg; ++s) {
					int n_BDDarcs = BDD_arcs[inst->no_seat_in_seg[s] - 1].size();
					for (int a = 0; a < n_BDDarcs; ++a) {
						//cout << "hey " << inst->subG_party[BDD_arcs[inst->no_seat_in_seg[s] - 1][a]->party] << endl;
						if (BDD_arcs[inst->no_seat_in_seg[s] - 1][a]->party >= 0) {
							if (inst->subG_party[BDD_arcs[inst->no_seat_in_seg[s] - 1][a]->party] == g) {
								LB_group_equ += BDD_arcs[inst->no_seat_in_seg[s] - 1][a]->value * x[s][a];
							}
						}
					}
				}
				
				double new_LB = 0;
				if (inst->LB_MT1[g] > new_LB) {
					new_LB = inst->LB_MT1[g];
				}
				if (inst->LB_MT2[g] > new_LB) {
					new_LB = inst->LB_MT2[g];
				}
				model.addConstr(LB_group_equ >= new_LB);
			}


			///////////////////////////
			//Adding objective function 
			GRBQuadExpr objective = 0.0;
			for (int s = 0; s < inst->n_seg; ++s) {
				int n_BDDarcs = BDD_arcs[inst->no_seat_in_seg[s] - 1].size();
				for (int a = 0; a < n_BDDarcs; ++a) {
					objective += BDD_arcs[inst->no_seat_in_seg[s] - 1][a]->value * x[s][a];
				}
			}

		

			model.update();
			model.setObjective(objective, GRB_MINIMIZE);

			model.getEnv().set(GRB_IntParam_Threads, 1);
			model.getEnv().set(GRB_DoubleParam_TimeLimit, time_limit);
			model.getEnv().set(GRB_IntParam_Method, 2);
			//model.getEnv().set(GRB_IntParam_OutputFlag, 0);
			//model.getEnv().set(GRB_IntParam_PreCrush, 1);
			//model.getEnv().set(GRB_IntParam_DualReductions, 0);
			//model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
			//model.getEnv().set(GRB_IntParam_Presolve, -1);
			//model.getEnv().set(GRB_IntParam_MIPFocus, 2);
			model.update();

			int counter_variables = model.get(GRB_IntAttr_NumVars);
			int counter_constraints = model.get(GRB_IntAttr_NumConstrs);


			model.optimize();
			cout << endl;
			cout << "Solution Status: " << model.get(GRB_IntAttr_Status) << endl;
			//file1 << model.get(GRB_IntAttr_Status) << ",";
			cout << "Time Taken: " << model.get(GRB_DoubleAttr_Runtime) << endl;
			file1 << model.get(GRB_DoubleAttr_Runtime) << ",";


			if (model.get(GRB_IntAttr_Status) == 2 || model.get(GRB_IntAttr_Status) == 9) {

				file1 << model.get(GRB_DoubleAttr_ObjBound) << ",";
				file1 << model.get(GRB_DoubleAttr_ObjVal) << ",";
				file1 << model.get(GRB_IntAttr_NumVars) << ",";
				file1 << model.get(GRB_IntAttr_NumConstrs) << ",";
				file1 << model.get(GRB_IntAttr_NumNZs) << ",";
				file1 << model.get(GRB_DoubleAttr_NodeCount) << endl;

				/*for (int s = 0; s < single_seg; s++) {
				int n_BDDarcs = BDD_arcs[inst->no_seat_in_seg[s] - 1].size();
				cout << "segment " << s << " :" << endl;
				for (int a = 0; a < n_BDDarcs; a++) {
				if (x[s][a].get(GRB_DoubleAttr_X) > 0.9) {
				print_arc(BDD_arcs[inst->no_seat_in_seg[s] - 1][a]);
				}
				}
				}*/

				double gap = -1.0;
				double ub = model.get(GRB_DoubleAttr_ObjVal);
				double lb = model.get(GRB_DoubleAttr_ObjBound);
				double runtime = model.get(GRB_DoubleAttr_Runtime);

				if (model.get(GRB_DoubleAttr_ObjVal) < 1e6) {
					gap = abs((ub - lb) / ub);
				}

				cout << "Results," << runtime << "," << lb << "," << ub << "," << gap << "\n";

			}

			delete x;

		}
		catch (GRBException e) {
			cout << "Error xmber: " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
			exit(1);
		}
		catch (...) {
			cout << "Other error ... " << endl;
			exit(1);
		}

	}

	else {
		file1 << " " << ",";
		file1 << " " << ",";
		file1 << 0 << ",";
	}

}

void BDDSolver::create_first_BDD_node(int segment_capacity) {

	BDDNode* first_BDD_node;
	first_BDD_node = new BDDNode();

	first_BDD_node->index = 0;
	first_BDD_node->parent_indexes.resize(1, -1);
	first_BDD_node->segment_capacity = segment_capacity;
	first_BDD_node->party = -1;

	first_BDD_node->capacity = segment_capacity;
	first_BDD_node->party_taken.push_back(false);
	first_BDD_node->value.resize(1, 0);


	BDD_nodes[segment_capacity - 1].push_back(first_BDD_node);
	BDD_nodes_left_to_branch[segment_capacity - 1].push_back(first_BDD_node);

}

void BDDSolver::merge_nodes(bool & merged, int seg_capacity, vector<vector<BDDNode* >> & BDD_nodes, vector<vector<BDDNode* >> & BDD_nodes_left_to_branch_on, BDDNode* new_BDD_node) {
	for (int n1 = 0; n1 < BDD_nodes_left_to_branch_on[seg_capacity - 1].size(); ++n1) {
		if (BDD_nodes_left_to_branch_on[seg_capacity - 1][n1]->capacity == new_BDD_node->capacity		 &&			BDD_nodes_left_to_branch_on[seg_capacity - 1][n1]->party == new_BDD_node->party && BDD_nodes_left_to_branch_on[seg_capacity - 1][n1]->if_party_from_same_group_selected == new_BDD_node->if_party_from_same_group_selected) {
			merged = true;
			int index = BDD_nodes_left_to_branch_on[seg_capacity - 1][n1]->index;
			BDD_nodes[seg_capacity - 1][index]->parent_indexes.push_back(new_BDD_node->parent_indexes[0]);
			BDD_nodes[seg_capacity - 1][index]->party_taken.push_back(new_BDD_node->party_taken[0]);
			BDD_nodes[seg_capacity - 1][index]->value.push_back(new_BDD_node->value[0]);

			BDD_nodes_left_to_branch_on[seg_capacity - 1][n1]->parent_indexes.push_back(new_BDD_node->parent_indexes[0]);
			BDD_nodes_left_to_branch_on[seg_capacity - 1][n1]->party_taken.push_back(new_BDD_node->party_taken[0]);
			BDD_nodes_left_to_branch_on[seg_capacity - 1][n1]->value.push_back(new_BDD_node->value[0]);
			break;
		}
	}
}

void BDDSolver::print_arc(BDDArc* arc) {

	cout << "*************BDD Arc! ************" << endl;
	cout << "\tArc_index: " << arc->index << endl;
	cout << "\tarc: (" << arc->i << "," << arc->j << ")" << endl;
	cout << "\tparty_taken: " << arc->party_taken << endl;
	cout << "\tparty: " << arc->party << endl;
	cout << "\tarc value: " << arc->value << endl;
	cout << "**************************************" << endl;

}



