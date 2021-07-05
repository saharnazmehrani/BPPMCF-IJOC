#include "MDD2.h"

void MDDSolver::create_MDDs(double & time_limit, ofstream& file1) {

	cout << "creating MDD" << endl;

	clock_t startTime = clock();
	double time_left;


	//creating the nodes in all MDDs (We have MDD per segemnt capacity)
	MDD_nodes.resize(inst->n_SeatRow);
	MDD_nodes_left_to_branch.resize(inst->n_SeatRow);
	T_conn_nodes_index.resize(inst->n_SeatRow);		//will include all the nodes which are connected to the terminal node
													//for (int cap = 0; cap < inst->n_SeatRow; ++cap) {
													//if (inst->no_seg_of_each_size[cap] > 0) {		
	int cap = inst->n_SeatRow - 1;

	create_first_MDD_node(cap + 1);

	while (MDD_nodes_left_to_branch[cap].size() > 0) {

		if (double(clock() - startTime) / (double)CLOCKS_PER_SEC >= time_limit) {
			cout << "time out and exit!!!!" << endl;
			break;
		}
		else {
			time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
		}

		MDDNode* current_node;
		current_node = MDD_nodes_left_to_branch[cap].front();
		MDD_nodes_left_to_branch[cap].erase(MDD_nodes_left_to_branch[cap].begin());

		MDDNode* new_MDD_node;

		//Branching the current node
		bool did_branching = false;
		if (current_node->size_group[0] < inst->max_party_size || current_node->size_group[1] < inst->n_subGroup - 1) {
			if (current_node->capacity > 1) {

				int party_size = 1;
				if (current_node->size_group[0] < inst->max_party_size) {
					party_size = current_node->size_group[0] + 1;
					while (inst->u[current_node->size_group[1]][party_size - 1] == 0 && party_size < inst->max_party_size) {
						party_size++;
					}
				}
				////branching on the same group but new party size:
				if (party_size > 1) {
					if (inst->u[current_node->size_group[1]][party_size - 1] > 0) {
						//// party is not selected
						new_MDD_node = new MDDNode();
						new_MDD_node->index = MDD_nodes[cap].size();
						new_MDD_node->parent_indexes.push_back(current_node->index);
						new_MDD_node->segment_capacity = cap + 1;
						new_MDD_node->size_group[0] = party_size;
						new_MDD_node->size_group[1] = current_node->size_group[1];
						new_MDD_node->capacity = current_node->capacity;
						new_MDD_node->n_party_selected.push_back(0);

						if (current_node->if_party_from_same_group_selected) {
							new_MDD_node->if_party_from_same_group_selected = true;
						}
						else {
							new_MDD_node->if_party_from_same_group_selected = false;
						}
						new_MDD_node->value.push_back(false);

						//merging this new node with any current one
						bool merged = false;
						merge_nodes(merged, cap + 1, MDD_nodes, MDD_nodes_left_to_branch, new_MDD_node);
						if (merged == false) {
							MDD_nodes[cap].push_back(new_MDD_node);
							MDD_nodes_left_to_branch[cap].push_back(new_MDD_node);
						}
						did_branching = true;

						int n_party = 1;
						while (current_node->capacity >= n_party*party_size && n_party <= inst->u[current_node->size_group[1]][party_size - 1]) {
							//// party is selected
							new_MDD_node = new MDDNode();
							new_MDD_node->index = MDD_nodes[cap].size();
							new_MDD_node->parent_indexes.push_back(current_node->index);
							new_MDD_node->segment_capacity = cap + 1;
							new_MDD_node->size_group[0] = party_size;
							new_MDD_node->size_group[1] = current_node->size_group[1];
							new_MDD_node->capacity = current_node->capacity - n_party*party_size;
							new_MDD_node->n_party_selected.push_back(n_party);
							new_MDD_node->if_party_from_same_group_selected = true;
							if (current_node->if_party_from_same_group_selected) {
								new_MDD_node->value.push_back(false);
							}
							else {
								new_MDD_node->value.push_back(true);
							}

							//merging this new node with any current one
							bool merged = false;
							merge_nodes(merged, cap + 1, MDD_nodes, MDD_nodes_left_to_branch, new_MDD_node);
							if (merged == false) {
								MDD_nodes[cap].push_back(new_MDD_node);
								MDD_nodes_left_to_branch[cap].push_back(new_MDD_node);
							}

							did_branching = true;
							++n_party;

						}

					}
				}

				////branching on the new group and new party size////
				///////// first scenario:
				if (party_size == 1 && current_node->size_group[1] + 1 < inst->n_subGroup) {
					while (inst->u[current_node->size_group[1] + 1][party_size - 1] == 0 && party_size < inst->max_party_size) {
						party_size++;
					}

					//// party is not selected
					new_MDD_node = new MDDNode();
					new_MDD_node->index = MDD_nodes[cap].size();
					new_MDD_node->parent_indexes.push_back(current_node->index);
					new_MDD_node->segment_capacity = cap + 1;
					new_MDD_node->size_group[0] = party_size;
					new_MDD_node->size_group[1] = current_node->size_group[1] + 1;
					new_MDD_node->capacity = current_node->capacity;
					new_MDD_node->n_party_selected.push_back(0);
					new_MDD_node->if_party_from_same_group_selected = false;
					new_MDD_node->value.push_back(false);

					//merging this new node with any current one
					bool merged = false;
					merge_nodes(merged, cap + 1, MDD_nodes, MDD_nodes_left_to_branch, new_MDD_node);
					if (merged == false) {
						MDD_nodes[cap].push_back(new_MDD_node);
						MDD_nodes_left_to_branch[cap].push_back(new_MDD_node);
					}
					did_branching = true;

					int n_party = 1;
					while (current_node->capacity >= n_party*party_size && n_party <= inst->u[current_node->size_group[1] + 1][party_size - 1]) {
						//// party is selected
						new_MDD_node = new MDDNode();
						new_MDD_node->index = MDD_nodes[cap].size();
						new_MDD_node->parent_indexes.push_back(current_node->index);
						new_MDD_node->segment_capacity = cap + 1;
						new_MDD_node->size_group[0] = party_size;
						new_MDD_node->size_group[1] = current_node->size_group[1] + 1;
						new_MDD_node->capacity = current_node->capacity - n_party*party_size;
						new_MDD_node->n_party_selected.push_back(n_party);
						new_MDD_node->if_party_from_same_group_selected = true;
						new_MDD_node->value.push_back(true);

						//merging this new node with any current one
						bool merged = false;
						merge_nodes(merged, cap + 1, MDD_nodes, MDD_nodes_left_to_branch, new_MDD_node);
						if (merged == false) {
							MDD_nodes[cap].push_back(new_MDD_node);
							MDD_nodes_left_to_branch[cap].push_back(new_MDD_node);
						}
						did_branching = true;
						++n_party;
					}

				}
				///////// second scenario:
				if (party_size == inst->max_party_size && inst->u[current_node->size_group[1]][party_size - 1] == 0 && current_node->size_group[1] + 1 < inst->n_subGroup) {
					party_size = 1;
					while (inst->u[current_node->size_group[1] + 1][party_size - 1] == 0 && party_size < inst->max_party_size) {
						party_size++;
					}

					//// party is not selected
					new_MDD_node = new MDDNode();
					new_MDD_node->index = MDD_nodes[cap].size();
					new_MDD_node->parent_indexes.push_back(current_node->index);
					new_MDD_node->segment_capacity = cap + 1;
					new_MDD_node->size_group[0] = party_size;
					new_MDD_node->size_group[1] = current_node->size_group[1] + 1;
					new_MDD_node->capacity = current_node->capacity;
					new_MDD_node->n_party_selected.push_back(0);
					new_MDD_node->if_party_from_same_group_selected = false;
					new_MDD_node->value.push_back(false);

					//merging this new node with any current one
					bool merged = false;
					merge_nodes(merged, cap + 1, MDD_nodes, MDD_nodes_left_to_branch, new_MDD_node);
					if (merged == false) {
						MDD_nodes[cap].push_back(new_MDD_node);
						MDD_nodes_left_to_branch[cap].push_back(new_MDD_node);
					}
					did_branching = true;

					int n_party = 1;
					while (current_node->capacity >= n_party*party_size && n_party <= inst->u[current_node->size_group[1] + 1][party_size - 1]) {
						//// party is selected
						new_MDD_node = new MDDNode();
						new_MDD_node->index = MDD_nodes[cap].size();
						new_MDD_node->parent_indexes.push_back(current_node->index);
						new_MDD_node->segment_capacity = cap + 1;
						new_MDD_node->size_group[0] = party_size;
						new_MDD_node->size_group[1] = current_node->size_group[1] + 1;
						new_MDD_node->capacity = current_node->capacity - n_party*party_size;
						new_MDD_node->n_party_selected.push_back(n_party);
						new_MDD_node->if_party_from_same_group_selected = true;
						new_MDD_node->value.push_back(true);

						//merging this new node with any current one
						bool merged = false;
						merge_nodes(merged, cap + 1, MDD_nodes, MDD_nodes_left_to_branch, new_MDD_node);
						if (merged == false) {
							MDD_nodes[cap].push_back(new_MDD_node);
							MDD_nodes_left_to_branch[cap].push_back(new_MDD_node);
						}
						did_branching = true;
						++n_party;
					}

				}
			}

		}

		if (did_branching == false) {
			T_conn_nodes_index[cap].push_back(current_node->index);
		}
		//end of Branching the current node

	}
	//}
	//}

	////Printing MDDs' nodes
	//for (int cap = 0; cap < inst->n_SeatRow; ++cap) {
	//	if (inst->no_seg_of_each_size[cap] > 0) {
	//		cout << "printing nodes in a segment with capacity " << cap + 1 << " :" << endl;
	//		for (int i = 0; i < (int)MDD_nodes[cap].size(); ++i) {
	//			print_node(MDD_nodes[cap][i]);
	//		}
	//	}
	//	cout << "number of connected nodes to the terminal: " << T_conn_nodes_index[cap].size() << endl;
	//}

	// Creating the Arcs in all MDDs
	if (double(clock() - startTime) / (double)CLOCKS_PER_SEC < time_limit) {
		time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);



		MDD_arcs.resize(inst->n_SeatRow);
		MDDArc* new_arc;
		for (int cap = 0; cap < inst->n_SeatRow; ++cap) {
			if (inst->no_seg_of_each_size[cap] > 0) {
				for (int jj = 1; jj < (int)MDD_nodes[cap].size(); jj++) {
					for (int aa = 0; aa < MDD_nodes[cap][jj]->parent_indexes.size(); ++aa) {
						new_arc = new MDDArc();
						new_arc->j = jj;
						new_arc->index = MDD_arcs[cap].size();
						new_arc->i = MDD_nodes[cap][jj]->parent_indexes[aa];
						new_arc->size_group[0] = MDD_nodes[cap][jj]->size_group[0];
						new_arc->size_group[1] = MDD_nodes[cap][jj]->size_group[1];
						new_arc->n_party_selected = MDD_nodes[cap][jj]->n_party_selected[aa];
						new_arc->value = MDD_nodes[cap][jj]->value[aa];
						MDD_nodes[cap][new_arc->j]->input_arcs.push_back(new_arc->index);
						MDD_nodes[cap][new_arc->i]->output_arcs.push_back(new_arc->index);

						MDD_arcs[cap].push_back(new_arc);
					}
				}

				//adding Arcs (i,terminal)
				for (int ii = 0; ii < (int)T_conn_nodes_index[cap].size(); ii++) {
					new_arc = new MDDArc();
					new_arc->index = MDD_arcs[cap].size();
					new_arc->j = MDD_nodes[cap].size();
					new_arc->i = T_conn_nodes_index[cap][ii];
					new_arc->size_group[0] = -1;
					new_arc->size_group[1] = -1;
					new_arc->n_party_selected = 0.0;
					new_arc->value = false;

					MDD_nodes[cap][new_arc->i]->output_arcs.push_back(new_arc->index);

					MDD_arcs[cap].push_back(new_arc);
				}

			}
		}

		////Printing MDDs' arcs
		//for (int cap = 0; cap < inst->n_SeatRow; ++cap) {
		//	if (inst->no_seg_of_each_size[cap] > 0) {
		//		cout << "printing arcs in a segment with capacity " << cap + 1 << " :" << endl;
		//		for (int a = 0; a < (int)MDD_arcs[cap].size(); ++a) {
		//			print_arc(MDD_arcs[cap][a]);
		//		}
		//	}
		//}

		cout << "Done creating MDD" << endl;
	}

	if (double(clock() - startTime) / (double)CLOCKS_PER_SEC < time_limit) {
		time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
	}
	else {
		time_left = 0;
	}
	time_limit = time_left;

}

void MDDSolver::solve_net_model(double & time_limit, ofstream& file1) {

	////creating MDD nodes and arcs
	//create_MDDs(time_limit, file1);


	if (time_limit > PRECISION) {

		cout << "=============================================" << endl << "Solving Arc Based Network Flow Model:" << endl << "=============================================";

		GRBEnv env;

		try {

			GRBModel model(env);

			//Defining and adding arc variables
			int count = 0;
			int n_const = 0;
			int n_arcs = 0;
			int n_nodes = 0;



			int n_MDDarcs = MDD_arcs[inst->n_SeatRow - 1].size();
			n_arcs = n_MDDarcs;
			n_nodes = MDD_nodes[inst->n_SeatRow - 1].size();



			GRBVar* x;																		//x[s][a]=1 : if arc "a" in the MDD of segment "s" is selected
			x = new GRBVar[n_MDDarcs];
			for (int a = 0; a < n_MDDarcs; ++a) {
				x[a] = model.addVar(0.0, inst->n_seg, 0.0, GRB_INTEGER);
				count++;
			}

			model.update();

			//adding constraints

			//sum_[(i,j) \in A(s)] x_s,(i,j) - sum_[(j,i) \in A(s)] x_s,(j,i) = 0
			for (int k = 1; k < n_nodes; ++k) {
				GRBLinExpr equ = 0.0;
				for (int in = 0; in < MDD_nodes[inst->n_SeatRow - 1][k]->input_arcs.size(); ++in) {
					equ += x[MDD_nodes[inst->n_SeatRow - 1][k]->input_arcs[in]];
				}
				for (int out = 0; out < MDD_nodes[inst->n_SeatRow - 1][k]->output_arcs.size(); ++out) {
					equ -= x[MDD_nodes[inst->n_SeatRow - 1][k]->output_arcs[out]];
				}
				model.addConstr(equ == 0.0);
				n_const++;
			}

			//x_s,(0,1) + x_s,(0,2) + ... = 1
			GRBLinExpr equ2 = 0.0;
			for (int out = 0; out < MDD_nodes[inst->n_SeatRow - 1][0]->output_arcs.size(); ++out) {
				equ2 += x[MDD_nodes[inst->n_SeatRow - 1][0]->output_arcs[out]];
			}
			model.addConstr(equ2 == inst->n_seg);
			n_const++;

			for (int cap = 0; cap < inst->n_SeatRow; ++cap) {
				if (inst->no_seg_of_each_size[cap] > 0) {

					GRBLinExpr equ3 = 0.0;
					for (int node = 0; node < T_conn_nodes_index[inst->n_SeatRow - 1].size(); ++node) {
						if ((inst->n_SeatRow - MDD_nodes[inst->n_SeatRow - 1][T_conn_nodes_index[inst->n_SeatRow - 1][node]]->capacity) <= cap + 1) {
							for (int aa = 0; aa < MDD_nodes[inst->n_SeatRow - 1][T_conn_nodes_index[inst->n_SeatRow - 1][node]]->input_arcs.size(); ++aa) {
								equ3 += x[MDD_nodes[inst->n_SeatRow - 1][T_conn_nodes_index[inst->n_SeatRow - 1][node]]->input_arcs[aa]];
							}
						}
					}

					int n_bin_of_size_cap_or_less = 0.0;
					for (int cap2 = 0; cap2 <= cap; ++cap2) {
						n_bin_of_size_cap_or_less += inst->no_seg_of_each_size[cap2];
					}

					//cout << " n_bin_of_size " << cap << " or less: " << n_bin_of_size_cap_or_less << endl;

					model.addConstr(equ3 >= n_bin_of_size_cap_or_less);

				}
			}


			////sum_[(i,T) \in A(s)] x_s,(i,T) = 1
			//for (int s = 0; s < inst->n_seg; s++) {
			//	int n_MDDnodes = MDD_nodes[inst->no_seat_in_seg[s] - 1].size();
			//	int n_MDDarcs = MDD_arcs[inst->no_seat_in_seg[s] - 1].size();
			//	GRBLinExpr equ = 0.0;
			//	for (int a = 0; a < n_MDDarcs; a++) {
			//		if (MDD_arcs[inst->no_seat_in_seg[s] - 1][a]->j == n_MDDnodes) {
			//			equ += x[s][a];
			//		}
			//	}
			//	model.addConstr(equ == 1.0);
			//}

			//total number of parties selected from each party size n from group g is equal to the total existing number of parties from g and of size n
			for (int g = 0; g < inst->n_subGroup; ++g) {
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					if (inst->u[g][nn] > 0) {
						GRBLinExpr party_equ = 0.0;

						for (int a = 0; a < n_MDDarcs; ++a) {
							if (MDD_arcs[inst->n_SeatRow - 1][a]->size_group[0] == nn + 1 && MDD_arcs[inst->n_SeatRow - 1][a]->size_group[1] == g) {
								party_equ += MDD_arcs[inst->n_SeatRow - 1][a]->n_party_selected * x[a];
							}
						}

						model.addConstr(party_equ == (int)inst->u[g][nn]);
						n_const++;
					}
				}
			}

			///////////////////////////
			//Adding objective function 
			GRBQuadExpr objective = 0.0;
			for (int a = 0; a < n_MDDarcs; ++a) {
				objective += MDD_arcs[inst->n_SeatRow - 1][a]->value * x[a];
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
			cout << "Time Taken: " << model.get(GRB_DoubleAttr_Runtime) << endl;

			//file1 << n_nodes << ",";			// number of nodes per diagram
			//file1 << n_arcs << ",";				// number of arcs per diagram
			//file1 << count << ",";				// number of variables 
			//file1 << n_const << ",";			// number of constraints
			//file1 << counter_variables << ",";				// number of variables 
			//file1 << counter_constraints << ",";			// number of constraints


			file1 << model.get(GRB_DoubleAttr_ObjBound) << ",";
			file1 << model.get(GRB_DoubleAttr_ObjVal) << ",";
			file1 << model.get(GRB_DoubleAttr_Runtime) << ",";



			// Output line for script
			double gap = -1.0;
			double ub = model.get(GRB_DoubleAttr_ObjVal);
			double lb = model.get(GRB_DoubleAttr_ObjBound);
			double runtime = model.get(GRB_DoubleAttr_Runtime);

			if (model.get(GRB_DoubleAttr_ObjVal) < 1e6)
				gap = abs((ub - lb) / ub);

			cout << "Results," <<
				runtime << "," <<
				lb << "," <<
				ub << "," <<
				gap << "\n";

			////// This gets you the relaxation bound
			//cout << "Relaxation Bound: " << model.get(GRB_DoubleAttr_ObjBound) << endl;
			////// This gets you the best known primal solution value
			//cout << "Best Known Objective : " << model.get(GRB_DoubleAttr_ObjVal) << endl;

			/*for (int s = 0; s < inst->n_seg; s++) {
			int n_MDDarcs = MDD_arcs[inst->no_seat_in_seg[s] - 1].size();
			cout << "segment " << s << " :" << endl;
			for (int a = 0; a < n_MDDarcs; a++) {
			if (x[s][a].get(GRB_DoubleAttr_X) > 0.9) {
			print_arc(MDD_arcs[inst->no_seat_in_seg[s] - 1][a]);
			}
			}
			}*/

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

void MDDSolver::merge_nodes(bool & merged, int seg_capacity, vector<vector<MDDNode* >> & MDD_nodes, vector<vector<MDDNode* >> & MDD_nodes_left_to_branch_on, MDDNode* new_MDD_node) {
	for (int n1 = 0; n1 < MDD_nodes_left_to_branch_on[seg_capacity - 1].size(); ++n1) {
		if (MDD_nodes_left_to_branch_on[seg_capacity - 1][n1]->capacity == new_MDD_node->capacity		 &&			MDD_nodes_left_to_branch_on[seg_capacity - 1][n1]->size_group[0] == new_MDD_node->size_group[0] && MDD_nodes_left_to_branch_on[seg_capacity - 1][n1]->size_group[1] == new_MDD_node->size_group[1] && MDD_nodes_left_to_branch_on[seg_capacity - 1][n1]->if_party_from_same_group_selected == new_MDD_node->if_party_from_same_group_selected) {
			merged = true;
			int index = MDD_nodes_left_to_branch_on[seg_capacity - 1][n1]->index;
			MDD_nodes[seg_capacity - 1][index]->parent_indexes.push_back(new_MDD_node->parent_indexes[0]);
			MDD_nodes[seg_capacity - 1][index]->n_party_selected.push_back(new_MDD_node->n_party_selected[0]);
			MDD_nodes[seg_capacity - 1][index]->value.push_back(new_MDD_node->value[0]);

			MDD_nodes_left_to_branch_on[seg_capacity - 1][n1]->parent_indexes.push_back(new_MDD_node->parent_indexes[0]);
			MDD_nodes_left_to_branch_on[seg_capacity - 1][n1]->n_party_selected.push_back(new_MDD_node->n_party_selected[0]);
			MDD_nodes_left_to_branch_on[seg_capacity - 1][n1]->value.push_back(new_MDD_node->value[0]);
			break;
		}
	}
}

void MDDSolver::print_node(MDDNode* node) {

	cout << "*************MDD NODE! ************" << endl;
	cout << "\tIndex of MDD node: " << node->index << endl;
	cout << "\tParty_size of MDD node: " << node->size_group[0] << endl;
	cout << "\tGroup of MDD node: " << node->size_group[1] << endl;
	cout << "\tCapacity of MDD node: " << node->capacity << endl;
	cout << "\tSegemnt capacity of MDD node: " << node->segment_capacity << endl;
	cout << "**************************************" << endl;

}

void MDDSolver::print_arc(MDDArc* arc) {

	cout << "*************MDD NODE! ************" << endl;
	cout << "\tArc_index: " << arc->index << endl;
	cout << "\tarc: (" << arc->i << "," << arc->j << ")" << endl;
	cout << "\tn_party_selected: " << arc->n_party_selected << endl;
	cout << "\t(party_size,party_group): (" << arc->size_group[0] << " , " << arc->size_group[1] << " )" << endl;
	cout << "\tarc value: " << arc->value << endl;
	cout << "**************************************" << endl;

}

void MDDSolver::create_first_MDD_node(int segment_capacity) {

	MDDNode* first_MDD_node;
	first_MDD_node = new MDDNode();

	first_MDD_node->index = 0;
	first_MDD_node->parent_indexes.resize(1, -1);
	first_MDD_node->segment_capacity = segment_capacity;
	first_MDD_node->size_group[0] = 0;
	first_MDD_node->size_group[1] = 0;

	first_MDD_node->capacity = segment_capacity;
	first_MDD_node->n_party_selected.resize(1, 0);
	first_MDD_node->value.resize(1, 0);


	MDD_nodes[segment_capacity - 1].push_back(first_MDD_node);
	MDD_nodes_left_to_branch[segment_capacity - 1].push_back(first_MDD_node);

}
