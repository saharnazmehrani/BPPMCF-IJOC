#include "BP_SS.h"
GRBEnv env;
bool compare_nodes(SearchNode* node1, SearchNode* node2) {
	return (node1->parent_lower_bound < node2->parent_lower_bound);
}

void BPSolver::pick_branch_variable(SearchNode* node, double & time_limit, ofstream& file1, int & branch_var_index, vector<int> & group_branch_index_and_value) {
	//cout << "=============================================" << endl << "Pick Branching Variable:" << endl << "=============================================";
	//GRBEnv env;
	char varname[256];

	try {
		GRBModel model(env);

		// Define variables

		GRBVar* x;
		x = new GRBVar[(int)variables.size()];
		for (int v = 0; v < (int)variables.size(); ++v) {
			//sprintf(varname, "x[%d]", v);
			x[v] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
		}

		model.update();         // required before you can use the variables

								// adding constraints

								// \sum_c m[c][g][nn] x_c == u[g][nn] for all group g and size nn 
		vector<vector<GRBConstr >>  group_size_constraints(inst->n_subGroup);
		for (int g = 0; g < inst->n_subGroup; ++g) {
			group_size_constraints[g].resize(inst->max_party_size);
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				GRBLinExpr group_size_con = 0.0;
				for (int c = 0; c < (int)variables.size(); ++c) {
					group_size_con += variables[c]->assigned_group_size[g][nn] * x[c];
				}
				group_size_constraints[g][nn] = model.addConstr(group_size_con >= inst->u[g][nn]);
			}
		}

		// number of selected columns from each capacity must be <= number of segemnets available of that capacity
		vector<GRBConstr>  capacity_constraints(inst->n_SeatRow);
		for (int cap = 0; cap < inst->n_SeatRow; ++cap) {
			GRBLinExpr capacity_con = 0.0;
			for (int c = 0; c < (int)variables.size(); ++c) {
				if (variables[c]->capacity == cap + 1) {
					capacity_con += x[c];
				}
			}
			capacity_constraints[cap] = model.addConstr(capacity_con <= inst->no_seg_of_each_size[cap]);
		}

		// LB and UB constraint for number of segments assigned to each group
		//vector<GRBConstr >  seg_LB_constraints(inst->n_subGroup);
		vector<GRBConstr >  seg_UB_constraints(inst->n_subGroup);
		vector<GRBConstr >  new_LB_constraints(inst->n_subGroup);
		for (int g = 0; g < inst->n_subGroup; ++g) {
			GRBLinExpr n_seg_assigned_to_g = 0.0;
			for (int c = 0; c < (int)variables.size(); ++c) {
				int n_party = 0.0;
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					if (variables[c]->assigned_group_size[g][nn] > 0.9) {
						n_party++;
					}
				}
				if (n_party > 0.9) {
					n_seg_assigned_to_g += x[c];
				}
			}
			//seg_LB_constraints[g] = model.addConstr(n_seg_assigned_to_g >= node->LB_seg_no_for_group[g]);
			seg_UB_constraints[g] = model.addConstr(n_seg_assigned_to_g <= node->UB_seg_no_for_group[g]);

			double LB_n_bin = inst->nn[g] / inst->n_SeatRow;
			if (LB_n_bin - floor(LB_n_bin) >= 0.001) {
				LB_n_bin = floor(LB_n_bin) + 1;
			}
			else {
				LB_n_bin = floor(LB_n_bin);
			}
			double new_LB = LB_n_bin;
			if (inst->LB_MT1[g] > new_LB) {
				new_LB = inst->LB_MT1[g];
			}
			if (inst->LB_MT2[g] > new_LB) {
				new_LB = inst->LB_MT2[g];
			}
			if (node->LB_seg_no_for_group[g] > new_LB) {
				new_LB = node->LB_seg_no_for_group[g];
			}
			new_LB_constraints[g] = model.addConstr(n_seg_assigned_to_g >= new_LB);
		}

		// for these two constraints, you might want to using vector iterators, and
		// find the position in the vector.  We will, however, make sure the index matches
		// the relateive position
		for (int i = 0; i < (int)node->variables_set_to_one.size(); ++i) {
			model.addConstr(x[node->variables_set_to_one[i]->index] == 1);
		}

		for (int i = 0; i < (int)node->variables_set_to_zero.size(); ++i) {
			model.addConstr(x[node->variables_set_to_zero[i]->index] == 0);
		}

		// objective function
		GRBLinExpr objective = 0.0;
		for (int c = 0; c < (int)variables.size(); ++c) {
			objective += variables[c]->value * x[c];
		}

		model.addConstr(objective >= global_lower_bound);
		model.addConstr(objective <= global_upper_bound);

		model.setObjective(objective, GRB_MINIMIZE);

		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_DoubleParam_TimeLimit, time_limit);
		model.getEnv().set(GRB_IntParam_OutputFlag, 0);
		// model.getEnv().set(GRB_IntParam_Cuts,0);
		//            model.getEnv().set(GRB_IntParam_PreCrush,1);
		//model.set(GRB_IntParam_PreCrush,1);
		//            model.getEnv().set(GRB_IntParam_DualReductions, 0);
		//            model.getEnv().set(GRB_IntParam_LazyConstraints, 1);

		model.update();
		model.optimize(); //solves the problem

		if (model.get(GRB_IntAttr_Status) == 2) {
			//cout << "Solution Status: " << model.get(GRB_IntAttr_Status) << endl;
			//cout << "Time Taken: " << model.get(GRB_DoubleAttr_Runtime) << endl;
			// This gets you the relaxation bound
			/*//cout << "Relaxation Bound: " << model.get(GRB_DoubleAttr_ObjBound) << endl;*/

			/*//cout << "x[c]: " << endl;
			for (int c = 0; c < variables.size(); ++c) {
			//cout << x[c].get(GRB_DoubleAttr_X) << "\t";
			}
			//cout << endl;*/


			double most_fractional_value = 0.0;
			//cout << "Final LP value: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
			//file1 << model.get(GRB_DoubleAttr_ObjVal) << endl;
			node->LP_lower_bound = model.get(GRB_DoubleAttr_ObjVal);
			if (node->LP_lower_bound < node->parent_lower_bound - PRECISION) {
				cout << "parent lower bound: " << node->parent_lower_bound << endl;
				cout << "current lower bound: " << node->LP_lower_bound << endl;
				cout << "ERROR IN BRANCH VARIABLE: " << endl;
				exit(1);
			}
			//cout << "Global upper bound: " << global_upper_bound << endl;
			// If the objective value is greater than the global upper bound, we don't have to branch
			// This way, branch_var_index will remain at -1
			if (model.get(GRB_DoubleAttr_ObjVal) > global_upper_bound - 0.99) {
				//cout << "DECIDED NOT TO CHECK " << endl;
				// don't pick variable to branch on
			}
			else {
				// pick variable that is most fractional
				// First checking number of segments assigned to each group
				double current_frac_value = 0.0;
				for (int g = 0; g < inst->n_subGroup; ++g) {
					double n_seg_assigned_to_g = 0.0;
					for (int c = 0; c < (int)variables.size(); ++c) {
						int n_party = 0;
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							if (variables[c]->assigned_group_size[g][nn] > 0.9) { //&& inst->u[g][nn] > 0.9
								n_party++;
								/*n_seg_assigned_to_g += x[c].get(GRB_DoubleAttr_X);
								break;*/
							}
						}
						if (n_party > 0.9) {
							n_seg_assigned_to_g += x[c].get(GRB_DoubleAttr_X);
						}
					}

					if (n_seg_assigned_to_g - floor(n_seg_assigned_to_g) <= 0.5) {
						current_frac_value = n_seg_assigned_to_g - floor(n_seg_assigned_to_g);
					}
					else {
						current_frac_value = 1.0 - (n_seg_assigned_to_g - floor(n_seg_assigned_to_g));
					}
					////cout << "Current frac value: " << current_frac_value << endl;
					if (current_frac_value > most_fractional_value + PRECISION) {
						group_branch_index_and_value[0] = g;
						group_branch_index_and_value[1] = floor(n_seg_assigned_to_g);
						most_fractional_value = current_frac_value;
					}

				}

				// Then Check fractional values for x_c
				if (most_fractional_value == 0.0) {
					for (int c = 0; c < (int)variables.size(); ++c) {
						double current_frac_value = 0.0;
						if (x[c].get(GRB_DoubleAttr_X) <= 0.5) {
							current_frac_value = x[c].get(GRB_DoubleAttr_X);
						}
						else {
							current_frac_value = 1.0 - x[c].get(GRB_DoubleAttr_X);
						}
						////cout << "Current frac value: " << current_frac_value << endl;
						if (current_frac_value > most_fractional_value + PRECISION) {
							branch_var_index = c;
							most_fractional_value = current_frac_value;
						}
					}
				}

				//cout << "most_fractional_value: " << most_fractional_value << endl;

			} // end of else
		} // end of feasibility check

		delete x;

	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		exit(1);
	}
	catch (...) {
		cout << "Other error ... " << endl;
		exit(1);
	}

}

class callback_pricing : public GRBCallback {

public:
	// put objects required

	GRBVar* y;
	GRBVar*** x_callback;
	vector <Variable*> variables;
	Stadium_Seating* inst;
	int capacity;
	//callback_vrpspd(int* _n_solutions_found,QuadKnap* _quad_knap, vector<vector<double > > * _quad_obj_matrix, GRBVar * _t, GRBVar* _x){
	callback_pricing(GRBVar* _y, GRBVar*** _x, vector <Variable*> _variables, Stadium_Seating* _inst, int _capacity) {

		////cout << "We are at a callback .. " << endl;
		y = _y;
		x_callback = _x;
		inst = _inst;
		variables = _variables;
		capacity = _capacity;
	};

protected:

	void callback() {
		try {

			if (where == GRB_CB_MIPSOL) {


				////// ADD CUTS BELOW THIS
				vector <vector<vector<int>>> X(inst->n_subGroup);
				for (int g = 0; g < inst->n_subGroup; ++g) {
					X[g].resize(inst->max_party_size);
					for (int nn = 0; nn < inst->max_party_size; ++nn) {
						X[g][nn].resize(inst->u[g][nn] + 1);
						for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
							if (getSolution(x_callback[g][nn][i]) >= 0.5) {
								X[g][nn][i] = 1.0;
							}
						}
					}
				}


				//Avoid constraint
				//cout << (int)variables.size() << endl;
				for (int c = 0; c < (int)variables.size(); ++c) {

					if (variables[c]->capacity == capacity) {
						int check = 0;
						int counter = 0;
						for (int g = 0; g < inst->n_subGroup; ++g) {
							for (int nn = 0; nn < inst->max_party_size; ++nn) {
								for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
									if (variables[c]->assigned_group_size[g][nn] == i) {
										check += X[g][nn][i];
										counter++;
									}
									else {
										check -= X[g][nn][i];
									}
								}
							}
						}
						if (check > counter - 1) {
							GRBLinExpr avoid_cut = 0.0;
							int count = 0;
							for (int g = 0; g < inst->n_subGroup; ++g) {
								for (int nn = 0; nn < inst->max_party_size; ++nn) {
									for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
										if (variables[c]->assigned_group_size[g][nn] == i) {
											avoid_cut += x_callback[g][nn][i];
											count++;
										}
										else {
											avoid_cut -= x_callback[g][nn][i];
										}
									}
								}
							}
							addLazy(avoid_cut, GRB_LESS_EQUAL, count - 1);
						}
					}
				}

				///// ADD CUTS ABOVE THIS
			}




		}
		catch (GRBException e) {
			cout << "Error number: " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
			exit(1);
		}
		catch (...) {
			cout << "Error during callback" << endl;
		}

	}
};

void BPSolver::node_feasibility_check(SearchNode* node, double & time_limit, bool & node_is_feasible, int n_node_explored) {
	//cout << "=============================================" << endl << "Checking feasibility of a search node:" << endl << "=============================================";

	///////checking whether the parent-LB for the node is greater than the global upper bound
	if (node->parent_lower_bound + PRECISION > global_upper_bound) {
		search_nodes.erase(search_nodes.begin());
		cout << "Node is pruned because of large parent-LB" << endl;
	}
	else {
		////////// sloving RMP without objective ////////////
		GRBEnv env;
		char varname[256];
		try {
			GRBModel model(env);

			// Define variables

			GRBVar* x;
			x = new GRBVar[(int)variables.size()];
			for (int v = 0; v < (int)variables.size(); ++v) {
				//sprintf(varname, "x[%d]", v);
				x[v] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
			}

			model.update();         // required before you can use the variables

									// adding constraints

									// \sum_c m[c][g][nn] x_c == u[g][nn] for all group g and size nn 
			vector<vector<GRBConstr >>  group_size_constraints(inst->n_subGroup);
			for (int g = 0; g < inst->n_subGroup; ++g) {
				group_size_constraints[g].resize(inst->max_party_size);
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					GRBLinExpr group_size_con = 0.0;
					for (int c = 0; c < (int)variables.size(); ++c) {
						group_size_con += variables[c]->assigned_group_size[g][nn] * x[c];
					}
					group_size_constraints[g][nn] = model.addConstr(group_size_con >= inst->u[g][nn]);
				}
			}


			// number of selected columns from each capacity must be <= number of segemnets available of that capacity
			vector<GRBConstr>  capacity_constraints(inst->n_SeatRow);
			for (int cap = 0; cap < inst->n_SeatRow; ++cap) {
				GRBLinExpr capacity_con = 0.0;
				for (int c = 0; c < (int)variables.size(); ++c) {
					if (variables[c]->capacity == cap + 1) {
						capacity_con += x[c];
					}
				}
				capacity_constraints[cap] = model.addConstr(capacity_con <= inst->no_seg_of_each_size[cap]);
			}
			// LB and UB constraint for number of segments assigned to each group
			//vector<GRBConstr >  seg_LB_constraints(inst->n_subGroup);
			vector<GRBConstr >  seg_UB_constraints(inst->n_subGroup);
			vector<GRBConstr >  new_LB_constraints(inst->n_subGroup);
			for (int g = 0; g < inst->n_subGroup; ++g) {
				GRBLinExpr n_seg_assigned_to_g = 0.0;
				for (int c = 0; c < (int)variables.size(); ++c) {
					int n_party = 0.0;
					for (int nn = 0; nn < inst->max_party_size; ++nn) {
						if (variables[c]->assigned_group_size[g][nn] > 0.9) {
							n_party++;
						}
					}
					if (n_party > 0.9) {
						n_seg_assigned_to_g += x[c];
					}
				}
				//seg_LB_constraints[g] = model.addConstr(n_seg_assigned_to_g >= node->LB_seg_no_for_group[g]);
				seg_UB_constraints[g] = model.addConstr(n_seg_assigned_to_g <= node->UB_seg_no_for_group[g]);
				double LB_n_bin = inst->nn[g] / inst->n_SeatRow;
				if (LB_n_bin - floor(LB_n_bin) >= 0.001) {
					LB_n_bin = floor(LB_n_bin) + 1;
				}
				else {
					LB_n_bin = floor(LB_n_bin);
				}
				double new_LB = LB_n_bin;
				if (inst->LB_MT1[g] > new_LB) {
					new_LB = inst->LB_MT1[g];
				}
				if (inst->LB_MT2[g] > new_LB) {
					new_LB = inst->LB_MT2[g];
				}
				if (node->LB_seg_no_for_group[g] > new_LB) {
					new_LB = node->LB_seg_no_for_group[g];
				}
				new_LB_constraints[g] = model.addConstr(n_seg_assigned_to_g >= new_LB);
			}

			// for these two constraints, you might want to using vector iterators, and
			// find the position in the vector.  We will, however, make sure the index matches
			// the relateive position
			for (int i = 0; i < (int)node->variables_set_to_one.size(); ++i) {
				model.addConstr(x[node->variables_set_to_one[i]->index] == 1);
			}

			for (int i = 0; i < (int)node->variables_set_to_zero.size(); ++i) {
				model.addConstr(x[node->variables_set_to_zero[i]->index] == 0);
			}

			// objective function
			GRBLinExpr objective = 0.0;

			model.setObjective(objective, GRB_MINIMIZE);

			model.getEnv().set(GRB_IntParam_Threads, 1);
			model.getEnv().set(GRB_DoubleParam_TimeLimit, time_limit);
			model.getEnv().set(GRB_IntParam_OutputFlag, 0);
			// model.getEnv().set(GRB_IntParam_Cuts,0);
			//            model.getEnv().set(GRB_IntParam_PreCrush,1);
			//model.set(GRB_IntParam_PreCrush,1);
			//            model.getEnv().set(GRB_IntParam_DualReductions, 0);
			//            model.getEnv().set(GRB_IntParam_LazyConstraints, 1);

			model.update();
			model.optimize(); //solves the problem
			if (model.get(GRB_IntAttr_Status) == 2) {
				node_is_feasible = true;
				//cout << "node is feasible" << endl;
				//search_nodes.erase(search_nodes.begin());
				////search_nodes.pop_back();
				//if (search_nodes.size() > 0) {
				//	//search_nodes.insert(search_nodes.begin(), node);
				//	search_nodes.push_back(node);
				//	//cout << "node number " << node->index << " is infeasible, and will be explored later" << endl;
				//}
				//else {
				//	//cout << "last node is cutted" << endl;
				//}
			}
			//if (model.get(GRB_IntAttr_Status) == 3) {
			//	//cout << "node is not feasible" << endl;
			//	search_nodes.erase(search_nodes.begin());
			//	search_nodes.push_back(node);
			//}

			delete x;

		}
		catch (GRBException e) {
			cout << "Feasibility check: Error number: " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
			exit(1);
		}
		catch (...) {
			cout << "Feasibility check: Other error ... " << endl;
			exit(1);
		}


		//////// If node is not feasible, check MIP for feasibility ///////////////////////////////
		if (!node_is_feasible) {
			char varname[256];

			GRBEnv env2;
			try {
				GRBModel model(env2);

				//Defining and adding variables
				GRBVar**** x;																		//x[g][n][i] : Number of party of size n from group g assigned to this column is equal to i
				x = new GRBVar***[inst->n_seg];
				for (int s = 0; s < inst->n_seg; s++) {
					x[s] = new GRBVar**[inst->n_subGroup];
					for (int g = 0; g < inst->n_subGroup; g++) {
						x[s][g] = new GRBVar*[inst->max_party_size];
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							x[s][g][nn] = new GRBVar[inst->u[g][nn] + 1];
							for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
								x[s][g][nn][i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
							}
						}
					}
				}

				GRBVar** y;																		//y[g][seg] is equal to one if group "g" has some party in segment "seg", and zero otherwise
				y = new GRBVar*[inst->n_subGroup];
				for (int g = 0; g < inst->n_subGroup; ++g) {
					y[g] = new GRBVar[inst->n_seg];
					for (int seg = 0; seg < inst->n_seg; ++seg) {
						y[g][seg] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					}
				}

				model.update();
				//adding constraints

				//sum_i x[s][g][nn][i] = 1 for each segment s, group g, and party size nn
				for (int seg = 0; seg < inst->n_seg; ++seg) {
					for (int g = 0; g < inst->n_subGroup; g++) {
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							GRBLinExpr equ = 0.0;
							for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
								equ += x[seg][g][nn][i];
							}
							model.addConstr(equ == 1);
						}
					}
				}

				//sum_seg,i i*x[g][nn][i] == number of parties of size nn+1 from group g
				for (int g = 0; g < inst->n_subGroup; g++) {
					for (int nn = 0; nn < inst->max_party_size; ++nn) {
						GRBLinExpr equ = 0.0;
						for (int seg = 0; seg < inst->n_seg; ++seg) {
							for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
								equ += i * x[seg][g][nn][i];
							}
						}
						model.addConstr(equ == inst->u[g][nn]);
					}
				}

				//x[g][nn]<= y[g] for all group g and size nn
				for (int s = 0; s < inst->n_seg; ++s) {
					for (int g = 0; g < inst->n_subGroup; ++g) {
						GRBLinExpr Yequ = 0.0;
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							model.addConstr(1 - x[s][g][nn][0] <= y[g][s]);
							Yequ += (1 - x[s][g][nn][0]);
						}
						model.addConstr(y[g][s] <= Yequ);
					}
				}


				//Knapsack constraint
				for (int s = 0; s < inst->n_seg; ++s) {
					GRBLinExpr cap_equ = 0.0;
					for (int g = 0; g < inst->n_subGroup; ++g) {
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
								cap_equ += (nn + 1) * i * x[s][g][nn][i];
							}
						}
					}
					model.addConstr(cap_equ <= inst->no_seat_in_seg[s]);
				}

				//constraint for variables set to one
				vector<bool> seg_fixed(inst->n_seg, false);		//if a segment is fixed to one in the current node, set_of_seg_fixed[seg] takes "true"
				for (int i = 0; i < (int)node->variables_set_to_one.size(); ++i) {
					for (int seg = 0; seg < inst->n_seg; ++seg) {
						if (!seg_fixed[seg] && inst->no_seat_in_seg[seg] == variables[node->variables_set_to_one[i]->index]->capacity) {
							for (int g = 0; g < inst->n_subGroup; ++g) {
								for (int nn = 0; nn < inst->max_party_size; ++nn) {
									if (variables[node->variables_set_to_one[i]->index]->assigned_group_size[g][nn] > 0.9) {
										model.addConstr(x[seg][g][nn][variables[node->variables_set_to_one[i]->index]->assigned_group_size[g][nn]] == 1);
									}
									else {
										model.addConstr(x[seg][g][nn][0] == 1);
									}
								}
							}
							seg_fixed[seg] = true;
						}
					}
					node->variables_set_to_one[i];
				}

				// LB and UB constraint for number of segments assigned to each group
				for (int g = 0; g < inst->n_subGroup; ++g) {
					GRBLinExpr n_seg_assigned_to_g = 0.0;
					for (int seg = 0; seg < inst->n_seg; ++seg) {
						n_seg_assigned_to_g += y[g][seg];
					}
					//model.addConstr(n_seg_assigned_to_g >= node->LB_seg_no_for_group[g]);
					model.addConstr(n_seg_assigned_to_g <= node->UB_seg_no_for_group[g]);
					double LB_n_bin = inst->nn[g] / inst->n_SeatRow;
					if (LB_n_bin - floor(LB_n_bin) >= 0.001) {
						LB_n_bin = floor(LB_n_bin) + 1;
					}
					else {
						LB_n_bin = floor(LB_n_bin);
					}
					double new_LB = LB_n_bin;
					if (inst->LB_MT1[g] > new_LB) {
						new_LB = inst->LB_MT1[g];
					}
					if (inst->LB_MT2[g] > new_LB) {
						new_LB = inst->LB_MT2[g];
					}
					if (node->LB_seg_no_for_group[g] > new_LB) {
						new_LB = node->LB_seg_no_for_group[g];
					}
					model.addConstr(n_seg_assigned_to_g >= new_LB);
				}


				///////////////////////////
				//Adding objective function 

				GRBLinExpr objective = 0.0;


				//adding objective[seg]
				//IloObjective obj;
				model.update();
				model.setObjective(objective, GRB_MINIMIZE);
				model.getEnv().set(GRB_IntParam_Threads, 1);
				model.getEnv().set(GRB_DoubleParam_TimeLimit, time_limit);
				model.getEnv().set(GRB_IntParam_OutputFlag, 0);
				model.getEnv().set(GRB_IntParam_PreCrush, 1);
				//model.getEnv().set(GRB_DoubleParam_Cutoff, PRECISION);
				model.getEnv().set(GRB_IntParam_DualReductions, 0);
				model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
				model.update();
				model.optimize();

				if (model.get(GRB_IntAttr_Status) == 3) {
					//if (n_node_explored < 25 || n_node_explored > 50) {
					search_nodes.erase(search_nodes.begin());
					//}
					/*else {
					search_nodes.pop_back();
					}*/
					cout << "node is pruned" << endl;
				}
				if (model.get(GRB_IntAttr_Status) == 2) {
					for (int seg = 0; seg < inst->n_seg; ++seg) {
						if (!seg_fixed[seg]) {
							Variable* p_var;
							p_var = new Variable();
							p_var->index = variables.size();
							p_var->assigned_group_size.resize(inst->n_subGroup);
							for (int g = 0; g < inst->n_subGroup; ++g) {
								p_var->assigned_group_size[g].resize(inst->max_party_size, 0);
								for (int nn = 0; nn < inst->max_party_size; ++nn) {
									for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
										if (x[seg][g][nn][i].get(GRB_DoubleAttr_X) > 0.9) {
											p_var->assigned_group_size[g][nn] = i;
										}
									}
								}
							}
							p_var->capacity = inst->no_seat_in_seg[seg];
							//computing value
							p_var->value = 0;
							for (int g = 0; g < inst->n_subGroup; ++g) {
								int n_party = 0;
								for (int nn = 0; nn < inst->max_party_size; ++nn) {
									if (x[seg][g][nn][0].get(GRB_DoubleAttr_X) < 0.1) {
										n_party++;
									}
								}
								if (n_party > 0.9) {
									p_var->value++;
								}
							}

							//adding the new variable		
							if (p_var->value > 0.5) {
								variables.push_back(p_var);
							}
						}
					}

					node_is_feasible = true;
				}

				delete x;
				delete y;

			} // end of try
			catch (GRBException e) {
				cout << "Feasibility check2: Error number: " << e.getErrorCode() << endl;
				cout << e.getMessage() << endl;
				cout << "time limit: " << time_limit << "\n";
				//exit(1);
			}
			catch (...) {
				cout << "Feasibility check2: Other error ... " << endl;
				exit(1);
			}
		}

	}



}

void BPSolver::Farkas_node_feasibility_check(SearchNode* node, double & time_limit, clock_t startTime, bool & node_is_feasible, int n_node_explored) {
	//cout << "=============================================" << endl << "Checking feasibility of a search node:" << endl << "=============================================";

	///////checking whether the parent-LB for the node is greater than the global upper bound
	if (node->parent_lower_bound + PRECISION > global_upper_bound) {
		search_nodes.erase(search_nodes.begin());
		cout << "Node is pruned because of large parent-LB" << endl;
	}
	else {

		bool RMP_infsbl = false; //indicating that RMP is infeasible

		while (!node_is_feasible && !RMP_infsbl) {

			////Setting the time_limit for the algorithm////
			double time_left;
			if (double(clock() - startTime) / (double)CLOCKS_PER_SEC >= time_limit) {
				cout << "time out and exit" << endl;
				break;
			}
			else {
				time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
			}

			////////// sloving RMP without objective ////////////

			vector<vector<double >> party_duals(inst->n_subGroup);	// one for every party
			for (int g = 0; g < inst->n_subGroup; ++g) {
				party_duals[g].resize(inst->max_party_size);
			}
			vector<double > capacity_duals(inst->n_SeatRow);	// one for every size from 1 to n_SeatRow
																//vector<double> seg_LB_duals(inst->n_subGroup);
			vector<double> seg_UB_duals(inst->n_subGroup);
			vector<double> new_LB_duals(inst->n_subGroup);

			//GRBEnv env;
			
			try {
				GRBModel model(env);

				// Define variables

				GRBVar* x;
				x = new GRBVar[(int)variables.size()];
				for (int v = 0; v < (int)variables.size(); ++v) {
					//sprintf(varname, "x[%d]", v);
					x[v] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
				}

				model.update();         // required before you can use the variables

										// adding constraints

										// \sum_c m[c][g][nn] x_c == u[g][nn] for all group g and size nn 
				vector<vector<GRBConstr >>  group_size_constraints(inst->n_subGroup);
				for (int g = 0; g < inst->n_subGroup; ++g) {
					group_size_constraints[g].resize(inst->max_party_size);
					for (int nn = 0; nn < inst->max_party_size; ++nn) {
						GRBLinExpr group_size_con = 0.0;
						for (int c = 0; c < (int)variables.size(); ++c) {
							group_size_con += variables[c]->assigned_group_size[g][nn] * x[c];
						}
						group_size_constraints[g][nn] = model.addConstr(group_size_con >= inst->u[g][nn]);
					}
				}


				// number of selected columns from each capacity must be <= number of segemnets available of that capacity
				vector<GRBConstr>  capacity_constraints(inst->n_SeatRow);
				for (int cap = 0; cap < inst->n_SeatRow; ++cap) {
					GRBLinExpr capacity_con = 0.0;
					for (int c = 0; c < (int)variables.size(); ++c) {
						if (variables[c]->capacity == cap + 1) {
							capacity_con += x[c];
						}
					}
					capacity_constraints[cap] = model.addConstr(capacity_con <= inst->no_seg_of_each_size[cap]);
				}
				// LB and UB constraint for number of segments assigned to each group
				//vector<GRBConstr >  seg_LB_constraints(inst->n_subGroup);
				vector<GRBConstr >  seg_UB_constraints(inst->n_subGroup);
				vector<GRBConstr >  new_LB_constraints(inst->n_subGroup);
				for (int g = 0; g < inst->n_subGroup; ++g) {
					GRBLinExpr n_seg_assigned_to_g = 0.0;
					for (int c = 0; c < (int)variables.size(); ++c) {
						int n_party = 0.0;
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							if (variables[c]->assigned_group_size[g][nn] > 0.9) {
								n_party++;
							}
						}
						if (n_party > 0.9) {
							n_seg_assigned_to_g += x[c];
						}
					}
					//seg_LB_constraints[g] = model.addConstr(n_seg_assigned_to_g >= node->LB_seg_no_for_group[g]);
					seg_UB_constraints[g] = model.addConstr(n_seg_assigned_to_g <= node->UB_seg_no_for_group[g]);
					double LB_n_bin = inst->nn[g] / inst->n_SeatRow;
					if (LB_n_bin - floor(LB_n_bin) >= 0.001) {
						LB_n_bin = floor(LB_n_bin) + 1;
					}
					else {
						LB_n_bin = floor(LB_n_bin);
					}
					double new_LB = LB_n_bin;
					if (inst->LB_MT1[g] > new_LB) {
						new_LB = inst->LB_MT1[g];
					}
					if (inst->LB_MT2[g] > new_LB) {
						new_LB = inst->LB_MT2[g];
					}
					if (node->LB_seg_no_for_group[g] > new_LB) {
						new_LB = node->LB_seg_no_for_group[g];
					}
					new_LB_constraints[g] = model.addConstr(n_seg_assigned_to_g >= new_LB);
				}



				// for these two constraints, you might want to using vector iterators, and
				// find the position in the vector.  We will, however, make sure the index matches
				// the relateive position
				for (int i = 0; i < (int)node->variables_set_to_one.size(); ++i) {
					model.addConstr(x[node->variables_set_to_one[i]->index] == 1);
				}

				for (int i = 0; i < (int)node->variables_set_to_zero.size(); ++i) {
					model.addConstr(x[node->variables_set_to_zero[i]->index] == 0);
				}

				// objective function
				GRBLinExpr objective = 0.0;

				model.setObjective(objective, GRB_MINIMIZE);

				model.getEnv().set(GRB_IntParam_Threads, 1);
				model.getEnv().set(GRB_DoubleParam_TimeLimit, time_limit);
				model.getEnv().set(GRB_IntParam_OutputFlag, 0);
				// model.getEnv().set(GRB_IntParam_Cuts,0);
				//            model.getEnv().set(GRB_IntParam_PreCrush,1);
				//model.set(GRB_IntParam_PreCrush,1);
				//            model.getEnv().set(GRB_IntParam_DualReductions, 0);
				//            model.getEnv().set(GRB_IntParam_LazyConstraints, 1);

				model.getEnv().set(GRB_IntParam_InfUnbdInfo, 1);


				model.update();
				model.optimize(); //solves the problem

				if (model.get(GRB_IntAttr_Status) == 2) {
					node_is_feasible = true;
					//cout << "node is feasible" << endl;
					//search_nodes.erase(search_nodes.begin());
					////search_nodes.pop_back();
					//if (search_nodes.size() > 0) {
					//	//search_nodes.insert(search_nodes.begin(), node);
					//	search_nodes.push_back(node);
					//	//cout << "node number " << node->index << " is infeasible, and will be explored later" << endl;
					//}
					//else {
					//	//cout << "last node is cutted" << endl;
					//}
				}
				//if (model.get(GRB_IntAttr_Status) == 3) {
				//	//cout << "node is not feasible" << endl;
				//	search_nodes.erase(search_nodes.begin());
				//	search_nodes.push_back(node);
				//}
				else {
					////cout << "party_duals:" << endl;
					for (int g = 0; g < inst->n_subGroup; ++g) {
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							party_duals[g][nn] = group_size_constraints[g][nn].get(GRB_DoubleAttr_FarkasDual);
						}
					}

					////cout << "capacity_duals:" << endl;
					for (int cap = 0; cap < inst->n_SeatRow; ++cap) {
						capacity_duals[cap] = capacity_constraints[cap].get(GRB_DoubleAttr_FarkasDual);
						////cout << capacity_duals[cap] << "\t";
					}
					////cout << endl;
					////cout << endl;

					////cout << "seg_LB_and_UB_duals" <<endl;
					for (int g = 0; g < inst->n_subGroup; ++g) {
						//seg_LB_duals[g] = seg_LB_constraints[g].get(GRB_DoubleAttr_FarkasDual);
						seg_UB_duals[g] = seg_UB_constraints[g].get(GRB_DoubleAttr_FarkasDual);
						new_LB_duals[g] = new_LB_constraints[g].get(GRB_DoubleAttr_FarkasDual);
					}
				}

				delete x;

			}
			catch (GRBException e) {
				cout << "Feasibility check: Error number: " << e.getErrorCode() << endl;
				cout << e.getMessage() << endl;
				exit(1);
			}
			catch (...) {
				cout << "Feasibility check: Other error ... " << endl;
				exit(1);
			}


			//////// If node is not feasible, check the pricing model to find new columns based on Farkas certificate to check infeasibility of the master problem ////////////
			if (!node_is_feasible) {
				//.... solve pricing problem with the new objective
				vector<vector<int>> xx; //xx[g][n] number of items of size n from group g, so far selected 
				xx.resize(inst->n_subGroup);
				for (int g = 0; g < inst->n_subGroup; ++g) {
					xx[g].resize(inst->max_party_size, 0);
				}

				int total_obj = 0;

				char varname[256];

				try {
					GRBModel model(env);

					//Defining and adding variables
					GRBVar*** x;																		//x[g][n][i] : Number of party of size n from group g assigned to this column is equal to i
					x = new GRBVar**[inst->n_subGroup];
					for (int g = 0; g < inst->n_subGroup; g++) {
						x[g] = new GRBVar*[inst->max_party_size];
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							x[g][nn] = new GRBVar[inst->u[g][nn] + 1];
							for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
								x[g][nn][i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
							}
						}
					}

					GRBVar* y;																		//y[g] is equal to one if group "g" has some party in this column, and zero otherwise
					y = new GRBVar[inst->n_subGroup];
					for (int g = 0; g < inst->n_subGroup; ++g) {
						y[g] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					}

					model.update();
					//adding constraints

					//sum_i x[g][nn][i] = 1
					for (int g = 0; g < inst->n_subGroup; g++) {
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							GRBLinExpr equ = 0.0;
							for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
								equ += x[g][nn][i];
							}
							model.addConstr(equ == 1);
						}
					}

					//x[g][nn]<= y[g] for all group g and size nn
					for (int g = 0; g < inst->n_subGroup; ++g) {
						//GRBLinExpr Yequ = 0.0;
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							model.addConstr(1 - x[g][nn][0] <= y[g]);
							//Yequ += (1 - x[g][nn][0]);
						}
						//model.addConstr(y[g] <= Yequ);
					}

					// \sum_{size} \sum_{i} x[g][nn]<= capacity * y[g] for all group g and size nn
					for (int g = 0; g < inst->n_subGroup; ++g) {
						GRBLinExpr equ_new = 0.0;
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
								equ_new += (nn + 1) * i * x[g][nn][i];
							}
						}
						model.addConstr(equ_new <= inst->n_SeatRow * y[g]);
					}

					//Knapsack constraint
					GRBLinExpr cap_equ = 0.0;
					for (int g = 0; g < inst->n_subGroup; ++g) {
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
								cap_equ += (nn + 1) * i * x[g][nn][i];
							}
						}
					}
					model.addConstr(cap_equ <= inst->n_SeatRow);


					///////////////////////////
					//Adding objective function 
					GRBLinExpr objective = 0.0;

					for (int g = 0; g < inst->n_subGroup; ++g) {
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							for (int i = 1; i < inst->u[g][nn] + 1; ++i) {
								objective += party_duals[g][nn] * i * x[g][nn][i];
							}
						}
					}
					objective += capacity_duals[inst->n_SeatRow - 1];
					for (int g = 0; g < inst->n_subGroup; ++g) {
						objective += (seg_UB_duals[g] + new_LB_duals[g]) * y[g];
					}

					// Using Callback to add the avoid constraints (i.e. to avoid from regenerating the existing columns)
					callback_pricing cb = callback_pricing(y, x, variables, inst, inst->n_SeatRow);
					model.setCallback(&cb);


					//calculating the total number of items selected
					GRBLinExpr Total_item = 0.0;
					for (int g = 0; g < inst->n_subGroup; ++g) {
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							for (int i = 1; i < inst->u[g][nn] + 1; ++i) {
								Total_item += (nn + 1) * i * x[g][nn][i];
							}
						}
					}

					double last_RC = 0.0;
					for (int iter = 0; iter < inst->n_seg; ++iter) { //iterations for solving the pricing problem
																	 //for (int iter = 0; iter < 1; ++iter) {


																	 //adding objective[seg]
																	 //IloObjective obj;
						model.update();
						model.setObjective(objective, GRB_MINIMIZE);
						model.getEnv().set(GRB_IntParam_Threads, 1);
						model.getEnv().set(GRB_DoubleParam_TimeLimit, time_limit);
						model.getEnv().set(GRB_IntParam_OutputFlag, 0);
						model.getEnv().set(GRB_IntParam_PreCrush, 1);
						//model.getEnv().set(GRB_DoubleParam_Cutoff, PRECISION);
						model.getEnv().set(GRB_IntParam_DualReductions, 0);
						model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
						//model.getEnv().set(GRB_IntParam_Presolve, -1);
						model.update();
						model.optimize();

						//cout << "Solution status: " << model.get(GRB_IntAttr_Status) << endl;
						//cout << ". Best objective: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
						//cout << "Pricing Time" << model.get(GRB_DoubleAttr_Runtime) << endl;

						//cout << "Time Taken for pricing iteration " << iter << ": " << model.get(GRB_DoubleAttr_Runtime) << endl;

						if (model.get(GRB_IntAttr_Status) == 2) {

							if (iter > 0 || model.get(GRB_DoubleAttr_ObjVal) < -PRECISION) {

								//cout << "Reduced Cost: " << model.get(GRB_DoubleAttr_ObjVal) << "\t";

								int ocup_cap = 0;

								Variable* p_var;
								p_var = new Variable();
								p_var->index = variables.size();
								p_var->assigned_group_size.resize(inst->n_subGroup);
								int column_repreat = 0; //number of times we want to repeat adding the current column because it may need to selected multiple times
								for (int g = 0; g < inst->n_subGroup; ++g) {
									p_var->assigned_group_size[g].resize(inst->max_party_size, 0);
									for (int nn = 0; nn < inst->max_party_size; ++nn) {
										if (inst->u[g][nn] > 0.5) {
											for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
												if (x[g][nn][i].get(GRB_DoubleAttr_X) > 0.9) {
													p_var->assigned_group_size[g][nn] = int(i);
													ocup_cap += i *nn;

													if (i > 0) {
														if (floor(inst->u[g][nn] / i) > column_repreat) {
															column_repreat = floor(inst->u[g][nn] / i);
														}
													}
												}
											}
										}
									}
								}
								p_var->capacity = inst->n_SeatRow;
								//computing value
								p_var->value = 0;
								for (int g = 0; g < inst->n_subGroup; ++g) {
									int n_party = 0;
									for (int nn = 0; nn < inst->max_party_size; ++nn) {
										if (x[g][nn][0].get(GRB_DoubleAttr_X) < 0.5) {
											n_party++;
										}
									}
									if (n_party > 0.9) {
										p_var->value++;
									}
								}
								int value1 = p_var->value;

								//adding the new variable
								//if (p_var->value > 0.5) {
								variables.push_back(p_var);
								for (int col = 0; col < column_repreat - 1; ++col) {
									variables.push_back(p_var);
								}
								//}

								//node_is_feasible = true;

								//cout << "occupied capacity: " << ocup_cap << endl;

								//updating the total objective
								total_obj += p_var->value;

								//updating xx[g][n]
								for (int g = 0; g < inst->n_subGroup; ++g) {
									for (int nn = 0; nn < inst->max_party_size; ++nn) {
										xx[g][nn] += p_var->assigned_group_size[g][nn];
										//adding new constraints to the model
										for (int i = inst->u[g][nn] - xx[g][nn] + 1; i < inst->u[g][nn] + 1; ++i) {
											model.addConstr(x[g][nn][i] == 0);
										}
									}
								}

								//if (model.get(GRB_DoubleAttr_ObjVal) > -PRECISION) {
								//	if (iter > 0 && last_RC > -PRECISION) {
								//		objective += Total_item * abs(last_RC) / inst->Tn_party;
								//	}
								//	objective -= Total_item * abs(model.get(GRB_DoubleAttr_ObjVal)) / inst->Tn_party;
								//	last_RC = model.get(GRB_DoubleAttr_ObjVal);
								//}

								objective = 0;
								for (int g = 0; g < inst->n_subGroup; ++g) {
									objective += y[g];
								}
								objective -= Total_item * 0.01;

							}

							else {
								iter = inst->n_seg;
								RMP_infsbl = true;
							}

						}


					}


					//bool feasible = true;
					//for (int g = 0; g < inst->n_subGroup; ++g) {
					//	for (int nn = 0; nn < inst->max_party_size; ++nn) {
					//		if (inst->u[g][nn] != xx[g][nn]) {
					//			feasible = false;
					//		}
					//	}
					//}

					//if (feasible) {
					//	//cout << "FEAS CHECK-FEASIBLE!!!!!!!!!!!!!" << endl;
					//	//cout << "Objective value: " << total_obj << endl;

					//	if (total_obj < global_upper_bound) {
					//		global_upper_bound = total_obj;
					//	}
					//}


					delete x;
					delete y;


				}
				catch (GRBException e) {
					cout << "Error number: " << e.getErrorCode() << endl;
					cout << e.getMessage() << endl;
					exit(1);
				}
				catch (...) {
					cout << "Other error ... " << endl;
					exit(1);
				}
			}
			///////////////////////////////////////////////////////////////////////////////////////////

		}
	}

	//cout << "node feasibility: " << node_is_feasible << endl;

}

void BPSolver::solve(double & time_limit, clock_t startTime, ofstream& file1) {

	//cout << "Start heuristics" << endl;

	//cout << "Creating initial set of variables .. " << endl;
	//create_first_set_of_variables();
	//create_first_set_of_variables_IP(file1);
	create_first_set_of_variables_bin_completion(file1);
	//create_first_set_of_variables2(file1);
	//cout << "Number of variables: " << (int)variables.size() << endl;

	cout << "End of heuristics" << endl;

	//creating the first search node:
	create_first_search_node();
	//print_search_node(search_nodes[0]);

	SearchNode* current_node;
	int n_nodes_explored = 0;
	int n_nodes_created = 1;
	int n_CG_iter = 0;
	int n_price_iter = 0;
	double total_price_time = 0;
	double total_root_node_time = 0;
	double root_node_LB = 0;
	double root_node_UB = 100000;
	int n_node_sofar = 1;

	while ((search_nodes.size() > 0) && (global_upper_bound > global_lower_bound + 0.99)) {


		////Setting the time_limit for the algorithm////
		double time_left;
		if (double(clock() - startTime) / (double)CLOCKS_PER_SEC >= time_limit) {
			cout << "time out and exit" << endl;
			break;
		}
		else {
			time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
		}

		//cout << "Number of nodes in the queue: " << search_nodes.size() << endl;
		//cout << "Number of nodes explored: " << n_nodes_explored << endl;
		//cout << "Global bounds: " << global_lower_bound << "\t" << global_upper_bound << endl;

		/*for (int i = 0; i < search_nodes.size(); ++i) {
		file1 << "(" << search_nodes[i]->index << "," << search_nodes[i]->parent_lower_bound << ")" << "\t";
		}
		file1 << endl;*/

		//check feasibility before picking a search node
		bool node_is_feasible = false;
		int no_checked_nodes = 0;
		while (!node_is_feasible && no_checked_nodes < search_nodes.size()) {
			//node_feasibility_check(search_nodes.front(), time_left, node_is_feasible, n_nodes_explored);
			Farkas_node_feasibility_check(search_nodes.front(), time_left, startTime, node_is_feasible, n_nodes_explored);
			no_checked_nodes++;
			time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
			if (time_left < 0.001) {
				cout << "time out and exit" << endl;
				break;
			}
		}

		time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
		if (time_left < 0.001) {
			cout << "time out and exit" << endl;
			break;
		}

		if (no_checked_nodes > search_nodes.size() && !node_is_feasible) {
			cout << "breaking (LB,UB): " << global_lower_bound << "," << global_upper_bound << endl;
			global_lower_bound = global_upper_bound;
			cout << "breaking (LB,UB): " << global_lower_bound << "," << global_upper_bound << endl;
			break;
		}

		// pick a search node to process
		//////Choosing the node with the smallest LB ////////
		current_node = search_nodes.front();
		search_nodes.erase(search_nodes.begin());
		//print_search_node(current_node);

		/*for (int g = 0; g < inst->n_subGroup; ++g) {
		//cout << "(" << current_node->LB_seg_no_for_group[g] << "," << current_node->UB_seg_no_for_group[g] << ")" << endl;
		}*/

		// solve the LP relaxation of that node (i.e., column generation)
		// FIND RELAXATION BOUND
		time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
		if (time_left < 0.001) {
			cout << "time out and exit" << endl;
			break;
		}
		double LP_start_time = time_left;

		solve_LP(current_node, time_left, time_limit, startTime, file1, n_CG_iter, n_price_iter, total_price_time);
		//solve_LP_dual_stabilization(current_node, time_left, time_limit, startTime, file1, n_CG_iter, n_price_iter, total_price_time);
		//solve_LP_with_pricing_and_DUALstabil(current_node, time_left, time_limit, startTime, file1, n_CG_iter, n_price_iter, total_price_time);
		//solve_LP_with_pricing_and_DUALstabil2(current_node, time_left, time_limit, startTime, file1, n_CG_iter, n_price_iter, total_price_time);
		if (n_node_sofar == 1) {
			total_root_node_time = LP_start_time - (time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC));
			root_node_LB = global_lower_bound;
			root_node_UB = global_upper_bound;
		}

		time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
		if (time_left < 0.001) {
			cout << "time out and exit" << endl;
			break;
		}

		// pick a variable to branch on
		// if branch_var_index still equals 1 after pick_branch_var, then I will know I don't have to branch
		int branch_var_index = -1;                            //index of the branch variable with fractional value (for the 1st branching strategy)
		vector<int> group_branch_index_and_value(2, -1);	  //includes: index of the group and value of the considered bound (for the 2nd branching strategy)

		time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
		if (time_left < 0.001) {
			cout << "time out and exit" << endl;
			break;
		}

		pick_branch_variable(current_node, time_left, file1, branch_var_index, group_branch_index_and_value);

		time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
		if (time_left < 0.001) {
			cout << "time out and exit" << endl;
			break;
		}

		//print_selected_node(current_node, file1);
		/*file1 << "(global_lower_bound,global_upper_bound): (" << global_lower_bound << "," << global_upper_bound << ")" << endl;
		file1 << "**************************************" << endl;*/

		// solve the IP for this node
		// FIND HEURISTIC SOLUTION
		solve_IP(current_node, time_left, file1);
		

		//cout << "(LB,UB): " << global_lower_bound << "," << global_upper_bound << endl;

		time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
		if (time_left < 0.001) {
			cout << "time out and exit" << endl;
			break;
		}

		//cout << "Branch variable: " << branch_var_index << endl;
		//cout << "group_branch_index_and_value[0]: " << group_branch_index_and_value[0] << endl;
		//cout << "group_branch_index_and_value[1]: " << group_branch_index_and_value[1] << endl;

		SearchNode* new_search_node;
		if (branch_var_index != -1 || group_branch_index_and_value[0] != -1) {


			// one branch node
			new_search_node = new SearchNode();
			new_search_node->index = n_nodes_created;
			n_nodes_created++;
			new_search_node->parent_lower_bound = current_node->LP_lower_bound;
			new_search_node->variables_set_to_one = current_node->variables_set_to_one;
			new_search_node->variables_set_to_zero = current_node->variables_set_to_zero;
			new_search_node->LB_seg_no_for_group = current_node->LB_seg_no_for_group;
			new_search_node->UB_seg_no_for_group = current_node->UB_seg_no_for_group;
			if (branch_var_index != -1) {
				new_search_node->variables_set_to_one.push_back(variables[branch_var_index]);
			}
			else {
				new_search_node->UB_seg_no_for_group[group_branch_index_and_value[0]] = group_branch_index_and_value[1];
			}
			search_nodes.push_back(new_search_node);

			// zero branch node
			new_search_node = new SearchNode();
			new_search_node->index = n_nodes_created;
			n_nodes_created++;
			new_search_node->parent_lower_bound = current_node->LP_lower_bound;
			new_search_node->variables_set_to_one = current_node->variables_set_to_one;
			new_search_node->variables_set_to_zero = current_node->variables_set_to_zero;
			new_search_node->LB_seg_no_for_group = current_node->LB_seg_no_for_group;
			new_search_node->UB_seg_no_for_group = current_node->UB_seg_no_for_group;
			if (branch_var_index != -1) {
				new_search_node->variables_set_to_zero.push_back(variables[branch_var_index]);
			}
			else {
				new_search_node->LB_seg_no_for_group[group_branch_index_and_value[0]] = group_branch_index_and_value[1] + 1;
			}
			search_nodes.push_back(new_search_node);

		}

		time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
		if (time_left < 0.001) {
			cout << "time out and exit" << endl;
			break;
		}

		//cout << "number of nodes in the queue: " << search_nodes.size() << endl;

		if (search_nodes.size() == 0) {
			delete current_node;
			cout << "breaking22222 (LB,UB): " << global_lower_bound << "," << global_upper_bound << endl;
			global_lower_bound = global_upper_bound;
			cout << "breaking22222 (LB,UB): " << global_lower_bound << "," << global_upper_bound << endl;
			break;
		}

		time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
		if (time_left < 0.001) {
			cout << "time out and exit" << endl;
			break;
		}

		// sorting nodes in the queue
		sort(search_nodes.begin(), search_nodes.end(), compare_nodes);

		time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
		if (time_left < 0.001) {
			cout << "time out and exit" << endl;
			break;
		}

		/*for (int s = 0; s < (int)search_nodes.size(); ++s) {
		print_search_node(search_nodes[s]);
		}*/

		if (global_lower_bound < search_nodes.front()->parent_lower_bound) {
			//cout << "************FOUND A NEW BETTER LOWER BOUND: " << search_nodes.front()->parent_lower_bound << endl;
			global_lower_bound = (search_nodes.front())->parent_lower_bound;
		}

		if (n_node_sofar == 1) {
			root_node_LB = global_lower_bound;
			root_node_UB = global_upper_bound;
		}
		n_node_sofar += 1;


		time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
		if (time_left < 0.001) {
			cout << "time out and exit" << endl;
			break;
		}

		/*file1 << n_nodes_explored << ",";
		file1 << variables.size() << ",";
		file1 << global_lower_bound << ",";
		file1 << global_upper_bound << endl;*/

		// delete this node
		delete current_node;
		n_nodes_explored++;
	}

	if (n_nodes_explored == 0) {
		root_node_LB = global_lower_bound;
		root_node_UB = global_upper_bound;
	}

	//	cout << "\nResults," << ceil(global_lower_bound-PRECISION) << "," << global_upper_bound <<","<<n_nodes_explored<<","<<variables.size()<< endl;
	//	cout << "BP-Time: " << double(clock() - startTime) / (double)CLOCKS_PER_SEC << endl;
	//print_best_solution();
	//cout << "BOUNDS: " << endl;
	//cout << "\tlower bound: " << global_lower_bound << endl;
	//cout << "\tupper bound: " << global_upper_bound << endl;

	//file1 << "(global_lower_bound,global_upper_bound): (" << global_lower_bound << "," << global_upper_bound << ")" << endl;
	//file1 << "**************************************" << endl;

	// Printing the results in the ofstream file
	//cout << "Global bounds: " << ceil(global_lower_bound) << "\t" << global_upper_bound << endl;
	file1 << double(clock() - startTime) / (double)CLOCKS_PER_SEC << ",";
	file1 << ceil(global_lower_bound - PRECISION) << ",";			    //lower bound
	file1 << global_upper_bound << ",";								    //upper bound
	file1 << total_root_node_time << ",";                               //total time to solve the root node
	file1 << root_node_LB << ",";                                       //root node relaxation solution
	file1 << root_node_UB << ",";                                       //best integer solution after solving the root node
	file1 << (n_nodes_explored + 1) << ",";								    //number of nodes explored
	file1 << variables.size() << ",";									//total number of added columns
	file1 << n_CG_iter << ",";										    //total number of CG iterations 
	file1 << n_price_iter << ",";									    //total number of pricing iterations
	file1 << total_price_time << ",";									//total pricing time
	file1 << total_price_time / n_price_iter << endl;					//average pricing time


	double gap = -1.0;
	double ub = global_upper_bound;
	double lb = ceil(global_lower_bound - PRECISION);
	double runtime = double(clock() - startTime) / (double)CLOCKS_PER_SEC;

	if (ub  < 1e6)
		gap = abs((ub - lb) / ub);

	cout << "Results," <<
		runtime << "," <<
		lb << "," <<
		ub << "," <<
		gap << "\n";


	// create a first search tree node
	// initialize queue of nodes

	// create while loop to search over the nodes
	// create a solver for the LP relaxations
	// create a solver for the IPs
	// create a branch function
	// implement pruning rules

}

void BPSolver::solve_RMP_LP(SearchNode* node, double & time_limit, ofstream& file1, vector<vector<double>> & party_duals, vector<double> & capacity_duals, vector<double> & seg_UB_duals, vector<double> & new_LB_duals) {

	//cout << "=============================================" << endl << "Solving Restriced Master Problem:" << endl << "=============================================";
	//GRBEnv env;
	char varname[256];

	try {
		GRBModel model(env);

		// Define variables

		GRBVar* x_rmp;
		x_rmp = new GRBVar[(int)variables.size()];
		for (int v = 0; v < (int)variables.size(); ++v) {
			//sprintf(varname, "x[%d]", v);
			x_rmp[v] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
		}

		model.update();         // required before you can use the variables

								// adding constraints

								// \sum_c m[c][g][nn] x_c == u[g][nn] for all group g and size nn 
		vector<vector<GRBConstr >>  group_size_constraints(inst->n_subGroup);
		for (int g = 0; g < inst->n_subGroup; ++g) {
			group_size_constraints[g].resize(inst->max_party_size);
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				GRBLinExpr group_size_con = 0.0;
				for (int c = 0; c < (int)variables.size(); ++c) {
					group_size_con += variables[c]->assigned_group_size[g][nn] * x_rmp[c];
				}
				group_size_constraints[g][nn] = model.addConstr(group_size_con >= inst->u[g][nn]);
			}
		}

		//// number of selected columns from each capacity must be <= number of segemnets available of that capacity
		vector<GRBConstr>  capacity_constraints(inst->n_SeatRow);
		for (int cap = 0; cap < inst->n_SeatRow; ++cap) {
			GRBLinExpr capacity_con = 0.0;
			for (int c = 0; c < (int)variables.size(); ++c) {
				if (variables[c]->capacity == cap + 1) {
					capacity_con += x_rmp[c];
				}
			}
			capacity_constraints[cap] = model.addConstr(capacity_con <= inst->no_seg_of_each_size[cap]);
		}

		// LB and UB constraint for number of segments assigned to each group
		// Also: similar LB suggested by IJOC-reviewer1 in revision1 for IP2
		//vector<GRBConstr >  seg_LB_constraints(inst->n_subGroup);
		vector<GRBConstr >  seg_UB_constraints(inst->n_subGroup);
		vector<GRBConstr >  new_LB_constraints(inst->n_subGroup);
		for (int g = 0; g < inst->n_subGroup; ++g) {
			GRBLinExpr n_seg_assigned_to_g = 0.0;
			for (int c = 0; c < (int)variables.size(); ++c) {
				int n_party = 0.0;
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					if (variables[c]->assigned_group_size[g][nn] > 0.9 && inst->u[g][nn] > 0.9) {
						n_party++;
					}
				}
				if (n_party > 0.9) {
					n_seg_assigned_to_g += x_rmp[c];
				}
			}
			//seg_LB_constraints[g] = model.addConstr(n_seg_assigned_to_g >= node->LB_seg_no_for_group[g]);
			seg_UB_constraints[g] = model.addConstr(n_seg_assigned_to_g <= node->UB_seg_no_for_group[g]);

			double LB_n_bin = inst->nn[g] / inst->n_SeatRow;
			if (LB_n_bin - floor(LB_n_bin) >= 0.001) {
				LB_n_bin = floor(LB_n_bin) + 1;
			}
			else {
				LB_n_bin = floor(LB_n_bin);
			}
			double new_LB = LB_n_bin;
			if (inst->LB_MT1[g] > new_LB) {
				new_LB = inst->LB_MT1[g];
			}
			if (inst->LB_MT2[g] > new_LB) {
				new_LB = inst->LB_MT2[g];
			}
			if (node->LB_seg_no_for_group[g] > new_LB) {
				new_LB = node->LB_seg_no_for_group[g];
			}
			new_LB_constraints[g] = model.addConstr(n_seg_assigned_to_g >= new_LB);
		}


		// for these two constraints, you might want to using vector iterators, and
		// find the position in the vector.  We will, however, make sure the index matches
		// the relateive position
		for (int i = 0; i < (int)node->variables_set_to_one.size(); ++i) {
			model.addConstr(x_rmp[node->variables_set_to_one[i]->index] == 1);
		}

		for (int i = 0; i < (int)node->variables_set_to_zero.size(); ++i) {
			model.addConstr(x_rmp[node->variables_set_to_zero[i]->index] == 0);
		}

		// objective function
		GRBLinExpr objective_RMP = 0.0;
		for (int c = 0; c < (int)variables.size(); ++c) {
			objective_RMP += variables[c]->value * x_rmp[c];
		}

		model.addConstr(objective_RMP >= global_lower_bound);
		model.addConstr(objective_RMP <= global_upper_bound);

		model.setObjective(objective_RMP, GRB_MINIMIZE);

		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_DoubleParam_TimeLimit, time_limit);
		model.getEnv().set(GRB_IntParam_OutputFlag, 0);
		// model.getEnv().set(GRB_IntParam_Cuts,0);
		//            model.getEnv().set(GRB_IntParam_PreCrush,1);
		//model.set(GRB_IntParam_PreCrush,1);
		//            model.getEnv().set(GRB_IntParam_DualReductions, 0);
		//            model.getEnv().set(GRB_IntParam_LazyConstraints, 1);

		model.update();
		model.optimize(); //solves the problem
						  //cout << "Solution Status: " << model.get(GRB_IntAttr_Status) << endl;
						  //cout << "Time Taken: " << model.get(GRB_DoubleAttr_Runtime) << endl;
		if (model.get(GRB_IntAttr_Status) == 2) {
			// This gets you the relaxation bound
			//cout << "Relaxation Bound: " << model.get(GRB_DoubleAttr_ObjBound) << endl;
			// This gets you the best known primal solution value

			//cout << "Current RMP LP value: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

			//cout << "number of columns *********:   " << variables.size() << endl;
			//cout << endl;

			/*cout << "value[c]: " << endl;
			for (int c = 0; c < variables.size(); ++c) {
			cout << variables[c]->value << "\t";
			}
			cout << endl;*/

			/*cout << "selected x[c]: " << endl;
			for (int c = 0; c < variables.size(); ++c) {
			if (x[c].get(GRB_DoubleAttr_X) > 0.01) {
			cout << "column " << c << ": " << x[c].get(GRB_DoubleAttr_X) << "\t";
			}
			if (x[c].get(GRB_DoubleAttr_X) > 0.99) {
			cout << endl;
			cout << "capacity: " << variables[c]->capacity << endl;
			cout << "parties in the variable: " << endl;
			for (int p = 0; p < inst->Tn_party; ++p) {
			if (variables[c]->assigned_party[p]) {
			cout << p << "\t";
			}
			}
			cout << endl;
			}
			}
			cout << endl;*/


			///////// CHECK /////////////
			//cout << "**** DUALS ***** " << endl;

			////cout << "party_duals:" << endl;
			for (int g = 0; g < inst->n_subGroup; ++g) {
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					party_duals[g][nn] = group_size_constraints[g][nn].get(GRB_DoubleAttr_Pi);
				}
			}

			////cout << "capacity_duals:" << endl;
			for (int cap = 0; cap < inst->n_SeatRow; ++cap) {
				capacity_duals[cap] = capacity_constraints[cap].get(GRB_DoubleAttr_Pi);
				////cout << capacity_duals[cap] << "\t";
			}
			////cout << endl;
			////cout << endl;

			////cout << "seg_LB_and_UB_duals" <<endl;
			for (int g = 0; g < inst->n_subGroup; ++g) {
				//seg_LB_duals[g] = seg_LB_constraints[g].get(GRB_DoubleAttr_Pi);
				seg_UB_duals[g] = seg_UB_constraints[g].get(GRB_DoubleAttr_Pi);
				new_LB_duals[g] = new_LB_constraints[g].get(GRB_DoubleAttr_Pi);
			}

			///////////////////Checking reduced costs//////////////////////////////
			////cout << "Given RC" << "\t" << "calculated RC" << endl;
			//for (int c = 0; c < variables.size(); ++c) {
			//	double calculated_rc = 0.0;
			//	calculated_rc += variables[c]->value;

			//	for (int g = 0; g < inst->n_subGroup; ++g) {
			//		for (int nn = 0; nn < inst->max_party_size; ++nn) {
			//			//if (variables[c]->assigned_group_size[g][nn] > 0.9) {
			//				calculated_rc -= variables[c]->assigned_group_size[g][nn] * party_duals[g][nn];
			//			//}
			//		}
			//	}

			//	calculated_rc -= capacity_duals[variables[c]->capacity - 1];

			//	for (int g = 0; g < inst->n_subGroup; ++g) {
			//		int n_party = 0;
			//		for (int nn = 0; nn < inst->max_party_size; ++nn) {
			//			if (variables[c]->assigned_group_size[g][nn] > 0.9) {
			//				n_party++;
			//			}
			//		}
			//		if (n_party > 0.9) {
			//			calculated_rc -= (seg_LB_duals[g] + seg_UB_duals[g]);
			//		}
			//	}

			//	//cout << x[c].get(GRB_DoubleAttr_RC) << "\t" << calculated_rc << endl;
			//	if ((x[c].get(GRB_DoubleAttr_RC) - calculated_rc) < -PRECISION || (x[c].get(GRB_DoubleAttr_RC) - calculated_rc) > PRECISION) {
			//		//cout << "ERROR" << "\t";
			//		//cout << x[c].get(GRB_DoubleAttr_RC) << "\t" << calculated_rc << endl; 
			//		system("pause");
			//	}	
			//}

		}

		delete x_rmp;
	}

	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		exit(1);
	}
	catch (...) {
		cout << "Other error ... " << endl;
		exit(1);
	}


}

void BPSolver::solve_pricing_problem(SearchNode* node, double & time_limit, ofstream& file1, vector<vector<double>> & party_duals, vector<double> & capacity_duals, vector<double> & seg_UB_duals, vector<double> & new_LB_duals, int capacity, bool & need_another_rmp_lp_solution, int& n_price_iter, double& total_price_time) {
	//cout << "=============================================" << endl << "Solving Pricing Problem:" << endl << "=============================================";

	vector<vector<int>> xx; //xx[g][n] number of items of size n from group g, so far selected 
	xx.resize(inst->n_subGroup);
	for (int g = 0; g < inst->n_subGroup; ++g) {
		xx[g].resize(inst->max_party_size, 0);
	}

	int total_obj = 0;

	char varname[256];

	//GRBEnv env;
	try {
		GRBModel model(env);

		//Defining and adding variables
		GRBVar*** x;																		//x[g][n][i] : Number of party of size n from group g assigned to this column is equal to i
		x = new GRBVar**[inst->n_subGroup];
		for (int g = 0; g < inst->n_subGroup; g++) {
			x[g] = new GRBVar*[inst->max_party_size];
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				x[g][nn] = new GRBVar[inst->u[g][nn] + 1];
				for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
					x[g][nn][i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				}
			}
		}

		GRBVar* y;																		//y[g] is equal to one if group "g" has some party in this column, and zero otherwise
		y = new GRBVar[inst->n_subGroup];
		for (int g = 0; g < inst->n_subGroup; ++g) {
			y[g] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}

		model.update();
		//adding constraints

		//sum_i x[g][nn][i] = 1
		for (int g = 0; g < inst->n_subGroup; g++) {
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				GRBLinExpr equ = 0.0;
				for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
					equ += x[g][nn][i];
				}
				model.addConstr(equ == 1);
			}
		}

		////i*x[g][nn][i] <= number of parties of size nn+1 from group g
		//for (int g = 0; g < inst->n_subGroup; g++) {
		//	for (int nn = 0; nn < inst->max_party_size; ++nn) {
		//		for (int i = 1; i < inst->u[g][nn] + 1; ++i) {
		//			model.addConstr(i * x[g][nn][i] <= inst->u[g][nn]);
		//		}
		//	}
		//}

		//////i*x[g][nn][i] <= number of parties of size nn+1 from group g
		////for (int g = 0; g < inst->n_subGroup; g++) {
		////	GRBLinExpr equ11 = 0.0;
		////	for (int nn = 0; nn < inst->max_party_size; ++nn) {
		////		equ11 += (x[g][nn][0] - 1);
		////	}
		////	model.addConstr(inst->max_party_size * y[g] + equ11 >= 0);
		////}

		//x[g][nn]<= y[g] for all group g and size nn
		for (int g = 0; g < inst->n_subGroup; ++g) {
			//GRBLinExpr Yequ = 0.0;
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				model.addConstr(1 - x[g][nn][0] <= y[g]);
				//Yequ += (1 - x[g][nn][0]);
			}
			//model.addConstr(y[g] <= Yequ);
		}


		// \sum_{size} \sum_{i} x[g][nn]<= capacity * y[g] for all group g and size nn
		for (int g = 0; g < inst->n_subGroup; ++g) {
			GRBLinExpr equ_new = 0.0;
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
					equ_new += (nn + 1) * i * x[g][nn][i];
				}
			}
			model.addConstr(equ_new <= capacity * y[g]);
		}

		//Knapsack constraint
		GRBLinExpr cap_equ = 0.0;
		for (int g = 0; g < inst->n_subGroup; ++g) {
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
					cap_equ += (nn + 1) * i * x[g][nn][i];
				}
			}
		}
		model.addConstr(cap_equ <= capacity);

		///////////////////////////
		//Adding objective function 
		GRBLinExpr objective = 0.0;
		for (int g = 0; g < inst->n_subGroup; ++g) {
			objective += y[g];
		}
		for (int g = 0; g < inst->n_subGroup; ++g) {
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				for (int i = 1; i < inst->u[g][nn] + 1; ++i) {
					objective -= party_duals[g][nn] * i * x[g][nn][i];
				}
			}
		}
		objective -= capacity_duals[capacity - 1];
		for (int g = 0; g < inst->n_subGroup; ++g) {
			objective -= (new_LB_duals[g] + seg_UB_duals[g]) * y[g];
		}



		// Using Callback to add the avoid constraints (i.e. to avoid from regenerating the existing columns)
		callback_pricing cb = callback_pricing(y, x, variables, inst, capacity);
		model.setCallback(&cb);

		//calculating the total number of items selected
		GRBLinExpr Total_item = 0.0;
		for (int g = 0; g < inst->n_subGroup; ++g) {
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				for (int i = 1; i < inst->u[g][nn] + 1; ++i) {
					Total_item += (nn + 1) * i * x[g][nn][i];
				}
			}
		}

		double last_RC = 0.0;


		for (int iter = 0; iter < inst->n_seg; ++iter) { //iterations for solving the pricing problem

			//if (iter >= 0) {
			//	// Using Callback to add the avoid constraints (i.e. to avoid from regenerating the existing columns)
			//	callback_pricing cb = callback_pricing(y, x, variables, inst, capacity);
			//	model.setCallback(&cb);
			//}

														 //adding objective[seg]
														 //IloObjective obj;
			model.update();
			model.setObjective(objective, GRB_MINIMIZE);
			model.getEnv().set(GRB_IntParam_Threads, 1);
			model.getEnv().set(GRB_DoubleParam_TimeLimit, time_limit);
			model.getEnv().set(GRB_IntParam_OutputFlag, 0);
			//model.getEnv().set(GRB_DoubleParam_Cutoff, PRECISION);
			//model.getEnv().set(GRB_IntParam_DualReductions, 0);
			if (iter >= 0) {
				model.getEnv().set(GRB_IntParam_PreCrush, 1);
				model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
			}
			model.getEnv().set(GRB_IntParam_Presolve, -1);
			model.update();
			model.optimize();

			n_price_iter += 1;
			total_price_time += model.get(GRB_DoubleAttr_Runtime);

			//cout << "Solution status: " << model.get(GRB_IntAttr_Status) << endl;
			//cout << ". Best objective: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
			//cout << "Pricing Time" << model.get(GRB_DoubleAttr_Runtime) << endl;

			//cout << "Time Taken for pricing iteration " << iter << ": " << model.get(GRB_DoubleAttr_Runtime) << endl;

			if (model.get(GRB_IntAttr_Status) == 2) {

				if (iter > 0 || model.get(GRB_DoubleAttr_ObjVal) < -PRECISION) {

					//cout << "Reduced Cost: " << model.get(GRB_DoubleAttr_ObjVal) << "\t";

					int ocup_cap = 0;

					Variable* p_var;
					p_var = new Variable();
					p_var->index = variables.size();
					p_var->assigned_group_size.resize(inst->n_subGroup);
					int column_repreat = 0; //number of times we want to repeat adding the current column because it may need to selected multiple times
					for (int g = 0; g < inst->n_subGroup; ++g) {
						p_var->assigned_group_size[g].resize(inst->max_party_size, 0);
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							if (inst->u[g][nn] > 0.5) {
								for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
									if (x[g][nn][i].get(GRB_DoubleAttr_X) > 0.9) {
										p_var->assigned_group_size[g][nn] = int(i);
										ocup_cap += i *nn;

										if (i > 0) {
											if (floor(inst->u[g][nn] / i) > column_repreat) {
												column_repreat = floor(inst->u[g][nn] / i);
											}
										}
									}
								}
							}
						}
					}
					p_var->capacity = capacity;
					//computing value
					p_var->value = 0;
					for (int g = 0; g < inst->n_subGroup; ++g) {
						int n_party = 0;
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							if (x[g][nn][0].get(GRB_DoubleAttr_X) < 0.5) {
								n_party++;
							}
						}
						if (n_party > 0.9) {
							p_var->value++;
						}
					}
					int value1 = p_var->value;

					//adding the new variable	
					if (p_var->value > 0.5) {
						variables.push_back(p_var);
						for (int col = 0; col < column_repreat - 1; ++col) {
							variables.push_back(p_var);
						}
					}

					need_another_rmp_lp_solution = true;

					//cout << "occupied capacity: " << ocup_cap << endl;

					//updating the total objective
					total_obj += p_var->value;

					//updating xx[g][n]
					for (int g = 0; g < inst->n_subGroup; ++g) {
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							xx[g][nn] += p_var->assigned_group_size[g][nn];
							//adding new constraints to the model
							for (int i = inst->u[g][nn] - xx[g][nn] + 1; i < inst->u[g][nn] + 1; ++i) {
								model.addConstr(x[g][nn][i] == 0);
							}
						}
					}

					//if (model.get(GRB_DoubleAttr_ObjVal) > -PRECISION) {
					//	if (iter > 0 && last_RC > -PRECISION) {
					//		objective += Total_item * abs(last_RC) / inst->Tn_party;
					//	}
					//	objective -= Total_item * abs(model.get(GRB_DoubleAttr_ObjVal)) / inst->Tn_party;
					//	last_RC = model.get(GRB_DoubleAttr_ObjVal);
					//}

					objective = 0;
					for (int g = 0; g < inst->n_subGroup; ++g) {
						objective += y[g];
					}
					objective -= Total_item / inst->no_seat_in_seg[0];

				}

				else {
					iter = inst->n_seg;
				}

			}


		}


		bool feasible = true;
		for (int g = 0; g < inst->n_subGroup; ++g) {
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				if (inst->u[g][nn] != xx[g][nn]) {
					feasible = false;
				}
			}
		}

		if (feasible) {
			//cout << "FEASIBLE!!!!!!!!!!!!!" << endl;
			//cout << "Objective value: " << total_obj << endl;

			if (total_obj < global_upper_bound) {
				global_upper_bound = total_obj;
			}
		}


		delete x;
		delete y;

	} // end of try
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		exit(1);
	}
	catch (...) {
		cout << "Other error ... " << endl;
		exit(1);
	}
}

void BPSolver::solve_LP(SearchNode* node, double & time_left, double time_limit, clock_t startTime, ofstream& file1, int& n_CG_iter, int& n_price_iter, double& total_price_time) {

	bool need_another_rmp_lp_solution = true;
	int n_lps_in_col_gen = 0;
	clock_t LP_startTime = clock();

	// the following vectors are defined to record all the duals for "dual price smoothing" method
	double alpha_weight = 0.5;     //The alpha in the Neame (2000) method
	vector<vector<vector<double >>> all_party_duals;
	vector<vector<double>> all_capacity_duals;
	//vector<vector<double>> all_seg_LB_duals;
	vector<vector<double>> all_seg_UB_duals;
	vector<vector<double>> all_new_LB_duals;
	bool use_dual_stabil = true;


	while (need_another_rmp_lp_solution) {
		n_CG_iter += 1;

		time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
		if (time_left < 0.001) {
			cout << "time out and exit" << endl;
			break;
		}

		n_lps_in_col_gen++;
		//cout << "LOOP IN COL GEN: " << n_lps_in_col_gen << endl;
		if (n_lps_in_col_gen > 1000) {
			cout << "solve_LP: Exceeded n_lps_in_col_gen \n";
			exit(1);
		}
		need_another_rmp_lp_solution = false;

		vector<vector<double >> party_duals(inst->n_subGroup);	// one for every party
		for (int g = 0; g < inst->n_subGroup; ++g) {
			party_duals[g].resize(inst->max_party_size);
		}
		vector<double > capacity_duals(inst->n_SeatRow);	// one for every size from 1 to n_SeatRow
		//vector<double> seg_LB_duals(inst->n_subGroup);
		vector<double> seg_UB_duals(inst->n_subGroup);
		vector<double> new_LB_duals(inst->n_subGroup);

		// solve restricted MP
		time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
		if (time_left < 0.001) {
			cout << "time out and exit" << endl;
			break;
		}
		solve_RMP_LP(node, time_left, file1, party_duals, capacity_duals, seg_UB_duals, new_LB_duals);


		//calculating the smoothed duals by Neame (2000) method
		//double alpha_weight = 0.3;     //The alpha in the Neame (2000) method
		all_party_duals.push_back(party_duals);
		all_capacity_duals.push_back(capacity_duals);
		//all_seg_LB_duals.push_back(seg_LB_duals);
		all_seg_UB_duals.push_back(seg_UB_duals);
		all_new_LB_duals.push_back(new_LB_duals);

		for (int g = 0; g < inst->n_subGroup; ++g) {
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				party_duals[g][nn] = 0;
				for (int iter = 0; iter < all_party_duals.size() - 1; ++iter) {
					party_duals[g][nn] += alpha_weight * (1 - alpha_weight) * pow(alpha_weight, (all_party_duals.size() - iter)) * all_party_duals[iter][g][nn];
				}
				party_duals[g][nn] += (1 - alpha_weight) * all_party_duals[all_party_duals.size() - 1][g][nn];
			}
		}

		for (int i = 0; i < inst->n_SeatRow; ++i) {
			capacity_duals[i] = 0;
			for (int iter = 0; iter < all_capacity_duals.size() - 1; ++iter) {
				capacity_duals[i] += alpha_weight * (1 - alpha_weight) * pow(alpha_weight, (all_capacity_duals.size() - iter)) * all_capacity_duals[iter][i];
			}
			capacity_duals[i] += (1 - alpha_weight) * all_capacity_duals[all_capacity_duals.size() - 1][i];
		}

		//for (int g = 0; g < inst->n_subGroup; ++g) {
		//	seg_LB_duals[g] = 0;
		//	for (int iter = 0; iter < all_seg_LB_duals.size() - 1; ++iter) {
		//		seg_LB_duals[g] += alpha_weight * (1 - alpha_weight) * pow(alpha_weight, all_seg_LB_duals.size() - iter) * all_seg_LB_duals[iter][g];
		//	}
		//	seg_LB_duals[g] += (1 - alpha_weight) * all_seg_LB_duals[all_seg_LB_duals.size() - 1][g];
		//}

		for (int g = 0; g < inst->n_subGroup; ++g) {
			seg_UB_duals[g] = 0;
			for (int iter = 0; iter < all_seg_UB_duals.size() - 1; ++iter) {
				seg_UB_duals[g] += alpha_weight * (1 - alpha_weight) * pow(alpha_weight, (all_seg_UB_duals.size() - iter)) * all_seg_UB_duals[iter][g];
			}
			seg_UB_duals[g] += (1 - alpha_weight) * all_seg_UB_duals[all_seg_UB_duals.size() - 1][g];
		}

		for (int g = 0; g < inst->n_subGroup; ++g) {
			new_LB_duals[g] = 0;
			for (int iter = 0; iter < all_new_LB_duals.size() - 1; ++iter) {
				new_LB_duals[g] += alpha_weight * (1 - alpha_weight) * pow(alpha_weight, all_new_LB_duals.size() - iter) * all_new_LB_duals[iter][g];
			}
			new_LB_duals[g] += (1 - alpha_weight) * all_new_LB_duals[all_new_LB_duals.size() - 1][g];
		}


		time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
		if (time_left < 0.001) {
			cout << "time out and exit" << endl;
			break;
		}

		// solve pricing problem for each available segemnt size
		if (use_dual_stabil) {
			for (int cap = 0; cap < inst->n_SeatRow; ++cap) {
				if (inst->no_seg_of_each_size[cap] > 0) {
					solve_pricing_problem(node, time_left, file1, party_duals, capacity_duals, seg_UB_duals, new_LB_duals, cap + 1, need_another_rmp_lp_solution, n_price_iter, total_price_time);
					//solve_MDD_pricing_problem(node, time_left, file1, party_duals, capacity_duals, seg_LB_duals, seg_UB_duals, cap + 1, need_another_rmp_lp_solution, n_price_iter, total_price_time);
				}
				time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
				if (time_left < 0.001) {
					cout << "time out and exit" << endl;
					break;
				}
			}

			if (!need_another_rmp_lp_solution) {
				need_another_rmp_lp_solution = true;
				use_dual_stabil = false;
			}
		}

		else {
			for (int cap = 0; cap < inst->n_SeatRow; ++cap) {
				if (inst->no_seg_of_each_size[cap] > 0) {
					solve_pricing_problem(node, time_left, file1, all_party_duals[all_party_duals.size() - 1], all_capacity_duals[all_capacity_duals.size() - 1], all_seg_UB_duals[all_seg_UB_duals.size() - 1], all_new_LB_duals[all_new_LB_duals.size() - 1], cap + 1, need_another_rmp_lp_solution, n_price_iter, total_price_time);
					//solve_MDD_pricing_problem(node, time_left, file1, party_duals, capacity_duals, seg_LB_duals, seg_UB_duals, cap + 1, need_another_rmp_lp_solution, n_price_iter, total_price_time);
				}
				time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
				if (time_left < 0.001) {
					cout << "time out and exit" << endl;
					break;
				}
			}
		}

		
		//cout << "Number of columns: " << variables.size() << endl;
	}

	double LP_Time = double(clock() - LP_startTime) / (double)CLOCKS_PER_SEC;
	//file1 << n_lps_in_col_gen << ",";
	//file1 << LP_Time << ",";
}

void BPSolver::solve_LP_with_pricing_and_DUALstabil(SearchNode* node, double & time_left, double time_limit, clock_t startTime, ofstream& file1, int& n_CG_iter, int& n_price_iter, double& total_price_time) {


	char varname[256];


	// creating the base of the pricing model
	GRBEnv env_price;
	try {
		GRBModel model_p(env_price);

		////////////////////////////////////////////////////////////
		//Defining and adding variables
		GRBVar*** x;																		//x[g][n][i] : Number of party of size n from group g assigned to this column is equal to i
		x = new GRBVar**[inst->n_subGroup];
		for (int g = 0; g < inst->n_subGroup; g++) {
			x[g] = new GRBVar*[inst->max_party_size];
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				x[g][nn] = new GRBVar[inst->u[g][nn] + 1];
				for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
					x[g][nn][i] = model_p.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				}
			}
		}

		GRBVar* y;																		//y[g] is equal to one if group "g" has some party in this column, and zero otherwise
		y = new GRBVar[inst->n_subGroup];
		for (int g = 0; g < inst->n_subGroup; ++g) {
			y[g] = model_p.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}

		model_p.update();

		//adding constraints

		//sum_i x[g][nn][i] = 1
		for (int g = 0; g < inst->n_subGroup; g++) {
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				GRBLinExpr equ = 0.0;
				for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
					equ += x[g][nn][i];
				}
				model_p.addConstr(equ == 1);
			}
		}

		////i*x[g][nn][i] <= number of parties of size nn+1 from group g
		//for (int g = 0; g < inst->n_subGroup; g++) {
		//	GRBLinExpr equ11 = 0.0;
		//	for (int nn = 0; nn < inst->max_party_size; ++nn) {
		//		equ11 += (x[g][nn][0] - 1);
		//	}
		//	model_p.addConstr(inst->max_party_size * y[g] + equ11 >= 0);
		//}

		//x[g][nn]<= y[g] for all group g and size nn
		for (int g = 0; g < inst->n_subGroup; ++g) {
			//GRBLinExpr Yequ = 0.0;
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				model_p.addConstr(1 - x[g][nn][0] <= y[g]);
				//Yequ += (1 - x[g][nn][0]);
			}
			//model_p.addConstr(y[g] <= Yequ);
		}


		// \sum_{size} \sum_{i} x[g][nn]<= capacity * y[g] for all group g and size nn
		for (int g = 0; g < inst->n_subGroup; ++g) {
			GRBLinExpr equ_new = 0.0;
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
					equ_new += (nn + 1) * i * x[g][nn][i];
				}
			}
			model_p.addConstr(equ_new <= inst->no_seat_in_seg[0] * y[g]);
		}

		//Knapsack constraint
		GRBLinExpr cap_equ = 0.0;
		for (int g = 0; g < inst->n_subGroup; ++g) {
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
					cap_equ += (nn + 1) * i * x[g][nn][i];
				}
			}
		}
		model_p.addConstr(cap_equ <= inst->no_seat_in_seg[0]);

		///////////////////////////////////////////////////////////////////

		bool need_another_rmp_lp_solution = true;
		int n_lps_in_col_gen = 0;
		clock_t LP_startTime = clock();

		// the following vectors are defined to record all the duals for "dual price smoothing" method
		double alpha_weight = 0.5;     //The alpha in the Neame (2000) method
		vector<vector<vector<double >>> all_party_duals;
		vector<vector<double>> all_capacity_duals;
		//vector<vector<double>> all_seg_LB_duals;
		vector<vector<double>> all_seg_UB_duals;
		vector<vector<double>> all_new_LB_duals;

		while (need_another_rmp_lp_solution) {

			if (global_upper_bound <= global_lower_bound + 0.99) {
				time_left = 0;
				break;
			}

			n_CG_iter += 1;

			time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
			if (time_left < 0.001) {
				//cout << "time out and exit" << endl;
				break;
			}

			n_lps_in_col_gen++;
			//cout << "LOOP IN COL GEN: " << n_lps_in_col_gen << endl;
			if (n_lps_in_col_gen > 1000) {
				cout << "solve_LP: Exceeded n_lps_in_col_gen \n";
				exit(1);
			}
			need_another_rmp_lp_solution = false;

			vector<vector<double >> party_duals(inst->n_subGroup);	// one for every party
			for (int g = 0; g < inst->n_subGroup; ++g) {
				party_duals[g].resize(inst->max_party_size);
			}
			vector<double > capacity_duals(inst->n_SeatRow);	// one for every size from 1 to n_SeatRow
																//vector<double> seg_LB_duals(inst->n_subGroup);
			vector<double> seg_UB_duals(inst->n_subGroup);
			vector<double> new_LB_duals(inst->n_subGroup);

			// solve restricted MP
			time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
			if (time_left < 0.001) {
				cout << "time out and exit" << endl;
				break;
			}
			solve_RMP_LP(node, time_left, file1, party_duals, capacity_duals, seg_UB_duals, new_LB_duals);

			//calculating the smoothed duals by Neame (2000) method
			//double alpha_weight = 0.3;     //The alpha in the Neame (2000) method
			all_party_duals.push_back(party_duals);
			all_capacity_duals.push_back(capacity_duals);
			//all_seg_LB_duals.push_back(seg_LB_duals);
			all_seg_UB_duals.push_back(seg_UB_duals);
			all_new_LB_duals.push_back(new_LB_duals);

			for (int g = 0; g < inst->n_subGroup; ++g) {
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					party_duals[g][nn] = 0;
					for (int iter = 0; iter < all_party_duals.size() - 1; ++iter) {
						party_duals[g][nn] += alpha_weight * (1 - alpha_weight) * pow(alpha_weight, all_party_duals.size() - iter) * all_party_duals[iter][g][nn];
					}
					party_duals[g][nn] += (1 - alpha_weight) * all_party_duals[all_party_duals.size() - 1][g][nn];
				}
			}

			for (int i = 0; i < inst->n_SeatRow; ++i) {
				capacity_duals[i] = 0;
				for (int iter = 0; iter < all_capacity_duals.size() - 1; ++iter) {
					capacity_duals[i] += alpha_weight * (1 - alpha_weight) * pow(alpha_weight, all_capacity_duals.size() - iter) * all_capacity_duals[iter][i];
				}
				capacity_duals[i] += (1 - alpha_weight) * all_capacity_duals[all_capacity_duals.size() - 1][i];
			}

			//for (int g = 0; g < inst->n_subGroup; ++g) {
			//	seg_LB_duals[g] = 0;
			//	for (int iter = 0; iter < all_seg_LB_duals.size() - 1; ++iter) {
			//		seg_LB_duals[g] += alpha_weight * (1 - alpha_weight) * pow(alpha_weight, all_seg_LB_duals.size() - iter) * all_seg_LB_duals[iter][g];
			//	}
			//	seg_LB_duals[g] += (1 - alpha_weight) * all_seg_LB_duals[all_seg_LB_duals.size() - 1][g];
			//}

			for (int g = 0; g < inst->n_subGroup; ++g) {
				seg_UB_duals[g] = 0;
				for (int iter = 0; iter < all_seg_UB_duals.size() - 1; ++iter) {
					seg_UB_duals[g] += alpha_weight * (1 - alpha_weight) * pow(alpha_weight, all_seg_UB_duals.size() - iter) * all_seg_UB_duals[iter][g];
				}
				seg_UB_duals[g] += (1 - alpha_weight) * all_seg_UB_duals[all_seg_UB_duals.size() - 1][g];
			}

			for (int g = 0; g < inst->n_subGroup; ++g) {
				new_LB_duals[g] = 0;
				for (int iter = 0; iter < all_new_LB_duals.size() - 1; ++iter) {
					new_LB_duals[g] += alpha_weight * (1 - alpha_weight) * pow(alpha_weight, all_new_LB_duals.size() - iter) * all_new_LB_duals[iter][g];
				}
				new_LB_duals[g] += (1 - alpha_weight) * all_new_LB_duals[all_new_LB_duals.size() - 1][g];
			}

			time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
			if (time_left < 0.001) {
				cout << "time out and exit" << endl;
				break;
			}

			// solve pricing problem for each available segemnt size
			///////////////////////////////////////

			vector<vector<int>> xx; //xx[g][n] number of items of size n from group g, so far selected 
			xx.resize(inst->n_subGroup);
			for (int g = 0; g < inst->n_subGroup; ++g) {
				xx[g].resize(inst->max_party_size, 0);
			}

			int total_obj = 0;

			//Adding objective function 
			GRBLinExpr objective = 0.0;
			for (int g = 0; g < inst->n_subGroup; ++g) {
				objective += y[g];
			}
			for (int g = 0; g < inst->n_subGroup; ++g) {
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					for (int i = 1; i < inst->u[g][nn] + 1; ++i) {
						objective -= party_duals[g][nn] * i * x[g][nn][i];
					}
				}
			}
			objective -= capacity_duals[inst->no_seat_in_seg[0] - 1];
			for (int g = 0; g < inst->n_subGroup; ++g) {
				objective -= (seg_UB_duals[g]) * y[g];
				objective -= new_LB_duals[g] * y[g];
			}
			// Using Callback to add the avoid constraints (i.e. to avoid from regenerating the existing columns)
			callback_pricing cb = callback_pricing(y, x, variables, inst, inst->no_seat_in_seg[0]);
			model_p.setCallback(&cb);

			//calculating the total number of items selected
			GRBLinExpr Total_item = 0.0;
			for (int g = 0; g < inst->n_subGroup; ++g) {
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					for (int i = 1; i < inst->u[g][nn] + 1; ++i) {
						Total_item += (nn + 1) * i * x[g][nn][i];
					}
				}
			}

			double last_RC = 0.0;


			vector<GRBConstr >  pricing_item_const;
			for (int iter = 0; iter < inst->n_seg; ++iter) { //iterations for solving the pricing problem



															 //adding objective[seg]
															 //IloObjective obj;
				model_p.update();
				model_p.setObjective(objective, GRB_MINIMIZE);
				model_p.getEnv().set(GRB_IntParam_Threads, 1);
				model_p.getEnv().set(GRB_DoubleParam_TimeLimit, time_limit);
				model_p.getEnv().set(GRB_IntParam_OutputFlag, 0);
				model_p.getEnv().set(GRB_IntParam_PreCrush, 1);
				//model.getEnv().set(GRB_DoubleParam_Cutoff, PRECISION);
				model_p.getEnv().set(GRB_IntParam_DualReductions, 0);
				model_p.getEnv().set(GRB_IntParam_LazyConstraints, 1);
				//model_p.getEnv().set(GRB_IntParam_Presolve, -1);
				model_p.update();
				model_p.optimize();

				n_price_iter += 1;
				total_price_time += model_p.get(GRB_DoubleAttr_Runtime);

				//cout << "Solution status: " << model.get(GRB_IntAttr_Status) << endl;
				//cout << ". Best objective: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
				//cout << "Pricing Time" << model.get(GRB_DoubleAttr_Runtime) << endl;

				//cout << "Time Taken for pricing iteration " << iter << ": " << model.get(GRB_DoubleAttr_Runtime) << endl;

				if (model_p.get(GRB_IntAttr_Status) == 2) {

					if (iter > 0 || model_p.get(GRB_DoubleAttr_ObjVal) < -PRECISION) {

						//cout << "Reduced Cost: " << model.get(GRB_DoubleAttr_ObjVal) << "\t";

						int ocup_cap = 0;

						Variable* p_var;
						p_var = new Variable();
						p_var->index = variables.size();
						p_var->assigned_group_size.resize(inst->n_subGroup);
						int column_repreat = 0; //number of times we want to repeat adding the current column because it may need to selected multiple times
						for (int g = 0; g < inst->n_subGroup; ++g) {
							p_var->assigned_group_size[g].resize(inst->max_party_size, 0);
							for (int nn = 0; nn < inst->max_party_size; ++nn) {
								if (inst->u[g][nn] > 0.5) {
									for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
										if (x[g][nn][i].get(GRB_DoubleAttr_X) > 0.9) {
											p_var->assigned_group_size[g][nn] = int(i);
											ocup_cap += i *nn;

											if (i > 0) {
												if (floor(inst->u[g][nn] / i) > column_repreat) {
													column_repreat = floor(inst->u[g][nn] / i);
												}
											}
										}
									}
								}
							}
						}
						p_var->capacity = inst->no_seat_in_seg[0];
						//computing value
						p_var->value = 0;
						for (int g = 0; g < inst->n_subGroup; ++g) {
							int n_party = 0;
							for (int nn = 0; nn < inst->max_party_size; ++nn) {
								if (x[g][nn][0].get(GRB_DoubleAttr_X) < 0.5) {
									n_party++;
								}
							}
							if (n_party > 0.9) {
								p_var->value++;
							}
						}
						int value1 = p_var->value;

						//adding the new variable	
						if (p_var->value > 0.5) {
							variables.push_back(p_var);
							for (int col = 0; col < column_repreat - 1; ++col) {
								variables.push_back(p_var);
							}
						}

						need_another_rmp_lp_solution = true;

						//cout << "occupied capacity: " << ocup_cap << endl;

						//updating the total objective
						total_obj += p_var->value;

						//updating xx[g][n]
						for (int g = 0; g < inst->n_subGroup; ++g) {
							for (int nn = 0; nn < inst->max_party_size; ++nn) {
								xx[g][nn] += p_var->assigned_group_size[g][nn];
								//adding new constraints to the model
								for (int i = inst->u[g][nn] - xx[g][nn] + 1; i < inst->u[g][nn] + 1; ++i) {
									pricing_item_const.push_back(model_p.addConstr(x[g][nn][i] == 0));
								}
							}
						}

						if (model_p.get(GRB_DoubleAttr_ObjVal) > -PRECISION) {
							if (iter > 0 && last_RC > -PRECISION) {
								objective += Total_item * abs(last_RC) / inst->Tn_party;
							}
							objective -= Total_item * abs(model_p.get(GRB_DoubleAttr_ObjVal)) / inst->Tn_party;
							last_RC = model_p.get(GRB_DoubleAttr_ObjVal);
						}

					}

					else {
						iter = inst->n_seg;
					}

				}


			}

			for (int con = 0; con < pricing_item_const.size(); con++) {
				model_p.remove(pricing_item_const[con]);
			}

			bool feasible = true;
			for (int g = 0; g < inst->n_subGroup; ++g) {
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					if (inst->u[g][nn] != xx[g][nn]) {
						feasible = false;
					}
				}
			}

			if (feasible) {
				//cout << "FEASIBLE!!!!!!!!!!!!!" << endl;
				//cout << "Objective value: " << total_obj << endl;

				if (total_obj < global_upper_bound) {
					global_upper_bound = total_obj;
					if (global_upper_bound <= global_lower_bound + 0.99) {
						break;
					}
				}
			}

			//////////////////////////////////////////

			time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
			if (time_left < 0.001) {
				cout << "time out and exit" << endl;
				break;
			}

			if (need_another_rmp_lp_solution == false) {

				party_duals = all_party_duals.back();
				capacity_duals = all_capacity_duals.back();
				//seg_LB_duals = all_seg_LB_duals.back();
				seg_UB_duals = all_seg_UB_duals.back();
				new_LB_duals = all_new_LB_duals.back();

				for (int cap = 0; cap < inst->n_SeatRow; ++cap) {
					if (inst->no_seg_of_each_size[cap] > 0) {
						solve_pricing_problem(node, time_left, file1, party_duals, capacity_duals, seg_UB_duals, new_LB_duals, cap + 1, need_another_rmp_lp_solution, n_price_iter, total_price_time);
						//solve_MDD_pricing_problem(node, time_left, file1, party_duals, capacity_duals, seg_LB_duals, seg_UB_duals, cap + 1, need_another_rmp_lp_solution, n_price_iter, total_price_time);
					}
					time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
					if (time_left < 0.001) {
						cout << "time out and exit" << endl;
						break;
					}
				}
			}

			//cout << "Number of columns: " << variables.size() << endl;


		} // end of one column generation iteration

		double LP_Time = double(clock() - LP_startTime) / (double)CLOCKS_PER_SEC;
		//file1 << n_lps_in_col_gen << ",";
		//file1 << LP_Time << ",";

		///////////////////////////////////////////////////////////////////

		delete x;
		delete y;

	} //end of try
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		exit(1);
	}
	catch (...) {
		cout << "Other error ... " << endl;
		exit(1);
	}
	/////////////////////////////////////////////////////////////

}

void BPSolver::solve_LP_with_pricing_and_DUALstabil2(SearchNode* node, double & time_left, double time_limit, clock_t startTime, ofstream& file1, int& n_CG_iter, int& n_price_iter, double& total_price_time) {


	char varname[256];


	// creating the base of the pricing model
	GRBEnv env_price;
	try {
		GRBModel model_p(env_price);

		////////////////////////////////////////////////////////////
		//Defining and adding variables
		GRBVar*** x;																		//x[g][n][i] : Number of party of size n from group g assigned to this column is equal to i
		x = new GRBVar**[inst->n_subGroup];
		for (int g = 0; g < inst->n_subGroup; g++) {
			x[g] = new GRBVar*[inst->max_party_size];
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				x[g][nn] = new GRBVar[inst->u[g][nn] + 1];
				for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
					x[g][nn][i] = model_p.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				}
			}
		}

		GRBVar* y;																		//y[g] is equal to one if group "g" has some party in this column, and zero otherwise
		y = new GRBVar[inst->n_subGroup];
		for (int g = 0; g < inst->n_subGroup; ++g) {
			y[g] = model_p.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}

		model_p.update();

		//adding constraints

		//sum_i x[g][nn][i] = 1
		for (int g = 0; g < inst->n_subGroup; g++) {
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				GRBLinExpr equ = 0.0;
				for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
					equ += x[g][nn][i];
				}
				model_p.addConstr(equ == 1);
			}
		}

		////i*x[g][nn][i] <= number of parties of size nn+1 from group g
		//for (int g = 0; g < inst->n_subGroup; g++) {
		//	GRBLinExpr equ11 = 0.0;
		//	for (int nn = 0; nn < inst->max_party_size; ++nn) {
		//		equ11 += (x[g][nn][0] - 1);
		//	}
		//	model_p.addConstr(inst->max_party_size * y[g] + equ11 >= 0);
		//}

		//x[g][nn]<= y[g] for all group g and size nn
		for (int g = 0; g < inst->n_subGroup; ++g) {
			//GRBLinExpr Yequ = 0.0;
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				model_p.addConstr(1 - x[g][nn][0] <= y[g]);
				//Yequ += (1 - x[g][nn][0]);
			}
			//model_p.addConstr(y[g] <= Yequ);
		}


		// \sum_{size} \sum_{i} x[g][nn]<= capacity * y[g] for all group g and size nn
		for (int g = 0; g < inst->n_subGroup; ++g) {
			GRBLinExpr equ_new = 0.0;
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
					equ_new += (nn + 1) * i * x[g][nn][i];
				}
			}
			model_p.addConstr(equ_new <= inst->no_seat_in_seg[0] * y[g]);
		}

		//Knapsack constraint
		GRBLinExpr cap_equ = 0.0;
		for (int g = 0; g < inst->n_subGroup; ++g) {
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
					cap_equ += (nn + 1) * i * x[g][nn][i];
				}
			}
		}
		model_p.addConstr(cap_equ <= inst->no_seat_in_seg[0]);

		///////////////////////////////////////////////////////////////////

		bool need_another_rmp_lp_solution = true;
		int n_lps_in_col_gen = 0;
		clock_t LP_startTime = clock();

		// the following vectors are defined to record all the duals for "dual price smoothing" method
		double alpha_weight = 0.5;     //The alpha in the Neame (2000) method
		vector<vector<vector<double >>> all_party_duals;
		vector<vector<double>> all_capacity_duals;
		//vector<vector<double>> all_seg_LB_duals;
		vector<vector<double>> all_seg_UB_duals;
		vector<vector<double>> all_new_LB_duals;

		while (need_another_rmp_lp_solution) {

			cout << "CG_iter " << n_CG_iter << endl;

			if (global_upper_bound <= global_lower_bound + 0.99) {
				time_left = 0;
				break;
			}

			n_CG_iter += 1;

			time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
			if (time_left < 0.001) {
				//cout << "time out and exit" << endl;
				break;
			}

			n_lps_in_col_gen++;
			//cout << "LOOP IN COL GEN: " << n_lps_in_col_gen << endl;
			//if (n_lps_in_col_gen > 1000) {
			//	cout << "solve_LP: Exceeded n_lps_in_col_gen \n";
			//	exit(1);
			//}
			need_another_rmp_lp_solution = false;

			vector<vector<double >> party_duals(inst->n_subGroup);	// one for every party
			for (int g = 0; g < inst->n_subGroup; ++g) {
				party_duals[g].resize(inst->max_party_size);
			}
			vector<double > capacity_duals(inst->n_SeatRow);	// one for every size from 1 to n_SeatRow
																//vector<double> seg_LB_duals(inst->n_subGroup);
			vector<double> seg_UB_duals(inst->n_subGroup);
			vector<double> new_LB_duals(inst->n_subGroup);

			// solve restricted MP
			time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
			if (time_left < 0.001) {
				cout << "time out and exit" << endl;
				break;
			}
			solve_RMP_LP(node, time_left, file1, party_duals, capacity_duals, seg_UB_duals, new_LB_duals);


			//calculating the smoothed duals by Neame (2000) method
			//double alpha_weight = 0.3;     //The alpha in the Neame (2000) method
			all_party_duals.push_back(party_duals);
			all_capacity_duals.push_back(capacity_duals);
			//all_seg_LB_duals.push_back(seg_LB_duals);
			all_seg_UB_duals.push_back(seg_UB_duals);
			all_new_LB_duals.push_back(new_LB_duals);

			for (int g = 0; g < inst->n_subGroup; ++g) {
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					party_duals[g][nn] = 0;
					for (int iter = 0; iter < all_party_duals.size() - 1; ++iter) {
						party_duals[g][nn] += alpha_weight * (1 - alpha_weight) * pow(alpha_weight, (all_party_duals.size() - iter)) * all_party_duals[iter][g][nn];
					}
					party_duals[g][nn] += (1 - alpha_weight) * all_party_duals[all_party_duals.size() - 1][g][nn];
				}
			}

			for (int i = 0; i < inst->n_SeatRow; ++i) {
				capacity_duals[i] = 0;
				for (int iter = 0; iter < all_capacity_duals.size() - 1; ++iter) {
					capacity_duals[i] += alpha_weight * (1 - alpha_weight) * pow(alpha_weight, (all_capacity_duals.size() - iter)) * all_capacity_duals[iter][i];
				}
				capacity_duals[i] += (1 - alpha_weight) * all_capacity_duals[all_capacity_duals.size() - 1][i];
			}

			//for (int g = 0; g < inst->n_subGroup; ++g) {
			//	seg_LB_duals[g] = 0;
			//	for (int iter = 0; iter < all_seg_LB_duals.size() - 1; ++iter) {
			//		seg_LB_duals[g] += alpha_weight * (1 - alpha_weight) * pow(alpha_weight, all_seg_LB_duals.size() - iter) * all_seg_LB_duals[iter][g];
			//	}
			//	seg_LB_duals[g] += (1 - alpha_weight) * all_seg_LB_duals[all_seg_LB_duals.size() - 1][g];
			//}

			for (int g = 0; g < inst->n_subGroup; ++g) {
				seg_UB_duals[g] = 0;
				for (int iter = 0; iter < all_seg_UB_duals.size() - 1; ++iter) {
					seg_UB_duals[g] += alpha_weight * (1 - alpha_weight) * pow(alpha_weight, (all_seg_UB_duals.size() - iter)) * all_seg_UB_duals[iter][g];
				}
				seg_UB_duals[g] += (1 - alpha_weight) * all_seg_UB_duals[all_seg_UB_duals.size() - 1][g];
			}

			for (int g = 0; g < inst->n_subGroup; ++g) {
				new_LB_duals[g] = 0;
				for (int iter = 0; iter < all_new_LB_duals.size() - 1; ++iter) {
					new_LB_duals[g] += alpha_weight * (1 - alpha_weight) * pow(alpha_weight, all_new_LB_duals.size() - iter) * all_new_LB_duals[iter][g];
				}
				new_LB_duals[g] += (1 - alpha_weight) * all_new_LB_duals[all_new_LB_duals.size() - 1][g];
			}

			time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
			if (time_left < 0.001) {
				cout << "time out and exit" << endl;
				break;
			}

			// solve pricing problem for each available segemnt size
			///////////////////////////////////////

			vector<vector<int>> xx; //xx[g][n] number of items of size n from group g, so far selected 
			xx.resize(inst->n_subGroup);
			for (int g = 0; g < inst->n_subGroup; ++g) {
				xx[g].resize(inst->max_party_size, 0);
			}

			int total_obj = 0;

			//Adding objective function 
			GRBLinExpr objective = 0.0;
			for (int g = 0; g < inst->n_subGroup; ++g) {
				objective += y[g];
			}
			for (int g = 0; g < inst->n_subGroup; ++g) {
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					for (int i = 1; i < inst->u[g][nn] + 1; ++i) {
						objective -= party_duals[g][nn] * i * x[g][nn][i];
					}
				}
			}
			objective -= capacity_duals[inst->no_seat_in_seg[0] - 1];
			for (int g = 0; g < inst->n_subGroup; ++g) {
				objective -= (seg_UB_duals[g]) * y[g];
				objective -= new_LB_duals[g] * y[g];
			}

			//calculating the total number of items selected
			GRBLinExpr Total_item = 0.0;
			for (int g = 0; g < inst->n_subGroup; ++g) {
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					for (int i = 1; i < inst->u[g][nn] + 1; ++i) {
						Total_item += (nn + 1) * i * x[g][nn][i];
					}
				}
			}

			double last_RC = 0.0;


			vector<GRBConstr >  pricing_item_const;
			for (int iter = 0; iter < inst->n_seg; ++iter) { //iterations for solving the pricing problem
				
															 // Using Callback to add the avoid constraints (i.e. to avoid from regenerating the existing columns)
				if (iter == 0) {
					callback_pricing cb = callback_pricing(y, x, variables, inst, inst->no_seat_in_seg[0]);
					model_p.setCallback(&cb);
				}


				//adding objective[seg]
				//IloObjective obj;
				model_p.update();
				model_p.setObjective(objective, GRB_MINIMIZE);
				model_p.getEnv().set(GRB_IntParam_Threads, 1);
				model_p.getEnv().set(GRB_DoubleParam_TimeLimit, time_limit);
				model_p.getEnv().set(GRB_IntParam_OutputFlag, 0);
				//model.getEnv().set(GRB_DoubleParam_Cutoff, PRECISION);
				//model_p.getEnv().set(GRB_IntParam_DualReductions, 0);
				if (iter >= 0) {
					model_p.getEnv().set(GRB_IntParam_PreCrush, 1);
					model_p.getEnv().set(GRB_IntParam_LazyConstraints, 1);
				}
				//model_p.getEnv().set(GRB_IntParam_Presolve, -1);
				model_p.update();
				model_p.optimize();


				n_price_iter += 1;
				total_price_time += model_p.get(GRB_DoubleAttr_Runtime);

				//cout << "Solution status: " << model.get(GRB_IntAttr_Status) << endl;
				//cout << ". Best objective: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
				//cout << "Pricing Time" << model.get(GRB_DoubleAttr_Runtime) << endl;

				//cout << "Time Taken for pricing iteration " << iter << ": " << model.get(GRB_DoubleAttr_Runtime) << endl;

				if (model_p.get(GRB_IntAttr_Status) == 2) {

					if (iter > 0 || model_p.get(GRB_DoubleAttr_ObjVal) < -PRECISION) {

						//cout << "Reduced Cost: " << model.get(GRB_DoubleAttr_ObjVal) << "\t";

						int ocup_cap = 0;

						Variable* p_var;
						p_var = new Variable();
						p_var->index = variables.size();
						p_var->assigned_group_size.resize(inst->n_subGroup);
						int column_repreat = 0; //number of times we want to repeat adding the current column because it may need to selected multiple times
						for (int g = 0; g < inst->n_subGroup; ++g) {
							p_var->assigned_group_size[g].resize(inst->max_party_size, 0);
							for (int nn = 0; nn < inst->max_party_size; ++nn) {
								if (inst->u[g][nn] > 0.5) {
									for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
										if (x[g][nn][i].get(GRB_DoubleAttr_X) > 0.9) {
											p_var->assigned_group_size[g][nn] = int(i);
											ocup_cap += i *nn;

											if (i > 0) {
												if (floor(inst->u[g][nn] / i) > column_repreat) {
													column_repreat = floor(inst->u[g][nn] / i);
												}
											}
										}
									}
								}
							}
						}
						p_var->capacity = inst->no_seat_in_seg[0];
						//computing value
						p_var->value = 0;
						for (int g = 0; g < inst->n_subGroup; ++g) {
							int n_party = 0;
							for (int nn = 0; nn < inst->max_party_size; ++nn) {
								if (x[g][nn][0].get(GRB_DoubleAttr_X) < 0.5) {
									n_party++;
								}
							}
							if (n_party > 0.9) {
								p_var->value++;
							}
						}
						int value1 = p_var->value;

						//adding the new variable	
						if (p_var->value > 0.5) {
							variables.push_back(p_var);
							for (int col = 0; col < column_repreat - 1; ++col) {
								variables.push_back(p_var);
							}
						}

						need_another_rmp_lp_solution = true;

						//cout << "occupied capacity: " << ocup_cap << endl;

						//updating the total objective
						total_obj += p_var->value;

						//updating xx[g][n]
						for (int g = 0; g < inst->n_subGroup; ++g) {
							for (int nn = 0; nn < inst->max_party_size; ++nn) {
								xx[g][nn] += p_var->assigned_group_size[g][nn];
								//adding new constraints to the model
								for (int i = inst->u[g][nn] - xx[g][nn] + 1; i < inst->u[g][nn] + 1; ++i) {
									pricing_item_const.push_back(model_p.addConstr(x[g][nn][i] == 0));
								}
							}
						}

						/*if (model_p.get(GRB_DoubleAttr_ObjVal) > -PRECISION) {
						if (iter > 0 && last_RC > -PRECISION) {
						objective += Total_item * abs(last_RC) / inst->Tn_party;
						}
						objective -= Total_item * abs(model_p.get(GRB_DoubleAttr_ObjVal)) / inst->Tn_party;
						last_RC = model_p.get(GRB_DoubleAttr_ObjVal);
						}*/
						objective = 0;
						for (int g = 0; g < inst->n_subGroup; ++g) {
							objective += y[g];
						}
						objective -= Total_item / inst->no_seat_in_seg[0];
					}

					else {
						iter = inst->n_seg;
					}

				}


			}

			for (int con = 0; con < pricing_item_const.size(); con++) {
				model_p.remove(pricing_item_const[con]);
			}

			bool feasible = true;
			for (int g = 0; g < inst->n_subGroup; ++g) {
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					if (inst->u[g][nn] != xx[g][nn]) {
						feasible = false;
					}
				}
			}

			if (feasible) {
				//cout << "FEASIBLE!!!!!!!!!!!!!" << endl;
				//cout << "Objective value: " << total_obj << endl;

				if (total_obj < global_upper_bound) {
					global_upper_bound = total_obj;
					cout << "LB , UB: " << global_lower_bound << " " << global_upper_bound;
					system("pause");
					if (global_upper_bound <= global_lower_bound + 0.99) {
						break;
					}
				}
			}

			//////////////////////////////////////////

			time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
			if (time_left < 0.001) {
				cout << "time out and exit" << endl;
				break;
			}

			if (need_another_rmp_lp_solution == false) {

				party_duals = all_party_duals.back();
				capacity_duals = all_capacity_duals.back();
				//seg_LB_duals = all_seg_LB_duals.back();
				seg_UB_duals = all_seg_UB_duals.back();
				new_LB_duals = all_new_LB_duals.back();

				for (int cap = 0; cap < inst->n_SeatRow; ++cap) {
					if (inst->no_seg_of_each_size[cap] > 0) {
						solve_pricing_problem(node, time_left, file1, party_duals, capacity_duals, seg_UB_duals, new_LB_duals, cap + 1, need_another_rmp_lp_solution, n_price_iter, total_price_time);
						//solve_MDD_pricing_problem(node, time_left, file1, party_duals, capacity_duals, seg_LB_duals, seg_UB_duals, cap + 1, need_another_rmp_lp_solution, n_price_iter, total_price_time);
					}
					time_left = time_limit - (double(clock() - startTime) / (double)CLOCKS_PER_SEC);
					if (time_left < 0.001) {
						cout << "time out and exit" << endl;
						break;
					}
				}
			}

			//cout << "Number of columns: " << variables.size() << endl;


		} // end of one column generation iteration

		double LP_Time = double(clock() - LP_startTime) / (double)CLOCKS_PER_SEC;
		//file1 << n_lps_in_col_gen << ",";
		//file1 << LP_Time << ",";

		///////////////////////////////////////////////////////////////////

		delete x;
		delete y;

		cout << "start of " << n_CG_iter -1  << endl;

	} //end of try
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		exit(1);
	}
	catch (...) {
		cout << "Other error ... " << endl;
		exit(1);
	}
	/////////////////////////////////////////////////////////////

}

void BPSolver::solve_IP(SearchNode* node, double & time_limit, ofstream& file1) {

	//cout << "=============================================" << endl << "Solving IP Heuristic Model:" << endl << "=============================================";
	//GRBEnv env;
	char varname[256];

	//system("pause");

	try {
		GRBModel model(env);

		// Define variables

		GRBVar* x;
		x = new GRBVar[(int)variables.size()];
		for (int v = 0; v < (int)variables.size(); ++v) {
			//sprintf(varname, "x[%d]", v);
			x[v] = model.addVar(0.0, 1, 0.0, GRB_BINARY);
		}

		model.update();         // required before you can use the variables

								// adding constraints

								// \sum_c m[c][g][nn] x_c == u[g][nn] for all group g and size nn 
		vector<vector<GRBConstr >>  group_size_constraints(inst->n_subGroup);
		for (int g = 0; g < inst->n_subGroup; ++g) {
			group_size_constraints[g].resize(inst->max_party_size);
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				GRBLinExpr group_size_con = 0.0;
				for (int c = 0; c < (int)variables.size(); ++c) {
					group_size_con += variables[c]->assigned_group_size[g][nn] * x[c];
				}
				group_size_constraints[g][nn] = model.addConstr(group_size_con >= inst->u[g][nn]);
			}
		}

		// number of selected columns from each capacity must be <= number of segemnets available of that capacity
		vector<GRBConstr>  capacity_constraints(inst->n_SeatRow);
		for (int cap = 0; cap < inst->n_SeatRow; ++cap) {
			GRBLinExpr capacity_con = 0.0;
			for (int c = 0; c < (int)variables.size(); ++c) {
				if (variables[c]->capacity == cap + 1) {
					capacity_con += x[c];
				}
			}
			capacity_constraints[cap] = model.addConstr(capacity_con <= inst->no_seg_of_each_size[cap]);
		}

		// LB and UB constraint for number of segments assigned to each group
		//vector<GRBConstr >  seg_LB_constraints(inst->n_subGroup);
		vector<GRBConstr >  seg_UB_constraints(inst->n_subGroup);
		vector<GRBConstr >  new_LB_constraints(inst->n_subGroup);
		for (int g = 0; g < inst->n_subGroup; ++g) {
			GRBLinExpr n_seg_assigned_to_g = 0.0;
			for (int c = 0; c < (int)variables.size(); ++c) {
				int n_party = 0.0;
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					if (variables[c]->assigned_group_size[g][nn] > 0.9) {
						n_party++;
					}
				}
				if (n_party > 0.9) {
					n_seg_assigned_to_g += x[c];
				}
			}
			//seg_LB_constraints[g] = model.addConstr(n_seg_assigned_to_g >= node->LB_seg_no_for_group[g]);
			seg_UB_constraints[g] = model.addConstr(n_seg_assigned_to_g <= node->UB_seg_no_for_group[g]);
			double LB_n_bin = inst->nn[g] / inst->n_SeatRow;
			if (LB_n_bin - floor(LB_n_bin) >= 0.001) {
				LB_n_bin = floor(LB_n_bin) + 1;
			}
			else {
				LB_n_bin = floor(LB_n_bin);
			}
			double new_LB = LB_n_bin;
			if (inst->LB_MT1[g] > new_LB) {
				new_LB = inst->LB_MT1[g];
			}
			if (inst->LB_MT2[g] > new_LB) {
				new_LB = inst->LB_MT2[g];
			}
			if (node->LB_seg_no_for_group[g] > new_LB) {
				new_LB = node->LB_seg_no_for_group[g];
			}
			new_LB_constraints[g] = model.addConstr(n_seg_assigned_to_g >= new_LB);
		}

		// for these two constraints, you might want to using vector iterators, and
		// find the position in the vector.  We will, however, make sure the index matches
		// the relateive position
		for (int i = 0; i < (int)node->variables_set_to_one.size(); ++i) {
			model.addConstr(x[node->variables_set_to_one[i]->index] == 1);
		}

		for (int i = 0; i < (int)node->variables_set_to_zero.size(); ++i) {
			model.addConstr(x[node->variables_set_to_zero[i]->index] == 0);
		}

		// objective function
		GRBLinExpr objective = 0.0;
		for (int c = 0; c < (int)variables.size(); ++c) {
			objective += variables[c]->value * x[c];
		}

		model.addConstr(objective >= global_lower_bound);
		model.addConstr(objective <= global_upper_bound - 1);

		model.setObjective(objective, GRB_MINIMIZE);

		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_DoubleParam_TimeLimit, time_limit);
		model.getEnv().set(GRB_IntParam_OutputFlag, 0);
		model.update();
		model.optimize();

		//cout << endl;
		//cout << "Solution Status: " << model.get(GRB_IntAttr_Status) << endl;

		if (model.get(GRB_IntAttr_Status) == 2) {
			//cout << "Time Taken: " << model.get(GRB_DoubleAttr_Runtime) << endl;
			// This gets you the relaxation bound
			//cout << "Relaxation Bound: " << model.get(GRB_DoubleAttr_ObjBound) << endl;
			// This gets you the best known primal solution value
			//cout << "Best Known Objective : " << model.get(GRB_DoubleAttr_ObjVal) << endl;
			//file1 << model.get(GRB_DoubleAttr_Runtime) << ",";
			//file1 << model.get(GRB_DoubleAttr_ObjVal) << ",";

			/*cout << "x[c]: " << endl;
			for (int c = 0; c < variables.size(); ++c) {
			cout << x[c].get(GRB_DoubleAttr_X) << "\t";
			}
			cout << endl;

			cout << "selected columns: " << endl;
			for (int c = 0; c < variables.size(); ++c) {
			if (x[c].get(GRB_DoubleAttr_X) > 0.5) {
			cout << "column[" << c << "]:" << endl;
			print_variable(variables[c]);
			}
			}
			cout << endl;*/

			if (model.get(GRB_DoubleAttr_ObjVal) < global_upper_bound) {
				global_upper_bound = model.get(GRB_DoubleAttr_ObjVal);
				//
				//cout << "******** found better solution: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
				best_solution.clear();
				for (int c = 0; c < (int)variables.size(); ++c) {
					if (x[c].get(GRB_DoubleAttr_X) > 0.9) {
						best_solution.push_back(variables[c]);
						////Printing the best solution
						//print_variable(variables[c]);
					}
				}

				//Printing the best solution


				////Printing features of the best solution 
				//file1 << "******** Found better solution: " << model.get(GRB_DoubleAttr_ObjVal) << " ********" << endl;
				//for (int g = 0; g < inst->n_subGroup; ++g) {
				//	file1 << inst->n_party_in_group[g] << " / " << inst->nn[g] << " / ";
				//	int n_assigned_seg = 0;
				//	int n_party_assigned_seg = 0;
				//	int n_people_assigned_seg = 0;
				//	for (int c = 0; c < (int)variables.size(); ++c) {
				//		n_party_assigned_seg = 0;
				//		n_people_assigned_seg = 0;
				//		if (x[c].get(GRB_DoubleAttr_X) > 0.9) {
				//			for (int p = 0; p < inst->Tn_party; ++p) {
				//				if (variables[c]->assigned_party[p] && inst->subG_party[p] == g) {
				//					n_party_assigned_seg++;
				//					n_people_assigned_seg += inst->n[p];
				//				}
				//			}
				//			if (n_party_assigned_seg > 0) {
				//				n_assigned_seg++;
				//			}
				//		}
				//	}
				//	file1 << n_assigned_seg;
				//	n_party_assigned_seg = 0;
				//	n_people_assigned_seg = 0;
				//	for (int c = 0; c < (int)variables.size(); ++c) {
				//		n_party_assigned_seg = 0;
				//		n_people_assigned_seg = 0;
				//		if (x[c].get(GRB_DoubleAttr_X) > 0.9) {
				//			for (int p = 0; p < inst->Tn_party; ++p) {
				//				if (variables[c]->assigned_party[p] && inst->subG_party[p] == g) {
				//					n_party_assigned_seg++;
				//					n_people_assigned_seg += inst->n[p];
				//				}
				//			}
				//			if (n_party_assigned_seg > 0) {
				//				file1 << " / " << "(" << n_party_assigned_seg << "," << n_people_assigned_seg << ")";
				//			}
				//		}
				//	}

				//	file1 << endl;
				//}
				//file1 << "*******************************************" << endl;
				////End of printing the best solution



			}
		}

		delete x;

	}
	catch (GRBException e) {
		cout << "solveIP: Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		exit(1);
	}
	catch (...) {
		cout << "solveIP: Other error ... " << endl;
		exit(1);
	}

}

void BPSolver::create_first_set_of_variables() {

	GRBEnv env;

	try {

		GRBModel model(env);

		//Defining and adding variables

		GRBVar**** x;																		//x[s][g][n][i]=1 : Number of party of size n from group g assigned to this column is equal to i
		x = new GRBVar***[inst->n_seg];
		for (int s = 0; s < inst->n_seg; s++) {
			x[s] = new GRBVar**[inst->n_subGroup];
			for (int g = 0; g < inst->n_subGroup; g++) {
				x[s][g] = new GRBVar*[inst->max_party_size];
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					x[s][g][nn] = new GRBVar[inst->u[g][nn] + 1];
					for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
						x[s][g][nn][i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					}
				}
			}
		}

		GRBVar** y;																		//y[g][seg] is equal to one if group "g" has some party in segment "seg", and zero otherwise
		y = new GRBVar*[inst->n_subGroup];
		for (int g = 0; g < inst->n_subGroup; ++g) {
			y[g] = new GRBVar[inst->n_seg];
			for (int seg = 0; seg < inst->n_seg; ++seg) {
				y[g][seg] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			}
		}

		model.update();
		//adding constraints

		//sum_i x[s][g][nn][i] = 1 for each segment s, group g, and party size nn
		for (int seg = 0; seg < inst->n_seg; ++seg) {
			for (int g = 0; g < inst->n_subGroup; g++) {
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					GRBLinExpr equ = 0.0;
					for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
						equ += x[seg][g][nn][i];
					}
					model.addConstr(equ == 1);
				}
			}
		}

		//sum_seg,i i*x[g][nn][i] == number of parties of size nn+1 from group g
		for (int g = 0; g < inst->n_subGroup; g++) {
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				GRBLinExpr equ = 0.0;
				for (int seg = 0; seg < inst->n_seg; ++seg) {
					for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
						equ += i * x[seg][g][nn][i];
					}
				}
				model.addConstr(equ == inst->u[g][nn]);
			}
		}

		//x[g][nn]<= y[g] for all group g and size nn
		for (int s = 0; s < inst->n_seg; ++s) {
			for (int g = 0; g < inst->n_subGroup; ++g) {
				GRBLinExpr Yequ = 0.0;
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					model.addConstr(1 - x[s][g][nn][0] <= y[g][s]);
					Yequ += (1 - x[s][g][nn][0]);
				}
				model.addConstr(y[g][s] <= Yequ);
			}
		}

		//Knapsack constraint
		for (int s = 0; s < inst->n_seg; ++s) {
			GRBLinExpr cap_equ = 0.0;
			for (int g = 0; g < inst->n_subGroup; ++g) {
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
						cap_equ += (nn + 1) * i * x[s][g][nn][i];
					}
				}
			}
			model.addConstr(cap_equ <= inst->no_seat_in_seg[s]);
		}

		///////////////////////////
		//Adding objective function 

		cout << "inst->n_subGroup " << inst->n_subGroup << endl;

		for (int group = 0; group < inst->n_subGroup; ++group) {

			cout << "color " << group << endl;

			GRBLinExpr objective = 0.0;
			for (int seg = 0; seg < inst->n_seg; ++seg) {
				objective += y[group][seg];
			}

			model.update();
			model.setObjective(objective, GRB_MINIMIZE);

			model.getEnv().set(GRB_IntParam_Threads, 1);
			model.getEnv().set(GRB_DoubleParam_TimeLimit, 10);
			model.getEnv().set(GRB_IntParam_OutputFlag, 0);
			model.optimize();

			//cout << endl;
			//cout << "Solution Status: " << model.get(GRB_IntAttr_Status) << endl;
			//cout << "Time Taken: " << model.get(GRB_DoubleAttr_Runtime) << endl;

			if (model.get(GRB_IntAttr_Status) == 2 || model.get(GRB_IntAttr_Status) == 9) {

				// This gets you the relaxation bound
				//cout << "Relaxation Bound: " << model.get(GRB_DoubleAttr_ObjBound) << endl;
				// This gets you the best known primal solution value
				//cout << "Best Known Objective : " << model.get(GRB_DoubleAttr_ObjVal) << endl;

				//cout << "Adding first set of variables:" << endl;
				// creating initial variables for column_generation:
				for (int seg = 0; seg < inst->n_seg; ++seg) {
					Variable* p_var;
					p_var = new Variable();
					p_var->index = variables.size();
					p_var->assigned_group_size.resize(inst->n_subGroup);
					for (int g = 0; g < inst->n_subGroup; ++g) {
						p_var->assigned_group_size[g].resize(inst->max_party_size, 0);
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
								if (x[seg][g][nn][i].get(GRB_DoubleAttr_X) > 0.9) {
									p_var->assigned_group_size[g][nn] = i;
								}
							}
						}
					}
					p_var->capacity = inst->no_seat_in_seg[seg];
					//computing value
					p_var->value = 0;
					for (int g = 0; g < inst->n_subGroup; ++g) {
						int n_party = 0;
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							if (x[seg][g][nn][0].get(GRB_DoubleAttr_X) < 0.1) {
								n_party++;
							}
						}
						if (n_party > 0.9) {
							p_var->value++;
						}
					}

					//adding the new variable	
					if (p_var->value > 0.5) {
						variables.push_back(p_var);
					}
				}

			}
		} // End of the loop for different groups

		delete x;
		delete y;

	} //end of try

	catch (GRBException e) {
		cout << "createFirstSet: Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		exit(1);
	}
	catch (...) {
		cout << "createFirstSet: Other error ... " << endl;
		exit(1);
	}

}

void BPSolver::create_first_set_of_variables_bin_completion(ofstream& file1) {

	//sum_i i*x[g][nn][i] <= number of parties of size nn+1 from group g
	vector<vector<int>> u_new = inst->u;

	GRBEnv env;

	try {

		GRBModel model(env);

		//Defining and adding variables
		GRBVar*** x;																		//x[g][n][i]=1 : Number of party of size n from group g assigned to one segment is equal to i
		x = new GRBVar**[inst->n_subGroup];
		for (int g = 0; g < inst->n_subGroup; g++) {
			x[g] = new GRBVar*[inst->max_party_size];
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				x[g][nn] = new GRBVar[inst->u[g][nn] + 1];
				for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
					x[g][nn][i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				}
			}
		}

		GRBVar* y;																		//y[g] is equal to one if group "g" has at least one party in the segment, and zero otherwise
		y = new GRBVar[inst->n_subGroup];
		for (int g = 0; g < inst->n_subGroup; ++g) {
			y[g] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}

		model.update();
		//adding constraints

		//sum_i x[g][nn][i] = 1 for each group g and party size nn
		for (int g = 0; g < inst->n_subGroup; g++) {
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				GRBLinExpr equ = 0.0;
				for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
					equ += x[g][nn][i];
				}
				model.addConstr(equ == 1);
			}
		}

		//sum_i i*x[g][nn][i] <= number of parties of size nn+1 from group g
		for (int g = 0; g < inst->n_subGroup; g++) {
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				GRBLinExpr equ = 0.0;
				for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
					equ += i * x[g][nn][i];
				}
				model.addConstr(equ <= inst->u[g][nn]);
			}
		}

		//x[g][nn]<= y[g] for all group g and size nn
		for (int g = 0; g < inst->n_subGroup; ++g) {
			GRBLinExpr Yequ = 0.0;
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				model.addConstr(1 - x[g][nn][0] <= y[g]);
				Yequ += (1 - x[g][nn][0]);
			}
			model.addConstr(y[g] <= Yequ);
		}

		//// \sum_{size} \sum_{i} x[g][nn]<= capacity * y[g] for all group g and size nn
		//for (int s = 0; s < inst->n_seg; ++s) {
		//	for (int g = 0; g < inst->n_subGroup; ++g) {
		//		GRBLinExpr equ_new = 0.0;
		//		for (int nn = 0; nn < inst->max_party_size; ++nn) {
		//			for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
		//				equ_new += (nn + 1) * i * x[s][g][nn][i];
		//			}
		//		}
		//		model.addConstr(equ_new <= inst->no_seat_in_seg[s] * y[g][s]);
		//	}
		//}

		//Knapsack constraint
		GRBLinExpr cap_equ = 0.0;
		for (int g = 0; g < inst->n_subGroup; ++g) {
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
					cap_equ += (nn + 1) * i * x[g][nn][i];
				}
			}
		}
		model.addConstr(cap_equ <= inst->n_SeatRow);

		///////////////////////////
		//Adding objective function 

		int total_obj = 0;

		for (int seg = 0; seg < inst->n_seg; ++seg) {

			GRBLinExpr objective = 0.0;
			for (int g = 0; g < inst->n_subGroup; ++g) {
				objective += y[g];
			}
			GRBLinExpr obj_equ = 0.0;
			for (int g = 0; g < inst->n_subGroup; g++) {
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
						obj_equ -= 0.001 * i * pow(nn, 3) * x[g][nn][i];
					}
				}
			}

			model.update();
			model.setObjective(objective, GRB_MINIMIZE);

			model.getEnv().set(GRB_IntParam_Threads, 1);
			model.getEnv().set(GRB_DoubleParam_TimeLimit, 5);
			model.getEnv().set(GRB_IntParam_OutputFlag, 0);
			model.optimize();

			//cout << endl;
			//cout << "Solution Status: " << model.get(GRB_IntAttr_Status) << endl;
			//cout << "Time Taken: " << model.get(GRB_DoubleAttr_Runtime) << endl;

			if (model.get(GRB_IntAttr_Status) == 2 || model.get(GRB_IntAttr_Status) == 9) {
				// This gets you the relaxation bound
				//cout << "Relaxation Bound: " << model.get(GRB_DoubleAttr_ObjBound) << endl;
				// This gets you the best known primal solution value
				//cout << "Best Known Objective : " << model.get(GRB_DoubleAttr_ObjVal) << endl;
				for (int g = 0; g < inst->n_subGroup; ++g) {
					total_obj += y[g].get(GRB_DoubleAttr_X);
				}

				//cout << "Adding first set of variables:" << endl;
				// creating initial variables for column_generation

				Variable* p_var;
				p_var = new Variable();
				p_var->index = variables.size();
				p_var->assigned_group_size.resize(inst->n_subGroup);
				int column_repeat = 0;
				for (int g = 0; g < inst->n_subGroup; ++g) {
					p_var->assigned_group_size[g].resize(inst->max_party_size, 0);
					for (int nn = 0; nn < inst->max_party_size; ++nn) {
						for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
							if (x[g][nn][i].get(GRB_DoubleAttr_X) > 0.9) {
								p_var->assigned_group_size[g][nn] = i;

								if (i > 0) {
									if (floor(inst->u[g][nn] / i) > column_repeat) {
										column_repeat = floor(inst->u[g][nn] / i);
									}
								}

							}
						}
					}
				}
				p_var->capacity = inst->no_seat_in_seg[seg];
				//computing value
				p_var->value = 0;
				for (int g = 0; g < inst->n_subGroup; ++g) {
					int n_party = 0;
					for (int nn = 0; nn < inst->max_party_size; ++nn) {
						if (x[g][nn][0].get(GRB_DoubleAttr_X) < 0.1) {
							n_party++;
						}
					}
					if (n_party > 0.9) {
						p_var->value++;
					}
				}
				//adding the new variable
				if (p_var->value > 0.5) {
					variables.push_back(p_var);
					for (int col = 0; col < column_repeat - 1; ++col) {
						variables.push_back(p_var);
					}
				}

				//Fixing the variables corresponding to segement "seg"
				for (int g = 0; g < inst->n_subGroup; ++g) {
					for (int nn = 0; nn < inst->max_party_size; ++nn) {
						int selected_i = 0;
						GRBLinExpr equ = 0.0;
						for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
							if (x[g][nn][i].get(GRB_DoubleAttr_X) > 0.9) {
								u_new[g][nn] -= int(i);
								selected_i = i;
							}
							equ += i * x[g][nn][i];
						}
						if (selected_i != 0) {
							model.addConstr(equ <= u_new[g][nn]);
						}
					}
				}

			}

		} // End of the loop for different segments

		int feasible_u_new = 0;
		for (int g = 0; g < inst->n_subGroup; g++) {
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
					feasible_u_new += u_new[g][nn];
				}
			}
		}
		if (feasible_u_new == 0) {
			if (total_obj < global_upper_bound) {
				global_upper_bound = total_obj;
				//file1 << model.get(GRB_DoubleAttr_ObjVal) << ",";
			}
		}

		delete x;
		delete y;

	} //end of try

	catch (GRBException e) {
		cout << "createFirstSet2:Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		exit(1);
	}
	catch (...) {
		cout << "createFirstSet2:Other error ... " << endl;
		exit(1);
	}

}

void BPSolver::create_first_set_of_variables2(ofstream& file1) {

	GRBEnv env;

	try {

		GRBModel model(env);

		//Defining and adding variables
		GRBVar**** x;																		//x[s][g][n][i]=1 : Number of party of size n from group g assigned to segment s is equal to i
		x = new GRBVar***[inst->n_seg];
		for (int s = 0; s < inst->n_seg; s++) {
			x[s] = new GRBVar**[inst->n_subGroup];
			for (int g = 0; g < inst->n_subGroup; g++) {
				x[s][g] = new GRBVar*[inst->max_party_size];
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					x[s][g][nn] = new GRBVar[inst->u[g][nn] + 1];
					for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
						x[s][g][nn][i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					}
				}
			}
		}

		//GRBVar** y;																		//y[g][seg] is equal to one if group "g" has at least one party in segment "seg", and zero otherwise
		//y = new GRBVar*[inst->n_subGroup];
		//for (int g = 0; g < inst->n_subGroup; ++g) {
		//	y[g] = new GRBVar[inst->n_seg];
		//	for (int seg = 0; seg < inst->n_seg; ++seg) {
		//		y[g][seg] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		//	}
		//}

		model.update();
		//adding constraints

		//sum_i x[s][g][nn][i] = 1 for each segment s, group g, and party size nn
		for (int seg = 0; seg < inst->n_seg; ++seg) {
			for (int g = 0; g < inst->n_subGroup; g++) {
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					GRBLinExpr equ = 0.0;
					for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
						equ += x[seg][g][nn][i];
					}
					model.addConstr(equ == 1);
				}
			}
		}

		//sum_seg,i i*x[g][nn][i] == number of parties of size nn+1 from group g
		for (int g = 0; g < inst->n_subGroup; g++) {
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				GRBLinExpr equ = 0.0;
				for (int seg = 0; seg < inst->n_seg; ++seg) {
					for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
						equ += i * x[seg][g][nn][i];
					}
				}
				model.addConstr(equ == inst->u[g][nn]);
			}
		}

		////x[g][nn]<= y[g] for all group g and size nn
		//for (int s = 0; s < inst->n_seg; ++s) {
		//	for (int g = 0; g < inst->n_subGroup; ++g) {
		//		GRBLinExpr Yequ = 0.0;
		//		for (int nn = 0; nn < inst->max_party_size; ++nn) {
		//			model.addConstr(1 - x[s][g][nn][0] <= y[g][s]);
		//			Yequ += (1 - x[s][g][nn][0]);
		//		}
		//		model.addConstr(y[g][s] <= Yequ);
		//	}
		//}



		//Knapsack constraint
		for (int s = 0; s < inst->n_seg; ++s) {
			GRBLinExpr cap_equ = 0.0;
			for (int g = 0; g < inst->n_subGroup; ++g) {
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
						cap_equ += (nn + 1) * i * x[s][g][nn][i];
					}
				}
			}
			model.addConstr(cap_equ <= inst->no_seat_in_seg[s]);
		}

		///////////////////////////
		//Adding objective function 

		int total_obj = 0;

		for (int seg = 0; seg < inst->n_seg; ++seg) {

			GRBLinExpr objective = 0.0;
			//for (int seg2 = 0; seg2 < seg + 1; seg2++) {
			//	for (int g = 0; g < inst->n_subGroup; ++g) {
			//		objective += y[g][seg2];
			//	}
			//}

			model.update();
			model.setObjective(objective, GRB_MINIMIZE);

			model.getEnv().set(GRB_IntParam_Threads, 1);
			model.getEnv().set(GRB_DoubleParam_TimeLimit, 10);
			model.getEnv().set(GRB_IntParam_MIPFocus, 1);
			model.getEnv().set(GRB_IntParam_OutputFlag, 0);
			model.optimize();

			//cout << endl;
			//cout << "Solution Status: " << model.get(GRB_IntAttr_Status) << endl;
			//cout << "Time Taken: " << model.get(GRB_DoubleAttr_Runtime) << endl;

			if (model.get(GRB_IntAttr_Status) == 2 || (model.get(GRB_IntAttr_Status) == 9 && model.get(GRB_DoubleAttr_ObjVal) > 0)) {
				// This gets you the relaxation bound
				//cout << "Relaxation Bound: " << model.get(GRB_DoubleAttr_ObjBound) << endl;
				// This gets you the best known primal solution value
				//cout << "Best Known Objective : " << model.get(GRB_DoubleAttr_ObjVal) << endl;

				total_obj += model.get(GRB_DoubleAttr_ObjVal);

				//cout << "Adding first set of variables:" << endl;
				// creating initial variables for column_generation

				Variable* p_var;
				p_var = new Variable();
				p_var->index = variables.size();
				p_var->assigned_group_size.resize(inst->n_subGroup);
				int column_repeat = 0;
				for (int g = 0; g < inst->n_subGroup; ++g) {
					p_var->assigned_group_size[g].resize(inst->max_party_size, 0);
					for (int nn = 0; nn < inst->max_party_size; ++nn) {
						for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
							if (x[seg][g][nn][i].get(GRB_DoubleAttr_X) > 0.9) {
								p_var->assigned_group_size[g][nn] = i;

								if (i > 0) {
									if (floor(inst->u[g][nn] / i) > column_repeat) {
										column_repeat = floor(inst->u[g][nn] / i);
									}
								}

							}
						}
					}
				}
				p_var->capacity = inst->no_seat_in_seg[seg];
				//computing value
				p_var->value = 0;
				for (int g = 0; g < inst->n_subGroup; ++g) {
					int n_party = 0;
					for (int nn = 0; nn < inst->max_party_size; ++nn) {
						if (x[seg][g][nn][0].get(GRB_DoubleAttr_X) < 0.1) {
							n_party++;
						}
					}
					if (n_party > 0.9) {
						p_var->value++;
					}
				}
				//adding the new variable	
				if (p_var->value > 0.5) {
					variables.push_back(p_var);
					for (int col = 0; col < column_repeat - 1; ++col) {
						variables.push_back(p_var);
					}
				}

				//Fixing the variables corresponding to segement "seg"
				for (int g = 0; g < inst->n_subGroup; ++g) {
					for (int nn = 0; nn < inst->max_party_size; ++nn) {
						for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
							if (x[seg][g][nn][i].get(GRB_DoubleAttr_X) > 0.9) {
								model.addConstr(x[seg][g][nn][i] == 1);
							}
							else {
								model.addConstr(x[seg][g][nn][i] == 0);
							}
						}
					}
				}

			}

		} // End of the loop for different segments
		  //if (model.get(GRB_IntAttr_Status) == 2 || model.get(GRB_IntAttr_Status) == 9) {
		  //	if (model.get(GRB_DoubleAttr_ObjVal) < global_upper_bound) {
		  //		global_upper_bound = model.get(GRB_DoubleAttr_ObjVal);
		  //		file1 << model.get(GRB_DoubleAttr_ObjVal) << ",";
		  //	}
		  //}

		delete x;
		//delete y;

	} //end of try

	catch (GRBException e) {
		cout << "createFirstSet2:Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		exit(1);
	}
	catch (...) {
		cout << "createFirstSet2:Other error ... " << endl;
		exit(1);
	}

}

void BPSolver::create_first_set_of_variables_IP(ofstream& file1) {

	GRBEnv env;

	try {

		GRBModel model(env);

		//Defining and adding variables
		GRBVar*** x;																		//x[s][g][n] : number of party of size n from group g is assigned to segement s
		x = new GRBVar**[inst->n_seg];
		for (int s = 0; s < inst->n_seg; s++) {
			x[s] = new GRBVar*[inst->n_subGroup];
			for (int g = 0; g < inst->n_subGroup; g++) {
				x[s][g] = new GRBVar[inst->max_party_size];
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					x[s][g][nn] = model.addVar(0.0, inst->u[g][nn], 0.0, GRB_INTEGER);
				}
			}
		}

		GRBVar** y;																		//y[g][seg] is equal to one if group "g" has at least one party in segment "seg", and zero otherwise
		y = new GRBVar*[inst->n_subGroup];
		for (int g = 0; g < inst->n_subGroup; ++g) {
			y[g] = new GRBVar[inst->n_seg];
			for (int seg = 0; seg < inst->n_seg; ++seg) {
				y[g][seg] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			}
		}

		model.update();
		//adding constraints

		//sum_seg x[seg][g][nn] == number of parties of size nn+1 from group g
		for (int g = 0; g < inst->n_subGroup; g++) {
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				GRBLinExpr equ = 0.0;
				for (int s = 0; s < inst->n_seg; ++s) {
					equ += x[s][g][nn];
				}
				model.addConstr(equ == inst->u[g][nn]);
			}
		}

		//x[p][seg]<= y[g][seg] for all seg,g,p \in P(g)
		for (int s = 0; s < inst->n_seg; ++s) {
			for (int g = 0; g < inst->n_subGroup; ++g) {
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					model.addConstr(x[s][g][nn] <= inst->u[g][nn] * y[g][s]);
				}
			}
		}

		//x[p][seg]<= y[g][seg] for all seg,g,p \in P(g)
		for (int s = 0; s < inst->n_seg; ++s) {
			for (int g = 0; g < inst->n_subGroup; ++g) {
				GRBLinExpr equ = 0.0;
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					equ += (nn + 1) * x[s][g][nn];
				}
				model.addConstr(equ <= inst->no_seat_in_seg[s] * y[g][s]);
			}
		}

		//Knapsack constraint
		for (int s = 0; s < inst->n_seg; ++s) {
			GRBLinExpr equ = 0.0;
			for (int g = 0; g < inst->n_subGroup; ++g) {
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					equ += (nn + 1) * x[s][g][nn];
				}
			}
			model.addConstr(equ <= inst->no_seat_in_seg[s]);
		}

		//new lower bound
		for (int g = 0; g < inst->n_subGroup; g++) {
			GRBLinExpr equ = 0.0;
			for (int s = 0; s < inst->n_seg; ++s) {
				equ += y[g][s];
			}
			double new_LB = 0;
			if (inst->LB_MT1[g] > new_LB) {
				new_LB = inst->LB_MT1[g];
			}
			if (inst->LB_MT2[g] > new_LB) {
				new_LB = inst->LB_MT2[g];
			}
			model.addConstr(equ >= new_LB);
		}

		///////////////////////////
		//Adding objective function 

		int total_obj = 0;

		GRBLinExpr objective = 0.0;
		for (int g = 0; g < inst->n_subGroup; ++g) {
			for (int seg = 0; seg < inst->n_seg; ++seg) {
				objective += y[g][seg];
			}
		}

		model.update();
		model.setObjective(objective, GRB_MINIMIZE);

		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_DoubleParam_TimeLimit, 30);
		//model.getEnv().set(GRB_DoubleParam_TimeLimit, 300);
		//model.getEnv().set(GRB_DoubleParam_MIPGap, 0.1);

		//model.getEnv().set(GRB_IntParam_MIPFocus, 1);

		//model.getEnv().set(GRB_IntParam_OutputFlag, 0);
		model.optimize();

		//cout << endl;
		//cout << "Solution Status: " << model.get(GRB_IntAttr_Status) << endl;
		//cout << "Time Taken: " << model.get(GRB_DoubleAttr_Runtime) << endl;

		if (model.get(GRB_IntAttr_Status) == 2 || model.get(GRB_IntAttr_Status) == 9) {
			// This gets you the relaxation bound
			//cout << "Relaxation Bound: " << model.get(GRB_DoubleAttr_ObjBound) << endl;
			// This gets you the best known primal solution value
			//cout << "Best Known Objective : " << model.get(GRB_DoubleAttr_ObjVal) << endl;

			if (model.get(GRB_DoubleAttr_ObjBound) > global_lower_bound) {
				global_lower_bound = model.get(GRB_DoubleAttr_ObjBound);
				//file1 << model.get(GRB_DoubleAttr_ObjBound) << ",";
			}
			if (model.get(GRB_DoubleAttr_ObjVal) < global_upper_bound) {
				global_upper_bound = model.get(GRB_DoubleAttr_ObjVal);
				//file1 << model.get(GRB_DoubleAttr_ObjVal) << ",";
			}

			//cout << "Adding first set of variables:" << endl;
			// creating initial variables for column_generation
			if (model.get(GRB_DoubleAttr_ObjVal) < 10000) {
				for (int s = 0; s < inst->n_seg; ++s) {
					Variable* p_var;
					p_var = new Variable();
					p_var->index = variables.size();
					p_var->assigned_group_size.resize(inst->n_subGroup);
					int column_repeat = 0;
					for (int g = 0; g < inst->n_subGroup; ++g) {
						p_var->assigned_group_size[g].resize(inst->max_party_size, 0);
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							if (x[s][g][nn].get(GRB_DoubleAttr_X) > 0.9) {

								if (x[s][g][nn].get(GRB_DoubleAttr_X) - floor(x[s][g][nn].get(GRB_DoubleAttr_X)) > 0.1) {
									p_var->assigned_group_size[g][nn] = floor(x[s][g][nn].get(GRB_DoubleAttr_X)) + 1;
								}
								else {
									p_var->assigned_group_size[g][nn] = floor(x[s][g][nn].get(GRB_DoubleAttr_X));
								}

								if (p_var->assigned_group_size[g][nn] > 0) {
									if (floor(inst->u[g][nn] / p_var->assigned_group_size[g][nn]) > column_repeat) {
										column_repeat = floor(inst->u[g][nn] / p_var->assigned_group_size[g][nn]);
									}
								}

							}
						}
					}
					p_var->capacity = inst->no_seat_in_seg[s];
					//computing value
					p_var->value = 0;
					for (int g = 0; g < inst->n_subGroup; ++g) {
						int n_party = 0;
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							if (x[s][g][nn].get(GRB_DoubleAttr_X) > 0.5) {
								n_party++;
							}
						}
						if (n_party > 0.9) {
							p_var->value++;
						}
					}
					//adding the new variable	
					if (p_var->value > 0.5) {
						variables.push_back(p_var);
						for (int col = 0; col < column_repeat - 1; ++col) {
							variables.push_back(p_var);
						}
					}
				}
			}
		}
		cout << "Hiii" << endl;
		//if (model.get(GRB_DoubleAttr_ObjVal) < global_upper_bound) {
		//	global_upper_bound = model.get(GRB_DoubleAttr_ObjVal);
		//	file1 << model.get(GRB_DoubleAttr_ObjVal) << ",";
		//}

		delete x;
		delete y;

	} //end of try

	catch (GRBException e) {
		cout << "createFirstSet2:Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		exit(1);
	}
	catch (...) {
		cout << "createFirstSet2:Other error ... " << endl;
		exit(1);
	}

}

void BPSolver::create_first_search_node() {

	SearchNode* first_search_node;
	first_search_node = new SearchNode();
	first_search_node->index = 0;
	first_search_node->LB_seg_no_for_group.resize(inst->n_subGroup);
	first_search_node->UB_seg_no_for_group.resize(inst->n_subGroup);
	for (int g = 0; g < inst->n_subGroup; ++g) {
		first_search_node->LB_seg_no_for_group[g] = 0; //min_seg[g];
	}
	for (int g = 0; g < inst->n_subGroup; ++g) {
		first_search_node->UB_seg_no_for_group[g] = inst->n_seg; //min(inst->n_seg, inst->n_party_in_group[g]);
	}
	search_nodes.push_back(first_search_node);

}

void BPSolver::print_search_node(SearchNode* node) {

	cout << "*************SEARCH NODE! ************" << endl;
	cout << "Index of search node: " << node->index << endl;
	cout << "\tparent lower bound: " << node->parent_lower_bound << endl;
	cout << "\tafter processing node lower bound: " << node->LP_lower_bound << endl;

	cout << "\tnumber of variables set to 1: " << node->variables_set_to_one.size() << endl;
	for (int v_index = 0; v_index < (int)node->variables_set_to_one.size(); ++v_index) {
		print_variable(node->variables_set_to_one[v_index]);
	}
	cout << "\tnumber of variables set to 0: " << node->variables_set_to_zero.size() << endl;
	for (int v_index = 0; v_index < (int)node->variables_set_to_zero.size(); ++v_index) {
		print_variable(node->variables_set_to_zero[v_index]);
	}
	cout << "**************************************" << endl;

}

void BPSolver::print_selected_node(SearchNode* node, ofstream& file1) {

	file1 << "*************SELECTED NODE! ************" << endl;
	file1 << "Index of search node: " << node->index << endl;
	file1 << "\tparent lower bound: " << node->parent_lower_bound << endl;
	file1 << "\tlower bound after CG: " << node->LP_lower_bound << endl;
	file1 << "\tnumber of variables set to one: " << node->variables_set_to_one.size() << endl;
	file1 << "\t\tindices of variables set to one:" << "\t";
	for (int v_index = 0; v_index < (int)node->variables_set_to_one.size(); ++v_index) {
		file1 << node->variables_set_to_one[v_index]->index << "\t";
	}
	file1 << endl;
	file1 << "\tnumber of variables set to zero: " << node->variables_set_to_zero.size() << endl;
	file1 << "\t\tindices of variables set to zero: " << "\t";
	for (int v_index = 0; v_index < (int)node->variables_set_to_zero.size(); ++v_index) {
		file1 << node->variables_set_to_zero[v_index]->index << "\t";
	}
	file1 << endl;
	file1 << "\t(LB,UB) for number of segemets assigned to each group:" << "\t";
	for (int g = 0; g < inst->n_subGroup; ++g) {
		file1 << "(" << node->LB_seg_no_for_group[g] << "," << node->UB_seg_no_for_group[g] << ")" << "\t";
	}
	file1 << endl;
	file1 << "**************************************" << endl;

}

void BPSolver::print_variable(Variable* var_to_print) {

	cout << "\t\tIndex of variable: " << var_to_print->index << endl;
	cout << "\t\t\tValue of variable: " << var_to_print->value << endl;
	cout << "\t\t\tCapacity of variable: " << var_to_print->capacity << endl;
	cout << "\t\t\t";
	cout << "assigned_group_size[g][nn]: " << endl;
	for (int g = 0; g < inst->n_subGroup; ++g) {
		for (int nn = 0; nn < inst->max_party_size; ++nn) {
			if (var_to_print->assigned_group_size[g][nn] > 0.9) {
				cout << "group: " << g << " size: " << nn + 1 << "\t";
			}
		}
	}
	cout << endl;

}

void BPSolver::print_best_solution() {
	//cout << "******* BEST SOLUTION ******** " << endl;
	for (int i = 0; i < (int)best_solution.size(); ++i) {
		print_variable(best_solution[i]);
	}
}

//// Codes related to using MDD ////

void BPSolver::solve_MDD_net_model(double & time_limit, ofstream& file1) {

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



			int n_MDDarcs = MDD->MDD_arcs[inst->n_SeatRow - 1].size();
			n_arcs = n_MDDarcs;
			n_nodes = MDD->MDD_nodes[inst->n_SeatRow - 1].size();



			GRBVar* x;																		//x[a] : number of times arc "a" in the MDD of segment "s" is selected
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
				for (int in = 0; in < MDD->MDD_nodes[inst->n_SeatRow - 1][k]->input_arcs.size(); ++in) {
					equ += x[MDD->MDD_nodes[inst->n_SeatRow - 1][k]->input_arcs[in]];
				}
				for (int out = 0; out < MDD->MDD_nodes[inst->n_SeatRow - 1][k]->output_arcs.size(); ++out) {
					equ -= x[MDD->MDD_nodes[inst->n_SeatRow - 1][k]->output_arcs[out]];
				}
				model.addConstr(equ == 0.0);
				n_const++;
			}

			//x_s,(0,1) + x_s,(0,2) + ... = 1
			GRBLinExpr equ2 = 0.0;
			for (int out = 0; out < MDD->MDD_nodes[inst->n_SeatRow - 1][0]->output_arcs.size(); ++out) {
				equ2 += x[MDD->MDD_nodes[inst->n_SeatRow - 1][0]->output_arcs[out]];
			}
			model.addConstr(equ2 == inst->n_seg);
			n_const++;

			for (int cap = 0; cap < inst->n_SeatRow; ++cap) {
				if (inst->no_seg_of_each_size[cap] > 0) {

					GRBLinExpr equ3 = 0.0;
					for (int node = 0; node < MDD->T_conn_nodes_index[inst->n_SeatRow - 1].size(); ++node) {
						if ((inst->n_SeatRow - MDD->MDD_nodes[inst->n_SeatRow - 1][MDD->T_conn_nodes_index[inst->n_SeatRow - 1][node]]->capacity) <= cap + 1) {
							for (int aa = 0; aa < MDD->MDD_nodes[inst->n_SeatRow - 1][MDD->T_conn_nodes_index[inst->n_SeatRow - 1][node]]->input_arcs.size(); ++aa) {
								equ3 += x[MDD->MDD_nodes[inst->n_SeatRow - 1][MDD->T_conn_nodes_index[inst->n_SeatRow - 1][node]]->input_arcs[aa]];
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
							if (MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[0] == nn + 1 && MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[1] == g) {
								party_equ += MDD->MDD_arcs[inst->n_SeatRow - 1][a]->n_party_selected * x[a];
							}
						}

						model.addConstr(party_equ == (int)inst->u[g][nn]);
						n_const++;
					}
				}
			}

			///////////////////////////
			//Adding objective function 
			GRBLinExpr objective = 0.0;
			for (int a = 0; a < n_MDDarcs; ++a) {
				objective += MDD->MDD_arcs[inst->n_SeatRow - 1][a]->value * x[a];
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

class callback_pricing_MDD : public GRBCallback {

public:
	// put objects required

	GRBVar* x_callback;
	vector <Variable*> variables;
	Stadium_Seating* inst;
	MDDSolver* MDD;
	int capacity;
	//callback_vrpspd(int* _n_solutions_found,QuadKnap* _quad_knap, vector<vector<double > > * _quad_obj_matrix, GRBVar * _t, GRBVar* _x){
	callback_pricing_MDD(GRBVar* _x, vector <Variable*> _variables, Stadium_Seating* _inst, MDDSolver* _MDD, int _capacity) {

		////cout << "We are at a callback .. " << endl;
		x_callback = _x;
		inst = _inst;
		MDD = _MDD;
		variables = _variables;
		capacity = _capacity;
	};

protected:

	void callback() {
		try {

			if (where == GRB_CB_MIPSOL) {


				////// ADD CUTS BELOW THIS
				vector<int> X(MDD->MDD_arcs[inst->n_SeatRow - 1].size());
				for (int a = 0; a < inst->n_subGroup; ++a) {
					if (getSolution(x_callback[a]) >= 0.5) {
						X[a] = 1.0;
					}
				}


				//Avoid constraint
				for (int c = 0; c < (int)variables.size(); ++c) {
					if (variables[c]->capacity == capacity) {
						int check = 0;
						int counter = 0;
						for (int a = 0; a < inst->n_subGroup; ++a) {
							if (variables[c]->assigned_group_size[MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[1]][MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[0] - 1] == MDD->MDD_arcs[inst->n_SeatRow - 1][a]->n_party_selected) {
								check += X[a];
								counter++;
							}
							else {
								check -= X[a];
							}
						}
						if (check > counter - 1) {
							GRBLinExpr avoid_cut = 0.0;
							int count = 0;
							for (int a = 0; a < inst->n_subGroup; ++a) {
								if (variables[c]->assigned_group_size[MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[1]][MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[0] - 1] == MDD->MDD_arcs[inst->n_SeatRow - 1][a]->n_party_selected) {
									avoid_cut += x_callback[a];
									count++;
								}
								else {
									avoid_cut -= x_callback[a];
								}
							}
							addLazy(avoid_cut, GRB_LESS_EQUAL, count - 1);
						}
					}
				}

				///// ADD CUTS ABOVE THIS
			}




		}
		catch (GRBException e) {
			cout << "Error number: " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
			exit(1);
		}
		catch (...) {
			cout << "Error during callback" << endl;
		}

	}
};

void BPSolver::solve_MDD_pricing_problem(SearchNode* node, double & time_limit, ofstream& file1, vector<vector<double>> & party_duals, vector<double> & capacity_duals, vector<double> & seg_LB_duals, vector<double> & seg_UB_duals, int capacity, bool & need_another_rmp_lp_solution, int& n_price_iter, double& total_price_time) {
	//cout << "=============================================" << endl << "Solving Pricing Problem using MDD and Net Flow formualtion:" << endl << "=============================================";

	vector<vector<int>> xx; //xx[g][n] number of items of size n from group g, so far selected 
	xx.resize(inst->n_subGroup);
	for (int g = 0; g < inst->n_subGroup; ++g) {
		xx[g].resize(inst->max_party_size, 0);
	}

	int total_obj = 0;


	char varname[256];

	GRBEnv env;
	try {

		GRBModel model(env);


		////Defining and adding variables
		//GRBVar*** x;																		//x[g][n][i] : Number of party of size n from group g assigned to this column is equal to i
		//x = new GRBVar**[inst->n_subGroup];
		//for (int g = 0; g < inst->n_subGroup; g++) {
		//	x[g] = new GRBVar*[inst->max_party_size];
		//	for (int nn = 0; nn < inst->max_party_size; ++nn) {
		//		x[g][nn] = new GRBVar[inst->u[g][nn] + 1];
		//		for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
		//			x[g][nn][i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		//		}
		//	}
		//}

		//GRBVar* y;																		//y[g] is equal to one if group "g" has some party in this column, and zero otherwise
		//y = new GRBVar[inst->n_subGroup];
		//for (int g = 0; g < inst->n_subGroup; ++g) {
		//	y[g] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		//}


		//Defining and adding arc variables
		int count = 0;
		int n_const = 0;
		int n_arcs = 0;
		int n_nodes = 0;



		int n_MDDarcs = MDD->MDD_arcs[inst->n_SeatRow - 1].size();
		n_arcs = n_MDDarcs;
		n_nodes = MDD->MDD_nodes[inst->n_SeatRow - 1].size();


		GRBVar* x;																		//x[a]=1 : if arc "a" in the MDD of segment "s" is selected
		x = new GRBVar[n_MDDarcs];
		for (int a = 0; a < n_MDDarcs; ++a) {
			x[a] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			count++;
		}

		model.update();


		//adding constraints


		//sum_[(i,j) \in A(s)] x_s,(i,j) - sum_[(j,i) \in A(s)] x_s,(j,i) = 0
		for (int k = 1; k < n_nodes; ++k) {
			GRBLinExpr equ = 0.0;
			for (int in = 0; in < MDD->MDD_nodes[inst->n_SeatRow - 1][k]->input_arcs.size(); ++in) {
				equ += x[MDD->MDD_nodes[inst->n_SeatRow - 1][k]->input_arcs[in]];
			}
			for (int out = 0; out < MDD->MDD_nodes[inst->n_SeatRow - 1][k]->output_arcs.size(); ++out) {
				equ -= x[MDD->MDD_nodes[inst->n_SeatRow - 1][k]->output_arcs[out]];
			}
			model.addConstr(equ == 0.0);
			n_const++;
		}


		//x_s,(0,1) + x_s,(0,2) + ... = 1
		GRBLinExpr equ2 = 0.0;
		for (int out = 0; out < MDD->MDD_nodes[inst->n_SeatRow - 1][0]->output_arcs.size(); ++out) {
			equ2 += x[MDD->MDD_nodes[inst->n_SeatRow - 1][0]->output_arcs[out]];
		}
		model.addConstr(equ2 == 1);
		n_const++;


		////total number of parties selected from each party size n from group g is equal to the total existing number of parties from g and of size n
		//for (int g = 0; g < inst->n_subGroup; ++g) {
		//	for (int nn = 0; nn < inst->max_party_size; ++nn) {
		//		if (inst->u[g][nn] > 0) {
		//			GRBLinExpr party_equ = 0.0;

		//			for (int a = 0; a < n_MDDarcs; ++a) {
		//				if (MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[0] == nn + 1 && MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[1] == g) {
		//					party_equ += MDD->MDD_arcs[inst->n_SeatRow - 1][a]->n_party_selected * x[a];
		//				}
		//			}

		//			model.addConstr(party_equ <= (int)inst->u[g][nn] - xx[g][nn]);
		//			n_const++;
		//		}
		//	}
		//}

		////sum_i x[g][nn][i] = 1
		//for (int g = 0; g < inst->n_subGroup; g++) {
		//	for (int nn = 0; nn < inst->max_party_size; ++nn) {
		//		GRBLinExpr equ = 0.0;
		//		for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
		//			equ += x[g][nn][i];
		//		}
		//		model.addConstr(equ == 1);
		//	}
		//}

		////i*x[g][nn][i] <= number of parties of size nn+1 from group g
		//for (int g = 0; g < inst->n_subGroup; g++) {
		//	GRBLinExpr equ11 = 0.0;
		//	for (int nn = 0; nn < inst->max_party_size; ++nn) {
		//		equ11 += (x[g][nn][0] - 1);
		//	}
		//	model.addConstr(inst->max_party_size * y[g] + equ11 >= 0);
		//}

		////x[g][nn]<= y[g] for all group g and size nn
		//for (int g = 0; g < inst->n_subGroup; ++g) {
		//	GRBLinExpr Yequ = 0.0;
		//	for (int nn = 0; nn < inst->max_party_size; ++nn) {
		//		model.addConstr(1 - x[g][nn][0] <= y[g]);
		//		Yequ += (1 - x[g][nn][0]);
		//	}
		//	model.addConstr(y[g] <= Yequ);
		//}

		////Knapsack constraint
		//GRBLinExpr cap_equ = 0.0;
		//for (int g = 0; g < inst->n_subGroup; ++g) {
		//	for (int nn = 0; nn < inst->max_party_size; ++nn) {
		//		for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
		//			cap_equ += (nn + 1) * i * x[g][nn][i];
		//		}
		//	}
		//}
		//model.addConstr(cap_equ <= capacity);

		///////////////////////////
		//Adding objective function 
		GRBLinExpr objective = 0.0;
		for (int a = 0; a < n_MDDarcs; ++a) {
			//cout << "arc " << a << " group " << MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[1] << " size " << MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[0] << endl;
			if (MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[1] != -1) {
				objective += MDD->MDD_arcs[inst->n_SeatRow - 1][a]->value * x[a];

				objective -= party_duals[MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[1]][MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[0] - 1] * MDD->MDD_arcs[inst->n_SeatRow - 1][a]->n_party_selected * x[a];

				objective -= MDD->MDD_arcs[inst->n_SeatRow - 1][a]->value * (seg_LB_duals[MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[1]] + seg_UB_duals[MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[1]]) * x[a];
			}
		}
		objective -= capacity_duals[capacity - 1];


		//for (int g = 0; g < inst->n_subGroup; ++g) {
		//	for (int nn = 0; nn < inst->max_party_size; ++nn) {
		//		for (int i = 1; i < inst->u[g][nn] + 1; ++i) {
		//			objective -= party_duals[g][nn] * i * x[g][nn][i];
		//		}
		//	}
		//}
		//objective -= capacity_duals[capacity - 1];
		//for (int g = 0; g < inst->n_subGroup; ++g) {
		//	objective -= (seg_LB_duals[g] + seg_UB_duals[g]) * y[g];
		//}

		// Using Callback to add the avoid constraints (i.e. to avoid from regenerating the existing columns)
		//callback_pricing_MDD cb = callback_pricing_MDD(x, variables, inst, MDD, capacity);
		//model.setCallback(&cb);

		//calculating the total number of items selected
		GRBLinExpr Total_item = 0.0;
		//for (int g = 0; g < inst->n_subGroup; ++g) {
		//	for (int nn = 0; nn < inst->max_party_size; ++nn) {
		//		for (int i = 1; i < inst->u[g][nn] + 1; ++i) {
		//			Total_item += i * x[g][nn][i];
		//		}
		//	}
		//}
		for (int a = 0; a < n_MDDarcs; ++a) {
			Total_item += MDD->MDD_arcs[inst->n_SeatRow - 1][a]->n_party_selected * x[a];
		}

		double last_RC = 0.0;


		for (int iter = 0; iter < inst->n_seg; ++iter) { //iterations for solving the pricing problem

														 //adding objective[seg]
														 //IloObjective obj;
			model.update();
			model.setObjective(objective, GRB_MINIMIZE);
			model.getEnv().set(GRB_IntParam_Threads, 1);
			model.getEnv().set(GRB_DoubleParam_TimeLimit, time_limit);
			model.getEnv().set(GRB_IntParam_OutputFlag, 0);
			model.getEnv().set(GRB_IntParam_PreCrush, 1);
			//model.getEnv().set(GRB_DoubleParam_Cutoff, PRECISION);
			model.getEnv().set(GRB_IntParam_DualReductions, 0);
			model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
			model.update();
			model.optimize();

			n_price_iter += 1;
			total_price_time += model.get(GRB_DoubleAttr_Runtime);

			//cout << "Solution status: " << model.get(GRB_IntAttr_Status) << endl;
			//cout << ". Best objective: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
			//cout << "Pricing Time" << model.get(GRB_DoubleAttr_Runtime) << endl;

			//cout << "Time Taken for pricing iteration " << iter << ": " << model.get(GRB_DoubleAttr_Runtime) << endl;

			if (model.get(GRB_IntAttr_Status) == 2) {

				if (iter > 0 || model.get(GRB_DoubleAttr_ObjVal) < -PRECISION) {

					//cout << "Reduced Cost: " << model.get(GRB_DoubleAttr_ObjVal) << "\t";

					int ocup_cap = 0;

					Variable* p_var;
					p_var = new Variable();
					p_var->index = variables.size();
					p_var->assigned_group_size.resize(inst->n_subGroup);
					int column_repreat = 0; //number of times we want to repeat adding the current column because it may need to selected multiple times
					for (int g = 0; g < inst->n_subGroup; ++g) {
						p_var->assigned_group_size[g].resize(inst->max_party_size, 0);
					}
					for (int a = 0; a < n_MDDarcs; ++a) {
						if (x[a].get(GRB_DoubleAttr_X) > 0.9 && MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[1] != -1) {
							p_var->assigned_group_size[MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[1]][MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[0] - 1] = int(MDD->MDD_arcs[inst->n_SeatRow - 1][a]->n_party_selected);
							//ocup_cap += MDD->MDD_arcs[inst->n_SeatRow - 1][a]->n_party_selected * MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[0];
							if (MDD->MDD_arcs[inst->n_SeatRow - 1][a]->n_party_selected > 0) {
								if (floor(inst->u[MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[1]][MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[0] - 1] / MDD->MDD_arcs[inst->n_SeatRow - 1][a]->n_party_selected) > column_repreat) {
									column_repreat = floor(inst->u[MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[1]][MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[0] - 1] / MDD->MDD_arcs[inst->n_SeatRow - 1][a]->n_party_selected);
								}
							}
						}
					}
					//for (int g = 0; g < inst->n_subGroup; ++g) {
					//	p_var->assigned_group_size[g].resize(inst->max_party_size, 0);
					//	for (int nn = 0; nn < inst->max_party_size; ++nn) {
					//		if (inst->u[g][nn] > 0.5) {
					//			for (int i = 0; i < inst->u[g][nn] + 1; ++i) {
					//				if (x[g][nn][i].get(GRB_DoubleAttr_X) > 0.9) {
					//					p_var->assigned_group_size[g][nn] = int(i);
					//					ocup_cap += i *nn;

					//					if (i > 0) {
					//						if (floor(inst->u[g][nn] / i) > column_repreat) {
					//							column_repreat = floor(inst->u[g][nn] / i);
					//						}
					//					}
					//				}
					//			}
					//		}
					//	}
					//}
					p_var->capacity = capacity;
					//computing value
					p_var->value = 0;
					for (int a = 0; a < n_MDDarcs; ++a) {
						if (x[a].get(GRB_DoubleAttr_X) > 0.9) {
							p_var->value += MDD->MDD_arcs[inst->n_SeatRow - 1][a]->value;
						}
					}
					//for (int g = 0; g < inst->n_subGroup; ++g) {
					//	int n_party = 0;
					//	for (int nn = 0; nn < inst->max_party_size; ++nn) {
					//		if (x[g][nn][0].get(GRB_DoubleAttr_X) < 0.5) {
					//			n_party++;
					//		}
					//	}
					//	if (n_party > 0.9) {
					//		p_var->value++;
					//	}
					//}
					//int value1 = p_var->value;

					//adding the new variable	
					if (p_var->value > 0.5) {
						variables.push_back(p_var);
						for (int col = 0; col < column_repreat - 1; ++col) {
							variables.push_back(p_var);
						}
					}

					need_another_rmp_lp_solution = true;

					//cout << "occupied capacity: " << ocup_cap << endl;

					//updating the total objective
					total_obj += p_var->value;

					//updating xx[g][n]
					for (int g = 0; g < inst->n_subGroup; ++g) {
						for (int nn = 0; nn < inst->max_party_size; ++nn) {
							xx[g][nn] += p_var->assigned_group_size[g][nn];
							//adding new constraints to the model
							for (int i = inst->u[g][nn] - xx[g][nn] + 1; i < inst->u[g][nn] + 1; ++i) {
								for (int a = 0; a < n_MDDarcs; ++a) {
									if (MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[1] == g && MDD->MDD_arcs[inst->n_SeatRow - 1][a]->size_group[0] == nn + 1 && MDD->MDD_arcs[inst->n_SeatRow - 1][a]->n_party_selected == int(i)) {
										model.addConstr(x[a] == 0);
									}
								}
							}
							//for (int i = inst->u[g][nn] - xx[g][nn] + 1; i < inst->u[g][nn] + 1; ++i) {
							//	model.addConstr(x[g][nn][i] == 0);
							//}
						}
					}

					if (model.get(GRB_DoubleAttr_ObjVal) > -PRECISION) {
						if (iter > 0 && last_RC > -PRECISION) {
							objective += Total_item * abs(last_RC) / inst->Tn_party;
						}
						objective -= Total_item * abs(model.get(GRB_DoubleAttr_ObjVal)) / inst->Tn_party;
						last_RC = model.get(GRB_DoubleAttr_ObjVal);
					}

				}

				else {
					iter = inst->n_seg;
				}

			}


		}


		bool feasible = true;
		for (int g = 0; g < inst->n_subGroup; ++g) {
			for (int nn = 0; nn < inst->max_party_size; ++nn) {
				if (inst->u[g][nn] != xx[g][nn]) {
					feasible = false;
				}
			}
		}

		if (feasible) {
			//cout << "FEASIBLE!!!!!!!!!!!!!" << endl;
			//cout << "Objective value: " << total_obj << endl;

			if (total_obj < global_upper_bound) {
				global_upper_bound = total_obj;
			}
		}


		delete x;

	} // end of try
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		exit(1);
	}
	catch (...) {
		cout << "Other error ... " << endl;
		exit(1);
	}
}
