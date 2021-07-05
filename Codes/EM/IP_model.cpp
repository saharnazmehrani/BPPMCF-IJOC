#include "IP_model.h"

void IP_Optimal_Model::solve_SS_IP(Stadium_Seating* inst, double time_limit, ofstream& file1) {
	cout << "=============================================" << endl << "Solving stadium seating problem:" << endl << "=============================================";

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
					x[s][g][nn] = model.addVar(0.0, INFINITY, 0.0, GRB_INTEGER);
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

		//Knapsack constraint
		for (int s = 0; s < inst->n_seg; ++s) {
			GRBLinExpr equ = 0.0;
			for (int g = 0; g < inst->n_subGroup; ++g) {
				for (int nn = 0; nn < inst->max_party_size; ++nn) {
					equ += (nn+1) * x[s][g][nn];
				}
			}
			model.addConstr(equ <= inst->no_seat_in_seg[s]);
		}


		///////////////////////////
		//Adding objective function 
		GRBQuadExpr objective = 0.0;
		for (int g = 0; g < inst->n_subGroup; ++g) {
			for (int seg = 0; seg < inst->n_seg; ++seg) {
				objective += y[g][seg];
			}
		}

		model.update();
		model.setObjective(objective, GRB_MINIMIZE);

		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_DoubleParam_TimeLimit, time_limit);
		//model.getEnv().set(GRB_IntParam_OutputFlag, 0);
		model.getEnv().set(GRB_IntParam_PreCrush, 1);
		model.getEnv().set(GRB_IntParam_DualReductions, 0);
		model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
		model.getEnv().set(GRB_IntParam_Presolve, -1);
		//model.getEnv().set(GRB_IntParam_MIPFocus, 2);
		model.update();
		model.optimize();
		cout << endl;
		cout << "Solution Status: " << model.get(GRB_IntAttr_Status) << endl;
		cout << "Time Taken: " << model.get(GRB_DoubleAttr_Runtime) << endl;
		file1 << model.get(GRB_DoubleAttr_Runtime) << ",";
		file1 << model.get(GRB_DoubleAttr_NodeCount) << ",";
		file1 << model.get(GRB_DoubleAttr_ObjBound) << ",";
		file1 << model.get(GRB_DoubleAttr_ObjVal) << endl;

		// This gets you the relaxation bound
		cout << "Relaxation Bound: " << model.get(GRB_DoubleAttr_ObjBound) << endl;
		// This gets you the best known primal solution value
		cout << "Best Known Objective : " << model.get(GRB_DoubleAttr_ObjVal) << endl;


		//printing the solution
		//for (int seg = 0; seg < inst->n_seg; ++seg) {
		//	cout << "Bin: " << seg << endl << "\t";
		//	for (int g = 0; g < inst->n_subGroup; ++g) {
		//		if (y[g][seg].get(GRB_DoubleAttr_X) > 0.9) {
		//			cout << "group: " << g << "\t";
		//		}
		//	}
		//	cout << endl;
		//}
		//cout << "*************************************************" << endl;
		//for (int s = 0; s < inst->n_seg; ++s) {
		//	cout << "Bin: " << s << endl << "\t";
		//	for (int g = 0; g < inst->n_subGroup; ++g) {
		//		for (int nn = 0; nn < inst->max_party_size; ++nn) {
		//			if (x[s][g][nn].get(GRB_DoubleAttr_X) > 0.9) {
		//				cout << "group: " << g << "\t";
		//			}
		//		}
		//	}
		//	cout << endl;
		//}


		delete x;
		delete y;

	} //end of try

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