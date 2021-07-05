#include "IP_model_new.h"

void IP_Optimal_Model_new::solve_SS_IP(Stadium_Seating* inst, double time_limit, ofstream& file1) {
	cout << "=============================================" << endl << "Solving stadium seating problem:" << endl << "=============================================";

	GRBEnv env;

	////// Creating LowerBound 1 by Silvano Martello and Paolo Toth ///////////
	vector<double> LB_MT1(inst->n_subGroup, 0);
	for (int g = 0; g < inst->n_subGroup; ++g) {
		for (int p = 0; p < inst->Tn_party; ++p) {
			if (inst->subG_party[p] == g) {
				LB_MT1[g] += inst->n[p];
			}
		}
		LB_MT1[g] = ceil(LB_MT1[g] / inst->n_SeatRow);
	}

	////// Creating LowerBound 2 by Silvano Martello and Paolo Toth ///////////
	vector<double> LB_MT2;
	LB_MT2.resize(inst->n_subGroup, 0);
	double K_MT = floor(inst->n_SeatRow / 2 + 0.001);
	for (int g = 0; g < inst->n_subGroup; ++g) {
		vector<bool> N_1;
		N_1.resize(inst->Tn_party, 0);
		vector<bool> N_2;
		N_2.resize(inst->Tn_party, 0);
		vector<bool> N_3;
		N_3.resize(inst->Tn_party, 0);
		for (int p = 0; p < inst->Tn_party; ++p) {
			if (inst->subG_party[p] == g) {
				if (inst->n[p] >(inst->n_SeatRow - K_MT)) {
					N_1[p] = 1;
				}
				if (inst->n[p] <= (inst->n_SeatRow - K_MT) && inst->n[p] >(inst->n_SeatRow / 2)) {
					N_2[p] = 1;
				}
				if (inst->n[p] >= K_MT && inst->n[p] <= (inst->n_SeatRow / 2)) {
					N_3[p] = 1;
				}
			}
		}
		double sum_term = 0;
		for (int p = 0; p < inst->Tn_party; ++p) {
			if (inst->subG_party[p] == g) {
				LB_MT2[g] += (N_1[p] + N_2[p]);
				if (N_3[p] == 1) {
					sum_term += inst->n[p];
				}
				if (N_2[p] == 1) {
					sum_term += (inst->n[p] - inst->n_SeatRow);
				}
			}
		}
		if (ceil(sum_term / inst->n_SeatRow - 0.001) > 0) {
			LB_MT2[g] += ceil(sum_term / inst->n_SeatRow - 0.001);
		}
	}

	try {

		GRBModel model(env);

		//Defining and adding variables
		GRBVar** x;																		//x[p][seg] is equal to one if party "p" is assigned segment "seg", and zero otherwise
		x = new GRBVar*[inst->Tn_party];
		for (int p = 0; p < inst->Tn_party; ++p) {
			x[p] = new GRBVar[inst->n_seg];
			for (int seg = 0; seg < inst->n_seg; ++seg) {
				x[p][seg] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
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

		//sum_r x[p][seg] == 1
		for (int p = 0; p < inst->Tn_party; ++p) {
			GRBLinExpr equ = 0.0;
			for (int seg = 0; seg < inst->n_seg; ++seg) {
				equ += x[p][seg];
			}
			model.addConstr(equ == 1);
		}

		//x[p][seg]<= y[g][seg] for all seg,g,p \in P(g)
		for (int seg = 0; seg < inst->n_seg; ++seg) {
			for (int g = 0; g < inst->n_subGroup; ++g) {
				for (int p = 0; p < inst->Tn_party; ++p) {
					if (inst->subG_party[p] == g) {
						model.addConstr(x[p][seg] <= y[g][seg]);
					}
				}
			}
		}

		for (int seg = 0; seg < inst->n_seg; ++seg) {
			for (int g = 0; g < inst->n_subGroup; ++g) {
				GRBLinExpr equ = 0.0;
				for (int p = 0; p < inst->Tn_party; ++p) {
					if (inst->subG_party[p] == g) {
						equ += inst->n[p] * x[p][seg];
					}
				}
				model.addConstr(equ <= inst->no_seat_in_seg[seg] * y[g][seg]);
			}
		}

		//Knapsack constraint
		vector<GRBConstr > cap_con(inst->n_seg);
		for (int seg = 0; seg < inst->n_seg; ++seg) {
			GRBLinExpr equ = 0.0;
			for (int p = 0; p < inst->Tn_party; ++p) {
				equ += inst->n[p] * x[p][seg];
			}
			cap_con[seg] = model.addConstr(equ <= inst->no_seat_in_seg[seg]);
		}

		for (int g = 0; g < inst->n_subGroup; g++) {
			GRBLinExpr equ = 0.0;
			for (int s = 0; s < inst->n_seg; ++s) {
				equ += y[g][s];
			}
			double new_LB = 0;
			if (LB_MT1[g] > new_LB) {
				new_LB = LB_MT1[g];
			}
			if (LB_MT2[g] > new_LB) {
				new_LB = LB_MT2[g];
			}
			model.addConstr(equ >= new_LB);
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
		//model.getEnv().set(GRB_IntParam_MIPFocus, 1);
		
		//model.getEnv().set(GRB_IntParam_PreCrush, 1);
		//model.getEnv().set(GRB_IntParam_DualReductions, 0);
		//model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
		//model.getEnv().set(GRB_IntParam_Presolve, -1);

		/////////////////////////////////////////////////////////
		//model.getEnv().set(GRB_IntParam_DualReductions, 1);
		//model.getEnv().set(GRB_IntParam_InfUnbdInfo, 1);
		/////////////////////////////////////////////////////////
		
		model.update();
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


			/////////////////////////////////////////////////////////
			//for (int seg = 0; seg < inst->n_seg; ++seg) {
				//cout << "FarkasDual for seg " << seg << " :" << cap_con[seg].get(GRB_DoubleAttr_FarkasDual) << endl;
				//cout << "FarkasProof for seg " << seg << " :" << cap_con[seg].get(GRB_DoubleAttr_FarkasProof) << endl;
			//}
			/////////////////////////////////////////////////////////

			// This gets you the relaxation bound
			cout << "Relaxation Bound: " << model.get(GRB_DoubleAttr_ObjBound) << endl;
			// This gets you the best known primal solution value
			cout << "Best Known Objective : " << model.get(GRB_DoubleAttr_ObjVal) << endl;

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
		}
		else {
			file1 << endl;
		}


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

}
