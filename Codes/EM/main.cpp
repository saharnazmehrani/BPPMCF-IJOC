#include "Instance_Stadium.h"
#include "IP_model.h"
#include "BP_SS.h"
#include "MDD2.h"

#include <iostream>
#include <fstream>
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

const char* filenameOut = "instances_ourModel.csv"; // creating the excel file to print the outputs
ofstream file1(filenameOut);

int main() {

	//for solve_IP_model : algo=1, for branch_and_price : algo=2
	int algo = 3;
	cout << "Select (for solve_IP_model : 1, for branch_and_price : 2) " << endl;
	//cin >> algo;

	//// what we are printing in the of stream file for the Branch-and-Price:
	/*file1 << "Instance" << ",";
	file1 << "Time" << ",";
	file1 << "global LB" << ",";
	file1 << "global UB" << ",";
	file1 << "# nodes explored" << ",";
	file1 << "# columns" << endl;*/

	//// what we are printing in the of stream file for the IP:
	/*file1 << "Instance" << ",";
	file1 << "Time" << ",";
	file1 << "# nodes explored" << ",";
	file1 << "LB" << ",";
	file1 << "UB" << endl;*/

	//// Headers for printing instances information:
	/*file1 << "Instance" << ",";
	file1 << "#secs" << ",";
	file1 << "#rows" << ",";
	file1 << "#Seat-row" << ",";
	file1 << "#Segment" << ",";
	file1 << "#empty seats" << ",";
	file1 << "#seat per seg" << ",";
	file1 << "#Groups" << ",";
	file1 << "#party" << ",";
	file1 << "# people" << ",";
	file1 << "Avg # people in party" << ",";
	file1 << "Min seg size" << ",";
	file1 << "Max seg size" << endl;*/

	
	for (int instance =1; instance < 4; ++instance) {


		string file;
		stringstream ss;
		ss << "10-10-" << instance << ".txt";
		file = ss.str();
		const char*     filename = file.c_str();
		double             time_limit = 1800;
		Stadium_Seating* inst = new Stadium_Seating(filename);

		cout << "Instance " << instance << endl;
		file1 << instance << ".txt" << ",";

		/*inst->print();*/

		////creating MDD
		MDDSolver* MDD_solver;
		MDD_solver = new MDDSolver(inst);
		//MDD_solver->create_MDDs(time_limit, file1);

		clock_t startTime = clock();

		if (algo == 1) {
			IP_Optimal_Model* IP_model;
			IP_model = new IP_Optimal_Model(inst);
			IP_model->solve_SS_IP(inst, time_limit, file1);
		}


		if (algo == 2) {
			//MDDSolver* MDD_solver;
			//MDD_solver = new MDDSolver(inst);
			//MDD_solver->create_MDDs(time_limit, file1);
			//clock_t endTime = clock();
			//time_limit -= (double(endTime - startTime) / (double)CLOCKS_PER_SEC);
			MDD_solver->solve_net_model(time_limit, file1);
		}

		if (algo == 3) {
			BPSolver* bp_solver;
			bp_solver = new BPSolver(inst, MDD_solver);
			bp_solver->solve(time_limit, startTime, file1);
			//bp_solver->solve_MDD_net_model(time_limit, file1);  // we solve the netflow model using MDD by this line too
		}


		//int n_diff_party_size = 0;
		//for (int nn = 0; nn < inst->max_party_size; ++nn) {
		//	int total_party = 0;
		//	for (int g = 0; g < inst->n_subGroup; ++g) {
		//		if (inst->u[g][nn] > 0) {
		//			total_party += inst->u[g][nn];
		//		}
		//	}
		//	if (total_party > 0) {
		//		n_diff_party_size++;
		//	}
		//}
		//cout << "n_diff_party_size: " << n_diff_party_size << endl;

		//int n_item_per_bin = 0;
		//int total_weight = 0;
		//int min_item_size = inst->max_party_size;
		//int max_item_size = 0;
		//for (int p = 0; p < inst->Tn_party; ++p) {
		//	total_weight += inst->n[p];
		//	if (inst->n[p] > max_item_size) {
		//		max_item_size = inst->n[p];
		//	}
		//	if (inst->n[p] < min_item_size) {
		//		min_item_size = inst->n[p];
		//	}
		//}
		//n_item_per_bin = inst->n_SeatRow * inst->Tn_party / total_weight;
		////cout << "avg. item per bin: " << n_item_per_bin << endl;
		//cout << "max size: " << max_item_size << endl;
		//cout << "min size: " << min_item_size << endl;

		//file1 << double(clock() - startTime) / (double)CLOCKS_PER_SEC << endl;
		//cout << "Total Time: " << double(clock() - startTime) / (double)CLOCKS_PER_SEC << endl;

		//// Printing instance information:
		//file1 << inst->n_sec << ",";
		//file1 << inst->n_row << ",";
		//file1 << inst->n_SeatRow << ",";
		//file1 << inst->n_seg <<",";
		//int T_ocp_seat = 0;
		//for (int s = 0; s < inst->n_seat; ++s) {
		//	T_ocp_seat += inst->OcSeat[s];
		//}
		//file1 << inst->n_seat - T_ocp_seat << ",";  // total # of empty seats in the stadium
		//file1 << (inst->n_seat - T_ocp_seat) / inst->n_seg << ",";
		//file1 << inst->n_subGroup << ",";
		//file1 << inst->Tn_party << ",";
		//file1 << inst->Tn_people << ",";
		//file1 << inst->Tn_people / inst->Tn_party << ",";
		//int max_seg_size = 0;
		//int min_seg_size = inst->n_seat;
		//for (int seg = 0; seg < inst->n_seg; ++seg) {
		//	if (inst->no_seat_in_seg[seg] > max_seg_size) {
		//		max_seg_size = inst->no_seat_in_seg[seg];
		//	}
		//	if (inst->no_seat_in_seg[seg] < min_seg_size) {
		//		min_seg_size = inst->no_seat_in_seg[seg];
		//	}
		//}
		//file1 << min_seg_size << ",";
		//file1 << max_seg_size << endl;
	}

	return 0;
}




