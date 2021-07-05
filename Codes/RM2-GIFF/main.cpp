#include "Instance_Stadium.h"
#include "IP_model.h"
#include "MDD.h"

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
#include <ilcplex/ilocplex.h>
#include <set>
#include "gurobi_c++.h"

using namespace std;

const char* filenameOut = "instances_ourModel.csv"; // creating the excel file to print the outputs
ofstream file1(filenameOut);

int main() {
	
	// algo = 1 is for solving the IP model, algo = 2 is for solving Net Flow model on MDD
	int algo=2;
	cout << "Select (for solve_IP_model : 1, for solve_Net_Flow_Model_on_MDD : 2) " << endl;


	////// What we are printing as outputs in ofstream file: ///////
	/*file1 << "Instance" << ",";
	file1 << "LB" << ",";
	file1 << "UB" << ",";
	file1 << "Solution Time" << ",";
	file1 << "Total Time" << endl;*/

	//// Headers for printing instance information:
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


	for (int instance = 1 ; instance < 11; ++instance) {
			

			string file;
			stringstream ss;
			ss << "60-300-3-" << instance << ".txt"; //<< "300-3-"
			file = ss.str();
			const char*     filename = file.c_str();
			double             time_limit = 1800;
			Stadium_Seating* inst = new Stadium_Seating(filename);

			cout << "Instance " << instance << endl;
			file1 << instance << ",";

			/*inst->print();*/
			clock_t startTime = clock();

			// Solving IP
			if (algo == 1) {
				IP_Optimal_Model* IP_model;
				IP_model = new IP_Optimal_Model(inst);
				IP_model->solve_SS_IP(inst, time_limit, file1);
			}

			// Solving Net Flow Model
			if (algo == 2) {
				MDDSolver* MDD_solver;
				MDD_solver = new MDDSolver(inst);
				MDD_solver->solve_net_model(time_limit, file1);
			}

			/////Printing the Total Time
			file1 << double(clock() - startTime) / (double)CLOCKS_PER_SEC << endl;
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




