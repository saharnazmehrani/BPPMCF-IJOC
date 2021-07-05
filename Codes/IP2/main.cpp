#include "Instance_Stadium.h"
#include "IP_model.h"
#include "IP_model_new.h"

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
	
	int algo = 2;
	//cout << "Select (for solve_IP_model : 1, for branch_and_price : 2) " << endl;
	//cin >> algo;

	/*file1 << "Instance" << ",";
	file1 << "# nodes explored" << ",";
	file1 << "# columns" << ",";
	file1 << "global LB" << ",";
	file1 << "global UB" << ",";
	file1 << "Time" << endl;*/

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


	for (int instance = 1 ; instance < 11 ; ++instance) {

			clock_t startTime = clock();

			string file;
			stringstream ss;
			ss << "50-10-" << instance << ".txt";
			file = ss.str();
			const char*     filename = file.c_str();
			int             time_limit = 1800;
			Stadium_Seating* inst = new Stadium_Seating(filename);

			cout << "Instance " << instance << endl;
			file1 << instance << ".txt" << ",";

			/*inst->print()*/;

			// running the old IP2, from first IJOC submission
			if (algo == 1) {
				IP_Optimal_Model* IP_model;
				IP_model = new IP_Optimal_Model(inst);
				IP_model->solve_SS_IP(inst, time_limit, file1);
			}

			// running the new IP2, afterning adding the constraint suggested by IJOC reviewer 1
			if (algo == 2) {
				IP_Optimal_Model_new* IP_model_new;
				IP_model_new = new IP_Optimal_Model_new(inst);
				IP_model_new->solve_SS_IP(inst, time_limit, file1);
			}

			//file1 << double(clock() - startTime) / (double)CLOCKS_PER_SEC << endl;

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




