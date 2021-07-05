#ifndef IP_model_HPP_
#define IP_model_HPP_

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
#include "Instance_Stadium.h"
#include <vector>
#include <algorithm>
#include <iomanip>
#include <set>
#include "gurobi_c++.h"

using namespace std;


struct IP_Optimal_Model{

	Stadium_Seating* inst;

	void solve_SS_IP(Stadium_Seating* inst, double time_limit, ofstream& file1); // Solving a two-stage IP for stadium seating

	IP_Optimal_Model(Stadium_Seating* _instance) {
		inst = _instance;
	}

};

#endif