#pragma once
#ifndef Stadium_Seating_
#define Stadium_Seating_

#include <vector>
#include <string>
#include <iostream>
#include "math.h"
#include <cstdlib>
#include <stdlib.h>
#include <fstream>
#include <algorithm>

using namespace std;

struct Stadium_Seating {


	// Internal elements
	const char*                 filename;                               // This is the name of the file in your system

	int n_sec;															 // Number of sections in the stadium
	int n_row;															 // Number of rows in each section
	int n_SeatRow;														 // Number of seats in each row
	int n_seat;			   									   			 // Total number seats = n_sec*n_row*n_SeatRow
	int n_seg;															 // Number of segments in the stadium
	vector<int> OcSeat;													 // binary parameter which is equal to one if the seat is occupied

	int n_subGroup;														 // Number of sub-groups;
	int Tn_party;														 // Total number of parties
	vector<int> subG_party;									             // sub-group number corresponding to each party
	vector<int> n;									         			 // Number of people in each party p
	int Tn_people;														 // \sum_p n  : total number of poeple					

	vector<int> a;											             // section number for each seat
	vector<int> r;														 // row number for each seat
	vector<int> sr;														 // seat number for each seat in the row
	vector<int> seg;													 // segment number for each seat

	vector<int> no_seat_in_seg;											 // number of seats in each segment
	vector<int> first_seat_seg;											 // index of the first seat in each segment
	vector<int> no_seg_of_each_size;									 // no_seg_of_each_size[sz]: number of segments available of size "sz";
																		 //Maximum possible size is the number of the seats in each row 
	vector<int> m;													     // number of available seats on the right side of each seat including the seat itself
	vector<int> o;														 // index of the seat in front of each seat

	vector<int> f;														 // binary parameter which is equal to one if seat s is the first seat in its corresponding segment
	vector<int> l;														 // binary parameter which is equal to one if seat s is the last seat in its corresponding segment

																		 // new parameters
	int max_party_size;													// maximum possible size for parties
	vector<vector<int>> u;												// number of parties of size n for each group
	vector<int> nn;														// number of people in each group
	vector<int> n_party_in_group;									    // number of parties in each group


	int n_diff_item_size;
	// Functions, implemented in .cpp file

	vector<double> LB_MT1;												// LowerBound 1 by Silvano Martello and Paolo Toth per group
	vector<double> LB_MT2;												// LowerBound 2 by Silvano Martello and Paolo Toth per group

	void read(const char* filename);                                     // Reads the instance
	void print();                                                        // Prints the instance


																		 // Constructor - this is run every time a new instnaces is created
	Stadium_Seating(const char* _filename) {

		filename = _filename;
		read(filename);

	};


};

#endif