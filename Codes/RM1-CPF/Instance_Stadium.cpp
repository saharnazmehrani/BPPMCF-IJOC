#include <ctime>
#include <algorithm>
#include "Instance_Stadium.h"


void Stadium_Seating::read(const char* filename) {

	cout << "Reading an instance .. " << endl;

	int reader_int;
	double reader_double;

	ifstream input;

	input.open(filename);

	// # of sections
	input >> reader_int;
	n_sec = reader_int;

	// # of rows in wach sections
	input >> reader_int;
	n_row = reader_int;

	// # of seat in each row
	input >> reader_int;
	n_SeatRow = reader_int;

	// total # of seats in the stadium
	n_seat = n_sec * n_row * n_SeatRow;

	// available seats in the stadium
	OcSeat.resize(n_seat);
	for (int s = 0; s < n_seat; ++s) {
		input >> reader_int;
		OcSeat[s] = reader_int;
	}

	// number of sub-groups
	input >> reader_int;
	n_subGroup = reader_int;

	//total number of parties
	input >> reader_int;
	Tn_party = reader_int;

	/// subG_party[p]: sub_group number for party p /// n[p]: number of people in party p
	vector<int> subG_party1(Tn_party);
	vector<int> n1(Tn_party);
	Tn_people = 0;
	for (int p = 0; p < Tn_party; ++p) {
		for (int i = 0; i < 2; ++i) {
			if (i == 0) {
				input >> reader_int;
				subG_party1[p] = reader_int;
			}
			else {
				input >> reader_int;
				n1[p] = reader_int;
				Tn_people += n1[p];
			}
		}
	}

	subG_party.resize(Tn_party);
	n.resize(Tn_party);
	int pp = 0;
	for (int g = 0; g < n_subGroup; ++g) {
		for (int p = 0; p < Tn_party; ++p) {
			if (subG_party1[p] == g) {
				n[pp] = n1[p];
				subG_party[pp] = subG_party1[p];
				pp++;
			}
		}
	}

	// a[s]: section number for seat s // r[s]: row number for seat s // sr[s]: seat number for seat s in its corresponding row
	a.resize(n_seat);
	r.resize(n_seat);
	sr.resize(n_seat);
	int sss = 0;
	for (int aa = 0; aa < n_sec; ++aa) {
		for (int rr = 0; rr < n_row; ++rr) {
			for (int ss = 0; ss < n_SeatRow; ++ss) {
				sr[sss] = ss;
				r[sss] = rr;
				a[sss] = aa;
				++sss;
			}
		}
	}

	// m[s]: number of seats available on the right of s, including seat s
	m.resize(n_seat);
	for (int s = 0; s < n_seat; ++s) {
		m[s] = 0;
		int ss = 0;		 //defined to be stopped by the first occupied seat in the following while loop
		sss = 0;	     //defined to count the number of empty seats on the right side of seat s
		while (ss == 0 && (sr[s] + sss) < n_SeatRow) {
			ss += OcSeat[s + sss];
			if (OcSeat[s + sss] == 0) {
				++sss;
			}
		}
		m[s] += sss;
	}

	// o[s]: index of the seat whichi is in front of s. o[s]= -1 if there is no seat in front of it
	o.resize(n_seat);
	for (int s = 0; s < n_seat; ++s) {
		o[s] = -1;
	}
	for (int a = 0; a < n_sec; ++a) {
		for (int s = (a*n_SeatRow*n_row); s< ((a + 1)*n_SeatRow*n_row - n_SeatRow); ++s) {
			o[s] = s + n_SeatRow;
		}
	}

	// f[s]
	f.resize(n_seat);
	for (int s = 0; s < n_seat; s++) {
		if (sr[s] == 0) {
			f[s] = 1;
		}
		else {
			if (OcSeat[s - 1] == 1) {
				f[s] = 1;
			}
			else {
				f[s] = 0;
			}
		}
	}

	// l[s]
	l.resize(n_seat);
	for (int s = 0; s < n_seat; s++) {
		if (sr[s] == n_SeatRow - 1) {
			l[s] = 1;
		}
		else {
			if (OcSeat[s + 1] == 1) {
				l[s] = 1;
			}
			else {
				l[s] = 0;
			}
		}
	}

	//n_seg: number of segments // seg[s]: segment number for seat s (if s is unavailable, seg[s]=-1)
	n_seg = 0;
	seg.resize(n_seat, -1);
	for (int s = 0; s < n_seat; s++) {
		if (OcSeat[s] == 0) {
			if (f[s] == 1) {
				n_seg++;
			}
			seg[s] = n_seg - 1;
		}
	}

	//number of seats in each segment "no_seat_in_seg[seg]"
	no_seat_in_seg.resize(n_seg, 0);
	for (int s = 0; s < n_seat; ++s) {
		if (OcSeat[s] == 0) {
			++no_seat_in_seg[seg[s]];
		}
	}

	//index of the first seat in each segment "first_seat_seg[seg]"
	int segment = 0;
	int seat = 0;
	first_seat_seg.resize(n_seg);
	while (segment < n_seg) {
		if (seg[seat] == segment) {
			first_seat_seg[segment] = seat;
			++segment;
		}
		++seat;
	}

	//finding the number of segment of each size {1,...,n_SeatRow}
	no_seg_of_each_size.resize(n_SeatRow, 0);
	for (int segm = 0; segm < n_seg; ++segm) {
		no_seg_of_each_size[no_seat_in_seg[segm] - 1]++;
	}

	//new parameters
	max_party_size = 0;
	for (int p = 0; p < Tn_party; ++p) {
		if (n[p] > max_party_size) {
			max_party_size = n[p];
		}
	}

	//u[g][n]: number of parties of size (n+1) in group g
	u.resize(n_subGroup);
	for (int g = 0; g < n_subGroup; ++g) {
		u[g].resize(max_party_size, 0);
		for (int nn = 1; nn < max_party_size + 1; ++nn) {
			for (int p = 0; p < Tn_party; ++p) {
				if (subG_party[p] == g && n[p] == nn) {
					u[g][nn - 1]++;
				}
			}
		}
	}

	//nn[g]: number of people in group g
	//n_party_in_group[g]
	nn.resize(n_subGroup, 0);
	n_party_in_group.resize(n_subGroup, 0);
	for (int g = 0; g < n_subGroup; ++g) {
		for (int p = 0; p < Tn_party; ++p) {
			if (subG_party[p] == g) {
				nn[g] += n[p];
				n_party_in_group[g] ++;
			}
		}
	}

	//number of different item sizes
	vector<bool> diff_item_size(n_SeatRow + 1);
	for (int p = 0; p < Tn_party; ++p) {
		diff_item_size[n[p]] = true;
	}
	n_diff_item_size = 0;
	for (int item_size = 0; item_size < n_SeatRow + 1; ++item_size) {
		n_diff_item_size += diff_item_size[item_size];
	}

	// LowerBound 1 by Silvano Martello and Paolo Toth per group
	LB_MT1.resize(n_subGroup, 0);
	for (int g = 0; g < n_subGroup; ++g) {
		for (int p = 0; p < Tn_party; ++p) {
			if (subG_party[p] == g) {
				LB_MT1[g] += n[p];
			}
		}
		LB_MT1[g] = ceil(LB_MT1[g] / n_SeatRow - 0.0001);
	}

	// LowerBound 2 by Silvano Martello and Paolo Toth per group
	LB_MT2.resize(n_subGroup, 0);
	double K_MT = floor(n_SeatRow / 2 + 0.001);
	for (int g = 0; g < n_subGroup; ++g) {
		vector<bool> N_1;
		N_1.resize(Tn_party, 0);
		vector<bool> N_2;
		N_2.resize(Tn_party, 0);
		vector<bool> N_3;
		N_3.resize(Tn_party, 0);
		for (int p = 0; p < Tn_party; ++p) {
			if (subG_party[p] == g) {
				if (n[p] >(n_SeatRow - K_MT)) {
					N_1[p] = 1;
				}
				if (n[p] <= (n_SeatRow - K_MT) && n[p] >(n_SeatRow / 2)) {
					N_2[p] = 1;
				}
				if (n[p] >= K_MT && n[p] <= (n_SeatRow / 2)) {
					N_3[p] = 1;
				}
			}
		}
		double sum_term = 0;
		for (int p = 0; p < Tn_party; ++p) {
			if (subG_party[p] == g) {
				LB_MT2[g] += (N_1[p] + N_2[p]);
				if (N_3[p] == 1) {
					sum_term += n[p];
				}
				if (N_2[p] == 1) {
					sum_term += (n[p] - n_SeatRow);
				}
			}
		}
		if (ceil(sum_term / n_SeatRow - 0.001) > 0) {
			LB_MT2[g] += ceil(sum_term / n_SeatRow - 0.001);
		}
	}

	input.close();

}

void Stadium_Seating::print() {
	cout << "\n*******************************************************" << endl;
	cout << "************** Data **********************" << endl;
	cout << "*******************************************************" << endl << endl;
	cout << "Instance file : " << filename << endl << endl;

	cout << "n_sec: " << n_sec << endl;
	cout << "n_row: " << n_row << endl;
	cout << "n_SeatRow: " << n_SeatRow << endl;
	cout << "n_seat:" << n_seat << endl;
	cout << "n_seg:" << n_seg << endl;
	cout << endl;

	cout << "Binary_OcSeat: " << endl;
	int sss = 0;
	for (int s = 0; s < n_sec; ++s) {
		for (int r = 0; r < n_row; ++r) {
			for (int ss = 0; ss < n_SeatRow; ++ss) {
				cout << OcSeat[sss] << "\t";
				sss++;
			}
			cout << endl;
		}
		cout << endl;
	}
	cout << endl;

	cout << "n_subGroup: " << n_subGroup << endl;
	cout << "Tn_party: " << Tn_party << endl;

	cout << "SubGroup corresponding to each party:" << endl;
	for (int p = 0; p < Tn_party; ++p) {
		cout << subG_party[p] << "\t";
	}
	cout << endl;
	cout << endl;

	cout << "Number of people in each party:" << endl;
	for (int p = 0; p < Tn_party; ++p) {
		cout << n[p] << "\t";
	}
	cout << endl;
	cout << endl;

	cout << "section number for each seat:" << endl;
	for (int s = 0; s < n_seat; ++s) {
		cout << a[s] << "\t";
	}
	cout << endl;
	cout << endl;

	cout << "row number for each seat:" << endl;
	for (int s = 0; s < n_seat; ++s) {
		cout << r[s] << "\t";
	}
	cout << endl;
	cout << endl;

	cout << "seat number for each seat:" << endl;
	for (int s = 0; s < n_seat; ++s) {
		cout << sr[s] << "\t";
	}
	cout << endl;
	cout << endl;

	cout << "segment number for each seat " << endl;
	sss = 0;
	for (int s = 0; s < n_sec; ++s) {
		for (int r = 0; r < n_row; ++r) {
			for (int ss = 0; ss < n_SeatRow; ++ss) {
				cout << seg[sss] << "\t";
				sss++;
			}
			cout << endl;
		}
		cout << endl;
	}
	cout << endl;

	cout << "number of seats in each segment" << endl;
	for (int seg = 0; seg < n_seg; ++seg) {
		cout << no_seat_in_seg[seg] << "\t";
	}
	cout << endl;
	cout << endl;

	cout << "number of segment of each size{ 1,...,n_SeatRow }" << endl;
	no_seg_of_each_size.resize(n_SeatRow, 0);
	for (int size = 0; size < n_SeatRow; ++size) {
		cout << no_seg_of_each_size[size] << "\t";
	}
	cout << endl;
	cout << endl;

	cout << "number of available seats on the right side of each seat:" << endl;
	for (int s = 0; s < n_seat; ++s) {
		cout << m[s] << "\t";
	}
	cout << endl;
	cout << endl;

	cout << "Index of the seat in front of each seat: o[s]: " << endl;
	sss = 0;
	for (int s = 0; s < n_sec; ++s) {
		for (int r = 0; r < n_row; ++r) {
			for (int ss = 0; ss < n_SeatRow; ++ss) {
				cout << o[sss] << "\t";
				sss++;
			}
			cout << endl;
		}
		cout << endl;
	}
	cout << endl;

	cout << "f[s]: " << endl;
	sss = 0;
	for (int s = 0; s < n_sec; ++s) {
		for (int r = 0; r < n_row; ++r) {
			for (int ss = 0; ss < n_SeatRow; ++ss) {
				cout << f[sss] << "\t";
				sss++;
			}
			cout << endl;
		}
		cout << endl;
	}
	cout << endl;

	cout << "l[s]: " << endl;
	sss = 0;
	for (int s = 0; s < n_sec; ++s) {
		for (int r = 0; r < n_row; ++r) {
			for (int ss = 0; ss < n_SeatRow; ++ss) {
				cout << l[sss] << "\t";
				sss++;
			}
			cout << endl;
		}
		cout << endl;
	}
	cout << endl;

	cout << "u[g][n]:" << endl;
	for (int g = 0; g < n_subGroup; ++g) {
		for (int nn = 1; nn < max_party_size + 1; ++nn) {
			cout << u[g][nn - 1] << "\t";
		}
		cout << endl;
	}
	cout << endl;

	/*cout << "nn[g]:" << endl;
	for (int g = 0; g < n_subGroup; ++g) {
	cout << nn[g] << "\t";
	}
	cout << endl;*/
}