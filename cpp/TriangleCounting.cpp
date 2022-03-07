#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <tuple>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mt19937ar.h"
#include "MemoryOperation.h"
#include "include/stats.hpp"

using namespace std;

string EdgeFile;
int NodeNum;
double Eps;
string Eps_s;
double Mu;
string Mu_s;
char *EpsMu_s[2];

double EpsNsDeg;
double EpsAllbutNsDeg;
double Eps1st, Eps2ndTrSt;

int NSType;
double EClip;
double TClip;
string Clip_s;
int ItrNum;
int Alg;
double Balloc[2];
char *Balloc_s[2];

// Initialization of statslib
stats::rand_engine_t engine(1776);

FILE *FileOpen(string filename, const char *mode) {
	FILE *fp;

	if ((fp = fopen(filename.c_str(), mode)) == NULL) {
		cout << "cannot open " << filename << endl;
		exit(-1);
	}
	return fp;
}

int compare_double(double *x, double *y) {
	if (*x > *y)       return(1);   /* return positive integer */
	else if (*x == *y) return(0);   /* return zero     integer */
	else              return(-1);  /* return negative integer */
}

bool checkFileExistence(const std::string& str) {
    std::ifstream ifs(str);
    return ifs.is_open();
}

// Randomly generate 0, 1, 2, ..., size-1, and store the first num values into rndperm
void MakeRndPerm(int *rndperm, int size, int num) {
	int rnd;
	int *ordperm;
	int i, j;

	// 0, 1, 2, ..., size-1 --> ordperm
	ordperm = (int *)malloc(size * sizeof(int));
	for (i = 0; i < size; i++) {
		ordperm[i] = i;
	}

	for (i = 0; i < num; i++) {
		rnd = genrand_int32() % (size - i);
		rndperm[i] = ordperm[rnd];
		for (j = rnd + 1; j < size - i; j++) {
			ordperm[j - 1] = ordperm[j];
		}
	}

	free(ordperm);
}

// Read edges from the edge file
void ReadEdges(map<int, int> *a_mat, int *node_order){
	int node1, node2;
	int i;
	char s[1025];
	char *tok;
	FILE *fp;

	fp = FileOpen(EdgeFile, "r");
	for(i=0;i<3;i++) fgets(s, 1024, fp);
	while(fgets(s, 1024, fp) != NULL){
		// 1st node --> node1
		tok = strtok(s, ",");
		node1 = atoi(tok);
		// 2nd node --> node2
		tok = strtok(NULL, ",");
		node2 = atoi(tok);
		if(node1 == node2) continue;
		// If both nodes exist, add the edge
		if(node_order[node1] < NodeNum && node_order[node2] < NodeNum){
			a_mat[node_order[node1]][node_order[node2]] = 1;
			a_mat[node_order[node2]][node_order[node1]] = 1;
		}
	}
	fclose(fp);
}

// Calculate #triangles in the non-interactive local model (RR)
void CalcNLocTri(map<int, int> *a_mat, string outfile, double &tri_num_ns, int emp){
	map<int, int> *a_mat_ns;			// noisy adjacency matrix
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	int *deg_ns;									// noisy degree
	long long tot_edge_num_ns;
	long long tri_num, st2_num, ed2_num, ed1_num, non_num;
	double q;
	double alp, alp_1_3, q_inv_11, q_inv_21, q_inv_31, q_inv_41;
	double rnd;
	int i, j, k;

	// Initialization
	a_mat_ns = new map<int, int>[NodeNum];
	malloc1D(&deg_ns, NodeNum);

	// Flip probability --> q
    q = 1.0 / (exp(Eps) + 1.0);

	// Flip 0/1 in a_mat with probability q --> a_mat_ns
	for(i=0;i<NodeNum;i++){
		for(j=i+1;j<NodeNum;j++){
			rnd = genrand_real2();
			// 0 --> 1 (flip)
			if(rnd < q && a_mat[i].count(j) == 0){
				a_mat_ns[i][j] = 1;
				a_mat_ns[j][i] = 1;
			}
			// 1 --> 1 (not flip)
			else if(rnd >= q && a_mat[i].count(j) == 1){
				a_mat_ns[i][j] = 1;
				a_mat_ns[j][i] = 1;
			}
		}
	}

	// Degree --> deg_ns
	for(i=0;i<NodeNum;i++){
		for (aitr = a_mat_ns[i].begin(); aitr != a_mat_ns[i].end(); aitr++) deg_ns[i] += 1;
	}

	// Total number of edges --> tot_edge_num_ns
	tot_edge_num_ns = 0;
	for(i=0;i<NodeNum;i++) tot_edge_num_ns += (long long)deg_ns[i];
	tot_edge_num_ns /= 2;

	// #triangles --> tri_num
	tri_num = 0;
	for(i=0;i<NodeNum;i++){
		for (aitr = a_mat_ns[i].begin(); aitr != a_mat_ns[i].end(); aitr++) {
			j = aitr->first;
			if (i >= j) continue;
			for (aitr2 = a_mat_ns[i].begin(); aitr2 != a_mat_ns[i].end(); aitr2++) {
				k = aitr2->first;
				if (j >= k) continue;
				if(a_mat_ns[j].count(k) > 0) tri_num++;
			}
		}
	}

	// With empirical estimation
	if(emp == 1){
		// #2-stars --> st2_num
		st2_num = 0;
		for(i=0;i<NodeNum;i++){
			st2_num += ((long long)deg_ns[i] * ((long long)deg_ns[i]-1)) / 2;
		}

		// #2-edges --> ed2_num
		ed2_num = st2_num - 3*tri_num;
		// #1-edge --> ed1_num
		ed1_num = (long long)tot_edge_num_ns*(NodeNum-2) - 2*ed2_num - 3*tri_num;
		// #none --> non_num
		non_num = (long long)NodeNum*(NodeNum-1)*(NodeNum-2)/6 - tri_num - ed2_num - ed1_num;

		alp = exp(Eps);
		alp_1_3 = (alp-1.0)*(alp-1.0)*(alp-1.0);
		q_inv_11 = (alp*alp*alp) / alp_1_3;
		q_inv_21 = - alp*alp / alp_1_3;
		q_inv_31 = alp / alp_1_3;
		q_inv_41 = - 1.0 / alp_1_3;

		tri_num_ns = (double)tri_num * q_inv_11 + (double)ed2_num * q_inv_21 + (double)ed1_num * q_inv_31 + (double)non_num * q_inv_41;
	}
	// Without empirical estimation
	else{
		tri_num_ns = (double)tri_num;
	}

	delete[] a_mat_ns;
	free1D(deg_ns);
}

// Calculate #triangles in the non-interactive local model (ARR)
void CalcNLocTriARR(map<int, int> *a_mat, string outfile, double &tri_num_ns, int emp){
	map<int, int> *a_mat_ns;			// noisy adjacency matrix
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	int *deg_ns;									// noisy degree
	long long tot_edge_num_ns;
	long long tri_num, st2_num, ed2_num, ed1_num, non_num;
	double tri_num_bs, ed2_num_bs, ed1_num_bs, non_num_bs;
	double q;
	double alp, alp_1_3, q_inv_11, q_inv_21, q_inv_31, q_inv_41;
	double rnd;
	int i, j, k;
	double murho, p1, p2, q2;

	// Initialization
	a_mat_ns = new map<int, int>[NodeNum];
	malloc1D(&deg_ns, NodeNum);

	// Parameters in asymmetric RR --> Mu (1 --> 1), murho (0 --> 1)
	murho = Mu / exp(Eps);
	// Sampling rate --> p2
	p1 = exp(Eps) / (exp(Eps) + 1.0);
	p2 = Mu / p1;

	// Flip 0/1 in a_mat with probability q --> a_mat_ns
	for(i=0;i<NodeNum;i++){
		for(j=i+1;j<NodeNum;j++){
			rnd = genrand_real2();
			// 0 --> 1 (flip)
			if(rnd < murho && a_mat[i].count(j) == 0){
				a_mat_ns[i][j] = 1;
				a_mat_ns[j][i] = 1;
			}
			// 1 --> 1 (not flip)
			else if(rnd < Mu && a_mat[i].count(j) == 1){
				a_mat_ns[i][j] = 1;
				a_mat_ns[j][i] = 1;
			}
		}
	}

	// Degree --> deg_ns
	for(i=0;i<NodeNum;i++){
		for (aitr = a_mat_ns[i].begin(); aitr != a_mat_ns[i].end(); aitr++) deg_ns[i] += 1;
	}

	// Total number of edges --> tot_edge_num_ns
	tot_edge_num_ns = 0;
	for(i=0;i<NodeNum;i++) tot_edge_num_ns += (long long)deg_ns[i];
	tot_edge_num_ns /= 2;

	// #triangles --> tri_num
	tri_num = 0;
	for(i=0;i<NodeNum;i++){
		for (aitr = a_mat_ns[i].begin(); aitr != a_mat_ns[i].end(); aitr++) {
			j = aitr->first;
			if (i >= j) continue;
			for (aitr2 = a_mat_ns[i].begin(); aitr2 != a_mat_ns[i].end(); aitr2++) {
				k = aitr2->first;
				if (j >= k) continue;
				if(a_mat_ns[j].count(k) > 0) tri_num++;
			}
		}
	}

	// With empirical estimation
	if(emp == 1){
		// #2-stars --> st2_num
		st2_num = 0;
		for(i=0;i<NodeNum;i++){
			st2_num += ((long long)deg_ns[i] * ((long long)deg_ns[i]-1)) / 2;
		}

		// #2-edges --> ed2_num
		ed2_num = st2_num - 3*tri_num;
		// #1-edge --> ed1_num
		ed1_num = (long long)tot_edge_num_ns*(NodeNum-2) - 2*ed2_num - 3*tri_num;

		// Calculate #triangles, #2-edges, #1-edge before sampling --> tri_num_bs, ed2_num_bs, ed1_num_bs
		q2 = 1.0 - p2;
		tri_num_bs = (double)tri_num / (p2 * p2 * p2);
		ed2_num_bs = (double)ed2_num / (p2 * p2) - 3.0 * q2 * tri_num_bs;
		ed1_num_bs = (double)ed1_num / p2 - 2.0 * q2 * ed2_num_bs - 3.0 * q2 * q2 * tri_num_bs;

		// #none --> non_num_bs
		non_num_bs = (double)NodeNum*(NodeNum-1)*(NodeNum-2)/6 - tri_num_bs - ed2_num_bs - ed1_num_bs;

		alp = exp(Eps);
		alp_1_3 = (alp-1.0)*(alp-1.0)*(alp-1.0);
		q_inv_11 = (alp*alp*alp) / alp_1_3;
		q_inv_21 = - alp*alp / alp_1_3;
		q_inv_31 = alp / alp_1_3;
		q_inv_41 = - 1.0 / alp_1_3;

		tri_num_ns = tri_num_bs * q_inv_11 + ed2_num_bs * q_inv_21 + ed1_num_bs * q_inv_31 + non_num_bs * q_inv_41;
	}
	// Without empirical estimation
	else{
		tri_num_ns = (double)tri_num;
	}

	delete[] a_mat_ns;
	free1D(deg_ns);
}

// Calculate #2-stars and #3-stars in the non-interactive local model
void CalcNLocSt(long long st2_num, long long st3_num, int *deg, string outfile, double &st2_num_ns, double &st3_num_ns, double &sen_st2, double &sen_st3){
	int max_deg;
	int i;
	FILE *fp;

	double sen;
	int deg_ns_floor;
	double *deg_ns;

    // Initialization
    st2_num_ns = st2_num;
    st3_num_ns = st3_num;
	malloc1D(&deg_ns, NodeNum);

	// no noise
	if(NSType == -1){
	}
	// Lap (max degree)
	else if(NSType == 0){
        // Sensitivity using max degree --> sen_st2, sen_st3
        sen_st2 = sen_st3 = 0;
		// max(deg) --> max_deg
		max_deg = 0;
		for(i=0;i<NodeNum;i++){
			if(max_deg < deg[i]) max_deg = deg[i];
		}
		sen_st2 = (double)max_deg;
		sen_st3 = (double)max_deg * ((double)max_deg - 1.0) / 2.0;

		// Add Lap for each user
		for(i=0;i<NodeNum;i++){
			// Add Lap(sen_st2/Eps) --> st2_num_ns
			st2_num_ns += stats::rlaplace(0.0, sen_st2/Eps, engine);
			// Add Lap(sen_st3/Eps) --> st3_num_ns
			st3_num_ns += stats::rlaplace(0.0, sen_st3/Eps, engine);
		}
	}
	// Lap (degree)
	else if(NSType == 1){
		// Average sensitivity using each user's degree --> sen_st2, sen_st3
        sen_st2 = sen_st3 = 0;
		// Add Lap for each user
		for(i=0;i<NodeNum;i++){
	        // Sensitivity using each user's degree --> sen
			sen = (double)deg[i];
			// Add Lap(sen/Eps) --> st2_num_ns
			st2_num_ns += stats::rlaplace(0.0, sen/Eps, engine);
			sen_st2 += sen;

	        // Sensitivity using each user's degree --> sen
			sen = (double)deg[i] * ((double)deg[i] - 1.0) / 2.0;
			// Add Lap(sen/Eps) --> st3_num_ns
			st3_num_ns += stats::rlaplace(0.0, sen/Eps, engine);
			sen_st3 += sen;
		}
		sen_st2 /= (double)NodeNum;
		sen_st3 /= (double)NodeNum;
	}
	// Lap (noisy degree)
	else if(NSType == 2 || NSType == 3){
		// Noisy degree --> deg_ns
		for(i=0;i<NodeNum;i++){
			deg_ns[i] = (double)deg[i] + stats::rlaplace(0.0, 1.0/EpsNsDeg, engine);
		}
		// Add positive bias (EClip) --> deg_ns
		if(EClip != -1){
			for(i=0;i<NodeNum;i++) deg_ns[i] += EClip;
		}
		for(i=0;i<NodeNum;i++) deg_ns[i] = max(deg_ns[i], 0.0);

		// Average sensitivity using each user's degree --> sen_st2, sen_st3
		sen_st2 = sen_st3 = 0.0;
		// Graph projection for each user
		st2_num_ns = st3_num_ns = 0;
		for(i=0;i<NodeNum;i++){
			deg_ns_floor = (int)floor(deg_ns[i]);
			// If deg[i] exceeds floor(deg_ns[i]), then perform graph projection
			if((double)deg[i] > deg_ns_floor){
				st2_num_ns += ((long long)deg_ns_floor * ((long long)deg_ns_floor-1)) / 2;
				st3_num_ns += ((long long)deg_ns_floor * ((long long)deg_ns_floor-1) * ((long long)deg_ns_floor-2)) / 6;
			}
			else{
				st2_num_ns += ((long long)deg[i] * ((long long)deg[i]-1)) / 2;
				st3_num_ns += ((long long)deg[i] * ((long long)deg[i]-1) * ((long long)deg[i]-2)) / 6;
			}

			// Sensitivity for #2-stars --> sen
			sen = (double)deg_ns_floor;
			// Add Lap(sen/Eps) --> st2_num_ns
			st2_num_ns += stats::rlaplace(0.0, sen/EpsAllbutNsDeg, engine);
			sen_st2 += sen;

			// Sensitivity for #3-stars --> sen
			sen = (double)deg_ns_floor * ((double)deg_ns_floor - 1.0) / 2.0;
			// Add Lap(sen/Eps) --> st3_num_ns
			st3_num_ns += stats::rlaplace(0.0, sen/EpsAllbutNsDeg, engine);
			sen_st3 += sen;
		}
		sen_st2 /= (double)NodeNum;
		sen_st3 /= (double)NodeNum;
	}

	free1D(deg_ns);
}

// Calculate #triangles by [Ye+, T-KDE]
void CalcTKDETri(map<int, int> *a_mat, int *deg, string outfile, double &tri_num_ns, double &st2_num_ns, double &sen_tri, double &sen_st2, int alg){
	double eps_opt, eps_ns, eps_rr, eps_lap;
	double sen_deg;
	double *deg_ns;
	double hd;
	double alp, alp_opt;
	double loss, loss_min;
	double exp_ae, exp_3ae, loss_left, loss_right, numerator, denominator;
	map<int, int> *a_mat_ns;			// noisy adjacency matrix
	double p, q;
	double gamma, gamma_ns;
	double r_tri_left, r_tri_right1, r_tri_right2, r_tri_right3, r_tri_right4;
	double *r_tri_num;
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	long long *tri_num;
	int i, j, k;

	// Initialization
	malloc1D(&deg_ns, NodeNum);
	malloc1D(&tri_num, NodeNum);
	malloc1D(&r_tri_num, NodeNum);
	a_mat_ns = new map<int, int>[NodeNum];

	// Epsilon for optimizing alpha (Theorem 5.1 in [Ye+, T-KDE]; [Ye+, T-KDE] used 10% of the privacy budget (see Section 7))
	eps_opt = Eps / 10;
	// Epsilon for adding noise (RR and Lap)
	eps_ns = Eps - eps_opt;

	// Global sensitivity of the degree (=1) --> sen_deg
	sen_deg = 1.0;

	// Add Lap((global sensitivity/eps_opt)) to each user's degree --> deg_ns
	for(i=0;i<NodeNum;i++){
		deg_ns[i] = (double)deg[i] + stats::rlaplace(0.0, sen_deg/eps_opt, engine);
	}

	// Calculate the mean of deg_ns --> hd
	if(alg == 7){
		hd = 0.0;
		for(i=0;i<NodeNum;i++) hd += deg_ns[i];
		hd /= (double)NodeNum;
	}
	// Calculate the median of deg_ns --> hd
	else if(alg == 8){
		double *deg_ns_sort;
		malloc1D(&deg_ns_sort, NodeNum);

		// Sort deg_ns --> deg_ns_sort
		for(i=0;i<NodeNum;i++) deg_ns_sort[i] = deg_ns[i];
		qsort(deg_ns_sort, NodeNum, sizeof(double), (int(*) (const void *, const void *))compare_double);

		// Median of deg_ns_sort --> hd
		hd = deg_ns_sort[int(NodeNum / 2)];

		free1D(deg_ns_sort);
	}
	// Calculate the most frequent degree in deg_ns --> hd
	else if(alg == 9){
		int *deg_hist;
		int deg_ns_int;
		int deg_hist_max;
		malloc1D(&deg_hist, NodeNum);
		// Histogram of deg_ns --> deg_hist
		for(i=0;i<NodeNum;i++) deg_hist[i] = 0;
		for(i=0;i<NodeNum;i++){
			deg_ns_int = int(round(deg_ns[i]));
			if(deg_ns_int >= 0 && deg_ns_int < NodeNum) deg_hist[deg_ns_int] += 1;
		}
		// Most frequent degree --> hd
		deg_hist_max = 0;
		for(i=0;i<NodeNum;i++){
			if(deg_hist_max < deg_hist[i]){
				deg_hist_max = deg_hist[i];
				hd = i;
			}
		}
		free1D(deg_hist);
	}

	// Optimization of alpha --> alp_opt
	numerator = 8.0 * (10.0 * hd * hd - 10.0 * hd + 3.0);
	alp_opt = 0.0;
	for(alp=0.0;alp<1.0;alp+=0.0001){
		// loss function (Theorem 5.1 in [Ye+, T-KDE]) --> loss
		exp_ae = exp(alp * eps_ns);
		exp_3ae = exp(3.0 * alp * eps_ns);
		loss_left = (exp_ae + 2.0) / (exp_3ae * (exp_ae - 1.0) * (exp_ae - 1.0));
		denominator = hd * hd * (hd - 1.0) * (hd - 1.0) * (1.0 - alp) * (1.0 - alp) * eps_ns * eps_ns;
		loss_right = 1.0 + numerator / denominator;
		loss = loss_left * loss_right;

		// Update the minimum loss and optimal alpha --> loss_min, alp_opt
		if(alp == 0.0 || loss_min > loss){
			loss_min = loss;
			alp_opt = alp;
		}
	}

	// epsilon for RR --> eps_rr
	eps_rr = alp_opt * eps_ns;
	// epsilon for Lap --> eps_lap
	eps_lap = (1.0 - alp_opt) * eps_ns;

	// Flip probability --> q
    q = 1.0 / (exp(eps_rr) + 1.0);
	// Non-flip probability --> p
	p = 1.0 - q;

	// Flip 0/1 in a_mat with probability q --> a_mat_ns
	// (NOTE: the following algorithm is exactly the same as RABV in [Ye+, T-KDE])
	for(i=0;i<NodeNum;i++){
		for(j=i+1;j<NodeNum;j++){
			// 0 --> 1 (flip)
			if(genrand_real2() < q && a_mat[i].count(j) == 0){
				a_mat_ns[i][j] = 1;
				a_mat_ns[j][i] = 1;
			}
			// 1 --> 1 (not flip)
			else if(genrand_real2() >= q && a_mat[i].count(j) == 1){
				a_mat_ns[i][j] = 1;
				a_mat_ns[j][i] = 1;
			}
		}
	}

	// Add Lap((global sensitivity/eps_lap)) to each user's degree --> deg_ns
	for(i=0;i<NodeNum;i++){
		deg_ns[i] = (double)deg[i] + stats::rlaplace(0.0, sen_deg/eps_lap, engine);
	}

	// Calculate gamma using deg_ns --> gamma
	gamma = 0.0;
	for(i=0;i<NodeNum;i++) gamma += deg_ns[i];
	gamma /= (NodeNum * (NodeNum - 1.0));
	// Noisy gamma --> gamma_ns
	gamma_ns = gamma * p + (1.0 - gamma) * (1.0 - p);

	// #triangles involving user i (j < k) --> tri_num[i]
	for(i=0;i<NodeNum;i++){
		tri_num[i] = 0;
		for (aitr = a_mat_ns[i].begin(); aitr != a_mat_ns[i].end(); aitr++) {
			j = aitr->first;
			if (i == j) continue;
			for (aitr2 = a_mat_ns[i].begin(); aitr2 != a_mat_ns[i].end(); aitr2++) {
				k = aitr2->first;
				if (i == k || j >= k) continue;
				if(a_mat_ns[j].count(k) > 0) tri_num[i]++;
			}
		}
	}

	// Bias-reduced estimate of #triangles for each user (Eq.(10) in [Ye+, T-KDE]) --> r_tri_num[i]
	for(i=0;i<NodeNum;i++){
		// Use p, tri_num[i], deg_ns[i], NodeNum, gamma_ns
		r_tri_left = p * p * (2.0 * p - 1.0);
		r_tri_right1 = (double)tri_num[i];
		r_tri_right2 = (deg_ns[i] * (deg_ns[i] - 1.0) * p * p * (1.0 - p)) / 2.0;
		r_tri_right3 = deg_ns[i] * ((double)NodeNum - deg_ns[i] - 1.0) * p * (1.0 - p) * gamma_ns;
		r_tri_right4 = (((double)NodeNum - deg_ns[i] - 1.0) * ((double)NodeNum - deg_ns[i] - 2.0) * (1.0 - p) * (1.0 - p) * gamma_ns) / 2.0;
		r_tri_num[i] = (r_tri_right1 - r_tri_right2 - r_tri_right3 - r_tri_right4) / r_tri_left;
	}

	// Bias-reduced estimate of #triangles --> tri_num_ns
	tri_num_ns = 0.0;
	for(i=0;i<NodeNum;i++) tri_num_ns += r_tri_num[i];
	// Sum of r_tri_num has three counts for each triangle (i, j, k) with i < j < k
	tri_num_ns /= 3.0;

	// Calculate the estimate of #2-stars using deg_ns --> st2_num_ns
	st2_num_ns = 0;
	for(i=0;i<NodeNum;i++){
		st2_num_ns += (deg_ns[i] * (deg_ns[i]-1)) / 2;
	}

	free1D(deg_ns);
	free1D(tri_num);
	free1D(r_tri_num);
	delete[] a_mat_ns;
}

// Calculate #2-stars ad #3-stars by [Ye+, T-KDE]
void CalcTKDESt(int *deg, string outfile, double &st2_num_ns, double &st3_num_ns, double &sen_st2, double &sen_st3){
	double sen_deg;
	double *deg_ns;
	int glb_max_deg;
	int max_deg;
	int sen_st2_ij, sen_st3_ij;
	double EpsLap;
	double max_deg_ns;
	int max_deg_ns_floor;
	int del_num;
	int *rndperm;
	int i, j, x;
	FILE *fp;

    // Initialization
	malloc1D(&deg_ns, NodeNum);

	// Global sensitivity of the degree (=1) --> sen_deg
	sen_deg = 1.0;

	// Add Lap((global sensitivity/Eps)) to each user's degree --> deg_ns
	for(i=0;i<NodeNum;i++){
		deg_ns[i] = (double)deg[i] + stats::rlaplace(0.0, sen_deg/Eps, engine);
	}

	// #2-stars, #3-stars --> st2_num_ns, st3_num_ns
	st2_num_ns = st3_num_ns = 0;
	for(i=0;i<NodeNum;i++){
		st2_num_ns += (deg_ns[i] * (deg_ns[i]-1)) / 2;
		st3_num_ns += (deg_ns[i] * (deg_ns[i]-1) * (deg_ns[i]-2)) / 6;
	}

	// sensitivity --> sen_st2, sen_st3
	sen_st2 = sen_deg;
    sen_st3 = sen_deg;

	free1D(deg_ns);
}

// Calculate #triangles and #2-stars in the interactive local model
void CalcILocTri(map<int, int> *a_mat, int *deg, string outfile, double &tri_num_ns, double &sen_tri){
	map<int, int> *a_mat_ns;			// noisy adjacency matrix
	map<int, int> *a_mat_del;			// deleted adjacency matrix after projection
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	map<tuple<int, int>, int> a2_mat;
	map<tuple<int, int>, int>::iterator a2_itr;
	double p, q;
	double *tri_num_u, *st2_num_u, *trist2_num_u;
	int max_deg;
	int sen_st2_ij;
	int del_num;
	int *rndperm;
	double rnd;
	int i, j, k, x;
	FILE *fp;

	double *deg_ns;
	double sen;
	int deg_ns_floor;

	// Initialization
	a_mat_ns = new map<int, int>[NodeNum];
	a_mat_del = new map<int, int>[NodeNum];
	malloc1D(&tri_num_u, NodeNum);
	malloc1D(&st2_num_u, NodeNum);
	malloc1D(&trist2_num_u, NodeNum);
	malloc1D(&deg_ns, NodeNum);

	// Flip probability --> q
    q = 1.0 / (exp(Eps1st) + 1.0);
    p = 1.0 - q;

	// RR --> a_mat_ns
    // Count #noisy triangles and #noisy 2-stars for each user --> tri_num_u, st2_num_u
	for(i=0;i<NodeNum;i++){
		for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
			j = aitr->first;
			for (aitr2 = a_mat[i].begin(); aitr2 != a_mat[i].end(); aitr2++) {
				k = aitr2->first;
				// it must be j < k < i.
				if (j >= k || k >= i) continue;

				st2_num_u[i] += 1.0;

				// If a_mat_ns[j][k] does not exist
				if(a_mat_ns[j].count(k) == 0){
					// Flip 0/1 in a_mat[j][k] with probability q --> a_mat_ns[j][k]
					rnd = genrand_real2();
					// 0 --> 1 (flip)
					if(rnd < q && a_mat[j].count(k) == 0){
						a_mat_ns[j][k] = 1;
					}
					// 1 --> 1 (not flip)
					else if(rnd >= q && a_mat[j].count(k) == 1){
						a_mat_ns[j][k] = 1;
					}
					// 1 --> 0 (flip) or 0 --> 0 (not flip)
					else{
						a_mat_ns[j][k] = 0;
					}
				}
				if(a_mat_ns[j][k] == 1) tri_num_u[i] += 1.0;
			}
		}
	}

	// #triangles - (1 - p) * #2-stars --> trist2_num_u
	for(i=0;i<NodeNum;i++) trist2_num_u[i] = tri_num_u[i] - (1 - p) * st2_num_u[i];

	// no noise
	if(NSType == -1){
	}
	// Lap (max degree)
	else if(NSType == 0){
		// Sensitivity using max degree --> sen_tri
		sen_tri = 0.0;
		// max(deg) --> max_deg
		max_deg = 0;
		for(i=0;i<NodeNum;i++){
			if(max_deg < deg[i]) max_deg = deg[i];
		}
		sen_tri = (double)max_deg;

		// Add Lap for each user
		for(i=0;i<NodeNum;i++){
			// Add Lap(sen_tri/Eps2ndTrSt) --> trist2_num_u
			trist2_num_u[i] += stats::rlaplace(0.0, sen_tri/Eps2ndTrSt, engine);
		}
	}
	// Lap (degree)
	else if(NSType == 1){
		// Average sensitivity using each user's degree --> sen_tri
		sen_tri = 0.0;
		// Add Lap for each user
		for(i=0;i<NodeNum;i++){
			// Sensitivity using each user's degree --> sen
			sen = (double)deg[i];
			// Add Lap(sen_tri/Eps2ndTrSt) --> trist2_num_u
			trist2_num_u[i] += stats::rlaplace(0.0, sen/Eps2ndTrSt, engine);
			sen_tri += sen;
		}
		sen_tri /= (double)NodeNum;
	}
	// Lap (noisy degree)
	else if(NSType == 2 || NSType == 3){
		// Noisy degree --> deg_ns
		for(i=0;i<NodeNum;i++){
			deg_ns[i] = (double)deg[i] + stats::rlaplace(0.0, 1.0/EpsNsDeg, engine);
		}
		// Add positive bias (EClip) --> deg_ns
		if(EClip != -1){
			for(i=0;i<NodeNum;i++) deg_ns[i] += EClip;
		}
		for(i=0;i<NodeNum;i++) deg_ns[i] = max(deg_ns[i], 0.0);

		// Graph projection for each user
		for(i=0;i<NodeNum;i++){
			// If deg[i] exceeds deg_ns[i], then perform graph projection
			if((double)deg[i] > deg_ns[i]){
				// Randomly generate 0, 1, ..., deg[i]-1 --> rndperm
				malloc1D(&rndperm, deg[i]);
				MakeRndPerm(rndperm, deg[i], deg[i]);

				// Randomly delete (deg[i] - deg_ns[i]) edges from a_mat[i]
				deg_ns_floor = (int)floor(deg_ns[i]);
				x = 0;
				for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
					if(rndperm[x] >= deg_ns_floor){
						j = aitr->first;
						// Deleted edge --> a_mat_del[i][j]
						a_mat_del[i][j] = 1;
					}
					x++;
				}
				free1D(rndperm);

				// Count #noisy triangles and #noisy 2-stars again --> tri_num_u, st2_num_u
				tri_num_u[i] = 0.0;
				st2_num_u[i] = 0.0;
				for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
					j = aitr->first;
					// Continue if the edge is deleted
					if(a_mat_del[i].count(j) == 1) continue;
					for (aitr2 = a_mat[i].begin(); aitr2 != a_mat[i].end(); aitr2++) {
						k = aitr2->first;
						// Continue if the edge is deleted
						if(a_mat_del[i].count(k) == 1) continue;
						if (j >= k || k >= i) continue;
						st2_num_u[i] += 1.0;
						if(a_mat_ns[j][k] == 1) tri_num_u[i] += 1.0;
					}
				}
				// #triangles - (1 - p) * #2-stars --> trist2_num_u
				trist2_num_u[i] = tri_num_u[i] - (1 - p) * st2_num_u[i];
			}
		}

		// Average sensitivity using each user's degree --> sen_tri
		sen_tri = 0.0;
		// Add Lap for each user
		for(i=0;i<NodeNum;i++){
			// Sensitivity using each user's degree --> sen
			sen = (double)deg[i];
			// Add Lap(sen_tri/Eps2ndTrSt) --> trist2_num_u
			trist2_num_u[i] += stats::rlaplace(0.0, sen/Eps2ndTrSt, engine);
			sen_tri += sen;
		}
		sen_tri /= (double)NodeNum;
	}

    // Empirical estimate --> tri_num_ns
    tri_num_ns = 0;
	for(i=0;i<NodeNum;i++) tri_num_ns += trist2_num_u[i];
	// Divide #triangles by 2p - 1 (= 1 - 2q)
	tri_num_ns /= (2.0 * p - 1.0);

	delete[] a_mat_ns;
	delete[] a_mat_del;
	free1D(tri_num_u);
	free1D(st2_num_u);
	free1D(trist2_num_u);
	free1D(deg_ns);
}

double CalcKL(double p1, double p2){
	return p1 * log(p1 / p2) + (1.0 - p1) * log((1.0 - p1) / (1.0 - p2));
}

double CalcRedSen(double deg_ns, int alg){
	double kappa, kappa_min, kappa_max;
	double mu_pow;
	double tclip_prob, tclip_prob_thr;
	double kappa_deg_ns, kl;

	if (alg <= 0 || alg >= 4){
		printf("Error: incorrect alg @ CalcRedSen\n");
		exit(-1);
	}

	if(deg_ns > 0.0){
		if (alg == 1) mu_pow = Mu;
		else mu_pow = Mu * Mu;

		if (alg == 1 || alg == 2) kappa_min = mu_pow * deg_ns;
		else kappa_min = Mu * Mu * Mu * deg_ns;
		kappa_max = deg_ns;
		tclip_prob_thr = pow(10, -TClip);

		// Find kappa s.t. the triangle clipping probability is small enough
		for(kappa = kappa_min; kappa <= kappa_max; kappa += kappa_min){
			// triangle clipping probability --> tclip_prob
			if (alg == 3 && kappa < Mu * Mu * deg_ns) tclip_prob = Mu;
			else{
				kappa_deg_ns = kappa / deg_ns;
				kl = CalcKL(kappa_deg_ns, mu_pow);
				tclip_prob = exp(- deg_ns * kl);
				if(alg == 3) tclip_prob = Mu * tclip_prob;
			}

			if(tclip_prob <= tclip_prob_thr) break;
		}
	}
	else{
		kappa = 0.0;
	}

	return kappa;
}

// Calculate #triangles and #2-stars in the interactive local model (efficient algorithm I)
void CalcILocTriE1(map<int, int> *a_mat, int *deg, string outfile, double &tri_num_ns, double &sen_tri, double &eclip_sum, double &tclip_sum, int &eclip_num, int &tclip_num){
	map<int, int> *a_mat_ns;			// noisy adjacency matrix
	map<int, int> *a_mat_del;			// deleted adjacency matrix after projection
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	map<tuple<int, int>, int> a2_mat;
	map<tuple<int, int>, int>::iterator a2_itr;
	double murho;
	double *tri_num_u, *st2_num_u, *trist2_num_u;
	int max_deg;
	int sen_st2_ij;
	int del_num;
	int *rndperm;
	double rnd;
	int i, j, k, x;
	FILE *fp;

	double *deg_ns;
	double *red_sen;
	double tri_num_u_ij;
	double sen;
	int deg_ns_floor;

	// Initialization
	a_mat_ns = new map<int, int>[NodeNum];
	a_mat_del = new map<int, int>[NodeNum];
	malloc1D(&tri_num_u, NodeNum);
	malloc1D(&st2_num_u, NodeNum);
	malloc1D(&trist2_num_u, NodeNum);
	malloc1D(&deg_ns, NodeNum);
	malloc1D(&red_sen, NodeNum);

	// Parameters in asymmetric RR --> Mu (1 --> 1), murho (0 --> 1)
	murho = Mu / exp(Eps1st);

	// no noise, Lap (max degree), Lap (degree)
	if(NSType == -1 || NSType == 0 || NSType == 1){
		// Reduced sensitivity --> red_sen
		if(TClip != -1){
			for(i=0;i<NodeNum;i++) red_sen[i] = CalcRedSen((double)deg[i], 1);
		}

		// asymmetric RR --> a_mat_ns
		// Count #noisy triangles and #noisy 2-stars for each user --> tri_num_u, st2_num_u
		for(i=0;i<NodeNum;i++){
			for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
				j = aitr->first;
				tri_num_u_ij = 0;
				for (aitr2 = a_mat[i].begin(); aitr2 != a_mat[i].end(); aitr2++) {
					k = aitr2->first;
					// it must be j < k < i.
					if (j >= k || k >= i) continue;

					st2_num_u[i] += 1.0;

					// If a_mat_ns[j][k] does not exist
					if(a_mat_ns[j].count(k) == 0){
						// Flip 0/1 in a_mat[j][k] --> a_mat_ns[j][k]
						rnd = genrand_real2();
						// 0 --> 1 (flip)
						if(rnd < murho && a_mat[j].count(k) == 0){
							a_mat_ns[j][k] = 1;
						}
						// 1 --> 1 (not flip)
						else if(rnd < Mu && a_mat[j].count(k) == 1){
							a_mat_ns[j][k] = 1;
						}
						// 1 --> 0 (flip) or 0 --> 0 (not flip)
						else{
							a_mat_ns[j][k] = 0;
						}
					}
					if(a_mat_ns[j][k] == 1){
						if(TClip != -1) tri_num_u_ij += 1.0;
						else tri_num_u[i] += 1.0;
					}
				}

				// Triangle clipping
				if(TClip != -1){
					if(tri_num_u_ij > red_sen[i]){
						tclip_sum += (tri_num_u_ij - red_sen[i]);
						tclip_num += 1;
						tri_num_u_ij = red_sen[i];
					}
					tri_num_u[i] += tri_num_u_ij;
				}
			}
		}

		// #triangles - murho * #2-stars --> trist2_num_u
		for(i=0;i<NodeNum;i++) trist2_num_u[i] = tri_num_u[i] - murho * st2_num_u[i];

		// no noise
		if(NSType == -1){
		}
		// Lap (max degree)
		else if(NSType == 0){
			// max(deg) --> max_deg
			max_deg = 0;
			for(i=0;i<NodeNum;i++){
				if(max_deg < deg[i]) max_deg = deg[i];
			}

			// Sensitivity using max degree --> sen_tri
			sen_tri = (double)max_deg;

			// Add Lap for each user
			for(i=0;i<NodeNum;i++){
				// Add Lap(sen_tri/Eps2ndTrSt) --> trist2_num_u
				trist2_num_u[i] += stats::rlaplace(0.0, sen_tri/Eps2ndTrSt, engine);
			}
		}
		// Lap (degree)
		else if(NSType == 1){
			// Reduced sensitivity
			if(TClip != -1){
				// Add Lap for each user
				sen_tri = 0.0;
				for(i=0;i<NodeNum;i++){
					sen_tri += red_sen[i];
					trist2_num_u[i] += stats::rlaplace(0.0, red_sen[i]/Eps2ndTrSt, engine);
				}
				sen_tri /= (double)NodeNum;
			}
			// Sensitivity using each user's degree
			else{
				// Average sensitivity using each user's degree --> sen_tri
				sen_tri = 0.0;
				// Add Lap for each user
				for(i=0;i<NodeNum;i++){
					// Sensitivity using each user's degree --> sen
					sen = (double)deg[i];
					// Add Lap(sen_tri/Eps2ndTrSt) --> trist2_num_u
					trist2_num_u[i] += stats::rlaplace(0.0, sen/Eps2ndTrSt, engine);
					sen_tri += sen;
				}
				sen_tri /= (double)NodeNum;
			}
		}
	}
	// Lap (noisy degree)
	else if(NSType == 2 || NSType == 3){
		// Noisy degree --> deg_ns
		for(i=0;i<NodeNum;i++){
			deg_ns[i] = (double)deg[i] + stats::rlaplace(0.0, 1.0/EpsNsDeg, engine);
		}
		// Add positive bias (EClip) --> deg_ns
		if(EClip != -1){
			for(i=0;i<NodeNum;i++) deg_ns[i] += EClip;
		}
		for(i=0;i<NodeNum;i++) deg_ns[i] = max(deg_ns[i], 0.0);

		// Reduced sensitivity --> red_sen
		if(TClip != -1){
			for(i=0;i<NodeNum;i++) red_sen[i] = CalcRedSen(deg_ns[i], 1);
		}

		// Graph projection (edge clipping) for each user
		for(i=0;i<NodeNum;i++){
			// If deg[i] exceeds deg_ns[i], then perform graph projection (edge clipping)
			if((double)deg[i] > deg_ns[i]){
				eclip_sum += ((double)deg[i] - floor(deg_ns[i]));
				eclip_num += 1;

				// Randomly generate 0, 1, ..., deg[i]-1 --> rndperm
				malloc1D(&rndperm, deg[i]);
				MakeRndPerm(rndperm, deg[i], deg[i]);

				// Randomly delete (deg[i] - floor(deg_ns[i])) edges from a_mat[i]
				deg_ns_floor = (int)floor(deg_ns[i]);
				x = 0;
				for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
					if(rndperm[x] >= deg_ns_floor){
						j = aitr->first;
						// Deleted edge --> a_mat_del[i][j]
						a_mat_del[i][j] = 1;
					}
					x++;
				}
				free1D(rndperm);
			}

			// Count #noisy triangles and #noisy 2-stars --> tri_num_u, st2_num_u
			tri_num_u[i] = 0.0;
			st2_num_u[i] = 0.0;
			for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
				j = aitr->first;
				// Continue if the edge is deleted
				if(a_mat_del[i].count(j) == 1) continue;
				tri_num_u_ij = 0;
				for (aitr2 = a_mat[i].begin(); aitr2 != a_mat[i].end(); aitr2++) {
					k = aitr2->first;
					// Continue if the edge is deleted
					if(a_mat_del[i].count(k) == 1) continue;
					if (j >= k || k >= i) continue;
					st2_num_u[i] += 1.0;

					// If a_mat_ns[j][k] does not exist
					if(a_mat_ns[j].count(k) == 0){
						// Flip 0/1 in a_mat[j][k] --> a_mat_ns[j][k]
						rnd = genrand_real2();
						// 0 --> 1 (flip)
						if(rnd < murho && a_mat[j].count(k) == 0){
							a_mat_ns[j][k] = 1;
						}
						// 1 --> 1 (not flip)
						else if(rnd < Mu && a_mat[j].count(k) == 1){
							a_mat_ns[j][k] = 1;
						}
						// 1 --> 0 (flip) or 0 --> 0 (not flip)
						else{
							a_mat_ns[j][k] = 0;
						}
					}

					if(a_mat_ns[j][k] == 1){
						if(TClip != -1) tri_num_u_ij += 1.0;
						else tri_num_u[i] += 1.0;
					}
				}
				// Triangle clipping
				if(TClip != -1){
					if(tri_num_u_ij > red_sen[i]){
						tclip_sum += (tri_num_u_ij - red_sen[i]);
						tclip_num += 1;
						tri_num_u_ij = red_sen[i];
					}
					tri_num_u[i] += tri_num_u_ij;
				}
			}
			// #triangles - murho * #2-stars --> trist2_num_u
			trist2_num_u[i] = tri_num_u[i] - murho * st2_num_u[i];
		}

		// Reduced sensitivity
		if(TClip != -1){
			// Add Lap for each user
			sen_tri = 0.0;
			for(i=0;i<NodeNum;i++){
				sen_tri += red_sen[i];
				trist2_num_u[i] += stats::rlaplace(0.0, red_sen[i]/Eps2ndTrSt, engine);
			}
			sen_tri /= (double)NodeNum;
		}
		// Sensitivity using each user's degree
		else{
			// Average sensitivity using each user's degree --> sen_tri
			sen_tri = 0.0;
			// Add Lap for each user
			for(i=0;i<NodeNum;i++){
				// Sensitivity using each user's degree --> sen
				sen = deg_ns[i];
				// Add Lap(sen_tri/Eps2ndTrSt) --> trist2_num_u
				trist2_num_u[i] += stats::rlaplace(0.0, sen/Eps2ndTrSt, engine);
				sen_tri += sen;
			}
			sen_tri /= (double)NodeNum;
		}
	}

	// Ignore users with few edges
	if(NSType == 3){
		for(i=0;i<NodeNum;i++){
			if(deg_ns[i] * Mu < 1.0){
				trist2_num_u[i] = 0;
				continue;
			}
		}
	}

    // Empirical estimate --> tri_num_ns
    tri_num_ns = 0;
	for(i=0;i<NodeNum;i++) tri_num_ns += trist2_num_u[i];
	// Divide #triangles by (Mu - murho)
	tri_num_ns /= (Mu - murho);

	delete[] a_mat_ns;
	delete[] a_mat_del;
	free1D(tri_num_u);
	free1D(st2_num_u);
	free1D(trist2_num_u);
	free1D(deg_ns);
	free1D(red_sen);
}

// Calculate #triangles and #2-stars in the interactive local model (efficient algorithm II)
void CalcILocTriE2(map<int, int> *a_mat, int *deg, string outfile, double &tri_num_ns, double &sen_tri, double &eclip_sum, double &tclip_sum, int &eclip_num, int &tclip_num){
	map<int, int> *a_mat_ns;			// noisy adjacency matrix
	map<int, int> *a_mat_del;			// deleted adjacency matrix after projection
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	map<tuple<int, int>, int> a2_mat;
	map<tuple<int, int>, int>::iterator a2_itr;
	double murho;
	double *tri_num_u, *st2_num_u, *trist2_num_u;
	int max_deg;
	int sen_st2_ij;
	int del_num;
	int *rndperm;
	double rnd;
	int i, j, k, x;
	FILE *fp;

	double *deg_ns;
	double *red_sen;
	double tri_num_u_ij;
	double sen;
	int deg_ns_floor;

	// Initialization
	a_mat_ns = new map<int, int>[NodeNum];
	a_mat_del = new map<int, int>[NodeNum];
	malloc1D(&tri_num_u, NodeNum);
	malloc1D(&st2_num_u, NodeNum);
	malloc1D(&trist2_num_u, NodeNum);
	malloc1D(&deg_ns, NodeNum);
	malloc1D(&red_sen, NodeNum);

	// Parameters in asymmetric RR --> Mu (1 --> 1), murho (0 --> 1)
	murho = Mu / exp(Eps1st);

	// no noise, Lap (max degree), Lap (degree)
	if(NSType == -1 || NSType == 0 || NSType == 1){
		// Reduced sensitivity --> red_sen
		if(TClip != -1){
			for(i=0;i<NodeNum;i++) red_sen[i] = CalcRedSen((double)deg[i], 2);
		}

		// asymmetric RR --> a_mat_ns
		// Count #noisy triangles and #noisy 2-stars for each user --> tri_num_u, st2_num_u
		for(i=0;i<NodeNum;i++){
			for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
				j = aitr->first;
				tri_num_u_ij = 0;
				for (aitr2 = a_mat[i].begin(); aitr2 != a_mat[i].end(); aitr2++) {
					k = aitr2->first;
					// it must be j < k < i.
					if (j >= k || k >= i) continue;

					st2_num_u[i] += 1.0;

					// If a_mat_ns[k][i] does not exist
					if(a_mat_ns[k].count(i) == 0){
						// Flip 0/1 in a_mat[k][i] --> a_mat_ns[k][i]
						rnd = genrand_real2();
						// 0 --> 1 (flip)
						if(rnd < murho && a_mat[k].count(i) == 0){
							a_mat_ns[k][i] = 1;
						}
						// 1 --> 1 (not flip)
						else if(rnd < Mu && a_mat[k].count(i) == 1){
							a_mat_ns[k][i] = 1;
						}
						// 1 --> 0 (flip) or 0 --> 0 (not flip)
						else{
							a_mat_ns[k][i] = 0;
						}
					}

					// If a_mat_ns[j][k] does not exist
					if(a_mat_ns[j].count(k) == 0){
						// Flip 0/1 in a_mat[j][k] --> a_mat_ns[j][k]
						rnd = genrand_real2();
						// 0 --> 1 (flip)
						if(rnd < murho && a_mat[j].count(k) == 0){
							a_mat_ns[j][k] = 1;
						}
						// 1 --> 1 (not flip)
						else if(rnd < Mu && a_mat[j].count(k) == 1){
							a_mat_ns[j][k] = 1;
						}
						// 1 --> 0 (flip) or 0 --> 0 (not flip)
						else{
							a_mat_ns[j][k] = 0;
						}
					}

					if(a_mat_ns[k][i] == 1 && a_mat_ns[j][k] == 1){
						if(TClip != -1) tri_num_u_ij += 1.0;
						else tri_num_u[i] += 1.0;
					}
				}

				// Triangle clipping
				if(TClip != -1){
					if(tri_num_u_ij > red_sen[i]){
						tclip_sum += (tri_num_u_ij - red_sen[i]);
						tclip_num += 1;
						tri_num_u_ij = red_sen[i];
					}
					tri_num_u[i] += tri_num_u_ij;
				}
			}
		}

		// #triangles - Mu * murho * #2-stars --> trist2_num_u
		for(i=0;i<NodeNum;i++) trist2_num_u[i] = tri_num_u[i] - Mu * murho * st2_num_u[i];

		// no noise
		if(NSType == -1){
		}
		// Lap (max degree)
		else if(NSType == 0){
			// max(deg) --> max_deg
			max_deg = 0;
			for(i=0;i<NodeNum;i++){
				if(max_deg < deg[i]) max_deg = deg[i];
			}

			// Sensitivity using max degree --> sen_tri
			sen_tri = (double)max_deg;

			// Add Lap for each user
			for(i=0;i<NodeNum;i++){
				// Add Lap(sen_tri/Eps2ndTrSt) --> trist2_num_u
				trist2_num_u[i] += stats::rlaplace(0.0, sen_tri/Eps2ndTrSt, engine);
			}
		}
		// Lap (degree)
		else if(NSType == 1){
			// Reduced sensitivity
			if(TClip != -1){
				// Add Lap for each user
				sen_tri = 0.0;
				for(i=0;i<NodeNum;i++){
					sen_tri += red_sen[i];
					trist2_num_u[i] += stats::rlaplace(0.0, red_sen[i]/Eps2ndTrSt, engine);
				}
				sen_tri /= (double)NodeNum;
			}
			// Sensitivity using max degree
			else{
				// Average sensitivity using each user's degree --> sen_tri
				sen_tri = 0.0;
				// Add Lap for each user
				for(i=0;i<NodeNum;i++){
					// Sensitivity using each user's degree --> sen
					sen = (double)deg[i];
					// Add Lap(sen_tri/Eps2ndTrSt) --> trist2_num_u
					trist2_num_u[i] += stats::rlaplace(0.0, sen/Eps2ndTrSt, engine);
					sen_tri += sen;
				}
				sen_tri /= (double)NodeNum;
			}
		}
	}
	// Lap (noisy degree)
	else if(NSType == 2 || NSType == 3){
		// Noisy degree --> deg_ns
		for(i=0;i<NodeNum;i++){
			deg_ns[i] = (double)deg[i] + stats::rlaplace(0.0, 1.0/EpsNsDeg, engine);
		}
		// Add positive bias (EClip) --> deg_ns
		if(EClip != -1){
			for(i=0;i<NodeNum;i++) deg_ns[i] += EClip;
		}
		for(i=0;i<NodeNum;i++) deg_ns[i] = max(deg_ns[i], 0.0);

		// Reduced sensitivity --> red_sen
		if(TClip != -1){
			for(i=0;i<NodeNum;i++) red_sen[i] = CalcRedSen(deg_ns[i], 2);
		}

		// Graph projection (edge clipping) for each user
		for(i=0;i<NodeNum;i++){
			// If deg[i] exceeds deg_ns[i], then perform graph projection (edge clipping)
			if((double)deg[i] > deg_ns[i]){
				eclip_sum += ((double)deg[i] - floor(deg_ns[i]));
				eclip_num += 1;

				// Randomly generate 0, 1, ..., deg[i]-1 --> rndperm
				malloc1D(&rndperm, deg[i]);
				MakeRndPerm(rndperm, deg[i], deg[i]);

				// Randomly delete (deg[i] - floor(deg_ns[i])) edges from a_mat[i]
				deg_ns_floor = (int)floor(deg_ns[i]);
				x = 0;
				for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
					if(rndperm[x] >= deg_ns_floor){
						j = aitr->first;
						// Deleted edge --> a_mat_del[i][j]
						a_mat_del[i][j] = 1;
					}
					x++;
				}
				free1D(rndperm);
			}

			// Count #noisy triangles and #noisy 2-stars --> tri_num_u, st2_num_u
			tri_num_u[i] = 0.0;
			st2_num_u[i] = 0.0;
			for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
				j = aitr->first;
				// Continue if the edge is deleted
				if(a_mat_del[i].count(j) == 1) continue;
				tri_num_u_ij = 0;
				for (aitr2 = a_mat[i].begin(); aitr2 != a_mat[i].end(); aitr2++) {
					k = aitr2->first;
					// Continue if the edge is deleted
					if(a_mat_del[i].count(k) == 1) continue;
					if (j >= k || k >= i) continue;
					st2_num_u[i] += 1.0;

					// If a_mat_ns[k][i] does not exist
					if(a_mat_ns[k].count(i) == 0){
						// Flip 0/1 in a_mat[k][i] --> a_mat_ns[k][i]
						rnd = genrand_real2();
						// 0 --> 1 (flip)
						if(rnd < murho && a_mat[k].count(i) == 0){
							a_mat_ns[k][i] = 1;
						}
						// 1 --> 1 (not flip)
						else if(rnd < Mu && a_mat[k].count(i) == 1){
							a_mat_ns[k][i] = 1;
						}
						// 1 --> 0 (flip) or 0 --> 0 (not flip)
						else{
							a_mat_ns[k][i] = 0;
						}
					}

					// If a_mat_ns[j][k] does not exist
					if(a_mat_ns[j].count(k) == 0){
						// Flip 0/1 in a_mat[j][k] --> a_mat_ns[j][k]
						rnd = genrand_real2();
						// 0 --> 1 (flip)
						if(rnd < murho && a_mat[j].count(k) == 0){
							a_mat_ns[j][k] = 1;
						}
						// 1 --> 1 (not flip)
						else if(rnd < Mu && a_mat[j].count(k) == 1){
							a_mat_ns[j][k] = 1;
						}
						// 1 --> 0 (flip) or 0 --> 0 (not flip)
						else{
							a_mat_ns[j][k] = 0;
						}
					}

					if(a_mat_ns[k][i] == 1 && a_mat_ns[j][k] == 1){
						if(TClip != -1) tri_num_u_ij += 1.0;
						else tri_num_u[i] += 1.0;
					}
				}
				// Triangle clipping
				if(TClip != -1){
					if(tri_num_u_ij > red_sen[i]){
						tclip_sum += (tri_num_u_ij - red_sen[i]);
						tclip_num += 1;
						tri_num_u_ij = red_sen[i];
					}
					tri_num_u[i] += tri_num_u_ij;
				}
			}
			// #triangles - Mu * murho * #2-stars --> trist2_num_u
			trist2_num_u[i] = tri_num_u[i] - Mu * murho * st2_num_u[i];
		}

		// Reduced sensitivity
		if(TClip != -1){
			// Add Lap for each user
			sen_tri = 0.0;
			for(i=0;i<NodeNum;i++){
				sen_tri += red_sen[i];
				trist2_num_u[i] += stats::rlaplace(0.0, red_sen[i]/Eps2ndTrSt, engine);
			}
			sen_tri /= (double)NodeNum;
		}
		// Sensitivity using each user's degree
		else{
			// Average sensitivity using each user's degree --> sen_tri
			sen_tri = 0.0;
			// Add Lap for each user
			for(i=0;i<NodeNum;i++){
				// Sensitivity using each user's degree --> sen
				sen = deg_ns[i];
				// Add Lap(sen_tri/Eps2ndTrSt) --> trist2_num_u
				trist2_num_u[i] += stats::rlaplace(0.0, sen/Eps2ndTrSt, engine);
				sen_tri += sen;
			}
			sen_tri /= (double)NodeNum;
		}
	}

	// Ignore users with few edges
	if(NSType == 3){
		for(i=0;i<NodeNum;i++){
			if(deg_ns[i] * Mu * Mu < 1.0){
				trist2_num_u[i] = 0;
				continue;
			}
		}
	}

    // Empirical estimate --> tri_num_ns
    tri_num_ns = 0;
	for(i=0;i<NodeNum;i++) tri_num_ns += trist2_num_u[i];
	// Divide #triangles by Mu * (Mu - murho)
	tri_num_ns /= (Mu * (Mu - murho));

	delete[] a_mat_ns;
	delete[] a_mat_del;
	free1D(tri_num_u);
	free1D(st2_num_u);
	free1D(trist2_num_u);
	free1D(deg_ns);
	free1D(red_sen);
}

// Calculate #triangles and #2-stars in the interactive local model (efficient algorithm III)
void CalcILocTriE3(map<int, int> *a_mat, int *deg, string outfile, double &tri_num_ns, double &sen_tri, double &eclip_sum, double &tclip_sum, int &eclip_num, int &tclip_num){
	map<int, int> *a_mat_ns;			// noisy adjacency matrix
	map<int, int> *a_mat_del;			// deleted adjacency matrix after projection
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	map<tuple<int, int>, int> a2_mat;
	map<tuple<int, int>, int>::iterator a2_itr;
	double murho;
	double *tri_num_u, *st2_num_u, *trist2_num_u;
	int max_deg;
	int sen_st2_ij;
	int del_num;
	int *rndperm;
	double rnd;
	int i, j, k, x;
	FILE *fp;

	double *deg_ns;
	double *red_sen;
	double tri_num_u_ij;
	double sen;
	int deg_ns_floor;

	// Initialization
	a_mat_ns = new map<int, int>[NodeNum];
	a_mat_del = new map<int, int>[NodeNum];
	malloc1D(&tri_num_u, NodeNum);
	malloc1D(&st2_num_u, NodeNum);
	malloc1D(&trist2_num_u, NodeNum);
	malloc1D(&deg_ns, NodeNum);
	malloc1D(&red_sen, NodeNum);

	// Parameters in asymmetric RR --> Mu (1 --> 1), murho (0 --> 1)
	murho = Mu / exp(Eps1st);

	// no noise, Lap (max degree), Lap (degree)
	if(NSType == -1 || NSType == 0 || NSType == 1){
		// Reduced sensitivity --> red_sen
		if(TClip != -1){
			for(i=0;i<NodeNum;i++) red_sen[i] = CalcRedSen((double)deg[i], 3);
		}

		// asymmetric RR --> a_mat_ns
		// Count #noisy triangles and #noisy 2-stars for each user --> tri_num_u, st2_num_u
		for(i=0;i<NodeNum;i++){
			for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
				j = aitr->first;
				tri_num_u_ij = 0;
				for (aitr2 = a_mat[i].begin(); aitr2 != a_mat[i].end(); aitr2++) {
					k = aitr2->first;
					// it must be j < k < i.
					if (j >= k || k >= i) continue;

					st2_num_u[i] += 1.0;

					// If a_mat_ns[j][i] does not exist
					if(a_mat_ns[j].count(i) == 0){
						// Flip 0/1 in a_mat[j][i] --> a_mat_ns[j][i]
						rnd = genrand_real2();
						// 0 --> 1 (flip)
						if(rnd < murho && a_mat[j].count(i) == 0){
							a_mat_ns[j][i] = 1;
						}
						// 1 --> 1 (not flip)
						else if(rnd < Mu && a_mat[j].count(i) == 1){
							a_mat_ns[j][i] = 1;
						}
						// 1 --> 0 (flip) or 0 --> 0 (not flip)
						else{
							a_mat_ns[j][i] = 0;
						}
					}

					// If a_mat_ns[k][i] does not exist
					if(a_mat_ns[k].count(i) == 0){
						// Flip 0/1 in a_mat[k][i] --> a_mat_ns[k][i]
						rnd = genrand_real2();
						// 0 --> 1 (flip)
						if(rnd < murho && a_mat[k].count(i) == 0){
							a_mat_ns[k][i] = 1;
						}
						// 1 --> 1 (not flip)
						else if(rnd < Mu && a_mat[k].count(i) == 1){
							a_mat_ns[k][i] = 1;
						}
						// 1 --> 0 (flip) or 0 --> 0 (not flip)
						else{
							a_mat_ns[k][i] = 0;
						}
					}

					// If a_mat_ns[j][k] does not exist
					if(a_mat_ns[j].count(k) == 0){
						// Flip 0/1 in a_mat[j][k] --> a_mat_ns[j][k]
						rnd = genrand_real2();
						// 0 --> 1 (flip)
						if(rnd < murho && a_mat[j].count(k) == 0){
							a_mat_ns[j][k] = 1;
						}
						// 1 --> 1 (not flip)
						else if(rnd < Mu && a_mat[j].count(k) == 1){
							a_mat_ns[j][k] = 1;
						}
						// 1 --> 0 (flip) or 0 --> 0 (not flip)
						else{
							a_mat_ns[j][k] = 0;
						}
					}

					if(a_mat_ns[j][i] == 1 && a_mat_ns[k][i] == 1 && a_mat_ns[j][k] == 1){
						if(TClip != -1) tri_num_u_ij += 1.0;
						else tri_num_u[i] += 1.0;
					}
				}

				// Triangle clipping
				if(TClip != -1){
					if(tri_num_u_ij > red_sen[i]){
						tclip_sum += (tri_num_u_ij - red_sen[i]);
						tclip_num += 1;
						tri_num_u_ij = red_sen[i];
					}
					tri_num_u[i] += tri_num_u_ij;
				}
			}
		}

		// #triangles - Mu * Mu * murho * #2-stars --> trist2_num_u
		for(i=0;i<NodeNum;i++) trist2_num_u[i] = tri_num_u[i] - Mu * Mu * murho * st2_num_u[i];

		// no noise
		if(NSType == -1){
		}
		// Lap (max degree)
		else if(NSType == 0){
			// max(deg) --> max_deg
			max_deg = 0;
			for(i=0;i<NodeNum;i++){
				if(max_deg < deg[i]) max_deg = deg[i];
			}

			// Sensitivity using max degree --> sen_tri
			sen_tri = (double)max_deg;

			// Add Lap for each user
			for(i=0;i<NodeNum;i++){
				// Add Lap(sen_tri/Eps2ndTrSt) --> trist2_num_u
				trist2_num_u[i] += stats::rlaplace(0.0, sen_tri/Eps2ndTrSt, engine);
			}
		}
		// Lap (degree)
		else if(NSType == 1){
			// Reduced sensitivity
			if(TClip != -1){
				// Add Lap for each user
				sen_tri = 0.0;
				for(i=0;i<NodeNum;i++){
					sen_tri += red_sen[i];
					trist2_num_u[i] += stats::rlaplace(0.0, red_sen[i]/Eps2ndTrSt, engine);
				}
				sen_tri /= (double)NodeNum;
			}
			// Sensitivity using max degree
			else{
				// Average sensitivity using each user's degree --> sen_tri
				sen_tri = 0.0;
				// Add Lap for each user
				for(i=0;i<NodeNum;i++){
					// Sensitivity using each user's degree --> sen
					sen = (double)deg[i];
					// Add Lap(sen_tri/Eps2ndTrSt) --> trist2_num_u
					trist2_num_u[i] += stats::rlaplace(0.0, sen/Eps2ndTrSt, engine);
					sen_tri += sen;
				}
				sen_tri /= (double)NodeNum;
			}
		}
	}
	// Lap (noisy degree)
	else if(NSType == 2 || NSType == 3){
		// Noisy degree --> deg_ns
		for(i=0;i<NodeNum;i++){
			deg_ns[i] = (double)deg[i] + stats::rlaplace(0.0, 1.0/EpsNsDeg, engine);
		}
		// Add positive bias (EClip) --> deg_ns
		if(EClip != -1){
			for(i=0;i<NodeNum;i++) deg_ns[i] += EClip;
		}
		for(i=0;i<NodeNum;i++) deg_ns[i] = max(deg_ns[i], 0.0);

		// Reduced sensitivity --> red_sen
		if(TClip != -1){
			for(i=0;i<NodeNum;i++) red_sen[i] = CalcRedSen(deg_ns[i], 3);
		}

		// Graph projection (edge clipping) for each user
		for(i=0;i<NodeNum;i++){
			// If deg[i] exceeds deg_ns[i], then perform graph projection (edge clipping)
			if((double)deg[i] > deg_ns[i]){
				eclip_sum += ((double)deg[i] - floor(deg_ns[i]));
				eclip_num += 1;

				// Randomly generate 0, 1, ..., deg[i]-1 --> rndperm
				malloc1D(&rndperm, deg[i]);
				MakeRndPerm(rndperm, deg[i], deg[i]);

				// Randomly delete (deg[i] - floor(deg_ns[i])) edges from a_mat[i]
				deg_ns_floor = (int)floor(deg_ns[i]);
				x = 0;
				for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
					if(rndperm[x] >= deg_ns_floor){
						j = aitr->first;
						// Deleted edge --> a_mat_del[i][j]
						a_mat_del[i][j] = 1;
					}
					x++;
				}
				free1D(rndperm);
			}

			// Count #noisy triangles and #noisy 2-stars --> tri_num_u, st2_num_u
			tri_num_u[i] = 0.0;
			st2_num_u[i] = 0.0;
			for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
				j = aitr->first;
				// Continue if the edge is deleted
				if(a_mat_del[i].count(j) == 1) continue;
				tri_num_u_ij = 0;
				for (aitr2 = a_mat[i].begin(); aitr2 != a_mat[i].end(); aitr2++) {
					k = aitr2->first;
					// Continue if the edge is deleted
					if(a_mat_del[i].count(k) == 1) continue;
					if (j >= k || k >= i) continue;
					st2_num_u[i] += 1.0;

					// If a_mat_ns[j][i] does not exist
					if(a_mat_ns[j].count(i) == 0){
						// Flip 0/1 in a_mat[j][i] --> a_mat_ns[j][i]
						rnd = genrand_real2();
						// 0 --> 1 (flip)
						if(rnd < murho && a_mat[j].count(i) == 0){
							a_mat_ns[j][i] = 1;
						}
						// 1 --> 1 (not flip)
						else if(rnd < Mu && a_mat[j].count(i) == 1){
							a_mat_ns[j][i] = 1;
						}
						// 1 --> 0 (flip) or 0 --> 0 (not flip)
						else{
							a_mat_ns[j][i] = 0;
						}
					}

					// If a_mat_ns[k][i] does not exist
					if(a_mat_ns[k].count(i) == 0){
						// Flip 0/1 in a_mat[k][i] --> a_mat_ns[k][i]
						rnd = genrand_real2();
						// 0 --> 1 (flip)
						if(rnd < murho && a_mat[k].count(i) == 0){
							a_mat_ns[k][i] = 1;
						}
						// 1 --> 1 (not flip)
						else if(rnd < Mu && a_mat[k].count(i) == 1){
							a_mat_ns[k][i] = 1;
						}
						// 1 --> 0 (flip) or 0 --> 0 (not flip)
						else{
							a_mat_ns[k][i] = 0;
						}
					}

					// If a_mat_ns[j][k] does not exist
					if(a_mat_ns[j].count(k) == 0){
						// Flip 0/1 in a_mat[j][k] --> a_mat_ns[j][k]
						rnd = genrand_real2();
						// 0 --> 1 (flip)
						if(rnd < murho && a_mat[j].count(k) == 0){
							a_mat_ns[j][k] = 1;
						}
						// 1 --> 1 (not flip)
						else if(rnd < Mu && a_mat[j].count(k) == 1){
							a_mat_ns[j][k] = 1;
						}
						// 1 --> 0 (flip) or 0 --> 0 (not flip)
						else{
							a_mat_ns[j][k] = 0;
						}
					}

					if(a_mat_ns[j][i] == 1 && a_mat_ns[k][i] == 1 && a_mat_ns[j][k] == 1){
						if(TClip != -1) tri_num_u_ij += 1.0;
						else tri_num_u[i] += 1.0;
					}
				}
				// Triangle clipping
				if(TClip != -1){
					if(tri_num_u_ij > red_sen[i]){
						tclip_sum += (tri_num_u_ij - red_sen[i]);
						tclip_num += 1;
						tri_num_u_ij = red_sen[i];
					}
					tri_num_u[i] += tri_num_u_ij;
				}
			}
			// #triangles - Mu * Mu * murho * #2-stars --> trist2_num_u
			trist2_num_u[i] = tri_num_u[i] - Mu * Mu * murho * st2_num_u[i];
		}

		// Reduced sensitivity
		if(TClip != -1){
			// Add Lap for each user
			sen_tri = 0.0;
			for(i=0;i<NodeNum;i++){
				sen_tri += red_sen[i];
				trist2_num_u[i] += stats::rlaplace(0.0, red_sen[i]/Eps2ndTrSt, engine);
			}
			sen_tri /= (double)NodeNum;
		}
		// Sensitivity using each user's degree
		else{
			// Average sensitivity using each user's degree --> sen_tri
			sen_tri = 0.0;
			// Add Lap for each user
			for(i=0;i<NodeNum;i++){
				// Sensitivity using each user's degree --> sen
				sen = deg_ns[i];
				// Add Lap(sen_tri/Eps2ndTrSt) --> trist2_num_u
				trist2_num_u[i] += stats::rlaplace(0.0, sen/Eps2ndTrSt, engine);
				sen_tri += sen;
			}
			sen_tri /= (double)NodeNum;
		}
	}

	// Ignore users with few edges
	if(NSType == 3){
		for(i=0;i<NodeNum;i++){
			if(deg_ns[i] * Mu * Mu * Mu < 1.0){
				trist2_num_u[i] = 0;
				continue;
			}
		}
	}

    // Empirical estimate --> tri_num_ns
    tri_num_ns = 0;
	for(i=0;i<NodeNum;i++) tri_num_ns += trist2_num_u[i];
	// Divide #triangles by Mu * Mu * (Mu - murho)
	tri_num_ns /= (Mu * Mu * (Mu - murho));

	delete[] a_mat_ns;
	delete[] a_mat_del;
	free1D(tri_num_u);
	free1D(st2_num_u);
	free1D(trist2_num_u);
	free1D(deg_ns);
	free1D(red_sen);
}

// Calculate the clustering-coefficient
double CalcClstCoef(double tri_num_ns, double st2_num_ns){
	double clst_ns;

    if(tri_num_ns < 0) tri_num_ns = 0;
    if(st2_num_ns < 0) st2_num_ns = 0;
    if(st2_num_ns == 0) clst_ns = 1.0;
    else clst_ns = 3.0 * tri_num_ns / st2_num_ns;
    if(clst_ns > 1.0) clst_ns = 1.0;
    else if(clst_ns < 0.0) clst_ns = 0.0;

	return clst_ns;
}

int main(int argc, char *argv[])
{
	int all_node_num;
	int triplet_num;
	int **node_order;
	map<int, int> *a_mat;			// adjacency matrix
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	int *deg;									// degree
	int *deg_lower;								// degree in the lower-triangular part of a_max
	int max_deg;
	double Balloc_sum = 0.0;
	long long tot_edge_num;
	long long tri_num, st2_num, st3_num, ed2_num, ed1_num, non_num;
	map<int, int> pa2_mat;			// 2-path matrix
	long long pa3_num, cy4_num, pa2_pow2;
	double clst;
	double tri_num_ns, sen_tri;
	double st2_num_ns, st3_num_ns, sen_st2, sen_st3;
	double st2_num_ns_2, sen_st_2;
	double clst_ns;
	double tri_re_ns, tri_l2_ns;
	double tri_re_ns_avg, tri_l2_ns_avg;
	double st2_re_ns, st2_l2_ns;
	double st2_re_ns_avg, st2_l2_ns_avg;
	double st3_re_ns, st3_l2_ns;
	double st3_re_ns_avg, st3_l2_ns_avg;
	double clst_re_ns, clst_l2_ns;
	double clst_re_ns_avg, clst_l2_ns_avg;
	double eclip_sum, tclip_sum;
	int eclip_num, tclip_num;
	int itr;
	int i, j, k, x;
	string outdir;
	string outfile;
	char s[1025], *str;
	char str_1[] = "1";
	char *tok;
	FILE *fp;

	int fix_perm;

	// Initialization of Mersennne Twister
	unsigned long init[4] = { 0x123, 0x234, 0x345, 0x456 }, length = 4;
	init_by_array(init, length);
//	init_genrand(1);

	if (argc < 2) {
	    printf("Usage: %s [EdgeFile (in)] ([#nodes (default: -1)] [epsilon-mu/mu2/mu3 (default: 1-1)] [NSType (default: -1)] [tclip-eclip (default: -1)] [#itr(-1) (default: 1)] [alg (default: 0)] [Balloc (default: 1-1)])\n\n", argv[0]);
		printf("[EdgeFile]: Edge file\n");
		printf("[#nodes]: Number of nodes (-1: all)\n");
		printf("[epsilon-mu/mu2/mu3]: Parameters epsilon and mu/mu2/mu3 (I/II/III) (mu/mu2/mu3 = exp(Eps1st)/(exp(Eps1st)+1) when it is set to 1)\n");
		printf("[NSType]: Noise type (-1: no noise, 0: Lap (max degree), 1: Lap (degree + clip), 2: Lap (noisy degree + clip), 3: Lap (noisy degree + clip + ignore users))\n");
		printf("[tclip(-eclip)]: Triangle and edge clipping parameters (-1: no clipping) (alg=2-4; set eclip when NSType=2)\n");
		printf("[#itr(-1)]: Number of iterations (set #itr-1 to fix the permutation of nodes)\n");
		printf("[alg]: Algorithm (1: interactive local, 2: efficient interactive local I, 3: efficient interactive local II, 4: efficient interactive local III, 5: non-interactive local (RR w/ emp), 6: non-interactive local (RR w/o emp), 7: [Ye+, T-KDE (mean)], 8: [Ye+, T-KDE (median)], 9: [Ye+, T-KDE (most frequent degree)], 10: non-interactive local (ARR w/ emp))\n");
		printf("[Balloc]: Privacy budget allocation (alg=1-3): Eps1st-Eps2ndTrSt\n");
		return -1;
	}

	EdgeFile = argv[1];

	NodeNum = -1;
	if (argc >= 3) NodeNum = atoi(argv[2]);

	Eps = 1.0;
	Eps_s = "1";
	Mu = 1.0;
	Mu_s = "1";

	if (argc >= 4){
		if((EpsMu_s[0] = strtok(argv[3], "-")) == NULL){
			printf("Error: incorrect [epsilon-mu]\n");
			exit(-1);
		}
		if((EpsMu_s[1] = strtok(NULL, "-")) == NULL){
			printf("Error: incorrect [epsilon-mu]\n");
			exit(-1);
		}
		Eps = atof(EpsMu_s[0]);
		Mu = atof(EpsMu_s[1]);
		Eps_s = EpsMu_s[0];
		Mu_s = EpsMu_s[1];
	}

	NSType = -1;
	if (argc >= 5) NSType = atoi(argv[4]);

	// Triangle and edge clipping parameters
	TClip = -1;
	EClip = -1;
	Clip_s = "-1";
	if (argc >= 6){
		Clip_s = argv[5];
		if(strcmp(argv[5], "-1") != 0){
			if(argv[5][0] !=  '-'){
				if((tok  = strtok(argv[5], "-")) == NULL){
					printf("Error: incorrect [tclip(-eclip)]\n");
					exit(-1);
				}
				TClip = atof(tok);
				if((tok  = strtok(NULL, "-")) != NULL) EClip = atof(tok);
			}
			else{
				if((tok  = strtok(argv[5], "-")) == NULL){
					printf("Error: incorrect [tclip(-eclip)]\n");
					exit(-1);
				}
				if((tok  = strtok(NULL, "-")) != NULL) EClip = atof(tok);
			}
		}

	}

	ItrNum = 1;
	fix_perm = 0;
	if (argc >= 7){
		tok  = strtok(argv[6], "-");
		ItrNum = atoi(tok);
		if((tok  = strtok(NULL, "-")) != NULL){
			if (strcmp(tok, "1") != 0){
				printf("Error: incorrect [#itr(-1)]\n");
				exit(-1);
			}
			else fix_perm = 1;
		}
	}

	Alg = 0;
	if (argc >= 8) Alg = atoi(argv[7]);
	if (Alg <= 0 || Alg > 10){
		printf("Error: incorrect [Alg]\n");
		exit(-1);
	}

	for(i=0;i<2;i++){
		Balloc[i] = 1.0;
		Balloc_s[i] = str_1;
	}
	if (argc >= 9){
		if((Balloc_s[0] = strtok(argv[8], "-")) == NULL){
			printf("Error: incorrect [Balloc]\n");
			exit(-1);
		}
		Balloc[0] = atof(Balloc_s[0]);
		if((Balloc_s[1] = strtok(NULL, "-")) == NULL){
			printf("Error: incorrect [Balloc]\n");
			exit(-1);
		}
		Balloc[1] = atof(Balloc_s[1]);
	}

	// Privacy budget allocation
	for(i=0;i<2;i++) Balloc_sum += Balloc[i];
	if(NSType == -1 || NSType == 0 || NSType == 1){
		EpsNsDeg = 0.0;				// Epsilon for calculating the noisy degree
		EpsAllbutNsDeg = Eps;
		Eps1st = Eps * Balloc[0] / Balloc_sum;
		Eps2ndTrSt = Eps * Balloc[1] / Balloc_sum;
	}
	else if(NSType == 2 || NSType == 3){
		EpsNsDeg = Eps / 10;		// Epsilon for calculating the noisy degree
		EpsAllbutNsDeg = Eps - EpsNsDeg;
		Eps1st = EpsAllbutNsDeg * Balloc[0] / Balloc_sum;
		Eps2ndTrSt = EpsAllbutNsDeg * Balloc[1]  / Balloc_sum;
	}

	// mu/mu^2/mu^3 = exp(Eps1st)/(exp(Eps1st)+1) when it is set to 1
	if(Mu == 1.0){
		Mu = exp(Eps1st) / (exp(Eps1st) + 1);
	}
	// When Alg == 3 (efficient interactive local II), calculate mu from mu^2
	if (Alg == 3){
		Mu = pow(Mu, 1.0/2.0);
	}
	// When Alg == 4 (efficient interactive local III), calculate mu from mu^3 
	if (Alg == 4){
		Mu = pow(Mu, 1.0/3.0);
	}

	// Total number of nodes --> all_node_num
	fp = FileOpen(EdgeFile, "r");
	for(i=0;i<2;i++) fgets(s, 1024, fp);
	all_node_num = atoi(s);
	fclose(fp);

	// malloc
	malloc2D(&node_order, ItrNum, all_node_num);

	// Use all nodes
	if (NodeNum == -1){
		NodeNum = all_node_num;
		for(j=0;j<NodeNum;j++) node_order[0][j] = j;
	}
	// Randomly generate the order of nodes --> node_order
	else{
		i = EdgeFile.find_last_of("/");
		outdir = EdgeFile.substr(0, i+1);
		outfile = outdir + "node-order_itr" + to_string(ItrNum) + ".csv";
		if(checkFileExistence(outfile)){
			fp = FileOpen(outfile, "r");
			for(j=0;j<all_node_num;j++){
				fgets(s, 1024, fp);
				strtok(s, ",");
				for(i=0;i<ItrNum;i++){
					node_order[i][j] = atoi(strtok(NULL, ","));
				}
			}
			fclose(fp);
		}
		else{
			for(i=0;i<ItrNum;i++){
				MakeRndPerm(node_order[i], all_node_num, all_node_num);
			}
			fp = FileOpen(outfile, "w");
			for(j=0;j<all_node_num;j++){
				fprintf(fp, "%d,", j);
				for(i=0;i<ItrNum;i++) fprintf(fp, "%d,", node_order[i][j]);
				fprintf(fp, "\n");
			}
			fclose(fp);
		}

		// Use only the first permutation
		if (fix_perm){
			for(j=0;j<all_node_num;j++){
				for(i=1;i<ItrNum;i++) node_order[i][j] = node_order[0][j];
			}
		}
	}

	// #triplet --> triplet_num
	triplet_num = NodeNum*(NodeNum-1)*(NodeNum-2)/6;

	// Initialization
	malloc1D(&deg, NodeNum);
	malloc1D(&deg_lower, NodeNum);
	tri_re_ns_avg = tri_l2_ns_avg = 0.0;
	st2_re_ns_avg = st2_l2_ns_avg = 0.0;
	st3_re_ns_avg = st3_l2_ns_avg = 0.0;
	clst_re_ns_avg = clst_l2_ns_avg = 0.0;

	// Output the header
	i = EdgeFile.find_last_of("/");
	outdir = EdgeFile.substr(0, i+1);
	for(i=0;i<3;i++){
		if(fix_perm) outfile = outdir + "res_n" + to_string(NodeNum) + "_alg" + to_string(Alg) + "_eps" + Eps_s + "-" + Mu_s + "_ns" + to_string(NSType) + "_cl" + Clip_s + "_ba" + Balloc_s[0] + "-" + Balloc_s[1] + "_itr" + to_string(ItrNum) + "-1.csv";
		else outfile = outdir + "res_n" + to_string(NodeNum) + "_alg" + to_string(Alg) + "_eps" + Eps_s + "-" + Mu_s + "_ns" + to_string(NSType) + "_cl" + Clip_s + "_ba" + Balloc_s[0] + "-" + Balloc_s[1] + "_itr" + to_string(ItrNum) + ".csv";
		fp = FileOpen(outfile, "w");
		fprintf(fp, "#tri(true),#tri(est),#tri(rel-err),#tri(l2-loss),#2st(true),#2st(est),#2st(rel-err),#2st(l2-loss),#3st(true),#3st(est),#3st(rel-err),#3st(l2-loss),clst(true),clst(est),clst(rel-err),clst(l2-loss),sen_tri,sen_2st,sen_3st,eclip_sum,tclip_sum,#eclip,#tclip,max_deg\n");
		fclose(fp);
	}

	// For each iteration
	for(itr=0;itr<ItrNum;itr++){
		// Read edges for each iteration when NodeNum < all_node_num
		if(NodeNum < all_node_num || itr == 0){
			// Initialization
			a_mat = new map<int, int>[NodeNum];

			// Read edges from the edge file --> a_mat
			ReadEdges(a_mat, node_order[itr]);

			// Degree --> deg
			for(i=0;i<NodeNum;i++) deg[i] = 0;
			for(i=0;i<NodeNum;i++){
				for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) deg[i] += 1;
			}

			// Degree --> deg_lower
			for(i=0;i<NodeNum;i++) deg_lower[i] = 0;
			for(i=0;i<NodeNum;i++){
				for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++){
					if(aitr->first < i) deg_lower[i] += 1;
				}
			}

			// max(deg) --> max_deg
			max_deg = 0;
			for(i=0;i<NodeNum;i++){
				if(max_deg < deg[i]) max_deg = deg[i];
			}

			// Total number of edges --> tot_edge_num
			tot_edge_num = 0;
			for(i=0;i<NodeNum;i++) tot_edge_num += (long long)deg[i];
			tot_edge_num /= 2;

			// #triangles --> tri_num
			tri_num = 0;
			for(i=0;i<NodeNum;i++){
				for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
					j = aitr->first;
					if (i >= j) continue;
					for (aitr2 = a_mat[i].begin(); aitr2 != a_mat[i].end(); aitr2++) {
						k = aitr2->first;
						if (j >= k) continue;
						if(a_mat[j].count(k) > 0) tri_num++;
					}
				}
			}

			// #2-stars, #3-stars --> st2_num, st3_num
			st2_num = st3_num = 0;
			for(i=0;i<NodeNum;i++){
				st2_num += ((long long)deg[i] * ((long long)deg[i]-1)) / 2;
				st3_num += ((long long)deg[i] * ((long long)deg[i]-1) * ((long long)deg[i]-2)) / 6;
			}

			/*
			// #2-path from i to k, #3-paths --> pa2_mat[i][k], pa3_num
			pa3_num = 0;
			cy4_num = 0;
			pa2_pow2 = 0;
			// i--j--k (2-path)
			for(i=0;i<NodeNum;i++){
				for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
					j = aitr->first;
					if (i==j) continue;
					for (aitr2 = a_mat[j].begin(); aitr2 != a_mat[j].end(); aitr2++) {
						k = aitr2->first;
						if (i==k || j==k) continue;
						// Update 2-path matrix --> pa2_mat
						if(pa2_mat.count(k) == 0) pa2_mat[k] = 1;
						else pa2_mat[k] += 1;
						// Add #3-paths (i--j--k--l) --> pa3_num
						if(a_mat[i].count(k) == 0){
							pa3_num += (long long)(deg[k] - 1); // Subtract 1 (i--j--k--j)
						}
						else{
							pa3_num += (long long)(deg[k] - 2); // Subtract 2 (i--j--k--i and i--j--k--j)
						}
					}
				}

				for (aitr = pa2_mat.begin(); aitr != pa2_mat.end(); aitr++) {
					j = aitr->first;
					if (i>=j) continue;
					pa2_pow2 += pa2_mat[j] * pa2_mat[j];
				}
				pa2_mat.clear();
			}
			// Divide pa3_num by 2 (because i--j--k--l and l--k--j--i are the same) --> pa3_num
			pa3_num /= 2;

			cy4_num = (pa2_pow2 - st2_num) / 2;

			if(fix_perm) outfile = outdir + "res_n" + to_string(NodeNum) + "_alg" + to_string(Alg) + "_eps" + Eps_s + "-" + Mu_s + "_ns" + to_string(NSType) + "_cl" + Clip_s + "_ba" + Balloc_s[0] + "-" + Balloc_s[1] + "_itr" + to_string(ItrNum) + "-1_4cycles.csv";
			else outfile = outdir + "res_n" + to_string(NodeNum) + "_alg" + to_string(Alg) + "_eps" + Eps_s + "-" + Mu_s + "_ns" + to_string(NSType) + "_cl" + Clip_s + "_ba" + Balloc_s[0] + "-" + Balloc_s[1] + "_itr" + to_string(ItrNum) + "_4cycles.csv";
			fp = FileOpen(outfile, "w");
			fprintf(fp, "#2-stars,#3-stars,#3-paths,#4-cycles\n");
			fprintf(fp, "%lld,%lld,%lld,%lld\n", st2_num, st3_num, pa3_num, cy4_num);
			fclose(fp);
			exit(0);
			*/

			// clustering coefficient --> clst
			if(st2_num != 0) clst = 3.0 * (double)tri_num / (double)st2_num;
			else clst = 1.0;

			// #2-edges --> ed2_num
			ed2_num = st2_num - 3*tri_num;
			// #1-edge --> ed1_num
			ed1_num = (long long)tot_edge_num*(NodeNum-2) - 2*ed2_num - 3*tri_num;
			// #none --> non_num
			non_num = (long long)NodeNum*(NodeNum-1)*(NodeNum-2)/6 - tri_num - ed2_num - ed1_num;
		}

		/************************ Calculate sub-graph counts ************************/
		eclip_sum = tclip_sum = 0.0;
		eclip_num = tclip_num = 0;
		// Interactive (2-rounds) local
		else if (Alg == 1){
			// Calculate #2-stars and #3-stars
        	CalcNLocSt(st2_num, st3_num, deg, outfile, st2_num_ns, st3_num_ns, sen_st2, sen_st3);
        	// Calculate #triangles
			CalcILocTri(a_mat, deg_lower, outfile, tri_num_ns, sen_tri);
		}
		// Efficient interactive (2-rounds) local I (asymmetric RR)
		else if (Alg == 2){
			// Calculate #2-stars and #3-stars
        	CalcNLocSt(st2_num, st3_num, deg, outfile, st2_num_ns, st3_num_ns, sen_st2, sen_st3);
        	// Calculate #triangles (efficient algorithm I)
			CalcILocTriE1(a_mat, deg_lower, outfile, tri_num_ns, sen_tri, eclip_sum, tclip_sum, eclip_num, tclip_num);
		}
		// Efficient interactive (2-rounds) local II (asymmetric RR + edge selection)
		else if (Alg == 3){
			// Calculate #2-stars and #3-stars
        	CalcNLocSt(st2_num, st3_num, deg, outfile, st2_num_ns, st3_num_ns, sen_st2, sen_st3);
        	// Calculate #triangles (efficient algorithm II)
			CalcILocTriE2(a_mat, deg_lower, outfile, tri_num_ns, sen_tri, eclip_sum, tclip_sum, eclip_num, tclip_num);
		}
		// Efficient interactive (2-rounds) local III (asymmetric RR + edge selection)
		else if (Alg == 4){
			// Calculate #2-stars and #3-stars
        	CalcNLocSt(st2_num, st3_num, deg, outfile, st2_num_ns, st3_num_ns, sen_st2, sen_st3);
        	// Calculate #triangles (efficient algorithm III)
			CalcILocTriE3(a_mat, deg_lower, outfile, tri_num_ns, sen_tri, eclip_sum, tclip_sum, eclip_num, tclip_num);
		}
		// Non-interactive (1-round) local (RR w/ emp)
		else if (Alg == 5){
			// Calculate #2-stars and #3-stars
        	CalcNLocSt(st2_num, st3_num, deg, outfile, st2_num_ns, st3_num_ns, sen_st2, sen_st3);
        	// Calculate #triangles
			if(NodeNum <= 10000) CalcNLocTri(a_mat, outfile, tri_num_ns, 1);
			else tri_num_ns = 0.0;
			sen_tri = 0.0;
		}
		// Non-interactive (1-round) local (RR w/o emp)
		else if (Alg == 6){
			// Calculate #2-stars and #3-stars
        	CalcNLocSt(st2_num, st3_num, deg, outfile, st2_num_ns, st3_num_ns, sen_st2, sen_st3);
        	// Calculate #triangles
			if(NodeNum <= 10000) CalcNLocTri(a_mat, outfile, tri_num_ns, 0);
			else tri_num_ns = 0.0;
			sen_tri = 0.0;
		}
		// [Ye+, T-KDE]
		else if (Alg >= 7 && Alg <= 9){
			// Calculate #2-stars ad #3-stars by [Ye+, T-KDE]
        	CalcTKDESt(deg, outfile, st2_num_ns, st3_num_ns, sen_st2, sen_st3);
        	// Calculate #triangles by [Ye+, T-KDE]
			if(NodeNum <= 10000) CalcTKDETri(a_mat, deg, outfile, tri_num_ns, st2_num_ns_2, sen_tri, sen_st_2, Alg);
			else tri_num_ns = st2_num_ns_2 = 0.0;
		}
		// Non-interactive (1-round) local (ARR w/ emp)
		else if (Alg == 10){
			// Calculate #2-stars and #3-stars
        	CalcNLocSt(st2_num, st3_num, deg, outfile, st2_num_ns, st3_num_ns, sen_st2, sen_st3);
        	// Calculate #triangles
			CalcNLocTriARR(a_mat, outfile, tri_num_ns, 1);
			sen_tri = 0.0;
		}

		/******************** Calculate the cluster coefficient *********************/
		clst_ns = CalcClstCoef(tri_num_ns, st2_num_ns);

		/**************************** Evaluate the loss *****************************/
		// relative error --> tri_re_ns
		tri_re_ns = fabs(tri_num_ns - (double)tri_num) / max((double)tri_num, 0.001 * NodeNum);
		tri_re_ns_avg += tri_re_ns;
		// l2_loss --> tri_l2_ns
		tri_l2_ns = (tri_num_ns - (double)tri_num)*(tri_num_ns - (double)tri_num);
		tri_l2_ns_avg += tri_l2_ns;

		// relative error --> st2_re_ns
		st2_re_ns = fabs(st2_num_ns - (double)st2_num) / max((double)st2_num, 0.001 * NodeNum);
		st2_re_ns_avg += st2_re_ns;
		// l2_loss --> st2_l2_ns
		st2_l2_ns = (st2_num_ns - (double)st2_num)*(st2_num_ns - (double)st2_num);
		st2_l2_ns_avg += st2_l2_ns;

		// relative error --> st3_re_ns
		st3_re_ns = fabs(st3_num_ns - (double)st3_num) / max((double)st3_num, 0.001 * NodeNum);
		st3_re_ns_avg+= st3_re_ns;
		// l2_loss --> st3_l2_ns
		st3_l2_ns = (st3_num_ns - (double)st3_num)*(st3_num_ns - (double)st3_num);
		st3_l2_ns_avg += st3_l2_ns;

		// relative error --> clst_re_ns
		clst_re_ns = fabs(clst_ns - clst) / (double)clst;
		clst_re_ns_avg += clst_re_ns;
		// l2_loss --> clst_l2_ns
		clst_l2_ns = (clst_ns - clst)*(clst_ns - clst);
		clst_l2_ns_avg += clst_l2_ns;

		/**************************** Output the results ****************************/
		fp = FileOpen(outfile, "a");
		fprintf(fp, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%d,%d\n", 
		(double)tri_num, tri_num_ns, tri_re_ns, tri_l2_ns, 
		(double)st2_num, st2_num_ns, st2_re_ns, st2_l2_ns, 
		(double)st3_num, st3_num_ns, st3_re_ns, st3_l2_ns, 
		clst, clst_ns, clst_re_ns, clst_l2_ns,
		sen_tri, sen_st2, sen_st3, 
		eclip_sum, tclip_sum, eclip_num, tclip_num, 
		max_deg);
		fclose(fp);

		if(NodeNum < all_node_num || itr == ItrNum - 1){
			delete[] a_mat;
		}
	}

	/************************* Output the results (AVG) *************************/
	tri_re_ns_avg /= (double)ItrNum;
	tri_l2_ns_avg /= (double)ItrNum;
	st2_re_ns_avg /= (double)ItrNum;
	st2_l2_ns_avg /= (double)ItrNum;
	st3_re_ns_avg /= (double)ItrNum;
	st3_l2_ns_avg /= (double)ItrNum;
	clst_re_ns_avg /= (double)ItrNum;
	clst_l2_ns_avg /= (double)ItrNum;

	fp = FileOpen(outfile, "a");
	fprintf(fp, "function,AVG(rel-err),AVG(l2-loss)\n");
	fprintf(fp, "Triangles,%e,%e\n", tri_re_ns_avg, tri_l2_ns_avg);
	fprintf(fp, "2-stars,%e,%e\n", st2_re_ns_avg, st2_l2_ns_avg);
	fprintf(fp, "3-stars,%e,%e\n", st3_re_ns_avg, st3_l2_ns_avg);
	fprintf(fp, "Clst,%e,%e\n", clst_re_ns_avg, clst_l2_ns_avg);
	fclose(fp);

	// free
	free2D(node_order, ItrNum);
	free1D(deg);
	free1D(deg_lower);

	return 0;
}
