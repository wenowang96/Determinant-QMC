#pragma once

#include <stdint.h>

struct params {
	int bps;
	int N, L;
	int *map_i, *map_ij;
	int *bonds, *bond2s, *bondws, *map_bs, *map_bb,*map_b2b,*map_bb2, *map_b2b2;
//	double *K, *U;
//	double dt;

	int n_matmul, n_delay;
	int n_sweep_warm, n_sweep_meas;
	int period_eqlt, period_uneqlt;
	int meas_bond_corr, meas_2bond_corr, meas_energy_corr, meas_nematic_corr, meas_correction;

	int num_i, num_ij;
	int num_b, num_b2, num_bs, num_bb, num_b2b, num_bb2, num_b2b2;
	int *degen_i, *degen_ij, *degen_bs, *degen_bb, *degen_b2b, *degen_bb2, *degen_b2b2;
	double *exp_K, *inv_exp_K;
	double *exp_halfK, *inv_exp_halfK;
	double *exp_lambda, *del;
	int F, n_sweep;
};

struct state {
	uint64_t rng[17];
	int sweep;
	int *hs;
};

struct meas_eqlt {
	int n_sample;
	double sign;

	double *density;
	double *double_occ;

	double *g00;
	double *nn;
	double *xx;
	double *zz;
	double *pair_sw;
	double *kk, *kv, *kn, *vv, *vn;
};

struct meas_uneqlt {
	int n_sample;
	double sign;

	double *gt0;
	double *nn;
	double *xx;
	double *zz;
	double *pair_sw;
	double *pair_bb;
	double *jj, *jsjs;
	double *kk, *ksks;
	double *pair_b2b2;
	double *j2j2, *js2js2;
	double *k2k2, *ks2ks2;
	double *kv, *kn, *vv, *vn;
	double *nem_nnnn, *nem_ssss;
        double *LLj1LLj1, *LLj1LLj2, *LLj2LLj2, *LLj1j2,*LLj2j2;
};

struct sim_data {
	struct params p;
	struct state s;
	struct meas_eqlt m_eq;
	struct meas_uneqlt m_ue;
};

int sim_data_read_alloc(struct sim_data *sim, const char *file);

int sim_data_save(const struct sim_data *sim, const char *file);

void sim_data_free(const struct sim_data *sim);
