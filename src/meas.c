#include "meas.h"
#include "data.h"

// number of bonds kept for nematic correlators. 2 by default
#define NEM_BONDS 2

void measure_eqlt(const struct params *const restrict p, const int sign,
		const double *const restrict gu,
		const double *const restrict gd,
		struct meas_eqlt *const restrict m)
{
	m->n_sample++;
	m->sign += sign;
	const int N = p->N, num_i = p->num_i, num_ij = p->num_ij;
	const int num_b = p->num_b, num_bs = p->num_bs, num_bb = p->num_bb;
	const int meas_energy_corr = p->meas_energy_corr;

	// 1 site measurements
	for (int i = 0; i < N; i++) {
		const int r = p->map_i[i];
		const double pre = (double)sign / p->degen_i[r];
		const double guii = gu[i + i*N], gdii = gd[i + i*N];
		m->density[r] += pre*(2. - guii - gdii);
		m->double_occ[r] += pre*(1. - guii)*(1. - gdii);
	}

	// 2 site measurements
	for (int j = 0; j < N; j++)
	for (int i = 0; i < N; i++) {
		const int delta = (i == j);
		const int r = p->map_ij[i + j*N];
		const double pre = (double)sign / p->degen_ij[r];
		const double guii = gu[i + i*N], gdii = gd[i + i*N];
		const double guij = gu[i + j*N], gdij = gd[i + j*N];
		const double guji = gu[j + i*N], gdji = gd[j + i*N];
		const double gujj = gu[j + j*N], gdjj = gd[j + j*N];
		m->g00[r] += 0.5*pre*(guij + gdij);
		const double x = delta*(guii + gdii) - (guji*guij + gdji*gdij);
		m->nn[r] += pre*((2. - guii - gdii)*(2. - gujj - gdjj) + x);
		m->xx[r] += 0.25*pre*(delta*(guii + gdii) - (guji*gdij + gdji*guij));
		m->zz[r] += 0.25*pre*((gdii - guii)*(gdjj - gujj) + x);
		m->pair_sw[r] += pre*guij*gdij;
		if (meas_energy_corr) {
			const double nuinuj = (1. - guii)*(1. - gujj) + (delta - guji)*guij;
			const double ndindj = (1. - gdii)*(1. - gdjj) + (delta - gdji)*gdij;
			m->vv[r] += pre*nuinuj*ndindj;
			m->vn[r] += pre*(nuinuj*(1. - gdii) + (1. - guii)*ndindj);
		}
	}

	if (!meas_energy_corr)
		return;

	// 1 bond 1 site measurements
	for (int j = 0; j < N; j++)
	for (int b = 0; b < num_b; b++) {
		const int i0 = p->bonds[b];
		const int i1 = p->bonds[b + num_b];
		const int bs = p->map_bs[b + num_b*j];
		const double pre = (double)sign / p->degen_bs[bs];
		const int delta_i0i1 = 0;
		const int delta_i0j = (i0 == j);
		const int delta_i1j = (i1 == j);
		const double gui0j = gu[i0 + N*j];
		const double guji0 = gu[j + N*i0];
		const double gdi0j = gd[i0 + N*j];
		const double gdji0 = gd[j + N*i0];
		const double gui1j = gu[i1 + N*j];
		const double guji1 = gu[j + N*i1];
		const double gdi1j = gd[i1 + N*j];
		const double gdji1 = gd[j + N*i1];
		const double gui0i1 = gu[i0 + N*i1];
		const double gui1i0 = gu[i1 + N*i0];
		const double gdi0i1 = gd[i0 + N*i1];
		const double gdi1i0 = gd[i1 + N*i0];
		const double gujj = gu[j + N*j];
		const double gdjj = gd[j + N*j];

		const double ku = 2.*delta_i0i1 - gui0i1 - gui1i0;
		const double kd = 2.*delta_i0i1 - gdi0i1 - gdi1i0;
		const double xu = (delta_i0j - guji0)*gui1j + (delta_i1j - guji1)*gui0j;
		const double xd = (delta_i0j - gdji0)*gdi1j + (delta_i1j - gdji1)*gdi0j;
		m->kv[bs] += pre*((ku*(1. - gujj) + xu)*(1. - gdjj)
		                + (kd*(1. - gdjj) + xd)*(1. - gujj));
		m->kn[bs] += pre*((ku + kd)*(2. - gujj - gdjj) + xu + xd);
	}

	// 2 bond measurements
	for (int c = 0; c < num_b; c++) {
		const int j0 = p->bonds[c];
		const int j1 = p->bonds[c + num_b];
	for (int b = 0; b < num_b; b++) {
		const int i0 = p->bonds[b];
		const int i1 = p->bonds[b + num_b];
		const int bb = p->map_bb[b + c*num_b];
		const double pre = (double)sign / p->degen_bb[bb];
		const int delta_i0j0 = (i0 == j0);
		const int delta_i1j0 = (i1 == j0);
		const int delta_i0j1 = (i0 == j1);
		const int delta_i1j1 = (i1 == j1);
		const double gui1i0 = gu[i1 + i0*N];
		const double gui0i1 = gu[i0 + i1*N];
		const double gui0j0 = gu[i0 + j0*N];
		const double gui1j0 = gu[i1 + j0*N];
		const double gui0j1 = gu[i0 + j1*N];
		const double gui1j1 = gu[i1 + j1*N];
		const double guj0i0 = gu[j0 + i0*N];
		const double guj1i0 = gu[j1 + i0*N];
		const double guj0i1 = gu[j0 + i1*N];
		const double guj1i1 = gu[j1 + i1*N];
		const double guj1j0 = gu[j1 + j0*N];
		const double guj0j1 = gu[j0 + j1*N];
		const double gdi1i0 = gd[i1 + i0*N];
		const double gdi0i1 = gd[i0 + i1*N];
		const double gdi0j0 = gd[i0 + j0*N];
		const double gdi1j0 = gd[i1 + j0*N];
		const double gdi0j1 = gd[i0 + j1*N];
		const double gdi1j1 = gd[i1 + j1*N];
		const double gdj0i0 = gd[j0 + i0*N];
		const double gdj1i0 = gd[j1 + i0*N];
		const double gdj0i1 = gd[j0 + i1*N];
		const double gdj1i1 = gd[j1 + i1*N];
		const double gdj1j0 = gd[j1 + j0*N];
		const double gdj0j1 = gd[j0 + j1*N];
		const double x = ((delta_i0j1 - guj1i0)*gui1j0 + (delta_i1j0 - guj0i1)*gui0j1
		                + (delta_i0j1 - gdj1i0)*gdi1j0 + (delta_i1j0 - gdj0i1)*gdi0j1);
		const double y = ((delta_i0j0 - guj0i0)*gui1j1 + (delta_i1j1 - guj1i1)*gui0j0
		                + (delta_i0j0 - gdj0i0)*gdi1j1 + (delta_i1j1 - gdj1i1)*gdi0j0);
		m->kk[bb] += pre*((gui0i1 + gui1i0 + gdi0i1 + gdi1i0)*(guj0j1 + guj1j0 + gdj0j1 + gdj1j0) + x + y);
	}
	}
}

void measure_uneqlt(const struct params *const restrict p, const int sign,
		const double *const Gu0t,
		const double *const Gutt,
		const double *const Gut0,
		const double *const Gd0t,
		const double *const Gdtt,
		const double *const Gdt0,
		struct meas_uneqlt *const restrict m)
{
	m->n_sample++;
	m->sign += sign;
	const int N = p->N, L = p->L, num_i = p->num_i, num_ij = p->num_ij;
	const int num_b = p->num_b, num_bs = p->num_bs, num_bb = p->num_bb;
        const int num_b2 = p->num_b2, num_b2b2 = p->num_b2b2, num_b2b = p->num_b2b, num_bb2 = p->num_bb2;
	const int meas_bond_corr = p->meas_bond_corr;
	const int meas_2bond_corr = p->meas_2bond_corr;
	const int meas_energy_corr = p->meas_energy_corr;
	const int meas_nematic_corr = p->meas_nematic_corr;
        const int meas_thermal = p->meas_thermal;

	const double *const restrict Gu00 = Gutt;
	const double *const restrict Gd00 = Gdtt;

	// 2 site measurements
	#pragma omp parallel for
	for (int t = 0; t < L; t++) {
		const int delta_t = (t == 0);
		const double *const restrict Gu0t_t = Gu0t + N*N*t;
		const double *const restrict Gutt_t = Gutt + N*N*t;
		const double *const restrict Gut0_t = Gut0 + N*N*t;
		const double *const restrict Gd0t_t = Gd0t + N*N*t;
		const double *const restrict Gdtt_t = Gdtt + N*N*t;
		const double *const restrict Gdt0_t = Gdt0 + N*N*t;
	for (int j = 0; j < N; j++)
	for (int i = 0; i < N; i++) {
		const int r = p->map_ij[i + j*N];
		const int delta_tij = delta_t * (i == j);
		const double pre = (double)sign / p->degen_ij[r];
		const double guii = Gutt_t[i + N*i];
		const double guij = Gut0_t[i + N*j];
		const double guji = Gu0t_t[j + N*i];
		const double gujj = Gu00[j + N*j];
		const double gdii = Gdtt_t[i + N*i];
		const double gdij = Gdt0_t[i + N*j];
		const double gdji = Gd0t_t[j + N*i];
		const double gdjj = Gd00[j + N*j];
		m->gt0[r + num_ij*t] += 0.5*pre*(guij + gdij);
		const double x = delta_tij*(guii + gdii) - (guji*guij + gdji*gdij);
		m->nn[r + num_ij*t] += pre*((2. - guii - gdii)*(2. - gujj - gdjj) + x);
		m->xx[r + num_ij*t] += 0.25*pre*(delta_tij*(guii + gdii) - (guji*gdij + gdji*guij));
		m->zz[r + num_ij*t] += 0.25*pre*((gdii - guii)*(gdjj - gujj) + x);
		m->pair_sw[r + num_ij*t] += pre*guij*gdij;
		if (meas_energy_corr) {
			const double nuinuj = (1. - guii)*(1. - gujj) + (delta_tij - guji)*guij;
			const double ndindj = (1. - gdii)*(1. - gdjj) + (delta_tij - gdji)*gdij;
			m->vv[r + num_ij*t] += pre*nuinuj*ndindj;
			m->vn[r + num_ij*t] += pre*(nuinuj*(1. - gdii) + (1. - guii)*ndindj);
		}
	}
	}

	// 1 bond 1 site measurements
	if (meas_energy_corr)
	#pragma omp parallel for
	for (int t = 0; t < L; t++) {
		const int delta_t = (t == 0);
		const double *const restrict Gu0t_t = Gu0t + N*N*t;
		const double *const restrict Gutt_t = Gutt + N*N*t;
		const double *const restrict Gut0_t = Gut0 + N*N*t;
		const double *const restrict Gd0t_t = Gd0t + N*N*t;
		const double *const restrict Gdtt_t = Gdtt + N*N*t;
		const double *const restrict Gdt0_t = Gdt0 + N*N*t;
	for (int j = 0; j < N; j++)
	for (int b = 0; b < num_b; b++) {
		const int i0 = p->bonds[b];
		const int i1 = p->bonds[b + num_b];
		const int bs = p->map_bs[b + num_b*j];
		const double pre = (double)sign / p->degen_bs[bs];
		const int delta_i0i1 = 0;
		const int delta_i0j = delta_t*(i0 == j);
		const int delta_i1j = delta_t*(i1 == j);
		const double gui0j = Gut0_t[i0 + N*j];
		const double guji0 = Gu0t_t[j + N*i0];
		const double gdi0j = Gdt0_t[i0 + N*j];
		const double gdji0 = Gd0t_t[j + N*i0];
		const double gui1j = Gut0_t[i1 + N*j];
		const double guji1 = Gu0t_t[j + N*i1];
		const double gdi1j = Gdt0_t[i1 + N*j];
		const double gdji1 = Gd0t_t[j + N*i1];
		const double gui0i1 = Gutt_t[i0 + N*i1];
		const double gui1i0 = Gutt_t[i1 + N*i0];
		const double gdi0i1 = Gdtt_t[i0 + N*i1];
		const double gdi1i0 = Gdtt_t[i1 + N*i0];
		const double gujj = Gu00[j + N*j];
		const double gdjj = Gd00[j + N*j];

		const double ku = 2.*delta_i0i1 - gui0i1 - gui1i0;
		const double kd = 2.*delta_i0i1 - gdi0i1 - gdi1i0;
		const double xu = (delta_i0j - guji0)*gui1j + (delta_i1j - guji1)*gui0j;
		const double xd = (delta_i0j - gdji0)*gdi1j + (delta_i1j - gdji1)*gdi0j;
		m->kv[bs + num_bs*t] += pre*((ku*(1. - gujj) + xu)*(1. - gdjj)
					   + (kd*(1. - gdjj) + xd)*(1. - gujj));
		m->kn[bs + num_bs*t] += pre*((ku + kd)*(2. - gujj - gdjj) + xu + xd);
	}
	}

	// 2 bond measurements
	// minor optimization: handle t = 0 separately, since there are no delta
	// functions for t > 0. not really needed in 2-site measurements above
	// as those are fast anyway
	if (meas_bond_corr || meas_thermal)
	for (int c = 0; c < num_b; c++) {
		const int j0 = p->bonds[c];
		const int j1 = p->bonds[c + num_b];
	for (int b = 0; b < num_b; b++) {
		const int i0 = p->bonds[b];
		const int i1 = p->bonds[b + num_b];
		const int bb = p->map_bb[b + c*num_b];
		const double pre = (double)sign / p->degen_bb[bb];
		const int delta_i0j0 = (i0 == j0);
		const int delta_i1j0 = (i1 == j0);
		const int delta_i0j1 = (i0 == j1);
		const int delta_i1j1 = (i1 == j1);
		const double gui1i0 = Gu00[i1 + i0*N];
		const double gui0i1 = Gu00[i0 + i1*N];
		const double gui0j0 = Gu00[i0 + j0*N];
		const double gui1j0 = Gu00[i1 + j0*N];
		const double gui0j1 = Gu00[i0 + j1*N];
		const double gui1j1 = Gu00[i1 + j1*N];
		const double guj0i0 = Gu00[j0 + i0*N];
		const double guj1i0 = Gu00[j1 + i0*N];
		const double guj0i1 = Gu00[j0 + i1*N];
		const double guj1i1 = Gu00[j1 + i1*N];
		const double guj1j0 = Gu00[j1 + j0*N];
		const double guj0j1 = Gu00[j0 + j1*N];
		const double gdi1i0 = Gd00[i1 + i0*N];
		const double gdi0i1 = Gd00[i0 + i1*N];
		const double gdi0j0 = Gd00[i0 + j0*N];
		const double gdi1j0 = Gd00[i1 + j0*N];
		const double gdi0j1 = Gd00[i0 + j1*N];
		const double gdi1j1 = Gd00[i1 + j1*N];
		const double gdj0i0 = Gd00[j0 + i0*N];
		const double gdj1i0 = Gd00[j1 + i0*N];
		const double gdj0i1 = Gd00[j0 + i1*N];
		const double gdj1i1 = Gd00[j1 + i1*N];
		const double gdj1j0 = Gd00[j1 + j0*N];
		const double gdj0j1 = Gd00[j0 + j1*N];
                
                if (meas_bond_corr){
		m->pair_bb[bb] += 0.5*pre*(gui0j0*gdi1j1 + gui1j0*gdi0j1 + gui0j1*gdi1j0 + gui1j1*gdi0j0);
		const double x = ((delta_i0j1 - guj1i0)*gui1j0 + (delta_i1j0 - guj0i1)*gui0j1
				+ (delta_i0j1 - gdj1i0)*gdi1j0 + (delta_i1j0 - gdj0i1)*gdi0j1);
		const double y = ((delta_i0j0 - guj0i0)*gui1j1 + (delta_i1j1 - guj1i1)*gui0j0
				+ (delta_i0j0 - gdj0i0)*gdi1j1 + (delta_i1j1 - gdj1i1)*gdi0j0);
		m->jj[bb]   += pre*((gui0i1 - gui1i0 + gdi0i1 - gdi1i0)*(guj0j1 - guj1j0 + gdj0j1 - gdj1j0) + x - y);
		m->jsjs[bb] += pre*((gui0i1 - gui1i0 - gdi0i1 + gdi1i0)*(guj0j1 - guj1j0 - gdj0j1 + gdj1j0) + x - y);
		m->kk[bb]   += pre*((gui0i1 + gui1i0 + gdi0i1 + gdi1i0)*(guj0j1 + guj1j0 + gdj0j1 + gdj1j0) + x + y);
		m->ksks[bb] += pre*((gui0i1 + gui1i0 - gdi0i1 - gdi1i0)*(guj0j1 + guj1j0 - gdj0j1 - gdj1j0) + x + y);
                }
                
                if (meas_thermal){
                const double gui0i0 = Gu00[i0 + i0*N];
                const double gui1i1 = Gu00[i1 + i1*N];
                const double guj0j0 = Gu00[j0 + j0*N];
                const double guj1j1 = Gu00[j1 + j1*N];
                const double gdi0i0 = Gd00[i0 + i0*N];
                const double gdi1i1 = Gd00[i1 + i1*N];
                const double gdj0j0 = Gd00[j0 + j0*N];
                const double gdj1j1 = Gd00[j1 + j1*N];      
                
		const double ujinj = (gui1j0*(delta_i0j0-guj0i0)-gui0j0*(delta_i1j0-guj0i1)+gui1j1*(delta_i0j1-guj1i0)-(gui0j1*(delta_i1j1-guj1i1))
                              +(gui0i1-gui1i0)*(2-guj0j0-guj1j1));
                const double djinj = (gdi1j0*(delta_i0j0-gdj0i0)-gdi0j0*(delta_i1j0-gdj0i1)+gdi1j1*(delta_i0j1-gdj1i0)-(gdi0j1*(delta_i1j1-gdj1i1))
                              +(gdi0i1-gdi1i0)*(2-gdj0j0-gdj1j1));
                const double ujijj = ((delta_i0j1-guj1i0)*gui1j0-(delta_i0j0-guj0i0)*gui1j1-(delta_i1j1-guj1i1)*gui0j0+(delta_i1j0-guj0i1)*gui0j1
                               +(gui1i0-gui0i1)*(guj1j0-guj0j1));
                const double djijj = ((delta_i0j1-gdj1i0)*gdi1j0-(delta_i0j0-gdj0i0)*gdi1j1-(delta_i1j1-gdj1i1)*gdi0j0+(delta_i1j0-gdj0i1)*gdi0j1
                               +(gdi1i0-gdi0i1)*(gdj1j0-gdj0j1));	
		const double unijj = (gui1j0*(delta_i1j1-guj1i1)-gui0j1*(delta_i0j0-guj0i0)+gui0j0*(delta_i0j1-guj1i0)-gui1j1*(delta_i1j0-guj0i1)
                              +(guj0j1-guj1j0)*(2-gui0i0-gui1i1));
                const double dnijj = (gdi1j0*(delta_i1j1-gdj1i1)-gdi0j1*(delta_i0j0-gdj0i0)+gdi0j0*(delta_i0j1-gdj1i0)-gdi1j1*(delta_i1j0-gdj0i1)
                              +(gdj0j1-gdj1j0)*(2-gdi0i0-gdi1i1));
		
                const double part1u = ujijj*((2-gdi0i0-gdi1i1)*(2-gdj0j0-gdj1j1)+gdi0j0*(delta_i0j0-gdj0i0)
                                      +gdi1j0*(delta_i1j0-gdj0i1) +gdi0j1*(delta_i0j1-gdj1i0) +gdi1j1*(delta_i1j1-gdj1i1) );
    
                const double part2u = dnijj*ujinj;
    
                const double part1d = djijj*((2-gui0i0-gui1i1)*(2-guj0j0-guj1j1)+gui0j0*(delta_i0j0-guj0i0)
                                      +gui1j0*(delta_i1j0-guj0i1) +gui0j1*(delta_i0j1-guj1i0) +gui1j1*(delta_i1j1-guj1i1) );
    
                const double part2d = unijj*djinj;
                
		m->jnjn[bb]   += pre*(part1u+part1d+part2u+part2d);
                }
	}
	}
	
        if (meas_thermal)
	for (int c = 0; c < num_b; c++) {
		const int j0 = p->bonds[c];
		const int j1 = p->bonds[c + num_b];
	for (int b = 0; b < num_b2; b++) {
		const int i0 = p->bond2s[b];
		const int i1 = p->bond2s[b + num_b2];
		const int bb2 = p->map_bb2[b + c*num_b2];
		const double pre = (double)sign / p->degen_bb2[bb2];
		const int delta_i0j0 = (i0 == j0);
		const int delta_i1j0 = (i1 == j0);
		const int delta_i0j1 = (i0 == j1);
		const int delta_i1j1 = (i1 == j1);
		const double gui1i0 = Gu00[i1 + i0*N];
		const double gui0i1 = Gu00[i0 + i1*N];
		const double gui0j0 = Gu00[i0 + j0*N];
		const double gui1j0 = Gu00[i1 + j0*N];
		const double gui0j1 = Gu00[i0 + j1*N];
		const double gui1j1 = Gu00[i1 + j1*N];
		const double guj0i0 = Gu00[j0 + i0*N];
		const double guj1i0 = Gu00[j1 + i0*N];
		const double guj0i1 = Gu00[j0 + i1*N];
		const double guj1i1 = Gu00[j1 + i1*N];
		const double guj1j0 = Gu00[j1 + j0*N];
		const double guj0j1 = Gu00[j0 + j1*N];
		const double gdi1i0 = Gd00[i1 + i0*N];
		const double gdi0i1 = Gd00[i0 + i1*N];
		const double gdi0j0 = Gd00[i0 + j0*N];
		const double gdi1j0 = Gd00[i1 + j0*N];
		const double gdi0j1 = Gd00[i0 + j1*N];
		const double gdi1j1 = Gd00[i1 + j1*N];
		const double gdj0i0 = Gd00[j0 + i0*N];
		const double gdj1i0 = Gd00[j1 + i0*N];
		const double gdj0i1 = Gd00[j0 + i1*N];
		const double gdj1i1 = Gd00[j1 + i1*N];
		const double gdj1j0 = Gd00[j1 + j0*N];
		const double gdj0j1 = Gd00[j0 + j1*N];
	
		const double gui0i0 = Gu00[i0 + i0*N];
                const double gui1i1 = Gu00[i1 + i1*N];
                const double guj0j0 = Gu00[j0 + j0*N];
                const double guj1j1 = Gu00[j1 + j1*N];
                const double gdi0i0 = Gd00[i0 + i0*N];
                const double gdi1i1 = Gd00[i1 + i1*N];
                const double gdj0j0 = Gd00[j0 + j0*N];
                const double gdj1j1 = Gd00[j1 + j1*N];     
		
		const double ujinj = (gui1j0*(delta_i0j0-guj0i0)-gui0j0*(delta_i1j0-guj0i1)+gui1j1*(delta_i0j1-guj1i0)-(gui0j1*(delta_i1j1-guj1i1))
                              +(gui0i1-gui1i0)*(2-guj0j0-guj1j1));
                const double djinj = (gdi1j0*(delta_i0j0-gdj0i0)-gdi0j0*(delta_i1j0-gdj0i1)+gdi1j1*(delta_i0j1-gdj1i0)-(gdi0j1*(delta_i1j1-gdj1i1))
                              +(gdi0i1-gdi1i0)*(2-gdj0j0-gdj1j1));
        
                const double part1u = (2-guj0j0-guj1j1)*((delta_i0j1-gdj1i0)*gdi1j0-(delta_i0j0-gdj0i0)*gdi1j1-(delta_i1j1-gdj1i1)*gdi0j0+(delta_i1j0-gdj0i1)*gdi0j1
                                      +(gdj1j0-gdj0j1)*(gdi1i0-gdi0i1));
    
                const double part2u = -(gdj1j0-gdj0j1)*ujinj;
    
                const double part1d = (2-gdj0j0-gdj1j1)*((delta_i0j1-guj1i0)*gui1j0-(delta_i0j0-guj0i0)*gui1j1-(delta_i1j1-guj1i1)*gui0j0+(delta_i1j0-guj0i1)*gui0j1
                                      +(guj1j0-guj0j1)*(gui1i0-gui0i1));
    
                const double part2d = -(guj1j0-guj0j1)*djinj;
                
		m->jjn[bb2]  += pre*(part1u+part1d+part2u+part2d);
	}
	}
	
	if (meas_thermal)
        for (int c = 0; c < num_b2; c++) {
		const int j0 = p->bond2s[c];
		const int j1 = p->bond2s[c + num_b2];
	for (int b = 0; b < num_b; b++) {
		const int i0 = p->bonds[b];
		const int i1 = p->bonds[b + num_b];
		const int b2b = p->map_b2b[b + c*num_b];
		const double pre = (double)sign / p->degen_b2b[b2b];
		const int delta_i0j0 = (i0 == j0);
		const int delta_i1j0 = (i1 == j0);
		const int delta_i0j1 = (i0 == j1);
		const int delta_i1j1 = (i1 == j1);
		const double gui1i0 = Gu00[i1 + i0*N];
		const double gui0i1 = Gu00[i0 + i1*N];
		const double gui0j0 = Gu00[i0 + j0*N];
		const double gui1j0 = Gu00[i1 + j0*N];
		const double gui0j1 = Gu00[i0 + j1*N];
		const double gui1j1 = Gu00[i1 + j1*N];
		const double guj0i0 = Gu00[j0 + i0*N];
		const double guj1i0 = Gu00[j1 + i0*N];
		const double guj0i1 = Gu00[j0 + i1*N];
		const double guj1i1 = Gu00[j1 + i1*N];
		const double guj1j0 = Gu00[j1 + j0*N];
		const double guj0j1 = Gu00[j0 + j1*N];
		const double gdi1i0 = Gd00[i1 + i0*N];
		const double gdi0i1 = Gd00[i0 + i1*N];
		const double gdi0j0 = Gd00[i0 + j0*N];
		const double gdi1j0 = Gd00[i1 + j0*N];
		const double gdi0j1 = Gd00[i0 + j1*N];
		const double gdi1j1 = Gd00[i1 + j1*N];
		const double gdj0i0 = Gd00[j0 + i0*N];
		const double gdj1i0 = Gd00[j1 + i0*N];
		const double gdj0i1 = Gd00[j0 + i1*N];
		const double gdj1i1 = Gd00[j1 + i1*N];
		const double gdj1j0 = Gd00[j1 + j0*N];
		const double gdj0j1 = Gd00[j0 + j1*N];
	
		const double gui0i0 = Gu00[i0 + i0*N];
                const double gui1i1 = Gu00[i1 + i1*N];
                const double guj0j0 = Gu00[j0 + j0*N];
                const double guj1j1 = Gu00[j1 + j1*N];
                const double gdi0i0 = Gd00[i0 + i0*N];
                const double gdi1i1 = Gd00[i1 + i1*N];
                const double gdj0j0 = Gd00[j0 + j0*N];
                const double gdj1j1 = Gd00[j1 + j1*N];     
		
                const double ujijj = ((delta_i0j1-guj1i0)*gui1j0-(delta_i0j0-guj0i0)*gui1j1-(delta_i1j1-guj1i1)*gui0j0+(delta_i1j0-guj0i1)*gui0j1
                               +(gui1i0-gui0i1)*(guj1j0-guj0j1));
                const double djijj = ((delta_i0j1-gdj1i0)*gdi1j0-(delta_i0j0-gdj0i0)*gdi1j1-(delta_i1j1-gdj1i1)*gdi0j0+(delta_i1j0-gdj0i1)*gdi0j1
                               +(gdi1i0-gdi0i1)*(gdj1j0-gdj0j1));
                const double unijj = (gui1j0*(delta_i1j1-guj1i1)-gui0j1*(delta_i0j0-guj0i0)+gui0j0*(delta_i0j1-guj1i0)-gui1j1*(delta_i1j0-guj0i1)
                              +(guj0j1-guj1j0)*(2-gui0i0-gui1i1));
                const double dnijj = (gdi1j0*(delta_i1j1-gdj1i1)-gdi0j1*(delta_i0j0-gdj0i0)+gdi0j0*(delta_i0j1-gdj1i0)-gdi1j1*(delta_i1j0-gdj0i1)
                              +(gdj0j1-gdj1j0)*(2-gdi0i0-gdi1i1));
                
                const double part1u = (2-gui0i0-gui1i1)*djijj;
    
                const double part2u = -(gdi1i0-gdi0i1)*unijj;
    
                const double part1d = (2-gdi0i0-gdi1i1)*ujijj;
    
                const double part2d = -(gui1i0-gui0i1)*dnijj;
                
		m->jnj[b2b] += pre*(part1u+part1d+part2u+part2d);
	}
	}
		
	if (meas_2bond_corr)
	for (int c = 0; c < num_b2; c++) {
		const int j0 = p->bond2s[c];
		const int j1 = p->bond2s[c + num_b2];
	for (int b = 0; b < num_b2; b++) {
		const int i0 = p->bond2s[b];
		const int i1 = p->bond2s[b + num_b2];
		const int bb = p->map_b2b2[b + c*num_b2];
		const double pre = (double)sign / p->degen_b2b2[bb];
		const int delta_i0j0 = (i0 == j0);
		const int delta_i1j0 = (i1 == j0);
		const int delta_i0j1 = (i0 == j1);
		const int delta_i1j1 = (i1 == j1);
		const double gui1i0 = Gu00[i1 + i0*N];
		const double gui0i1 = Gu00[i0 + i1*N];
		const double gui0j0 = Gu00[i0 + j0*N];
		const double gui1j0 = Gu00[i1 + j0*N];
		const double gui0j1 = Gu00[i0 + j1*N];
		const double gui1j1 = Gu00[i1 + j1*N];
		const double guj0i0 = Gu00[j0 + i0*N];
		const double guj1i0 = Gu00[j1 + i0*N];
		const double guj0i1 = Gu00[j0 + i1*N];
		const double guj1i1 = Gu00[j1 + i1*N];
		const double guj1j0 = Gu00[j1 + j0*N];
		const double guj0j1 = Gu00[j0 + j1*N];
		const double gdi1i0 = Gd00[i1 + i0*N];
		const double gdi0i1 = Gd00[i0 + i1*N];
		const double gdi0j0 = Gd00[i0 + j0*N];
		const double gdi1j0 = Gd00[i1 + j0*N];
		const double gdi0j1 = Gd00[i0 + j1*N];
		const double gdi1j1 = Gd00[i1 + j1*N];
		const double gdj0i0 = Gd00[j0 + i0*N];
		const double gdj1i0 = Gd00[j1 + i0*N];
		const double gdj0i1 = Gd00[j0 + i1*N];
		const double gdj1i1 = Gd00[j1 + i1*N];
		const double gdj1j0 = Gd00[j1 + j0*N];
		const double gdj0j1 = Gd00[j0 + j1*N];
		m->pair_b2b2[bb] += 0.5*pre*(gui0j0*gdi1j1 + gui1j0*gdi0j1 + gui0j1*gdi1j0 + gui1j1*gdi0j0);
		const double x = ((delta_i0j1 - guj1i0)*gui1j0 + (delta_i1j0 - guj0i1)*gui0j1
				+ (delta_i0j1 - gdj1i0)*gdi1j0 + (delta_i1j0 - gdj0i1)*gdi0j1);
		const double y = ((delta_i0j0 - guj0i0)*gui1j1 + (delta_i1j1 - guj1i1)*gui0j0
				+ (delta_i0j0 - gdj0i0)*gdi1j1 + (delta_i1j1 - gdj1i1)*gdi0j0);
		m->j2j2[bb]   += pre*((gui0i1 - gui1i0 + gdi0i1 - gdi1i0)*(guj0j1 - guj1j0 + gdj0j1 - gdj1j0) + x - y);
		m->js2js2[bb] += pre*((gui0i1 - gui1i0 - gdi0i1 + gdi1i0)*(guj0j1 - guj1j0 - gdj0j1 + gdj1j0) + x - y);
		m->k2k2[bb]   += pre*((gui0i1 + gui1i0 + gdi0i1 + gdi1i0)*(guj0j1 + guj1j0 + gdj0j1 + gdj1j0) + x + y);
		m->ks2ks2[bb] += pre*((gui0i1 + gui1i0 - gdi0i1 - gdi1i0)*(guj0j1 + guj1j0 - gdj0j1 - gdj1j0) + x + y);
	}
	}

	if (meas_nematic_corr)
	for (int c = 0; c < NEM_BONDS*N; c++) {
		const int j0 = p->bonds[c];
		const int j1 = p->bonds[c + num_b];
	for (int b = 0; b < NEM_BONDS*N; b++) {
		const int i0 = p->bonds[b];
		const int i1 = p->bonds[b + num_b];
		const int bb = p->map_bb[b + c*num_b];
		const double pre = (double)sign / p->degen_bb[bb];
		const int delta_i0j0 = (i0 == j0);
		const int delta_i1j0 = (i1 == j0);
		const int delta_i0j1 = (i0 == j1);
		const int delta_i1j1 = (i1 == j1);
		const double gui0i0 = Gu00[i0 + i0*N];
		const double gui1i0 = Gu00[i1 + i0*N];
		const double gui0i1 = Gu00[i0 + i1*N];
		const double gui1i1 = Gu00[i1 + i1*N];
		const double gui0j0 = Gu00[i0 + j0*N];
		const double gui1j0 = Gu00[i1 + j0*N];
		const double gui0j1 = Gu00[i0 + j1*N];
		const double gui1j1 = Gu00[i1 + j1*N];
		const double guj0i0 = Gu00[j0 + i0*N];
		const double guj1i0 = Gu00[j1 + i0*N];
		const double guj0i1 = Gu00[j0 + i1*N];
		const double guj1i1 = Gu00[j1 + i1*N];
		const double guj0j0 = Gu00[j0 + j0*N];
		const double guj1j0 = Gu00[j1 + j0*N];
		const double guj0j1 = Gu00[j0 + j1*N];
		const double guj1j1 = Gu00[j1 + j1*N];
		const double gdi0i0 = Gd00[i0 + i0*N];
		const double gdi1i0 = Gd00[i1 + i0*N];
		const double gdi0i1 = Gd00[i0 + i1*N];
		const double gdi1i1 = Gd00[i1 + i1*N];
		const double gdi0j0 = Gd00[i0 + j0*N];
		const double gdi1j0 = Gd00[i1 + j0*N];
		const double gdi0j1 = Gd00[i0 + j1*N];
		const double gdi1j1 = Gd00[i1 + j1*N];
		const double gdj0i0 = Gd00[j0 + i0*N];
		const double gdj1i0 = Gd00[j1 + i0*N];
		const double gdj0i1 = Gd00[j0 + i1*N];
		const double gdj1i1 = Gd00[j1 + i1*N];
		const double gdj0j0 = Gd00[j0 + j0*N];
		const double gdj1j0 = Gd00[j1 + j0*N];
		const double gdj0j1 = Gd00[j0 + j1*N];
		const double gdj1j1 = Gd00[j1 + j1*N];
		const int delta_i0i1 = 0;
		const int delta_j0j1 = 0;
		const double uuuu = +(1.-gui0i0)*(1.-gui1i1)*(1.-guj0j0)*(1.-guj1j1)+(1.-gui0i0)*(1.-gui1i1)*(delta_j0j1-guj1j0)*guj0j1+(1.-gui0i0)*(delta_i1j0-guj0i1)*gui1j0*(1.-guj1j1)-(1.-gui0i0)*(delta_i1j0-guj0i1)*gui1j1*(delta_j0j1-guj1j0)+(1.-gui0i0)*(delta_i1j1-guj1i1)*gui1j0*guj0j1+(1.-gui0i0)*(delta_i1j1-guj1i1)*gui1j1*(1.-guj0j0)+(delta_i0i1-gui1i0)*gui0i1*(1.-guj0j0)*(1.-guj1j1)+(delta_i0i1-gui1i0)*gui0i1*(delta_j0j1-guj1j0)*guj0j1-(delta_i0i1-gui1i0)*gui0j0*(delta_i1j0-guj0i1)*(1.-guj1j1)-(delta_i0i1-gui1i0)*gui0j0*(delta_i1j1-guj1i1)*guj0j1+(delta_i0i1-gui1i0)*gui0j1*(delta_i1j0-guj0i1)*(delta_j0j1-guj1j0)-(delta_i0i1-gui1i0)*gui0j1*(delta_i1j1-guj1i1)*(1.-guj0j0)+(delta_i0j0-guj0i0)*gui0i1*gui1j0*(1.-guj1j1)-(delta_i0j0-guj0i0)*gui0i1*gui1j1*(delta_j0j1-guj1j0)+(delta_i0j0-guj0i0)*gui0j0*(1.-gui1i1)*(1.-guj1j1)+(delta_i0j0-guj0i0)*gui0j0*(delta_i1j1-guj1i1)*gui1j1-(delta_i0j0-guj0i0)*gui0j1*(1.-gui1i1)*(delta_j0j1-guj1j0)-(delta_i0j0-guj0i0)*gui0j1*(delta_i1j1-guj1i1)*gui1j0+(delta_i0j1-guj1i0)*gui0i1*gui1j0*guj0j1+(delta_i0j1-guj1i0)*gui0i1*gui1j1*(1.-guj0j0)+(delta_i0j1-guj1i0)*gui0j0*(1.-gui1i1)*guj0j1-(delta_i0j1-guj1i0)*gui0j0*(delta_i1j0-guj0i1)*gui1j1+(delta_i0j1-guj1i0)*gui0j1*(1.-gui1i1)*(1.-guj0j0)+(delta_i0j1-guj1i0)*gui0j1*(delta_i1j0-guj0i1)*gui1j0;
		const double uuud = +(1.-gui0i0)*(1.-gui1i1)*(1.-guj0j0)*(1.-gdj1j1)+(1.-gui0i0)*(delta_i1j0-guj0i1)*gui1j0*(1.-gdj1j1)+(delta_i0i1-gui1i0)*gui0i1*(1.-guj0j0)*(1.-gdj1j1)-(delta_i0i1-gui1i0)*gui0j0*(delta_i1j0-guj0i1)*(1.-gdj1j1)+(delta_i0j0-guj0i0)*gui0i1*gui1j0*(1.-gdj1j1)+(delta_i0j0-guj0i0)*gui0j0*(1.-gui1i1)*(1.-gdj1j1);
		const double uudu = +(1.-gui0i0)*(1.-gui1i1)*(1.-gdj0j0)*(1.-guj1j1)+(1.-gui0i0)*(delta_i1j1-guj1i1)*gui1j1*(1.-gdj0j0)+(delta_i0i1-gui1i0)*gui0i1*(1.-gdj0j0)*(1.-guj1j1)-(delta_i0i1-gui1i0)*gui0j1*(delta_i1j1-guj1i1)*(1.-gdj0j0)+(delta_i0j1-guj1i0)*gui0i1*gui1j1*(1.-gdj0j0)+(delta_i0j1-guj1i0)*gui0j1*(1.-gui1i1)*(1.-gdj0j0);
		const double uudd = +(1.-gui0i0)*(1.-gui1i1)*(1.-gdj0j0)*(1.-gdj1j1)+(1.-gui0i0)*(1.-gui1i1)*(delta_j0j1-gdj1j0)*gdj0j1+(delta_i0i1-gui1i0)*gui0i1*(1.-gdj0j0)*(1.-gdj1j1)+(delta_i0i1-gui1i0)*gui0i1*(delta_j0j1-gdj1j0)*gdj0j1;
		const double uduu = +(1.-gui0i0)*(1.-gdi1i1)*(1.-guj0j0)*(1.-guj1j1)+(1.-gui0i0)*(1.-gdi1i1)*(delta_j0j1-guj1j0)*guj0j1+(delta_i0j0-guj0i0)*gui0j0*(1.-gdi1i1)*(1.-guj1j1)-(delta_i0j0-guj0i0)*gui0j1*(1.-gdi1i1)*(delta_j0j1-guj1j0)+(delta_i0j1-guj1i0)*gui0j0*(1.-gdi1i1)*guj0j1+(delta_i0j1-guj1i0)*gui0j1*(1.-gdi1i1)*(1.-guj0j0);
		const double udud = +(1.-gui0i0)*(1.-gdi1i1)*(1.-guj0j0)*(1.-gdj1j1)+(1.-gui0i0)*(delta_i1j1-gdj1i1)*gdi1j1*(1.-guj0j0)+(delta_i0j0-guj0i0)*gui0j0*(1.-gdi1i1)*(1.-gdj1j1)+(delta_i0j0-guj0i0)*gui0j0*(delta_i1j1-gdj1i1)*gdi1j1;
		const double uddu = +(1.-gui0i0)*(1.-gdi1i1)*(1.-gdj0j0)*(1.-guj1j1)+(1.-gui0i0)*(delta_i1j0-gdj0i1)*gdi1j0*(1.-guj1j1)+(delta_i0j1-guj1i0)*gui0j1*(1.-gdi1i1)*(1.-gdj0j0)+(delta_i0j1-guj1i0)*gui0j1*(delta_i1j0-gdj0i1)*gdi1j0;
		const double uddd = +(1.-gui0i0)*(1.-gdi1i1)*(1.-gdj0j0)*(1.-gdj1j1)+(1.-gui0i0)*(1.-gdi1i1)*(delta_j0j1-gdj1j0)*gdj0j1+(1.-gui0i0)*(delta_i1j0-gdj0i1)*gdi1j0*(1.-gdj1j1)-(1.-gui0i0)*(delta_i1j0-gdj0i1)*gdi1j1*(delta_j0j1-gdj1j0)+(1.-gui0i0)*(delta_i1j1-gdj1i1)*gdi1j0*gdj0j1+(1.-gui0i0)*(delta_i1j1-gdj1i1)*gdi1j1*(1.-gdj0j0);
		const double duuu = +(1.-gdi0i0)*(1.-gui1i1)*(1.-guj0j0)*(1.-guj1j1)+(1.-gdi0i0)*(1.-gui1i1)*(delta_j0j1-guj1j0)*guj0j1+(1.-gdi0i0)*(delta_i1j0-guj0i1)*gui1j0*(1.-guj1j1)-(1.-gdi0i0)*(delta_i1j0-guj0i1)*gui1j1*(delta_j0j1-guj1j0)+(1.-gdi0i0)*(delta_i1j1-guj1i1)*gui1j0*guj0j1+(1.-gdi0i0)*(delta_i1j1-guj1i1)*gui1j1*(1.-guj0j0);
		const double duud = +(1.-gdi0i0)*(1.-gui1i1)*(1.-guj0j0)*(1.-gdj1j1)+(1.-gdi0i0)*(delta_i1j0-guj0i1)*gui1j0*(1.-gdj1j1)+(delta_i0j1-gdj1i0)*gdi0j1*(1.-gui1i1)*(1.-guj0j0)+(delta_i0j1-gdj1i0)*gdi0j1*(delta_i1j0-guj0i1)*gui1j0;
		const double dudu = +(1.-gdi0i0)*(1.-gui1i1)*(1.-gdj0j0)*(1.-guj1j1)+(1.-gdi0i0)*(delta_i1j1-guj1i1)*gui1j1*(1.-gdj0j0)+(delta_i0j0-gdj0i0)*gdi0j0*(1.-gui1i1)*(1.-guj1j1)+(delta_i0j0-gdj0i0)*gdi0j0*(delta_i1j1-guj1i1)*gui1j1;
		const double dudd = +(1.-gdi0i0)*(1.-gui1i1)*(1.-gdj0j0)*(1.-gdj1j1)+(1.-gdi0i0)*(1.-gui1i1)*(delta_j0j1-gdj1j0)*gdj0j1+(delta_i0j0-gdj0i0)*gdi0j0*(1.-gui1i1)*(1.-gdj1j1)-(delta_i0j0-gdj0i0)*gdi0j1*(1.-gui1i1)*(delta_j0j1-gdj1j0)+(delta_i0j1-gdj1i0)*gdi0j0*(1.-gui1i1)*gdj0j1+(delta_i0j1-gdj1i0)*gdi0j1*(1.-gui1i1)*(1.-gdj0j0);
		const double dduu = +(1.-gdi0i0)*(1.-gdi1i1)*(1.-guj0j0)*(1.-guj1j1)+(1.-gdi0i0)*(1.-gdi1i1)*(delta_j0j1-guj1j0)*guj0j1+(delta_i0i1-gdi1i0)*gdi0i1*(1.-guj0j0)*(1.-guj1j1)+(delta_i0i1-gdi1i0)*gdi0i1*(delta_j0j1-guj1j0)*guj0j1;
		const double ddud = +(1.-gdi0i0)*(1.-gdi1i1)*(1.-guj0j0)*(1.-gdj1j1)+(1.-gdi0i0)*(delta_i1j1-gdj1i1)*gdi1j1*(1.-guj0j0)+(delta_i0i1-gdi1i0)*gdi0i1*(1.-guj0j0)*(1.-gdj1j1)-(delta_i0i1-gdi1i0)*gdi0j1*(delta_i1j1-gdj1i1)*(1.-guj0j0)+(delta_i0j1-gdj1i0)*gdi0i1*gdi1j1*(1.-guj0j0)+(delta_i0j1-gdj1i0)*gdi0j1*(1.-gdi1i1)*(1.-guj0j0);
		const double dddu = +(1.-gdi0i0)*(1.-gdi1i1)*(1.-gdj0j0)*(1.-guj1j1)+(1.-gdi0i0)*(delta_i1j0-gdj0i1)*gdi1j0*(1.-guj1j1)+(delta_i0i1-gdi1i0)*gdi0i1*(1.-gdj0j0)*(1.-guj1j1)-(delta_i0i1-gdi1i0)*gdi0j0*(delta_i1j0-gdj0i1)*(1.-guj1j1)+(delta_i0j0-gdj0i0)*gdi0i1*gdi1j0*(1.-guj1j1)+(delta_i0j0-gdj0i0)*gdi0j0*(1.-gdi1i1)*(1.-guj1j1);
		const double dddd = +(1.-gdi0i0)*(1.-gdi1i1)*(1.-gdj0j0)*(1.-gdj1j1)+(1.-gdi0i0)*(1.-gdi1i1)*(delta_j0j1-gdj1j0)*gdj0j1+(1.-gdi0i0)*(delta_i1j0-gdj0i1)*gdi1j0*(1.-gdj1j1)-(1.-gdi0i0)*(delta_i1j0-gdj0i1)*gdi1j1*(delta_j0j1-gdj1j0)+(1.-gdi0i0)*(delta_i1j1-gdj1i1)*gdi1j0*gdj0j1+(1.-gdi0i0)*(delta_i1j1-gdj1i1)*gdi1j1*(1.-gdj0j0)+(delta_i0i1-gdi1i0)*gdi0i1*(1.-gdj0j0)*(1.-gdj1j1)+(delta_i0i1-gdi1i0)*gdi0i1*(delta_j0j1-gdj1j0)*gdj0j1-(delta_i0i1-gdi1i0)*gdi0j0*(delta_i1j0-gdj0i1)*(1.-gdj1j1)-(delta_i0i1-gdi1i0)*gdi0j0*(delta_i1j1-gdj1i1)*gdj0j1+(delta_i0i1-gdi1i0)*gdi0j1*(delta_i1j0-gdj0i1)*(delta_j0j1-gdj1j0)-(delta_i0i1-gdi1i0)*gdi0j1*(delta_i1j1-gdj1i1)*(1.-gdj0j0)+(delta_i0j0-gdj0i0)*gdi0i1*gdi1j0*(1.-gdj1j1)-(delta_i0j0-gdj0i0)*gdi0i1*gdi1j1*(delta_j0j1-gdj1j0)+(delta_i0j0-gdj0i0)*gdi0j0*(1.-gdi1i1)*(1.-gdj1j1)+(delta_i0j0-gdj0i0)*gdi0j0*(delta_i1j1-gdj1i1)*gdi1j1-(delta_i0j0-gdj0i0)*gdi0j1*(1.-gdi1i1)*(delta_j0j1-gdj1j0)-(delta_i0j0-gdj0i0)*gdi0j1*(delta_i1j1-gdj1i1)*gdi1j0+(delta_i0j1-gdj1i0)*gdi0i1*gdi1j0*gdj0j1+(delta_i0j1-gdj1i0)*gdi0i1*gdi1j1*(1.-gdj0j0)+(delta_i0j1-gdj1i0)*gdi0j0*(1.-gdi1i1)*gdj0j1-(delta_i0j1-gdj1i0)*gdi0j0*(delta_i1j0-gdj0i1)*gdi1j1+(delta_i0j1-gdj1i0)*gdi0j1*(1.-gdi1i1)*(1.-gdj0j0)+(delta_i0j1-gdj1i0)*gdi0j1*(delta_i1j0-gdj0i1)*gdi1j0;
		m->nem_nnnn[bb] += pre*(uuuu + uuud + uudu + uudd
				      + uduu + udud + uddu + uddd
				      + duuu + duud + dudu + dudd
				      + dduu + ddud + dddu + dddd);
		m->nem_ssss[bb] += pre*(uuuu - uuud - uudu + uudd
				      - uduu + udud + uddu - uddd
				      - duuu + duud + dudu - dudd
				      + dduu - ddud - dddu + dddd);
	}
	}

	// no delta functions here.
	if (meas_bond_corr || meas_thermal)
	#pragma omp parallel for
	for (int t = 1; t < L; t++) {
		const double *const restrict Gu0t_t = Gu0t + N*N*t;
		const double *const restrict Gutt_t = Gutt + N*N*t;
		const double *const restrict Gut0_t = Gut0 + N*N*t;
		const double *const restrict Gd0t_t = Gd0t + N*N*t;
		const double *const restrict Gdtt_t = Gdtt + N*N*t;
		const double *const restrict Gdt0_t = Gdt0 + N*N*t;
	for (int c = 0; c < num_b; c++) {
		const int j0 = p->bonds[c];
		const int j1 = p->bonds[c + num_b];
	for (int b = 0; b < num_b; b++) {
		const int i0 = p->bonds[b];
		const int i1 = p->bonds[b + num_b];
		const int bb = p->map_bb[b + c*num_b];
		const double pre = (double)sign / p->degen_bb[bb];
		const double gui0i0 = Gutt_t[i0 + i0*N];
		const double gui1i0 = Gutt_t[i1 + i0*N];
		const double gui0i1 = Gutt_t[i0 + i1*N];
		const double gui1i1 = Gutt_t[i1 + i1*N];
		const double gui0j0 = Gut0_t[i0 + j0*N];
		const double gui1j0 = Gut0_t[i1 + j0*N];
		const double gui0j1 = Gut0_t[i0 + j1*N];
		const double gui1j1 = Gut0_t[i1 + j1*N];
		const double guj0i0 = Gu0t_t[j0 + i0*N];
		const double guj1i0 = Gu0t_t[j1 + i0*N];
		const double guj0i1 = Gu0t_t[j0 + i1*N];
		const double guj1i1 = Gu0t_t[j1 + i1*N];
		const double guj0j0 = Gu00[j0 + j0*N];
		const double guj1j0 = Gu00[j1 + j0*N];
		const double guj0j1 = Gu00[j0 + j1*N];
		const double guj1j1 = Gu00[j1 + j1*N];
		const double gdi0i0 = Gdtt_t[i0 + i0*N];
		const double gdi1i0 = Gdtt_t[i1 + i0*N];
		const double gdi0i1 = Gdtt_t[i0 + i1*N];
		const double gdi1i1 = Gdtt_t[i1 + i1*N];
		const double gdi0j0 = Gdt0_t[i0 + j0*N];
		const double gdi1j0 = Gdt0_t[i1 + j0*N];
		const double gdi0j1 = Gdt0_t[i0 + j1*N];
		const double gdi1j1 = Gdt0_t[i1 + j1*N];
		const double gdj0i0 = Gd0t_t[j0 + i0*N];
		const double gdj1i0 = Gd0t_t[j1 + i0*N];
		const double gdj0i1 = Gd0t_t[j0 + i1*N];
		const double gdj1i1 = Gd0t_t[j1 + i1*N];
		const double gdj0j0 = Gd00[j0 + j0*N];
		const double gdj1j0 = Gd00[j1 + j0*N];
		const double gdj0j1 = Gd00[j0 + j1*N];
		const double gdj1j1 = Gd00[j1 + j1*N];
                if (meas_bond_corr){
		m->pair_bb[bb + num_bb*t] += 0.5*pre*(gui0j0*gdi1j1 + gui1j0*gdi0j1 + gui0j1*gdi1j0 + gui1j1*gdi0j0);
		const double x = -guj1i0*gui1j0 - guj0i1*gui0j1 - gdj1i0*gdi1j0 - gdj0i1*gdi0j1;
		const double y = -guj0i0*gui1j1 - guj1i1*gui0j0 - gdj0i0*gdi1j1 - gdj1i1*gdi0j0;
		m->jj[bb + num_bb*t]   += pre*((gui0i1 - gui1i0 + gdi0i1 - gdi1i0)*(guj0j1 - guj1j0 + gdj0j1 - gdj1j0) + x - y);
		m->jsjs[bb + num_bb*t] += pre*((gui0i1 - gui1i0 - gdi0i1 + gdi1i0)*(guj0j1 - guj1j0 - gdj0j1 + gdj1j0) + x - y);
		m->kk[bb + num_bb*t]   += pre*((gui0i1 + gui1i0 + gdi0i1 + gdi1i0)*(guj0j1 + guj1j0 + gdj0j1 + gdj1j0) + x + y);
		m->ksks[bb + num_bb*t] += pre*((gui0i1 + gui1i0 - gdi0i1 - gdi1i0)*(guj0j1 + guj1j0 - gdj0j1 - gdj1j0) + x + y);
                }
                if (meas_thermal){
                const double ujinj = (gui1j0*(-guj0i0)-gui0j0*(-guj0i1)+gui1j1*(-guj1i0)-(gui0j1*(-guj1i1))
                              +(gui0i1-gui1i0)*(2-guj0j0-guj1j1));
                const double djinj = (gdi1j0*(-gdj0i0)-gdi0j0*(-gdj0i1)+gdi1j1*(-gdj1i0)-(gdi0j1*(-gdj1i1))
                              +(gdi0i1-gdi1i0)*(2-gdj0j0-gdj1j1));
                const double ujijj = ((-guj1i0)*gui1j0-(-guj0i0)*gui1j1-(-guj1i1)*gui0j0+(-guj0i1)*gui0j1
                               +(gui1i0-gui0i1)*(guj1j0-guj0j1));
                const double djijj = ((-gdj1i0)*gdi1j0-(-gdj0i0)*gdi1j1-(-gdj1i1)*gdi0j0+(-gdj0i1)*gdi0j1
                               +(gdi1i0-gdi0i1)*(gdj1j0-gdj0j1));
                const double unijj = (gui1j0*(-guj1i1)-gui0j1*(-guj0i0)+gui0j0*(-guj1i0)-gui1j1*(-guj0i1)
                              +(guj0j1-guj1j0)*(2-gui0i0-gui1i1));
                const double dnijj = (gdi1j0*(-gdj1i1)-gdi0j1*(-gdj0i0)+gdi0j0*(-gdj1i0)-gdi1j1*(-gdj0i1)
                              +(gdj0j1-gdj1j0)*(2-gdi0i0-gdi1i1));
		
 		const double part1u = ujijj*((2-gdi0i0-gdi1i1)*(2-gdj0j0-gdj1j1)+gdi0j0*(-gdj0i0)
                         +gdi1j0*(-gdj0i1) +gdi0j1*(-gdj1i0) +gdi1j1*(-gdj1i1) );
    
                const double part2u = dnijj*ujinj;
    
                const double part1d = djijj*((2-gui0i0-gui1i1)*(2-guj0j0-guj1j1)+gui0j0*(-guj0i0)
                         +gui1j0*(-guj0i1) +gui0j1*(-guj1i0) +gui1j1*(-guj1i1) );
    
                const double part2d = unijj*djinj;
                        
                m->jnjn[bb + num_bb*t]   += pre*(part1u+part1d+part2u+part2d);
                }
	}
	}
	}
	
        if (meas_thermal)
	#pragma omp parallel for
	for (int t = 1; t < L; t++) {
		const double *const restrict Gu0t_t = Gu0t + N*N*t;
		const double *const restrict Gutt_t = Gutt + N*N*t;
		const double *const restrict Gut0_t = Gut0 + N*N*t;
		const double *const restrict Gd0t_t = Gd0t + N*N*t;
		const double *const restrict Gdtt_t = Gdtt + N*N*t;
		const double *const restrict Gdt0_t = Gdt0 + N*N*t;
	for (int c = 0; c < num_b; c++) {
		const int j0 = p->bonds[c];
		const int j1 = p->bonds[c + num_b];
	for (int b = 0; b < num_b2; b++) {
		const int i0 = p->bond2s[b];
		const int i1 = p->bond2s[b + num_b2];
		const int bb2 = p->map_bb2[b + c*num_b2];
		const double pre = (double)sign / p->degen_bb2[bb2];
		const double gui0i0 = Gutt_t[i0 + i0*N];
		const double gui1i0 = Gutt_t[i1 + i0*N];
		const double gui0i1 = Gutt_t[i0 + i1*N];
		const double gui1i1 = Gutt_t[i1 + i1*N];
		const double gui0j0 = Gut0_t[i0 + j0*N];
		const double gui1j0 = Gut0_t[i1 + j0*N];
		const double gui0j1 = Gut0_t[i0 + j1*N];
		const double gui1j1 = Gut0_t[i1 + j1*N];
		const double guj0i0 = Gu0t_t[j0 + i0*N];
		const double guj1i0 = Gu0t_t[j1 + i0*N];
		const double guj0i1 = Gu0t_t[j0 + i1*N];
		const double guj1i1 = Gu0t_t[j1 + i1*N];
		const double guj0j0 = Gu00[j0 + j0*N];
		const double guj1j0 = Gu00[j1 + j0*N];
		const double guj0j1 = Gu00[j0 + j1*N];
		const double guj1j1 = Gu00[j1 + j1*N];
		const double gdi0i0 = Gdtt_t[i0 + i0*N];
		const double gdi1i0 = Gdtt_t[i1 + i0*N];
		const double gdi0i1 = Gdtt_t[i0 + i1*N];
		const double gdi1i1 = Gdtt_t[i1 + i1*N];
		const double gdi0j0 = Gdt0_t[i0 + j0*N];
		const double gdi1j0 = Gdt0_t[i1 + j0*N];
		const double gdi0j1 = Gdt0_t[i0 + j1*N];
		const double gdi1j1 = Gdt0_t[i1 + j1*N];
		const double gdj0i0 = Gd0t_t[j0 + i0*N];
		const double gdj1i0 = Gd0t_t[j1 + i0*N];
		const double gdj0i1 = Gd0t_t[j0 + i1*N];
		const double gdj1i1 = Gd0t_t[j1 + i1*N];
		const double gdj0j0 = Gd00[j0 + j0*N];
		const double gdj1j0 = Gd00[j1 + j0*N];
		const double gdj0j1 = Gd00[j0 + j1*N];
		const double gdj1j1 = Gd00[j1 + j1*N];

                const double ujinj = (gui1j0*(-guj0i0)-gui0j0*(-guj0i1)+gui1j1*(-guj1i0)-(gui0j1*(-guj1i1))
                              +(gui0i1-gui1i0)*(2-guj0j0-guj1j1));
                const double djinj = (gdi1j0*(-gdj0i0)-gdi0j0*(-gdj0i1)+gdi1j1*(-gdj1i0)-(gdi0j1*(-gdj1i1))
                              +(gdi0i1-gdi1i0)*(2-gdj0j0-gdj1j1));

             
                const double part1u = (2-guj0j0-guj1j1)*((-gdj1i0)*gdi1j0-(-gdj0i0)*gdi1j1-(-gdj1i1)*gdi0j0+(-gdj0i1)*gdi0j1
                                      +(gdj1j0-gdj0j1)*(gdi1i0-gdi0i1));
    
		const double part2u = -(gdj1j0-gdj0j1)*ujinj;
    
		const double part1d = (2-gdj0j0-gdj1j1)*((-guj1i0)*gui1j0-(-guj0i0)*gui1j1-(-guj1i1)*gui0j0+(-guj0i1)*gui0j1
                                      +(guj1j0-guj0j1)*(gui1i0-gui0i1));
    
 		const double part2d = -(guj1j0-guj0j1)*djinj;

		m->jjn[bb2 + num_bb2*t]  += pre*(part1u+part1d+part2u+part2d);
	}
	}
	}

	if (meas_thermal)
	#pragma omp parallel for
	for (int t = 1; t < L; t++) {
		const double *const restrict Gu0t_t = Gu0t + N*N*t;
		const double *const restrict Gutt_t = Gutt + N*N*t;
		const double *const restrict Gut0_t = Gut0 + N*N*t;
		const double *const restrict Gd0t_t = Gd0t + N*N*t;
		const double *const restrict Gdtt_t = Gdtt + N*N*t;
		const double *const restrict Gdt0_t = Gdt0 + N*N*t;
	for (int c = 0; c < num_b2; c++) {
		const int j0 = p->bond2s[c];
		const int j1 = p->bond2s[c + num_b2];
	for (int b = 0; b < num_b; b++) {
		const int i0 = p->bonds[b];
		const int i1 = p->bonds[b + num_b];
		const int b2b = p->map_b2b[b + c*num_b];
		const double pre = (double)sign / p->degen_b2b[b2b];
		const double gui0i0 = Gutt_t[i0 + i0*N];
		const double gui1i0 = Gutt_t[i1 + i0*N];
		const double gui0i1 = Gutt_t[i0 + i1*N];
		const double gui1i1 = Gutt_t[i1 + i1*N];
		const double gui0j0 = Gut0_t[i0 + j0*N];
		const double gui1j0 = Gut0_t[i1 + j0*N];
		const double gui0j1 = Gut0_t[i0 + j1*N];
		const double gui1j1 = Gut0_t[i1 + j1*N];
		const double guj0i0 = Gu0t_t[j0 + i0*N];
		const double guj1i0 = Gu0t_t[j1 + i0*N];
		const double guj0i1 = Gu0t_t[j0 + i1*N];
		const double guj1i1 = Gu0t_t[j1 + i1*N];
		const double guj0j0 = Gu00[j0 + j0*N];
		const double guj1j0 = Gu00[j1 + j0*N];
		const double guj0j1 = Gu00[j0 + j1*N];
		const double guj1j1 = Gu00[j1 + j1*N];
		const double gdi0i0 = Gdtt_t[i0 + i0*N];
		const double gdi1i0 = Gdtt_t[i1 + i0*N];
		const double gdi0i1 = Gdtt_t[i0 + i1*N];
		const double gdi1i1 = Gdtt_t[i1 + i1*N];
		const double gdi0j0 = Gdt0_t[i0 + j0*N];
		const double gdi1j0 = Gdt0_t[i1 + j0*N];
		const double gdi0j1 = Gdt0_t[i0 + j1*N];
		const double gdi1j1 = Gdt0_t[i1 + j1*N];
		const double gdj0i0 = Gd0t_t[j0 + i0*N];
		const double gdj1i0 = Gd0t_t[j1 + i0*N];
		const double gdj0i1 = Gd0t_t[j0 + i1*N];
		const double gdj1i1 = Gd0t_t[j1 + i1*N];
		const double gdj0j0 = Gd00[j0 + j0*N];
		const double gdj1j0 = Gd00[j1 + j0*N];
		const double gdj0j1 = Gd00[j0 + j1*N];
		const double gdj1j1 = Gd00[j1 + j1*N];

                const double ujijj = ((-guj1i0)*gui1j0-(-guj0i0)*gui1j1-(-guj1i1)*gui0j0+(-guj0i1)*gui0j1
                               +(gui1i0-gui0i1)*(guj1j0-guj0j1));
                const double djijj = ((-gdj1i0)*gdi1j0-(-gdj0i0)*gdi1j1-(-gdj1i1)*gdi0j0+(-gdj0i1)*gdi0j1
                               +(gdi1i0-gdi0i1)*(gdj1j0-gdj0j1));
                const double unijj = (gui1j0*(-guj1i1)-gui0j1*(-guj0i0)+gui0j0*(-guj1i0)-gui1j1*(-guj0i1)
                              +(guj0j1-guj1j0)*(2-gui0i0-gui1i1));
                const double dnijj = (gdi1j0*(-gdj1i1)-gdi0j1*(-gdj0i0)+gdi0j0*(-gdj1i0)-gdi1j1*(-gdj0i1)
                              +(gdj0j1-gdj1j0)*(2-gdi0i0-gdi1i1));
                
 		const double part1u = (2-gui0i0-gui1i1)*djijj;
    
                const double part2u = -(gdi1i0-gdi0i1)*unijj;
    
                const double part1d = (2-gdi0i0-gdi1i1)*ujijj;
    
                const double part2d = -(gui1i0-gui0i1)*dnijj;

		m->jnj[b2b + num_b2b*t] += pre*(part1u+part1d+part2u+part2d);
	}
	}
	}
	
    // no delta functions here.
	if (meas_2bond_corr)
	#pragma omp parallel for
	for (int t = 1; t < L; t++) {
		const double *const restrict Gu0t_t = Gu0t + N*N*t;
		const double *const restrict Gutt_t = Gutt + N*N*t;
		const double *const restrict Gut0_t = Gut0 + N*N*t;
		const double *const restrict Gd0t_t = Gd0t + N*N*t;
		const double *const restrict Gdtt_t = Gdtt + N*N*t;
		const double *const restrict Gdt0_t = Gdt0 + N*N*t;
	for (int c = 0; c < num_b2; c++) {
		const int j0 = p->bond2s[c];
		const int j1 = p->bond2s[c + num_b2];
	for (int b = 0; b < num_b2; b++) {
		const int i0 = p->bond2s[b];
		const int i1 = p->bond2s[b + num_b2];
		const int bb = p->map_b2b2[b + c*num_b2];
		const double pre = (double)sign / p->degen_b2b2[bb];
		const double gui0i0 = Gutt_t[i0 + i0*N];
		const double gui1i0 = Gutt_t[i1 + i0*N];
		const double gui0i1 = Gutt_t[i0 + i1*N];
		const double gui1i1 = Gutt_t[i1 + i1*N];
		const double gui0j0 = Gut0_t[i0 + j0*N];
		const double gui1j0 = Gut0_t[i1 + j0*N];
		const double gui0j1 = Gut0_t[i0 + j1*N];
		const double gui1j1 = Gut0_t[i1 + j1*N];
		const double guj0i0 = Gu0t_t[j0 + i0*N];
		const double guj1i0 = Gu0t_t[j1 + i0*N];
		const double guj0i1 = Gu0t_t[j0 + i1*N];
		const double guj1i1 = Gu0t_t[j1 + i1*N];
		const double guj0j0 = Gu00[j0 + j0*N];
		const double guj1j0 = Gu00[j1 + j0*N];
		const double guj0j1 = Gu00[j0 + j1*N];
		const double guj1j1 = Gu00[j1 + j1*N];
		const double gdi0i0 = Gdtt_t[i0 + i0*N];
		const double gdi1i0 = Gdtt_t[i1 + i0*N];
		const double gdi0i1 = Gdtt_t[i0 + i1*N];
		const double gdi1i1 = Gdtt_t[i1 + i1*N];
		const double gdi0j0 = Gdt0_t[i0 + j0*N];
		const double gdi1j0 = Gdt0_t[i1 + j0*N];
		const double gdi0j1 = Gdt0_t[i0 + j1*N];
		const double gdi1j1 = Gdt0_t[i1 + j1*N];
		const double gdj0i0 = Gd0t_t[j0 + i0*N];
		const double gdj1i0 = Gd0t_t[j1 + i0*N];
		const double gdj0i1 = Gd0t_t[j0 + i1*N];
		const double gdj1i1 = Gd0t_t[j1 + i1*N];
		const double gdj0j0 = Gd00[j0 + j0*N];
		const double gdj1j0 = Gd00[j1 + j0*N];
		const double gdj0j1 = Gd00[j0 + j1*N];
		const double gdj1j1 = Gd00[j1 + j1*N];
		m->pair_b2b2[bb + num_b2b2*t] += 0.5*pre*(gui0j0*gdi1j1 + gui1j0*gdi0j1 + gui0j1*gdi1j0 + gui1j1*gdi0j0);
		const double x = -guj1i0*gui1j0 - guj0i1*gui0j1 - gdj1i0*gdi1j0 - gdj0i1*gdi0j1;
		const double y = -guj0i0*gui1j1 - guj1i1*gui0j0 - gdj0i0*gdi1j1 - gdj1i1*gdi0j0;
		m->j2j2[bb + num_b2b2*t]   += pre*((gui0i1 - gui1i0 + gdi0i1 - gdi1i0)*(guj0j1 - guj1j0 + gdj0j1 - gdj1j0) + x - y);
		m->js2js2[bb + num_b2b2*t] += pre*((gui0i1 - gui1i0 - gdi0i1 + gdi1i0)*(guj0j1 - guj1j0 - gdj0j1 + gdj1j0) + x - y);
		m->k2k2[bb + num_b2b2*t]   += pre*((gui0i1 + gui1i0 + gdi0i1 + gdi1i0)*(guj0j1 + guj1j0 + gdj0j1 + gdj1j0) + x + y);
		m->ks2ks2[bb + num_b2b2*t] += pre*((gui0i1 + gui1i0 - gdi0i1 - gdi1i0)*(guj0j1 + guj1j0 - gdj0j1 - gdj1j0) + x + y);
	}
	}
	}


	if (meas_nematic_corr)
	#pragma omp parallel for
	for (int t = 1; t < L; t++) {
		const double *const restrict Gu0t_t = Gu0t + N*N*t;
		const double *const restrict Gutt_t = Gutt + N*N*t;
		const double *const restrict Gut0_t = Gut0 + N*N*t;
		const double *const restrict Gd0t_t = Gd0t + N*N*t;
		const double *const restrict Gdtt_t = Gdtt + N*N*t;
		const double *const restrict Gdt0_t = Gdt0 + N*N*t;
	for (int c = 0; c < NEM_BONDS*N; c++) {
		const int j0 = p->bonds[c];
		const int j1 = p->bonds[c + num_b];
	for (int b = 0; b < NEM_BONDS*N; b++) {
		const int i0 = p->bonds[b];
		const int i1 = p->bonds[b + num_b];
		const int bb = p->map_bb[b + c*num_b];
		const double pre = (double)sign / p->degen_bb[bb];
		const double gui0i0 = Gutt_t[i0 + i0*N];
		const double gui1i0 = Gutt_t[i1 + i0*N];
		const double gui0i1 = Gutt_t[i0 + i1*N];
		const double gui1i1 = Gutt_t[i1 + i1*N];
		const double gui0j0 = Gut0_t[i0 + j0*N];
		const double gui1j0 = Gut0_t[i1 + j0*N];
		const double gui0j1 = Gut0_t[i0 + j1*N];
		const double gui1j1 = Gut0_t[i1 + j1*N];
		const double guj0i0 = Gu0t_t[j0 + i0*N];
		const double guj1i0 = Gu0t_t[j1 + i0*N];
		const double guj0i1 = Gu0t_t[j0 + i1*N];
		const double guj1i1 = Gu0t_t[j1 + i1*N];
		const double guj0j0 = Gu00[j0 + j0*N];
		const double guj1j0 = Gu00[j1 + j0*N];
		const double guj0j1 = Gu00[j0 + j1*N];
		const double guj1j1 = Gu00[j1 + j1*N];
		const double gdi0i0 = Gdtt_t[i0 + i0*N];
		const double gdi1i0 = Gdtt_t[i1 + i0*N];
		const double gdi0i1 = Gdtt_t[i0 + i1*N];
		const double gdi1i1 = Gdtt_t[i1 + i1*N];
		const double gdi0j0 = Gdt0_t[i0 + j0*N];
		const double gdi1j0 = Gdt0_t[i1 + j0*N];
		const double gdi0j1 = Gdt0_t[i0 + j1*N];
		const double gdi1j1 = Gdt0_t[i1 + j1*N];
		const double gdj0i0 = Gd0t_t[j0 + i0*N];
		const double gdj1i0 = Gd0t_t[j1 + i0*N];
		const double gdj0i1 = Gd0t_t[j0 + i1*N];
		const double gdj1i1 = Gd0t_t[j1 + i1*N];
		const double gdj0j0 = Gd00[j0 + j0*N];
		const double gdj1j0 = Gd00[j1 + j0*N];
		const double gdj0j1 = Gd00[j0 + j1*N];
		const double gdj1j1 = Gd00[j1 + j1*N];
		const int delta_i0i1 = 0;
		const int delta_j0j1 = 0;
		const int delta_i0j0 = 0;
		const int delta_i1j0 = 0;
		const int delta_i0j1 = 0;
		const int delta_i1j1 = 0;
		const double uuuu = +(1.-gui0i0)*(1.-gui1i1)*(1.-guj0j0)*(1.-guj1j1)+(1.-gui0i0)*(1.-gui1i1)*(delta_j0j1-guj1j0)*guj0j1+(1.-gui0i0)*(delta_i1j0-guj0i1)*gui1j0*(1.-guj1j1)-(1.-gui0i0)*(delta_i1j0-guj0i1)*gui1j1*(delta_j0j1-guj1j0)+(1.-gui0i0)*(delta_i1j1-guj1i1)*gui1j0*guj0j1+(1.-gui0i0)*(delta_i1j1-guj1i1)*gui1j1*(1.-guj0j0)+(delta_i0i1-gui1i0)*gui0i1*(1.-guj0j0)*(1.-guj1j1)+(delta_i0i1-gui1i0)*gui0i1*(delta_j0j1-guj1j0)*guj0j1-(delta_i0i1-gui1i0)*gui0j0*(delta_i1j0-guj0i1)*(1.-guj1j1)-(delta_i0i1-gui1i0)*gui0j0*(delta_i1j1-guj1i1)*guj0j1+(delta_i0i1-gui1i0)*gui0j1*(delta_i1j0-guj0i1)*(delta_j0j1-guj1j0)-(delta_i0i1-gui1i0)*gui0j1*(delta_i1j1-guj1i1)*(1.-guj0j0)+(delta_i0j0-guj0i0)*gui0i1*gui1j0*(1.-guj1j1)-(delta_i0j0-guj0i0)*gui0i1*gui1j1*(delta_j0j1-guj1j0)+(delta_i0j0-guj0i0)*gui0j0*(1.-gui1i1)*(1.-guj1j1)+(delta_i0j0-guj0i0)*gui0j0*(delta_i1j1-guj1i1)*gui1j1-(delta_i0j0-guj0i0)*gui0j1*(1.-gui1i1)*(delta_j0j1-guj1j0)-(delta_i0j0-guj0i0)*gui0j1*(delta_i1j1-guj1i1)*gui1j0+(delta_i0j1-guj1i0)*gui0i1*gui1j0*guj0j1+(delta_i0j1-guj1i0)*gui0i1*gui1j1*(1.-guj0j0)+(delta_i0j1-guj1i0)*gui0j0*(1.-gui1i1)*guj0j1-(delta_i0j1-guj1i0)*gui0j0*(delta_i1j0-guj0i1)*gui1j1+(delta_i0j1-guj1i0)*gui0j1*(1.-gui1i1)*(1.-guj0j0)+(delta_i0j1-guj1i0)*gui0j1*(delta_i1j0-guj0i1)*gui1j0;
		const double uuud = +(1.-gui0i0)*(1.-gui1i1)*(1.-guj0j0)*(1.-gdj1j1)+(1.-gui0i0)*(delta_i1j0-guj0i1)*gui1j0*(1.-gdj1j1)+(delta_i0i1-gui1i0)*gui0i1*(1.-guj0j0)*(1.-gdj1j1)-(delta_i0i1-gui1i0)*gui0j0*(delta_i1j0-guj0i1)*(1.-gdj1j1)+(delta_i0j0-guj0i0)*gui0i1*gui1j0*(1.-gdj1j1)+(delta_i0j0-guj0i0)*gui0j0*(1.-gui1i1)*(1.-gdj1j1);
		const double uudu = +(1.-gui0i0)*(1.-gui1i1)*(1.-gdj0j0)*(1.-guj1j1)+(1.-gui0i0)*(delta_i1j1-guj1i1)*gui1j1*(1.-gdj0j0)+(delta_i0i1-gui1i0)*gui0i1*(1.-gdj0j0)*(1.-guj1j1)-(delta_i0i1-gui1i0)*gui0j1*(delta_i1j1-guj1i1)*(1.-gdj0j0)+(delta_i0j1-guj1i0)*gui0i1*gui1j1*(1.-gdj0j0)+(delta_i0j1-guj1i0)*gui0j1*(1.-gui1i1)*(1.-gdj0j0);
		const double uudd = +(1.-gui0i0)*(1.-gui1i1)*(1.-gdj0j0)*(1.-gdj1j1)+(1.-gui0i0)*(1.-gui1i1)*(delta_j0j1-gdj1j0)*gdj0j1+(delta_i0i1-gui1i0)*gui0i1*(1.-gdj0j0)*(1.-gdj1j1)+(delta_i0i1-gui1i0)*gui0i1*(delta_j0j1-gdj1j0)*gdj0j1;
		const double uduu = +(1.-gui0i0)*(1.-gdi1i1)*(1.-guj0j0)*(1.-guj1j1)+(1.-gui0i0)*(1.-gdi1i1)*(delta_j0j1-guj1j0)*guj0j1+(delta_i0j0-guj0i0)*gui0j0*(1.-gdi1i1)*(1.-guj1j1)-(delta_i0j0-guj0i0)*gui0j1*(1.-gdi1i1)*(delta_j0j1-guj1j0)+(delta_i0j1-guj1i0)*gui0j0*(1.-gdi1i1)*guj0j1+(delta_i0j1-guj1i0)*gui0j1*(1.-gdi1i1)*(1.-guj0j0);
		const double udud = +(1.-gui0i0)*(1.-gdi1i1)*(1.-guj0j0)*(1.-gdj1j1)+(1.-gui0i0)*(delta_i1j1-gdj1i1)*gdi1j1*(1.-guj0j0)+(delta_i0j0-guj0i0)*gui0j0*(1.-gdi1i1)*(1.-gdj1j1)+(delta_i0j0-guj0i0)*gui0j0*(delta_i1j1-gdj1i1)*gdi1j1;
		const double uddu = +(1.-gui0i0)*(1.-gdi1i1)*(1.-gdj0j0)*(1.-guj1j1)+(1.-gui0i0)*(delta_i1j0-gdj0i1)*gdi1j0*(1.-guj1j1)+(delta_i0j1-guj1i0)*gui0j1*(1.-gdi1i1)*(1.-gdj0j0)+(delta_i0j1-guj1i0)*gui0j1*(delta_i1j0-gdj0i1)*gdi1j0;
		const double uddd = +(1.-gui0i0)*(1.-gdi1i1)*(1.-gdj0j0)*(1.-gdj1j1)+(1.-gui0i0)*(1.-gdi1i1)*(delta_j0j1-gdj1j0)*gdj0j1+(1.-gui0i0)*(delta_i1j0-gdj0i1)*gdi1j0*(1.-gdj1j1)-(1.-gui0i0)*(delta_i1j0-gdj0i1)*gdi1j1*(delta_j0j1-gdj1j0)+(1.-gui0i0)*(delta_i1j1-gdj1i1)*gdi1j0*gdj0j1+(1.-gui0i0)*(delta_i1j1-gdj1i1)*gdi1j1*(1.-gdj0j0);
		const double duuu = +(1.-gdi0i0)*(1.-gui1i1)*(1.-guj0j0)*(1.-guj1j1)+(1.-gdi0i0)*(1.-gui1i1)*(delta_j0j1-guj1j0)*guj0j1+(1.-gdi0i0)*(delta_i1j0-guj0i1)*gui1j0*(1.-guj1j1)-(1.-gdi0i0)*(delta_i1j0-guj0i1)*gui1j1*(delta_j0j1-guj1j0)+(1.-gdi0i0)*(delta_i1j1-guj1i1)*gui1j0*guj0j1+(1.-gdi0i0)*(delta_i1j1-guj1i1)*gui1j1*(1.-guj0j0);
		const double duud = +(1.-gdi0i0)*(1.-gui1i1)*(1.-guj0j0)*(1.-gdj1j1)+(1.-gdi0i0)*(delta_i1j0-guj0i1)*gui1j0*(1.-gdj1j1)+(delta_i0j1-gdj1i0)*gdi0j1*(1.-gui1i1)*(1.-guj0j0)+(delta_i0j1-gdj1i0)*gdi0j1*(delta_i1j0-guj0i1)*gui1j0;
		const double dudu = +(1.-gdi0i0)*(1.-gui1i1)*(1.-gdj0j0)*(1.-guj1j1)+(1.-gdi0i0)*(delta_i1j1-guj1i1)*gui1j1*(1.-gdj0j0)+(delta_i0j0-gdj0i0)*gdi0j0*(1.-gui1i1)*(1.-guj1j1)+(delta_i0j0-gdj0i0)*gdi0j0*(delta_i1j1-guj1i1)*gui1j1;
		const double dudd = +(1.-gdi0i0)*(1.-gui1i1)*(1.-gdj0j0)*(1.-gdj1j1)+(1.-gdi0i0)*(1.-gui1i1)*(delta_j0j1-gdj1j0)*gdj0j1+(delta_i0j0-gdj0i0)*gdi0j0*(1.-gui1i1)*(1.-gdj1j1)-(delta_i0j0-gdj0i0)*gdi0j1*(1.-gui1i1)*(delta_j0j1-gdj1j0)+(delta_i0j1-gdj1i0)*gdi0j0*(1.-gui1i1)*gdj0j1+(delta_i0j1-gdj1i0)*gdi0j1*(1.-gui1i1)*(1.-gdj0j0);
		const double dduu = +(1.-gdi0i0)*(1.-gdi1i1)*(1.-guj0j0)*(1.-guj1j1)+(1.-gdi0i0)*(1.-gdi1i1)*(delta_j0j1-guj1j0)*guj0j1+(delta_i0i1-gdi1i0)*gdi0i1*(1.-guj0j0)*(1.-guj1j1)+(delta_i0i1-gdi1i0)*gdi0i1*(delta_j0j1-guj1j0)*guj0j1;
		const double ddud = +(1.-gdi0i0)*(1.-gdi1i1)*(1.-guj0j0)*(1.-gdj1j1)+(1.-gdi0i0)*(delta_i1j1-gdj1i1)*gdi1j1*(1.-guj0j0)+(delta_i0i1-gdi1i0)*gdi0i1*(1.-guj0j0)*(1.-gdj1j1)-(delta_i0i1-gdi1i0)*gdi0j1*(delta_i1j1-gdj1i1)*(1.-guj0j0)+(delta_i0j1-gdj1i0)*gdi0i1*gdi1j1*(1.-guj0j0)+(delta_i0j1-gdj1i0)*gdi0j1*(1.-gdi1i1)*(1.-guj0j0);
		const double dddu = +(1.-gdi0i0)*(1.-gdi1i1)*(1.-gdj0j0)*(1.-guj1j1)+(1.-gdi0i0)*(delta_i1j0-gdj0i1)*gdi1j0*(1.-guj1j1)+(delta_i0i1-gdi1i0)*gdi0i1*(1.-gdj0j0)*(1.-guj1j1)-(delta_i0i1-gdi1i0)*gdi0j0*(delta_i1j0-gdj0i1)*(1.-guj1j1)+(delta_i0j0-gdj0i0)*gdi0i1*gdi1j0*(1.-guj1j1)+(delta_i0j0-gdj0i0)*gdi0j0*(1.-gdi1i1)*(1.-guj1j1);
		const double dddd = +(1.-gdi0i0)*(1.-gdi1i1)*(1.-gdj0j0)*(1.-gdj1j1)+(1.-gdi0i0)*(1.-gdi1i1)*(delta_j0j1-gdj1j0)*gdj0j1+(1.-gdi0i0)*(delta_i1j0-gdj0i1)*gdi1j0*(1.-gdj1j1)-(1.-gdi0i0)*(delta_i1j0-gdj0i1)*gdi1j1*(delta_j0j1-gdj1j0)+(1.-gdi0i0)*(delta_i1j1-gdj1i1)*gdi1j0*gdj0j1+(1.-gdi0i0)*(delta_i1j1-gdj1i1)*gdi1j1*(1.-gdj0j0)+(delta_i0i1-gdi1i0)*gdi0i1*(1.-gdj0j0)*(1.-gdj1j1)+(delta_i0i1-gdi1i0)*gdi0i1*(delta_j0j1-gdj1j0)*gdj0j1-(delta_i0i1-gdi1i0)*gdi0j0*(delta_i1j0-gdj0i1)*(1.-gdj1j1)-(delta_i0i1-gdi1i0)*gdi0j0*(delta_i1j1-gdj1i1)*gdj0j1+(delta_i0i1-gdi1i0)*gdi0j1*(delta_i1j0-gdj0i1)*(delta_j0j1-gdj1j0)-(delta_i0i1-gdi1i0)*gdi0j1*(delta_i1j1-gdj1i1)*(1.-gdj0j0)+(delta_i0j0-gdj0i0)*gdi0i1*gdi1j0*(1.-gdj1j1)-(delta_i0j0-gdj0i0)*gdi0i1*gdi1j1*(delta_j0j1-gdj1j0)+(delta_i0j0-gdj0i0)*gdi0j0*(1.-gdi1i1)*(1.-gdj1j1)+(delta_i0j0-gdj0i0)*gdi0j0*(delta_i1j1-gdj1i1)*gdi1j1-(delta_i0j0-gdj0i0)*gdi0j1*(1.-gdi1i1)*(delta_j0j1-gdj1j0)-(delta_i0j0-gdj0i0)*gdi0j1*(delta_i1j1-gdj1i1)*gdi1j0+(delta_i0j1-gdj1i0)*gdi0i1*gdi1j0*gdj0j1+(delta_i0j1-gdj1i0)*gdi0i1*gdi1j1*(1.-gdj0j0)+(delta_i0j1-gdj1i0)*gdi0j0*(1.-gdi1i1)*gdj0j1-(delta_i0j1-gdj1i0)*gdi0j0*(delta_i1j0-gdj0i1)*gdi1j1+(delta_i0j1-gdj1i0)*gdi0j1*(1.-gdi1i1)*(1.-gdj0j0)+(delta_i0j1-gdj1i0)*gdi0j1*(delta_i1j0-gdj0i1)*gdi1j0;
		m->nem_nnnn[bb + num_bb*t] += pre*(uuuu + uuud + uudu + uudd
			                         + uduu + udud + uddu + uddd
			                         + duuu + duud + dudu + dudd
			                         + dduu + ddud + dddu + dddd);
		m->nem_ssss[bb + num_bb*t] += pre*(uuuu - uuud - uudu + uudd
			                         - uduu + udud + uddu - uddd
			                         - duuu + duud + dudu - dudd
			                         + dduu - ddud - dddu + dddd);
	}
	}
	}
}

// unmaintained, outdated
void measure_uneqlt_full(const struct params *const restrict p, const int sign,
		const double *const restrict Gu,
		const double *const restrict Gd,
		struct meas_uneqlt *const restrict m)
{
	m->n_sample++;
	m->sign += sign;
	const int N = p->N, L = p->L, num_i = p->num_i, num_ij = p->num_ij;
	const int num_b = p->num_b, num_bb = p->num_bb;

	// 2-site measurements
	#pragma omp parallel for
	for (int t = 0; t < L; t++)
	for (int l = 0; l < L; l++) {
		const int k = (l + t) % L;
		const int T_sign = ((k >= l) ? 1.0 : -1.0); // for fermionic
		const int delta_t = (t == 0);
	for (int j = 0; j < N; j++)
	for (int i = 0; i < N; i++) {
		const int r = p->map_ij[i + j*N];
		const int delta_tij = delta_t * (i == j);
		const double pre = (double)sign / p->degen_ij[r] / L;
		const double guii = Gu[(i + N*i) + N*N*(k + L*k)];
		const double guij = Gu[(i + N*j) + N*N*(k + L*l)];
		const double guji = Gu[(j + N*i) + N*N*(l + L*k)];
		const double gujj = Gu[(j + N*j) + N*N*(l + L*l)];
		const double gdii = Gd[(i + N*i) + N*N*(k + L*k)];
		const double gdij = Gd[(i + N*j) + N*N*(k + L*l)];
		const double gdji = Gd[(j + N*i) + N*N*(l + L*k)];
		const double gdjj = Gd[(j + N*j) + N*N*(l + L*l)];
		m->gt0[r + num_ij*t] += 0.5*T_sign*pre*(guij + gdij);
		const double x = delta_tij*(guii + gdii) - (guji*guij + gdji*gdij);
		m->nn[r + num_ij*t] += pre*((2. - guii - gdii)*(2. - gujj - gdjj) + x);
		m->xx[r + num_ij*t] += 0.25*pre*(delta_tij*(guii + gdii) - (guji*gdij + gdji*guij));
		m->zz[r + num_ij*t] += 0.25*pre*((gdii - guii)*(gdjj - gujj) + x);
		m->pair_sw[r + num_ij*t] += pre*guij*gdij;
	}
	}

	// 2-bond measurements
	#pragma omp parallel for
	for (int t = 0; t < L; t++)
	for (int l = 0; l < L; l++) {
		const int k = (l + t) % L;
		const int delta_t = (t == 0);
	for (int c = 0; c < num_b; c++) {
		const int j0 = p->bonds[c];
		const int j1 = p->bonds[c + num_b];
	for (int b = 0; b < num_b; b++) {
		const int i0 = p->bonds[b];
		const int i1 = p->bonds[b + num_b];
		const int bb = p->map_bb[b + c*num_b];
		const double pre = (double)sign / p->degen_bb[bb] / L;
		const int delta_ti0j0 = delta_t * (i0 == j0);
		const int delta_ti1j0 = delta_t * (i1 == j0);
		const int delta_ti0j1 = delta_t * (i0 == j1);
		const int delta_ti1j1 = delta_t * (i1 == j1);
		const double gui0i0 = Gu[i0 + i0*N + N*N*(k + L*k)];
		const double gui1i0 = Gu[i1 + i0*N + N*N*(k + L*k)];
		const double gui0i1 = Gu[i0 + i1*N + N*N*(k + L*k)];
		const double gui1i1 = Gu[i1 + i1*N + N*N*(k + L*k)];
		const double gui0j0 = Gu[i0 + j0*N + N*N*(k + L*l)];
		const double gui1j0 = Gu[i1 + j0*N + N*N*(k + L*l)];
		const double gui0j1 = Gu[i0 + j1*N + N*N*(k + L*l)];
		const double gui1j1 = Gu[i1 + j1*N + N*N*(k + L*l)];
		const double guj0i0 = Gu[j0 + i0*N + N*N*(l + L*k)];
		const double guj1i0 = Gu[j1 + i0*N + N*N*(l + L*k)];
		const double guj0i1 = Gu[j0 + i1*N + N*N*(l + L*k)];
		const double guj1i1 = Gu[j1 + i1*N + N*N*(l + L*k)];
		const double guj0j0 = Gu[j0 + j0*N + N*N*(l + L*l)];
		const double guj1j0 = Gu[j1 + j0*N + N*N*(l + L*l)];
		const double guj0j1 = Gu[j0 + j1*N + N*N*(l + L*l)];
		const double guj1j1 = Gu[j1 + j1*N + N*N*(l + L*l)];
		const double gdi0i0 = Gd[i0 + i0*N + N*N*(k + L*k)];
		const double gdi1i0 = Gd[i1 + i0*N + N*N*(k + L*k)];
		const double gdi0i1 = Gd[i0 + i1*N + N*N*(k + L*k)];
		const double gdi1i1 = Gd[i1 + i1*N + N*N*(k + L*k)];
		const double gdi0j0 = Gd[i0 + j0*N + N*N*(k + L*l)];
		const double gdi1j0 = Gd[i1 + j0*N + N*N*(k + L*l)];
		const double gdi0j1 = Gd[i0 + j1*N + N*N*(k + L*l)];
		const double gdi1j1 = Gd[i1 + j1*N + N*N*(k + L*l)];
		const double gdj0i0 = Gd[j0 + i0*N + N*N*(l + L*k)];
		const double gdj1i0 = Gd[j1 + i0*N + N*N*(l + L*k)];
		const double gdj0i1 = Gd[j0 + i1*N + N*N*(l + L*k)];
		const double gdj1i1 = Gd[j1 + i1*N + N*N*(l + L*k)];
		const double gdj0j0 = Gd[j0 + j0*N + N*N*(l + L*l)];
		const double gdj1j0 = Gd[j1 + j0*N + N*N*(l + L*l)];
		const double gdj0j1 = Gd[j0 + j1*N + N*N*(l + L*l)];
		const double gdj1j1 = Gd[j1 + j1*N + N*N*(l + L*l)];
		m->pair_bb[bb + num_bb*t] += 0.5*pre*(gui0j0*gdi1j1 + gui1j0*gdi0j1 + gui0j1*gdi1j0 + gui1j1*gdi0j0);
		m->jj[bb + num_bb*t] += pre*((gui0i1 - gui1i0 + gdi0i1 - gdi1i0)*(guj0j1 - guj1j0 + gdj0j1 - gdj1j0)
		                             + (delta_ti0j1 - guj1i0)*gui1j0 - (delta_ti0j0 - guj0i0)*gui1j1
		                             - (delta_ti1j1 - guj1i1)*gui0j0 + (delta_ti1j0 - guj0i1)*gui0j1
		                             + (delta_ti0j1 - gdj1i0)*gdi1j0 - (delta_ti0j0 - gdj0i0)*gdi1j1
		                             - (delta_ti1j1 - gdj1i1)*gdi0j0 + (delta_ti1j0 - gdj0i1)*gdi0j1);
		m->kk[bb + num_bb*t] += pre*((gui0i1 + gui1i0 + gdi0i1 + gdi1i0)*(guj0j1 + guj1j0 + gdj0j1 + gdj1j0)
		                             + (delta_ti0j1 - guj1i0)*gui1j0 + (delta_ti0j0 - guj0i0)*gui1j1
		                             + (delta_ti1j1 - guj1i1)*gui0j0 + (delta_ti1j0 - guj0i1)*gui0j1
		                             + (delta_ti0j1 - gdj1i0)*gdi1j0 + (delta_ti0j0 - gdj0i0)*gdi1j1
		                             + (delta_ti1j1 - gdj1i1)*gdi0j0 + (delta_ti1j0 - gdj0i1)*gdi0j1);
	}
	}
	}
}
