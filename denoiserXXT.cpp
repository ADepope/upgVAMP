#include <vector>
#include <algorithm>
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <numeric> // contains std::accumulate
#include <random>
#include <omp.h>
#include <fstream>
#include <cfloat>
#include "na_lut.hpp"
#include "vamp.hpp"
#include "utilities.hpp"

std::vector<double> vamp::lmmse_multAAT(std::vector<double> u, double tau, data* dataset){ // multiplying with (tau^(-1)*u + gam2^(-1)*AA^T)

    size_t phen_size = 4 * (*dataset).get_mbytes();
    if (u == std::vector<double>(N, 0.0) || u == std::vector<double>(phen_size, 0.0))
        return std::vector<double>(phen_size, 0.0);
    std::vector<double> res(phen_size, 0.0);
    std::vector<double> res_temp(M, 0.0);
    res_temp = (*dataset).ATx(u.data());
    res = (*dataset).Ax(res_temp.data());
    //for (int i = 0; i < N; i++){
    //    res[i] /= gam2;
    //    res[i] += u[i] / tau;
    //}
    for (int i = 0; i < N; i++){
        res[i] *= tau;
        res[i] += u[i] * gam2;
    }
    return res;
}

// returns A^T*u_A
std::vector<double> vamp::lmmse_denoiserAAT(std::vector<double> r2, std::vector<double> mu_CG_last, std::vector<double> p_last, double *correction, double *new_vars, data* dataset){
    
    std::vector<double> v = std::vector<double> (4*(*dataset).get_mbytes(), 0.0);
    for (int i=0; i<N; i++)
        v[i] = y[i] - z2[i];

    //if (rank == 0)
    //    std::cout << "v[0] = " << v[0] << std::endl;
    std::vector<double> u = CG_solverAAT(v, mu_CG_last, p_last, gamw, correction, new_vars, 1, dataset);

    //if (rank == 0)
    //    std::cout << "u[0] = " << u[0] << std::endl;
    std::vector<double> res = (*dataset).ATx(u.data());

    *new_vars = l2_norm2(res, 1) / Mt / *correction / *correction  - 1.0 / gam2;

    if (rank == 0)
        std::cout << "new_vars = " << *new_vars << std::endl;

    return res;
}

std::vector<double> vamp::CG_solverAAT(std::vector<double> v, std::vector<double> mu_start, std::vector<double> & p_start, double tau, double *correction, double *new_vars, int save, data* dataset){
    
    int mbytes4 = 4*(*dataset).get_mbytes();

    double scalar_fact = gam2 * gamw;
    //scalar_fact = 1;

    std::vector<double> noise(N, 0.0);
    noise = (*dataset).Ax(true_signal.data());
    for (int i=0; i<N; i++)
        noise[i] = y[i] - noise[i] / c_scale;

    std::vector<double> mu = mu_start;
    // change from CG to WS-CG
    mu = std::vector<double> (4*(*dataset).get_mbytes(), 0.0);
    std::vector<double> p = p_start; 
    std::vector<double> d; // d = Ap
    std::vector<double> r = lmmse_multAAT(mu, tau, dataset); // r = residual
    for (int i0 = 0; i0 < N; i0++)
        r[i0] = v[i0] - r[i0];

    // change from CG to WS-CG
    if (p == std::vector<double> (mbytes4, 0.0) || 1)
        for (int i0 = 0; i0 < N; i0++)
            p[i0] = v[i0];

    std::vector<double> alphas;
    std::vector<double> betas;
    std::vector<double> etas;
    double eta = (double) N / (double) Mt / gamw;
    if (rank == 0)
        std::cout << "eta = " << eta << std::endl;
    etas.push_back(eta);
    std::vector<double> ksis;
    double ksi = 0;
    ksis.push_back(ksi);
    std::vector<double> Apalpha(mbytes4, 0.0);
    std::vector<double> palpha(mbytes4, 0.0);
    double alpha, beta;

    for (int i = 0; i < CG_max_iter; i++){

        d = lmmse_multAAT(p, tau, dataset);

        //if (rank == 0)
        //    std::cout << "d[0] = " << d[0] << std::endl;

        alpha = l2_norm2(r, 0) / inner_prod(d, p, 0);

        if (rank == 0)
            std::cout << "alpha = " << alpha << std::endl;

        alphas.push_back(alpha * scalar_fact);
        
        for (int j = 0; j < N; j++)
            palpha[j] = alpha * p[j];

        std::transform (mu.begin(), mu.end(), palpha.begin(), mu.begin(), std::plus<double>());

        for (int j = 0; j < N; j++)
            Apalpha[j] = d[j] * alpha;

        beta = 1.0 / l2_norm2(r, 0); // optimize further!

        ksi = ksi + alphas.back() * l2_norm2(r, 0) / Mt;

        std::transform(r.begin(), r.end(), Apalpha.begin(), r.begin(), std::minus<double>());

        double inner_prod_rr = inner_prod(r, r, 0);
        beta *= inner_prod_rr;
        
        if (rank == 0)
            std::cout << "beta = " << beta << std::endl;

        betas.push_back(beta);

        for (int j = 0; j < N; j++)
            p[j] = r[j] + beta * p[j];

        eta = 1.0/gamw * ((double) N / Mt - scalar_fact * inner_prod(v, mu, 0) / Mt) + etas.back() * betas.back();

        if (rank == 0)
            std::cout << "eta = " << eta << std::endl;

        etas.push_back(eta);

        ksis.push_back(ksi);

        // stopping criteria
        double rel_err = sqrt( inner_prod_rr / l2_norm2(v, 0) );
        double norm_mu = sqrt( l2_norm2(mu, 0) ) * scalar_fact;
        double err_tol = 1e-4;
        err_tol = 1e-2;
        if (rank == 0)
            std::cout << "[CG] it = " << i << ": ||r_it|| / ||RHS|| = " << rel_err << ", ||x_it|| = " << norm_mu << std::endl;
        if (rel_err < err_tol) 
            break;

        if (rank == 0)
            std::cout << "inner_prod(alphas, etas, 0) = " << inner_prod(alphas, etas, 0) << ", <w, mu> / Mt = " << scalar_fact * inner_prod(noise, mu, 0) / Mt << std::endl;
    
        if (rank == 0){
            std::cout << "alphas elems = ";
            for (int j=0; j<alphas.size(); j++)
                std::cout << alphas[j] << ", ";
            std::cout << std::endl;
            
            std::cout << "etas elems = ";
            for (int j=0; j<etas.size(); j++)
                std::cout << etas[j] << ", ";
            std::cout << std::endl;
        }

        if (rank == 0)
            std::cout << "inner_prod(v, mu, 0) / Mt = " << scalar_fact * inner_prod(v, mu, 0) / Mt << std::endl;

        double yWy = inner_prod(y, lmmse_multAAT(y, tau, dataset), 0);

        if (rank == 0)
            std::cout << "a0 = " << gam2 * gamw * l2_norm2(y, 0) / yWy << std::endl;

        if (rank == 0)
            std::cout << "eta0 = " << (double) N / (double) Mt / gamw << std::endl;

        if (rank == 0)
            std::cout << "zero iter aeta = " << gam2 * l2_norm2(y, 0) * (double) N / (double) Mt / yWy << std::endl;
    }

    double alphas_etas = inner_prod(alphas, etas, 0);
    // alphas_etas = inner_prod(noise, mu, 0) / Mt;

    *correction = gam2 * ( alphas_etas - scalar_fact * inner_prod(v, mu, 0) / Mt ) * (-1);

    if (rank == 0)
        std::cout << "correction = " << *correction << std::endl;

    if (*correction < 0)
        *correction = - *correction;

    if (rank == 0)
        std::cout << "ksis.back() = " << ksis.back() << ", other part = " << scalar_fact * scalar_fact * l2_norm2(mu, 0) / gamw / Mt << std::endl;
    
    *new_vars = (ksis.back() - scalar_fact * scalar_fact * l2_norm2(mu, 0) / gamw / Mt) / (*correction) / (*correction) * gam2  - 1.0 / gam2; 

    if (rank == 0)
        std::cout << "new_vars = " << *new_vars << std::endl;

    for (int j = 0; j < N; j++)
        p_start[j] = p[j] - r[j];

    if (save == 1)
        mu_CG_last = mu; // we save the non-scaled version

    std::vector<double> out_mu = mu;
    for (int i=0; i<N; i++)
        out_mu[i] *= scalar_fact;
    return out_mu;
 }

 

void vamp::updateNoisePrecAAT(data* dataset){

    std::vector<double> temp = (*dataset).Ax(x2_hat.data());
    // std::vector<double> y = (*dataset).get_phen();

    for (int i = 0; i < N; i++)  // because length(y) = N
        temp[i] -= y[i];
    
    double temp_norm2 = l2_norm2(temp, 0); 
    double trace_corr = Mt * (1-alpha2) / gamw;
    if (rank == 0){
        std::cout << "l2_norm2(temp) / N = " << temp_norm2 / N << std::endl;
        std::cout << "trace_correction / N = " << trace_corr / N << std::endl;
    }
    gamw = (double) N / (temp_norm2 + trace_corr);
     
}


double vamp::g2d_onsagerAAT(double gam2, double tau, data* dataset) { // shared between linear and binary classification model
    
    std::random_device rd;
    std::bernoulli_distribution bern(0.5);

    bern_vec = std::vector<double> (M, 0.0);
    for (int i = 0; i < M; i++)
        bern_vec[i] = (2*bern(rd) - 1) / sqrt(Mt); // Bernoulli variables are sampled independently
    std::vector<double> res = (*dataset).Ax(bern_vec.data());
    double correct_empty = 0;
    double new_vars_empty = 0;    
    std::vector<double> p_start = std::vector<double> (4*(*dataset).get_mbytes(), 0.0);
    invQ_bern_vec = CG_solverAAT(res, std::vector<double> (4*(*dataset).get_mbytes(), 0.0), p_start, tau, &correct_empty, &new_vars_empty, 0, dataset); // precond_change
    res = (*dataset).ATx(invQ_bern_vec.data());
    double onsager = inner_prod(bern_vec, res, 1) * gamw; 
    return gam2 * (1 + onsager);    
}