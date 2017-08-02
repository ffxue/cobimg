/**
 * COBIMG, Constrained Optimization-based BIM Generator
 * Copyright (c) 2016-2017 The University of Hong Kong
 * Author: Fan Xue <xuef@hku.hk; fanxue@outlook.com>
 *
 * This file is part of cobimg.
 *
 * COBIMG is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * COBIMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with COBIMG.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "lib/cmaespp/cmaes.h"

#include "lib/libcmaes/cmaes.h"
#include "lib/libcmaes/esoptimizer.h"
#include "lib/libcmaes/cmastrategy.h"
#include "lib/libcmaes/ipopcmastrategy.h"
#include "lib/libcmaes/bipopcmastrategy.h"

#include "lib/ga/ga.h"
#include "lib/ga/GAGenome.h"

#include <ruby.h>
 
#include "cobimg.h"
#include "photomeasure.h"
 
#include <limits>
#include <chrono>

namespace cobimg
{
    void Cobimg::load_ref_measurement(char* filename, int metric)
    {
        reference = new PhotoMeasurement(metric);
        reference->readFile(filename, false);
    };
    
    double Cobimg::evaluate_pop_representation(char* filename)
    {
        return reference->evaluate(filename);
    };
    void Cobimg::clear()
    {
        reference->clear();
    };
    VALUE Cobimg::ruby_get_best_solution()
    {
        VALUE ruby_sol = rb_ary_new2(dimension);
        for (std::size_t j = 0; j < dimension; j++)
            ruby_sol = rb_ary_push(  ruby_sol, rb_float_new(bestSol[j]) );
        return ruby_sol;
    };
    double Cobimg::get_best_fitness()
    {
        return bestVal;
    };
    double Cobimg::evolve_cmaespp(VALUE obj, int dim, int max_iterations, char* pop_filename, char* measure_func, char* progress_func)
    {
        // allocate mem
        dimension = dim;
        bestVal = std::numeric_limits<double>::max();
        bestSol = new double[dim];
        algorithm_time = 0.0;
        client_time = 0.0;
        evaluation_time = 0.0;
        overall_time = 0.0;
        ruby_obj = obj;
        ruby_measure_ID = rb_intern(measure_func);
        ruby_progress_ID = rb_intern(progress_func);
        
        // set up cmaespp
        CMAES<double> evo;
        double *fitvals, *const*pop;
        
        Parameters<double> parameters;
        double *xstart;
        double *stddev;
        double *stepev;
        xstart = new double[dim];
        stddev = new double[dim];
        stepev = new double[dim];
        for(int i=0; i<dim; i++) xstart[i] = 0.0;
        for(int i=0; i<dim; i++) stddev[i] = 0.6;
        for(int i=0; i<dim; i++) stepev[i] = 0.001;
        parameters.init(dim, xstart, stddev);
        parameters.stopMaxFunEvals = max_iterations;
        parameters.rgDiffMinChange = stepev;
        fitvals = evo.init(parameters); // allocs fitvals
        
        
        auto overall_start = std::chrono::high_resolution_clock::now();
        // each population has Lambda solutions
        size_t pop_size = evo.get(CMAES<double>::Lambda);
        while(!evo.testForTermination())
        {
            auto sample_start = std::chrono::high_resolution_clock::now();
            pop = evo.samplePopulation();
            algorithm_time += get_time_span ( sample_start );
            for (std::size_t i = 0; i < pop_size; i++)
            {
                // get measurement from client plugin
                auto client_start = std::chrono::high_resolution_clock::now();
                VALUE ruby_sol = rb_ary_new2(dim);
                for (std::size_t j = 0; j < dim; j++)
                    ruby_sol = rb_ary_push(  ruby_sol, rb_float_new(pop[i][j]) );
                rb_funcall(ruby_obj, ruby_measure_ID, 1, ruby_sol);
                client_time += get_time_span ( client_start );
                // evaluate the measurement
                auto eval_start = std::chrono::high_resolution_clock::now();
                double v = evaluate_pop_representation(pop_filename);
                evaluation_time += get_time_span ( eval_start );
                fitvals[i] = v;
                // tell client plugin current progress
                rb_funcall(ruby_obj, ruby_progress_ID, 1, rb_float_new(v));
                if (v < bestVal)
                {
                    // record when new best solution found
                    bestVal = v;
                    std::memcpy(bestSol, pop[i], sizeof(double)*dim);
                }
            }
            // cmaes update
            auto update_start = std::chrono::high_resolution_clock::now();
            evo.updateDistribution(fitvals);
            algorithm_time += get_time_span ( update_start );
        }
        overall_time = get_time_span ( overall_start );
        return bestVal;
    };
    
    double Cobimg::get_time_span(std::chrono::high_resolution_clock::time_point t0)
    {
        using namespace std::chrono;
        auto time_span = duration_cast<duration<double>>(high_resolution_clock::now() - t0);
        return time_span.count();
    }
    double Cobimg::get_algorithm_time()
    {
        return algorithm_time;
    };
    double Cobimg::get_client_time()
    {
        return client_time;
    };
    double Cobimg::get_evaluation_time()
    {
        return evaluation_time;
    };
    double Cobimg::get_overall_time()
    {
        return overall_time;
    };
    
    double Cobimg::evolve_libcmaes(VALUE obj, int dim, int max_iterations, char* pop_filename, char* measure_func, char* progress_func, int algorithm)
    {
        // allocate mem
        dimension = dim;
        bestVal = std::numeric_limits<double>::max();
        bestSol = new double[dim];
        algorithm_time = 0.0;
        client_time = 0.0;
        evaluation_time = 0.0;
        overall_time = 0.0;
        ruby_obj = obj;
        ruby_temp_filename = pop_filename;
        ruby_measure_ID = rb_intern(measure_func);
        ruby_progress_ID = rb_intern(progress_func);
        
        st_max_iter = max_iterations;
        st_iter = 0;
        // Lambda function of evaluation of a solution at client software
        libcmaes::FitFunc fitfunc = [&](const double *x, const int N)
        {
            // skip if exceeding the limit
            if (++st_iter > st_max_iter)
                return 0.0;
            
            // get measurement from client plugin
            auto client_start = std::chrono::high_resolution_clock::now();
            VALUE ruby_sol = rb_ary_new2(dimension);
            for (std::size_t j = 0; j < dimension; j++)
                ruby_sol = rb_ary_push(  ruby_sol, rb_float_new(x[j]) );
            rb_funcall(ruby_obj, ruby_measure_ID, 1, ruby_sol);
            client_time += Cobimg::get_time_span ( client_start );
            
            // evaluate the measurement
            auto eval_start = std::chrono::high_resolution_clock::now();
            double v = reference->evaluate(ruby_temp_filename);
            evaluation_time += Cobimg::get_time_span ( eval_start );
            
            // tell client plugin current progress
            rb_funcall(ruby_obj, ruby_progress_ID, 1, rb_float_new(v));
            if (v < bestVal)
            {
                // record when new best solution found
                bestVal = v;
                std::memcpy(bestSol, x, sizeof(double)*dimension);
            }
            return v;
        };
        
        // set up libcmaes
		std::vector<double> x0(dim, 0.5);
		double lbounds[dim],ubounds[dim];
		for (int i=0;i<dim;i++)
		{
			lbounds[i] = 0.0;
			ubounds[i] = 1.0;
		}
		libcmaes::GenoPheno<libcmaes::pwqBoundStrategy> gp(lbounds,ubounds,dim); // genotype / phenotype transform associated to bounds.
        
        libcmaes::CMAParameters<libcmaes::GenoPheno<libcmaes::pwqBoundStrategy>> cmaparams(x0, 0, -1, 0, gp); // -1 for automatically decided lambda
        
        cmaparams.set_seed(RANDOM_SEED);
        cmaparams.set_max_fevals(max_iterations);
        cmaparams.set_noisy();
        cmaparams.set_sep();
        cmaparams.set_algo(algorithm); 
        // Available algorithm : CMAES_DEFAULT, IPOP_CMAES, BIPOP_CMAES, aCMAES, aIPOP_CMAES, aBIPOP_CMAES, sepCMAES, sepIPOP_CMAES, sepBIPOP_CMAES, sepaCMAES, sepaIPOP_CMAES, sepaBIPOP_CMAES, VD_CMAES, VD_IPOP_CMAES, VD_BIPOP_CMAES 
        
        libcmaes::ESOptimizer<libcmaes::BIPOPCMAStrategy<libcmaes::ACovarianceUpdate,libcmaes::GenoPheno<libcmaes::pwqBoundStrategy>>,libcmaes::CMAParameters<libcmaes::GenoPheno<libcmaes::pwqBoundStrategy>>,libcmaes::CMASolutions> optim(fitfunc, cmaparams);
        optim.set_progress_func(libcmaes::CMAStrategy<libcmaes::CovarianceUpdate,libcmaes::GenoPheno<libcmaes::pwqBoundStrategy>>::_defaultPFunc);
        
        auto overall_start = std::chrono::high_resolution_clock::now();
        while(!optim.stop())
        {
            // generate population
            auto sample_start = std::chrono::high_resolution_clock::now();
            dMat candidates = optim.ask();
            algorithm_time += get_time_span ( sample_start );
            // evaluation 
            auto client_start = std::chrono::high_resolution_clock::now();
            optim.eval(candidates);
            // update distribution
            auto update_start = std::chrono::high_resolution_clock::now();
            optim.tell();
            optim.inc_iter(); // important step: signals next iteration.
            algorithm_time += get_time_span ( update_start );
        }
        overall_time = get_time_span ( overall_start );
        return optim.get_solutions().get_best_seen_candidate().get_fvalue();
    };
    
    
    float galib_fitfunc (GAGenome & c)
    {
        if (++st_iter > st_max_iter)
            return 1.0f;
        GABin2DecGenome & genome = (GABin2DecGenome &)c;
        // get measurement from client plugin
        auto client_start = std::chrono::high_resolution_clock::now();
        VALUE ruby_sol = rb_ary_new2(st_dim);
        for (std::size_t j = 0; j < st_dim; j++)
            ruby_sol = rb_ary_push(  ruby_sol, rb_float_new(genome.phenotype(j)) );
        rb_funcall(st_ruby_obj, st_ruby_measure_ID, 1, ruby_sol);
        
        auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - client_start);
        client_time += time_span.count();
        
        // evaluate the measurement
        auto eval_start = std::chrono::high_resolution_clock::now();
        double v = st_reference->evaluate(st_ruby_temp_filename);
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - eval_start);
        evaluation_time += time_span.count();
        
        // tell client plugin current progress
        rb_funcall(st_ruby_obj, st_ruby_progress_ID, 1, rb_float_new(v));
        return (float)(v);
    }
    
    double Cobimg::evolve_galib(VALUE obj, int dim, int max_iterations, char* pop_filename, char* measure_func, char* progress_func, float mutation, float crossover)
    {
        // allocate mem
        dimension = dim;
        bestVal = std::numeric_limits<double>::max();
        bestSol = new double[dim];
        algorithm_time = 0.0;
        client_time = 0.0;
        evaluation_time = 0.0;
        overall_time = 0.0;
        ruby_obj = obj;
        ruby_temp_filename = pop_filename;
        ruby_measure_ID = rb_intern(measure_func);
        ruby_progress_ID = rb_intern(progress_func);
        
        st_ruby_obj = ruby_obj;
        st_ruby_temp_filename = ruby_temp_filename;
        st_ruby_measure_ID = ruby_measure_ID;
        st_ruby_progress_ID = ruby_progress_ID;
        st_reference = reference;
        st_dim = dim;
        st_max_iter = max_iterations;
        st_iter = 0;
        
        // set up GALib
        
        // seed
        GARandomSeed((unsigned int)RANDOM_SEED);
        
        // Create a phenotype for two variables. We use 16 bits to represent a floating point number whose value
        // can range from 0 to 1, inclusive.  The bounds on x1 and x2 can be applied
        // here and/or in the objective function.

        GABin2DecPhenotype map;
        for (std::size_t j = 0; j < dimension; j++)
            map.add(16, 0, 1);

        // Create the genome using the phenotype map.
        // default: max(fitfunc)
        GABin2DecGenome genome(map, galib_fitfunc);

        auto overall_start = std::chrono::high_resolution_clock::now();
        // simple GA
        int popsize = (int)(sqrt(max_iterations)/1.732)+1;
        int ngen = (int)(max_iterations/popsize)+1;
        GASimpleGA ga(genome);
        GASigmaTruncationScaling scaling;
        ga.minimize();
        ga.populationSize(popsize);
        ga.nGenerations(ngen);
        ga.pMutation(mutation);
        ga.pCrossover(crossover);
        ga.scaling(scaling);
        // seed = 0
        ga.evolve(RANDOM_SEED);
        
        overall_time = get_time_span ( overall_start );
        genome = ga.statistics().bestIndividual();
        for (std::size_t j = 0; j < dimension; j++)
            bestSol[j] = (double)(genome.phenotype(j));
        bestVal = (double)(- ga.statistics().maxEver());
        
        // return best value
        return bestVal;
    };
}