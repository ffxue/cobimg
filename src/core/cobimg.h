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
 
#ifndef COBIMG_COBIMG_H
#define COBIMG_COBIMG_H

#include <chrono>
#include <ruby.h>

#include "measurement.h"
#include "../../lib/ga/GAGenome.h"

namespace cobimg
{
    int st_dim, st_max_iter, st_iter;
    char *st_ruby_temp_filename;
    double algorithm_time, client_time, evaluation_time, overall_time;
    
    VALUE st_ruby_obj;
    ID st_ruby_measure_ID , st_ruby_progress_ID;
    Measurement *st_reference;
    
    float galib_fitfunc (GAGenome & c);
        
    class Cobimg
    {
    private :
        size_t dimension;
        Measurement *reference;
        double bestVal, *bestSol;
        
        // for Ruby (SketchUp Pro)
        VALUE ruby_obj;
        char *ruby_temp_filename;
        ID ruby_measure_ID, ruby_progress_ID;
        
        const int RANDOM_SEED = 0;
        
    public :
        void load_ref_measurement(char* filename, int metric = 0);
        double evaluate_pop_representation(char* filename);
        void clear();
        VALUE ruby_get_best_solution();
        double get_best_fitness();
        double evolve_cmaespp(VALUE obj, int dim, int max_iterations, char* pop_filename, char* measure_func, char* progress_func);
        double get_time_span(std::chrono::high_resolution_clock::time_point t0);
        double get_algorithm_time();
        double get_client_time();
        double get_evaluation_time();
        double get_overall_time();
        double evolve_libcmaes(VALUE obj, int dim, int max_iterations, char* pop_filename, char* measure_func, char* progress_func, int algorithm = 11);
        double evolve_galib(VALUE obj, int dim, int max_iterations, char* pop_filename, char* measure_func, char* progress_func, float mutation = 0.01, float crossover = 0.6);
    };
}

#endif