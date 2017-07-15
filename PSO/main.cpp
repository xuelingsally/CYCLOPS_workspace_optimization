/*
 * FILE: main.cpp, v.1.7.1, 4/28/2014
 * Author: Tomas V. Arredondo
 *
 * Example of a C++ application using SimGALib and SimPSOLib.
 *
 * This is NOT a tutorial on GA as there are many tutorials and books on this
 * subject see http://en.wikipedia.org/wiki/Genetic_algorithm for more information.
 *
 * To change the function and range being optimized go to ga.h and change
 * FUNCTION_LOWER_RANGE, FUNCTION_UPPER_RANGE and FUNCTION definitions.
 *
 * DISCLAIMER: No liability is assumed by the author for any use made
 * of this program.
 * DISTRIBUTION: Any use may be made of this program, as long as the
 * clear acknowledgment is made to the author in code and runtime executables
 */
#include <iostream>
#include <ctime>
#include <cmath>
#include <queue>
#include "simgalib.h"
#include "simpsolib.h"
#include "simtstlib.h"

#include "Eigen/Eigen"
#include "cyclops.h"

using namespace simgalib;
using namespace simpsolib;
using namespace simtstlib;

int main()
{
    /*
    int number_runs=50;
    //pso parms
    int pso_pop=20;
    int pso_iters=100;
    float phi_p=1.49445;
    float phi_g=1.49445;
    float omega=.729;
    bool rand_update=false;

   
    // perform PSO experiment
    vector<double> lower_range(DEJONG2_FN_NUM_VARS);
    vector<double> upper_range(DEJONG2_FN_NUM_VARS);

    for (int i=0; i< DEJONG2_FN_NUM_VARS; i++)
    {
        lower_range[i]=DEJONG2_FN_LOWER_RANGE;
        upper_range[i]=DEJONG2_FN_UPPER_RANGE;
    }

    simpsolib::EvalFN pso_eval_fn((char *)"DEJONG2", DEJONG2_FN_NUM_VARS,lower_range, upper_range, dejong2_test_fn_real);
    run_pso(pso_eval_fn, number_runs, pso_pop, pso_iters, phi_p, phi_g, omega, rand_update);
    */
    
    cout << "Enter Position" << endl;
    double a, b, c;
    std::cin >> a;
    std::cin >> b;
    std::cin >> c;

    Eigen::Vector3d position;
    position << a, b, c; 
    
    std::cout << "First Input is " << position[0] << " Second Input is " << position[1] << "Third Input is " << position[2] << std::endl;

    Eigen::Matrix3d Output = cyclops::test_function(int(1), position);
    
    // WIP
    //Eigen::Map<Eigen::VectorXf,0,Eigen::OuterStride<> > Output2(Output.data, Output.rows(), Eigen::InnerStride<>(3));

    std::cout << "The output is: " << std::endl;
    std::cout << Output << std::endl;
    std::cout << "Output 2 is " << std::endl << Output2 << std::endl;

    return 0;
}
