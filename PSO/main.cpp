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
    
    Eigen::Matrix3d Output_temp;
    Output_temp << 1,2,3,4,5,6,7,8,9;

    // WIP
    Eigen::Vector3d Output2 = Output_temp.block<3,1>(0,2);

    std::cout << "The output is: " << std::endl;
    std::cout << Output << std::endl;
    std::cout << "Output_temp is " << std::endl << Output_temp << std::endl;
    std::cout << "Output 2 is " << std::endl << Output2 << std::endl;

    Eigen::Matrix<double,3,6> a1, B;
    a1 << 30, 30, 30, -30, -30, -30,
         0.5, -1, 0.5, 0.5, -1, 0.5,
         0.866, 0, -0.866, 0.866, 0, -0.866;
    B << 0, 0, 0, -100, -100, -100,
         15, -30, 15, 15, -30, 15,
         25.98, 0, -25.98, 25.98, 0, -25.98;
    Eigen::Vector3d f_ee, r_ee;
    f_ee << 0, 0, 0;
    r_ee << 70, 0, 0;

    Eigen::Matrix<double,5,1> P;
    P << -46, 0, 0, 0, 0;
   // P << -200, 0, 0, 0, 0;

    Eigen::Matrix<double, 6, 1> t_min, t_max;
    t_min << 10,10,10,10,10,10;
    t_max << 60,60,60,60,60,60;

    bool result = cyclops::feasible_pose(P/1000, a1/1000, B/1000, f_ee, r_ee/1000, t_min, t_max);

    std::cout << "Feasibility of Given Pose: " << result << std::endl;
    return 0;
}
