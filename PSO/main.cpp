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

#define PI 3.141592654

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
    
    Eigen::Matrix<double,3,6> a, B;
    a << 30, 30, 30, -30, -30, -30,
         0.5, -1, 0.5, 0.5, -1, 0.5,
         0.866, 0, -0.866, 0.866, 0, -0.866;
    B << 0, 0, 0, -100, -100, -100,
         15, -30, 15, 15, -30, 15,
         25.98, 0, -25.98, 25.98, 0, -25.98;

    Eigen::Vector3d f_ee, r_ee;
    f_ee << 0, 0, 0;
    r_ee << 70, 0, 0;
    Eigen::Matrix<double,6,1> W;
    W << 0,0,0,0,0,0;
    std::vector<Eigen::Vector3d> f_ee_vec;
    f_ee_vec.push_back(f_ee);


    Eigen::Matrix<double,5,1> P;
    P << -50, 0, 0, 0, 0;
    //P << -200, 0, 0, 0, 0;
    //P << -0.0660, -0.0060, -0.0000, 0, 0;

    Eigen::Matrix<double, 6, 1> t_min, t_max;
    t_min << 5,5,5,5,5,5;
    t_max << 60,60,60,60,60,60;
    
    Eigen::Vector2d phi_min, phi_max;
    phi_min << -10/180 * PI, -10/180 * PI;
    phi_max << 10/180 * PI, 10/180 * PI;

    bool result = cyclops::feasible_pose(P/1000, a/1000, B/1000, W, f_ee, r_ee/1000, t_min, t_max);

    std::cout << "Feasibility of Given Pose: " << result << std::endl;

    cyclops::dw_result test_dw = cyclops::dex_workspace(a/1000, B/1000, W, f_ee_vec, r_ee/1000, phi_min, phi_max, t_min, t_max);
    

    std::cout << "No. of Feasible pts is: " << (test_dw.feasible).size() << std::endl;
    std::cout << "No. of UnFeasible pts is: " << (test_dw.unfeasible).size() << std::endl;
    std::cout << "Workspace Size: " << (test_dw.size) << std::endl;

    std::vector<Eigen::Vector3d> feasible = test_dw.feasible;

    std::vector<Eigen::Vector3d>::iterator feasible_iter;

    std::cout << "Positions:" << std::endl;
    for(feasible_iter=feasible.begin(); feasible_iter!=feasible.end(); ++feasible_iter)
    {
        std::cout << (*feasible_iter)(0) << ", " << (*feasible_iter)(1) << ", " << (*feasible_iter)(2) << ";" << std::endl;
    }

    return 0;
}
