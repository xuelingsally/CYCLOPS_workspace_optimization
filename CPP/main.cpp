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



    //---- Defining CYCLOPS parameters ----
    cyclops::fnInputs fnInputs;
    // Radius of Overtube
    double radius_tool = 1.0;
    fnInputs.radius_tool = radius_tool;
    // Radius of Scaffold
    double radius_scaffold = 30.0;
    fnInputs.radius_scaffold = radius_scaffold;
    // Lenght of Scaffold
    double length_scaffold = 100.0;
    fnInputs.length_scaffold = length_scaffold;
    // Length of Overtube
    double length_overtube = 60.0;
    // Distance of back of tool to CG of tool
    double dist_tool_b_cg = 30.0;



    // Tensions
    Eigen::Matrix<double,6,1> t_min;
    t_min << 5.0,5.0,5.0,5.0,5.0,5.0;
    fnInputs.t_min = t_min;
    Eigen::Matrix<double,6,1> t_max;
    t_max << 60.0,60.0,60.0,60.0,60.0,60.0;
    fnInputs.t_max = t_max;

    // Orientation Limits
    Eigen::Vector2d phi_min;
    phi_min << -10.0/180.0 * PI, -10.0/180.0 * PI;
    fnInputs.phi_min = phi_min;
    Eigen::Vector2d phi_max;
    phi_max << 10.0/180.0 * PI, 10.0/180.0 * PI;
    fnInputs.phi_max = phi_max;


    // Wrench Vector
    Eigen::Matrix<double,6,1> W;
    W << 0.0,0.0,0.0,0.0,0.0,0.0;
    fnInputs.W = W;

    // end effector force
    Eigen::Vector3d f_ee;
    f_ee << -1.0,0.0,0.0;
    fnInputs.f_ee_vec.push_back(f_ee);
    f_ee(0) = 1.0;
    fnInputs.f_ee_vec.push_back(f_ee);
    f_ee(0) = 0.0;
    f_ee(1) = -1.0;
    fnInputs.f_ee_vec.push_back(f_ee);
    f_ee(1) = 1.0;
    fnInputs.f_ee_vec.push_back(f_ee);
    f_ee(1) = 0.0;
    f_ee(2) = -1.0;
    fnInputs.f_ee_vec.push_back(f_ee);
    f_ee(2) = 1.0;
    fnInputs.f_ee_vec.push_back(f_ee);


    // Taskspace Definition
    Eigen::Vector3d tp_temp;
    tp_temp(0) = 20; tp_temp(1) = 5, tp_temp(2) = 5;
    fnInputs.taskspace.push_back(tp_temp/1000.0);
    tp_temp(0) = 20; tp_temp(1) = -5, tp_temp(2) = 5;
    fnInputs.taskspace.push_back(tp_temp/1000.0);
    tp_temp(0) = 20; tp_temp(1) = -5, tp_temp(2) = -5;
    fnInputs.taskspace.push_back(tp_temp/1000.0);
    tp_temp(0) = 20; tp_temp(1) = 5, tp_temp(2) = -5;
    fnInputs.taskspace.push_back(tp_temp/1000.0);

    tp_temp(0) = 25; tp_temp(1) = 5, tp_temp(2) = -6;
    fnInputs.taskspace.push_back(tp_temp/1000.0);
    tp_temp(0) = 25; tp_temp(1) = 5, tp_temp(2) = 6;
    fnInputs.taskspace.push_back(tp_temp/1000.0);
    tp_temp(0) = 25; tp_temp(1) = -5, tp_temp(2) = -6;
    fnInputs.taskspace.push_back(tp_temp/1000.0);
    tp_temp(0) = 25; tp_temp(1) = -5, tp_temp(2) = 6;
    fnInputs.taskspace.push_back(tp_temp/1000.0);

    tp_temp(0) = 30; tp_temp(1) = 5, tp_temp(2) = -5;
    fnInputs.taskspace.push_back(tp_temp/1000.0);
    tp_temp(0) = 30; tp_temp(1) = 5, tp_temp(2) = 5;
    fnInputs.taskspace.push_back(tp_temp/1000.0);
    tp_temp(0) = 30; tp_temp(1) = -5, tp_temp(2) = -5;
    fnInputs.taskspace.push_back(tp_temp/1000.0);
    tp_temp(0) = 30; tp_temp(1) = -5, tp_temp(2) = 5;
    fnInputs.taskspace.push_back(tp_temp/1000.0);
    tp_temp(0) = 35; tp_temp(1) = 5, tp_temp(2) = -10;
    fnInputs.taskspace.push_back(tp_temp/1000.0);


    // Design Vector and limits
    vector<double> lower_range(15);
    vector<double> upper_range(15);

    // Limits for e (angles)
    for (int i=0; i<6; i++)
    {
        lower_range[i]=0.0;
        upper_range[i]= 2.0 * PI;
    }
    // Limits for a
    // Front 3 tendons
    for (int i=6; i<9; i++)
    {
        lower_range[i] = -dist_tool_b_cg;
        upper_range[i] = 0.0;
    }
    // Back 3 Tendons
    for (int i=9; i<12; i++)
    {
        lower_range[i] = 0.0;
        upper_range[i] = (length_overtube - dist_tool_b_cg);
    }
    // Limits for B
    for (int i=12; i<14; i++)
    {
        lower_range[i] = -length_scaffold;
        upper_range[i] = 0.0;
    }
    // Limits for tooltip
    lower_range[14] = 0.0;
    upper_range[14] = 70.0;


    // Parameters for PSO
    int number_runs=1;
    int pso_pop=500;
    int pso_iters=80;
    float phi_p=0.50;
    float phi_g=0.50;;
    double omega_initial=0.90;
    double omega_final=0.50;
    bool rand_update=false; 

    // perform PSO experiment
    //simpsolib::EvalFN pso_eval_fn((char *)"DEJONG2", DEJONG2_FN_NUM_VARS,lower_range, upper_range, dejong2_test_fn_real);
    simpsolib::EvalFN pso_eval_fn((char *)"CYCLOPS_OPT", 15, lower_range, upper_range, dejong2_test_fn_real);
    pso_eval_fn.Input = fnInputs;
    run_pso(pso_eval_fn, number_runs, pso_pop, pso_iters, phi_p, phi_g, omega_initial, omega_final, rand_update);
    
 
/*

    Eigen::Matrix<double,3,6> a, B;
    
//    a << 30.0, 30.0, 30.0, -30.0, -30.0, -30.0,
//         0.5, -1.0, 0.5, 0.5, -1.0, 0.5,
//         0.866, 0.0, -0.866, 0.866, 0.0, -0.866;
//    B << 0.0, 0.0, 0.0, -100.0, -100.0, -100.0,
//         15.0, -30.0, 15.0, 15.0, -30.0, 15.0,
//         25.98, 0.0, -25.98, 25.98, 0.0, -25.98; 

    a << -7.2716,  -11.3936,  -21.7060,   11.5591,  16.7089,   19.7650,
   -0.6269,    0.9821,   -0.2583,   -0.3502,    1.0000,   -0.4169,
   -0.7791,   -0.1882,    0.9661,   -0.9367,   -0.0000,    0.9090;

    B << -85.1078, -85.1078, -85.1078, -4.9541, -4.9541, -4.9541,
    -18.8056, 29.4638, -7.7491, -10.5072, 30.0000, -12.5071,
    -23.3741, -5.6468, 28.9819, -28.0998, -0.0000, 27.2685;

    Eigen::Vector3d f_ee;
    Eigen::Vector3d r_ee;
    f_ee << 0.0, 0.0, 0.0;
    r_ee << 70.0, 0.0, 0.0;
    Eigen::Matrix<double,6,1> W;
    W << 0.0,0.0,0.0,0.0,0.0,0.0;
    std::vector<Eigen::Vector3d> f_ee_vec;
    f_ee_vec.push_back(f_ee);
    
    std::cout << "f_ee is: " << f_ee << std::endl << std::endl;


    Eigen::Matrix<double,5,1> P;
    P << -50.0, 0.0, 0.0, 0.0, 0.0;
    //P << -200, 0, 0, 0, 0;
    //P << -0.0660, -0.0060, -0.0000, 0, 0;

    Eigen::Matrix<double, 6, 1> t_min, t_max;
    t_min << 5.0,5.0,5.0,5.0,5.0,5.0;
    t_max << 60.0,60.0,60.0,60.0,60.0,60.0;
    
    Eigen::Vector2d phi_min, phi_max;
    phi_min << -10.0/180.0 * PI, -10.0/180.0 * PI;
    phi_max << 10.0/180.0 * PI, 10.0/180.0 * PI;

    bool result = cyclops::feasible_pose(P/1000.0, a/1000.0, B/1000.0, W, f_ee, r_ee/1000.0, t_min, t_max);

    std::cout << "Feasibility of Given Pose: " << result << std::endl;

    cyclops::dw_result test_dw = cyclops::dex_workspace(a/1000.0, B/1000.0, W, f_ee_vec, r_ee/1000.0, phi_min, phi_max, t_min, t_max);
    

    std::cout << "No. of Feasible pts is: " << (test_dw.feasible).size() << std::endl;
    std::cout << "No. of UnFeasible pts is: " << (test_dw.unfeasible).size() << std::endl;
    std::cout << "Workspace Size: " << (test_dw.size) << std::endl;

    std::vector<Eigen::Vector3d> feasible = test_dw.feasible;
    std::vector<Eigen::Vector3d> unfeasible = test_dw.unfeasible;

    std::vector<Eigen::Vector3d>::iterator feasible_iter;
    std::vector<Eigen::Vector3d>::iterator unfeasible_iter;

    std::cout << "Positions:" << std::endl;
    for(feasible_iter=feasible.begin(); feasible_iter!=feasible.end(); ++feasible_iter)
    {
        std::cout << (*feasible_iter)(0) << ", " << (*feasible_iter)(1) << ", " << (*feasible_iter)(2) << ";" << std::endl;
    }

    Eigen::Matrix<double, 15, 1> eaB;
    eaB << 4.0349, 6.0938, 1.8321, 4.3546, 6.2832, 2.0008,
    -7.2716, -11.3936, -21.7060, 11.5591, 16.7089, 19.7650,
    -85.1078, -4.9541,
    70.0000;

    double val = cyclops::objective_function(eaB, W, f_ee_vec, phi_min, phi_max,
                          t_min, t_max, fnInputs.taskspace, radius_tool, 
                          radius_scaffold);

    std::cout << "Value of Objective Function is: " << val << std::endl;
*/
    return 0;
}
