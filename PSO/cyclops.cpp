/*
This contains various functions required for calculating the workspace of the CYCLOPS.
*/

#include "cyclops.h"
#include "Eigen/Eigen"
#include <math.h>

using namespace Eigen;

Matrix3d cyclops::test_function(int num_vars, Vector3d position)
{
	Matrix3d Position_mat;
	Position_mat = position * position.transpose();
	return Position_mat;
}

bool cyclops::feasible_pose(Matrix<double, 5,1> P, Matrix<double,3,6> a, 
	                        Matrix<double,3,6> B, Vector3d f_ee,
	                        Vector3d r_ee, VectorXd t_min, VectorXd t_max)
{

	Vector3d p;
	p << P[0], P[1], P[2]; //The position of the tool
	double phi_y = P[4]; //Orientation of the tool (y)
	double phi_z = P[5]; //Orientation of the tool (z)

    // Define rotation Matrix
    Matrix3d R_y, R_z, T_r;
    R_y << cos(phi_y), 0, sin(phi_y),
           0, 1, 0,
           -sin(phi_y), 0, cos(phi_y);
    R_z << cos(phi_z), -sin(phi_z), 0,
           sin(phi_z), cos(phi_z), 0,
           0, 0, 1;
    T_r = R_y * R_z;

    Matrix<double, 6, Dynamic> A;

    // Create Structural Matrix, A
    for (int i=0; i<a.cols(); i++)
    {
    	// Getting Attachment points in base frame
    	Vector3d a_temp = a.block<3,1>(0,i);
    	Vector3d a_b = (T_r * a_temp) + p;

    	// Find the vector, l for each tendon
    	Vector3d B_temp = B.block<3,1>(0,i);
    	Vector3d l = B_temp - a_b;

    	//Get the normalised unit vector
    	Vector3d l_hat = l/l.squaredNorm();
        
        // Calculate Torque Component of Structural Matrix
        Vector3d tau = a_temp.cross(l);

        // Adding to the structural Matrix
        A.block<3,1>(0,i) = l_hat;
        A.block<2,1>(3,i) = tau.block<2,1>(1,0);    
    }

    // Compute overall wrench
    // We ignore torque in the x direction
    Vector3d tau_f_ee = f_ee.cross(r_ee);
    Matrix<double,5,1> f;
    f << f_ee(0), f_ee(1), f_ee(2), tau_f_ee(1), tau_f_ee(2);

    // Obtaining Tension Solution
    // Analytical method with L1-norm Solution

    Matrix<double,5,5> Partition_A = A.block<5,5>(0,0);
    Matrix<double,5,1> Partition_B = A.block<5,1>(0,6);

    Matrix<double,5,1> M = - Partition_A.inverse() * f;
    Matrix<double,5,1> N = - Partition_A.inverse() * Partition_B;

    Matrix<double, 6, 1> t_low;
    Matrix<double, 6, 1> t_high;

    for (int i=0; i<5; i++)
    {
        if (N(i,0) > 0)
        {
            t_low(i,0) = (t_min(i) - M(i,0))/N(i,0);
            t_high(i,0) = (t_max(i) - M(i))/N(i);
        }
        else
        {
            t_low(i,0) = (t_min(i) - M(i,0))/N(i,0);
            t_high(i,0) = (t_min(i) - M(i,0))/N(i,0);
        }
    }

    t_low(6,0) = t_min(6);
    t_high(6,0) = t_max(6);

    double t_B_min, t_B_max;

    MatrixXf::Index tempRow, tempCol;

    t_B_min = t_low.maxCoeff(&tempRow, &tempCol);
    t_B_max = t_high.minCoeff(&tempRow, &tempCol);

    bool feasible;

    if (t_B_min <= t_B_max)
    {
        feasible = true;
    }
    else
    {
        feasible = false;
    }


	return feasible;
}