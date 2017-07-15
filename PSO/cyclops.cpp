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
	                        Matrix<double,3,6> B, Matrix<double,Dynamic,3> f_ee,
	                        double r_ee, VectorXd t_min, VectorXd t_max)
{

	Vector3d p << P[0], P[1], P[2]; //The position of the tool
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

    // Create Structural Matrix, A
    for (int i=0; i<a.cols; i++)
    {
    	// Getting Attachment points in base frame WIP

    }

	return true;
}