/*
This contains various functions required for calculating the workspace of the CYCLOPS.
*/

#include "cyclops.h"
#include "Eigen/Eigen"
#include <math.h>

using namespace Eigen;
using namespace std;

Matrix3d cyclops::test_function(int num_vars, Vector3d position)
{
	Matrix3d Position_mat;
	Position_mat = position * position.transpose();
	return Position_mat;
}

bool cyclops::feasible_pose(Matrix<double, 5,1> P, Matrix<double,3,6> a, 
	                        Matrix<double,3,6> B, Matrix<double,6,1> W,
                            Vector3d f_ee, Vector3d r_ee,
                            VectorXd t_min, VectorXd t_max)
{

	Vector3d p;
    
	p << P(0,0), P(1,0), P(2,0); //The position of the tool
	double phi_y = P(3,0); //Orientation of the tool (y)
	double phi_z = P(4,0); //Orientation of the tool (z)

    // Define rotation Matrix
    Matrix3d R_y, R_z, T_r;
    R_y << cos(phi_y), 0, sin(phi_y),
           0, 1, 0,
           -sin(phi_y), 0, cos(phi_y);
    R_z << cos(phi_z), -sin(phi_z), 0,
           sin(phi_z), cos(phi_z), 0,
           0, 0, 1;
    T_r = R_y * R_z;

    Matrix<double, 5, 6> A;


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
    	Vector3d l_hat = l/l.norm();
        
        // Calculate Torque Component of Structural Matrix
        Vector3d tau = a_temp.cross(l);

        // Adding to the structural Matrix
        A.block<3,1>(0,i) = l_hat;
        A.block<2,1>(3,i) = tau.segment(1,2);    
    }

    // Compute overall wrench
    // We ignore torque in the x direction
    Vector3d tau_f_ee = f_ee.cross(r_ee);
    Matrix<double,6,1> f_temp;
    f_temp << f_ee(0), f_ee(1), f_ee(2), tau_f_ee(0), tau_f_ee(1), tau_f_ee(2);
    Matrix<double,6,1> f_temp2 = W + f_temp;
    Matrix<double,5,1> f;
    f << f_temp2(0), f_temp2(1), f_temp2(2), f_temp2(4), f_temp2(5);

    // Obtaining Tension Solution
    // Analytical method with L1-norm Solution

    Matrix<double,5,5> Partition_A = A.block<5,5>(0,0);
    Matrix<double,5,1> Partition_B = A.block<5,1>(0,5);

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

    t_low(5,0) = t_min(5);
    t_high(5,0) = t_max(5);

    double t_B_min, t_B_max;

    MatrixXf::Index tempRow, tempCol;

    t_B_min = t_low.maxCoeff(&tempRow, &tempCol);
    t_B_max = t_high.minCoeff(&tempRow, &tempCol);

    bool feasible = false;

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

cyclops::dw_result cyclops::dex_workspace(Matrix<double,3,6> a, Matrix<double,3,6> B,
                        Matrix<double,6,1> W, vector<Vector3d> f_ee, 
                        Vector3d r_ee, Vector2d phi_min, Vector2d phi_max,
                        VectorXd t_min, VectorXd t_max)
{
    int x_res = 10;
    int y_res = 10;
    int z_res = 10;

    // Determining Search Volume;
    Matrix<double,1,6> B_x = B.block<1,6>(0,0);
    MatrixXf::Index tempRow, tempCol;
    double x_middle = (B_x.maxCoeff(&tempRow, &tempCol) + B_x.minCoeff(&tempRow, &tempCol))/2;
    
    Vector3d p_arbitrary;
    p_arbitrary << x_middle, 0 , 0;

    double x_space_length1 = x_space_length(a, B, p_arbitrary);
    double x_space_length2 = x_space_length(a, B, -p_arbitrary);

    double radius = (B.block<2,1>(1,0)).norm();

    double x_step = (x_space_length1 + x_space_length2)/x_res;
    double y_step = 2 * radius /y_res;
    double z_step = 2 * radius /z_res;

    double x_temp = x_middle - x_space_length1;
    double y_temp = -radius;
    double z_temp = -radius;

    vector<Vector3d> vol_grid;
    for (int i=0; i<x_res+1; i++)
    {
        for (int j=0; j<y_res+1; j++)
        {
            for (int k=0; k<z_res+1; k++)
            {
                if ((z_temp * z_temp + y_temp * y_temp) <= (radius * radius))
                {
                    Vector3d temp_vec;
                    temp_vec << x_temp, y_temp, z_temp;
                    vol_grid.push_back(temp_vec);
                }
                z_temp+= z_step;
            }
            z_temp = -radius;
            y_temp+= y_step;
        }
        y_temp = -radius;
        x_temp+= x_step;
    }

    dw_result Temp;
    Temp.feasible = vol_grid;
    Temp.unfeasible = vol_grid;
    Temp.size = 0;
    return Temp;
}

double cyclops::x_space_length(Matrix<double,3,6> a, Matrix<double,3,6> B, Vector3d p)
{
    double max_length = 0;

    for (int i=0; i<a.cols(); i++)
    {
        Vector3d a_temp = a.block<3,1>(0,i);
        Vector3d a_b = a_temp + p;

        Vector3d B_temp = B.block<3,1>(0,i);

        double x_length = a_b(0) - B_temp(0);

        if (max_length < x_length)
            max_length = x_length;
    }
    return max_length;
}

