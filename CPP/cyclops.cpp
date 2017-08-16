/*
This contains various functions required for calculating the workspace of the CYCLOPS.
*/

#include "cyclops.h"

using namespace Eigen;
using namespace std;

Matrix3d cyclops::test_function(int num_vars, Vector3d position)
{
	Matrix3d Position_mat;
	Position_mat = position * position.transpose();
	return Position_mat;
}

bool cyclops::feasible_pose(Matrix<double, 5,1> P, Matrix<double,3,Dynamic> a, 
	                        Matrix<double,3,Dynamic> B, Matrix<double,6,1> W,
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

    Matrix<double, 5, Dynamic> A;

    A.resize(5,a.cols());

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
        Vector3d tau = a_temp.cross(l_hat);

        // Adding to the structural Matrix
        A.block<3,1>(0,i) = l_hat;
        A.block<2,1>(3,i) = tau.block<2,1>(1,0);    
    }


    // Compute overall wrench
    // We ignore torque in the x direction
    Vector3d tau_f_ee = f_ee.cross(r_ee);
    Matrix<double,6,1> f_temp;
    f_temp << f_ee(0), f_ee(1), f_ee(2), tau_f_ee(0), tau_f_ee(1), tau_f_ee(2);
    Matrix<double,6,1> f_temp2 = W + f_temp;
    Matrix<double,5,1> f;
    f << f_temp2(0), f_temp2(1), f_temp2(2), f_temp2(4), f_temp2(5);

    //cout << "f = " << endl << f << endl << endl;

    // Obtaining Tension Solution
    // Analytical method with L1-norm Solution

    int redundant = a.cols() - 5;


    bool feasible = false;


    Matrix<double,5,5> Partition_A = A.block<5,5>(0,0);
    //Matrix<double,5,1> Partition_B = A.block<5,1>(0,5);
    Matrix<double,5,Dynamic> Partition_B;
    Partition_B.resize(5, redundant);
    for (int i=0; i<redundant; i++)
    {
        Partition_B.block<5,1>(0,i) = A.block<5,1>(0,i+5);
    }

    //cout << "A = " << endl << A << endl << endl;

    Matrix<double,5,1> M = -Partition_A.inverse() * f;
    //Matrix<double,5,1> N = -Partition_A.inverse() * Partition_B;

    Matrix<double,5,Dynamic> N;
    N.resize(5,redundant);
    N = -Partition_A.inverse() * Partition_B;

    //cout << "M = " << endl << M << endl << endl;
    //cout << "N = " << endl << N << endl << endl;


    // Using the analytical method to test if a pose is feasible.
    // Only can be used if there is only one redundant cable.
    if (redundant == 1)
    {

        Matrix<double, 6, 1> t_low;
        Matrix<double, 6, 1> t_high;

        for (int i=0; i<5; i++)
        {
            if (N(i,0) > 0)
            {
                t_low(i,0) = (t_min(i) - M(i,0))/N(i,0);
                t_high(i,0) = (t_max(i) - M(i,0))/N(i,0);
            }
            else
            {
                t_low(i,0) = (t_max(i) - M(i,0))/N(i,0);
                t_high(i,0) = (t_min(i) - M(i,0))/N(i,0);
            }
        }


        t_low(5,0) = t_min(5);
        t_high(5,0) = t_max(5);

        double t_B_min, t_B_max;

        MatrixXf::Index tempRow, tempCol;

        t_B_min = t_low.maxCoeff(&tempRow, &tempCol);
        t_B_max = t_high.minCoeff(&tempRow, &tempCol);


        if (t_B_min <= t_B_max)
        {
            feasible = true;
        }
        else
        {
            feasible = false;
        }

    } 
    

    // Using the L1 Norm solution
    // Used if there is more than 1 redundant Cable
    // As per Hay and Snyman (2005)

    else
    {
        Matrix<double,1,5> ones_5;
        Matrix<double,1,5> a_;
        for (int i=0; i<5; i++)
        {
            ones_5(0,i) = 1;
        }

        a_ = ones_5 * (-Partition_A.inverse());


        Matrix<double, 1, Dynamic> ones_redunt;
        ones_redunt.resize(1, redundant);
        for (int i=0; i<redundant; i++)
        {
            ones_redunt(0,i) = 1;
        }

        Matrix<double, 1, Dynamic> c;
        c.resize(1,redundant);
        c = a_ * Partition_B + ones_redunt;
        Matrix<double, Dynamic, 1> c_;
        c_.resize(redundant, 1);
        c_ = c.transpose();


        Matrix<double, Dynamic, Dynamic> A_;
        A_.resize(10+redundant, redundant);
        for (int i=0; i<redundant; i++)
        {
            A_.block<5,1>(0,i) = -N.block<5,1>(0,i);
            A_.block<5,1>(5,i) = N.block<5,1>(0,i);
        }

        Matrix<double, Dynamic, 1> b_;
        b_.resize(10+redundant, 1);

        b_.block<5,1>(0,0) = - t_min.block<5,1>(0,0) + M;
        b_.block<5,1>(5,0) = t_max.block<5,1>(0,0) - M;


        for (int i=0; i<redundant; i++)
        {
            for (int j=0; j<redundant; j++)
            {
                if (i==j)
                {
                    A_(i+10,j) = 1;
                }
                else
                {
                    A_(i+10,j) = 0; 
                }
            }

            b_(10+i) = t_max(0,0) - t_min(0,0);
        } 

        Matrix<double, Dynamic, 1> y;
        y.resize(redundant, 1);


        feasible = igl::linprog(c_, A_, b_, 10+redundant, y);
    }

	return feasible;
}

cyclops::dw_result cyclops::dex_workspace(Matrix<double,3,Dynamic> a, Matrix<double,3,Dynamic> B,
                        Matrix<double,6,1> W, vector<Vector3d> f_ee_vec, 
                        Vector3d r_ee, Vector2d phi_min, Vector2d phi_max,
                        VectorXd t_min, VectorXd t_max,
                        double length_scaffold)
{
    int x_res = 10;
    int y_res = 10;
    int z_res = 10;

    // Determining Search Volume;
    Matrix<double,1,Dynamic> B_x; 
    B_x.resize(1,B.cols());
    for (int i=0; i<B.cols(); i++)
    {
        B_x(0,i) = B(0,i);
    }
    //B_x = B.block<1,B.cols()>(0,0);

    MatrixXf::Index tempRow, tempCol;
    double x_middle = (B_x.maxCoeff(&tempRow, &tempCol) + B_x.minCoeff(&tempRow, &tempCol))/2;

    Vector3d p_arbitrary;
    p_arbitrary << x_middle, 0 , 0;

    double x_space_length1 = x_space_length(a, B, p_arbitrary);
    double x_space_length2 = x_space_length(-a, -B, -p_arbitrary);

    double radius = (B.block<2,1>(1,0)).norm();

    //double x_search = (x_space_length1 + x_space_length2);
    double x_end = x_middle + x_space_length2;
    double x_step = length_scaffold/x_res;
    double y_step = 2 * radius /y_res;
    double z_step = 2 * radius /z_res;

    double x_temp = x_middle - x_space_length1;
    double y_temp = -radius;
    double z_temp = -radius;

    vector<Vector3d> vol_grid;
    //for (int i=0; i<x_res+1; i++)
    while (x_temp <= x_end)
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

    // Checking if the points in the search volume are feasible across the given forces on the end effector
    vector<Vector3d> feasible, unfeasible;
    vector<Vector3d>::iterator vol_grid_iter;
    vector<Vector3d>::iterator f_ee_iter;


    for (vol_grid_iter = vol_grid.begin(); vol_grid_iter!=vol_grid.end(); ++vol_grid_iter)
    {

        unsigned int feasible_counter = 0;
        for (f_ee_iter = f_ee_vec.begin(); f_ee_iter!=f_ee_vec.end(); ++f_ee_iter)
        {
            bool feasible_temp = false;
            Matrix<double,5,1> P;
            
            P << (*vol_grid_iter)(0), (*vol_grid_iter)(1), (*vol_grid_iter)(2), 0, 0;
            feasible_temp = feasible_pose(P, a, B, W, *f_ee_iter, r_ee, t_min, t_max);
            if (feasible_temp)
                feasible_counter++;

            P(3,0) = phi_min(0);
            feasible_temp = feasible_pose(P, a, B, W, *f_ee_iter, r_ee, t_min, t_max);
            if (feasible_temp)
                feasible_counter++;

            P(3,0) = phi_max(0);
            feasible_temp = feasible_pose(P, a, B, W, *f_ee_iter, r_ee, t_min, t_max);
            if (feasible_temp)
                feasible_counter++;

            P(3,0) = 0;
            P(4,0) = phi_min(1);
            feasible_temp = feasible_pose(P, a, B, W, *f_ee_iter, r_ee, t_min, t_max);
            if (feasible_temp)
                feasible_counter++;

            P(4,0) = phi_max(1);
            feasible_temp = feasible_pose(P, a, B, W, *f_ee_iter, r_ee, t_min, t_max);
            if (feasible_temp)
                feasible_counter++;
        }

        if (feasible_counter == 5 * f_ee_vec.size())
            feasible.push_back(*vol_grid_iter);
        else
            unfeasible.push_back(*vol_grid_iter);
    }

    //double yz_search = PI * radius * radius;
    //double search_vol = x_search * yz_search;
    //double wp_size = double(feasible.size())/double(vol_grid.size()) * search_vol;
    double wp_size = double(feasible.size())/(x_res*y_res*z_res);

    dw_result Result;
    Result.feasible = feasible;
    Result.unfeasible = unfeasible;
    Result.size = wp_size;
    return Result;
}

double cyclops::x_space_length(Matrix<double,3,Dynamic> a, Matrix<double,3,Dynamic> B, Vector3d p)
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

double cyclops::objective_function(Matrix<double,Dynamic,1> eaB, Matrix<double,6,1> W,
	                      vector<Vector3d> f_ee_vec,
	                      Vector2d phi_min, Vector2d phi_max,
	                      VectorXd t_min, VectorXd t_max,
	                      vector<VectorXd> taskspace,
	                      double radius_tool, double radius_scaffold,
                          double length_scaffold)
{
	double val = 0.0;
    //cout << eaB << endl;
	// Finding the attachment and feeding points based on design vector
	Matrix<double,3,Dynamic> a, B;

    int num_tendons = (eaB.rows() - 15)/3 + 6;
    a.resize(3,num_tendons);
    B.resize(3,num_tendons);

    // Feeding and attachment points for the 6 main tendons.
    for(int i=0; i<6; i++)
    {
        double cos_mul = cos(eaB(i,0));
        double sin_mul = sin(eaB(i,0));

        a(0,i) = eaB(i+6,0);
        a(1,i) = radius_tool * sin_mul; //y-component
        a(2,i) = radius_tool * cos_mul; //z-component

        if(i<3)
            B(0,i) = eaB(12,0);
        else
            B(0,i) = eaB(13,0);

        B(1,i) = radius_scaffold * sin_mul; //y-component
        B(2,i) = radius_scaffold * cos_mul; //z-component
    }

    // Feeding and attachment points for the 'extra' tendons
    for (int i=0; i<num_tendons - 6; i++)
    {
        double cos_mul = cos(eaB(i*3+15, 0));
        double sin_mul = sin(eaB(i*3+15, 0));

        a(0,i+6) = eaB(i*3+15+1,0);
        a(1,i+6) = radius_tool * sin_mul; 
        a(2,i+6) = radius_tool * cos_mul;

        B(0,i+6) = eaB(i*3+15+2,0);
        B(1,i+6) = radius_scaffold * sin_mul;
        B(2,i+6) = radius_scaffold * cos_mul;
    }

/*
    cout << "a = [";
    for (int i=0; i<a.rows(); i++){
    	for (int j=0; j<a.cols(); j++) {
    		cout << a(i,j);
    		if (j < a.cols()-1)
    			cout << ",";
    	}
    	cout << ";" << endl;
    }
    cout << "];" << endl;
	
	cout << "B = [";
	for (int i=0; i<B.rows(); i++){
    	for (int j=0; j<B.cols(); j++) {
    		cout << B(i,j);
    		if (j < B.cols()-1)
    			cout << ",";
    	}
    	cout << ";" << endl;
    }
    cout << "];" << endl;

*/


	// Checking for crossing of cables
	for (int i=0;i<a.cols();i++)
	{
		for (int j=0;j<a.cols();j++)
		{
			if(a(0,i) < a(0,j) && B(0,i) > B(0,j))
			{
				val = -1.5;
			    //cout << "Returned0 " << val << endl;
				return val;
			}
		}
        if (a(0,i) >= eaB(14))
        {
            val = -1.5;
            return val;
        }
	}

	// Based on the tooltip, find the poses that the CG of the tool has to reach
	double dist_tooltip = eaB(14);
	Vector3d r_ee;
	r_ee << dist_tooltip, 0, 0;


	// Checking if the points in the taskspace are feasible across the given forces on the end effector
    vector<VectorXd>::iterator taskspace_iter;
    vector<Vector3d>::iterator f_ee_iter;


    for (taskspace_iter = taskspace.begin(); taskspace_iter!=taskspace.end(); ++taskspace_iter)
    {
        unsigned int feasible_counter = 0;

        Vector3d taskspace_temp = (*taskspace_iter) - r_ee/1000;

        for (f_ee_iter = f_ee_vec.begin(); f_ee_iter!=f_ee_vec.end(); ++f_ee_iter)
        {

            bool feasible_temp = false;
            Matrix<double,5,1> P;
            
            P << (taskspace_temp)(0), (taskspace_temp)(1), (taskspace_temp)(2), 0, 0;
            feasible_temp = feasible_pose(P, a/1000.0, B/1000.0, W, *f_ee_iter, r_ee/1000.0, t_min, t_max);
            if (feasible_temp)
            {
                feasible_counter++;
                feasible_temp = false;
            }

            P(3,0) = phi_min(0);
            feasible_temp = feasible_pose(P, a/1000.0, B/1000.0, W, *f_ee_iter, r_ee/1000.0, t_min, t_max);
            if (feasible_temp)
            {
                feasible_counter++;
                feasible_temp = false;
            }

            P(3,0) = phi_max(0);
            feasible_temp = feasible_pose(P, a/1000.0, B/1000.0, W, *f_ee_iter, r_ee/1000.0, t_min, t_max);
            if (feasible_temp)
            {
                feasible_counter++;
                feasible_temp = false;
            }

            P(3,0) = 0;
            P(4,0) = phi_min(1);
            feasible_temp = feasible_pose(P, a/1000.0, B/1000.0, W, *f_ee_iter, r_ee/1000.0, t_min, t_max);
            if (feasible_temp)
            {
                feasible_counter++;
                feasible_temp = false;
            }

            P(4,0) = phi_max(1);
            feasible_temp = feasible_pose(P, a/1000.0, B/1000.0, W, *f_ee_iter, r_ee/1000.0, t_min, t_max);
            if (feasible_temp)
            {
                feasible_counter++;
                feasible_temp = false;
            }

            if (feasible_counter != 5)
                val = val - 1.0;
            feasible_counter = 0;
        }
    }

    if (val < 0.0)
    {
    	//cout << val << endl;
    	double val2 = double(val)/(double(taskspace.size()) * double(f_ee_vec.size()));
    	//cout << "Returned1 " << val2 << endl;
    	return val2;
    }


    // Calculate and return the zero wrench workspace if the taskspace condition is fulfilled.

    Vector3d zero_f_ee;
    zero_f_ee << 0,0,0;

    vector<Vector3d> zero_f_ee_vec;
    zero_f_ee_vec.push_back(zero_f_ee);

    dw_result dex_wp = dex_workspace(a/1000.0, B/1000.0, W, zero_f_ee_vec, r_ee/1000.0, phi_min, phi_max, t_min, t_max, length_scaffold/1000.0);

    val = dex_wp.size;
    //cout << "Returned2 " << val << endl;
    return val;
}

double cyclops::objective_function2(Matrix<double,Dynamic,1> eaB, Matrix<double,6,1> W,
                          vector<Vector3d> f_ee_vec,
                          Vector2d phi_min, Vector2d phi_max,
                          VectorXd t_min, VectorXd t_max,
                          vector<VectorXd> taskspace,
                          double radius_tool, double radius_scaffold,
                          double length_scaffold)
{
    double val = 0.0;
    //cout << eaB << endl;
    // Finding the attachment and feeding points based on design vector
    Matrix<double,3,Dynamic> a, B;

    int num_tendons = (eaB.rows() - 15)/3 + 6;
    a.resize(3,num_tendons);
    B.resize(3,num_tendons);


    // Feeding and attachment points for the 6 main tendons.
    for(int i=0; i<6; i++)
    {
        double cos_mul = cos(eaB(i,0));
        double sin_mul = sin(eaB(i,0));

        a(0,i) = eaB(i+6,0);
        a(1,i) = radius_tool * sin_mul; //y-component
        a(2,i) = radius_tool * cos_mul; //z-component

        if(i<3)
            B(0,i) = eaB(12,0);
        else
            B(0,i) = eaB(13,0);

        B(1,i) = radius_scaffold * sin_mul; //y-component
        B(2,i) = radius_scaffold * cos_mul; //z-component
    }

    // Feeding and attachment points for the 'extra' tendons
    for (int i=0; i<num_tendons - 6; i++)
    {
        double cos_mul = cos(eaB(i*3+15, 0));
        double sin_mul = sin(eaB(i*3+15, 0));

        a(0,i+6) = eaB(i*3+15+1,0);
        a(1,i+6) = radius_tool * sin_mul; 
        a(2,i+6) = radius_tool * cos_mul;

        B(0,i+6) = eaB(i*3+15+2,0);
        B(1,i+6) = radius_scaffold * sin_mul;
        B(2,i+6) = radius_scaffold * cos_mul;
    }

/*
    cout << "a = [";
    for (int i=0; i<a.rows(); i++){
        for (int j=0; j<a.cols(); j++) {
            cout << a(i,j);
            if (j < a.cols()-1)
                cout << ",";
        }
        cout << ";" << endl;
    }
    cout << "];" << endl;
    
    cout << "B = [";
    for (int i=0; i<B.rows(); i++){
        for (int j=0; j<B.cols(); j++) {
            cout << B(i,j);
            if (j < B.cols()-1)
                cout << ",";
        }
        cout << ";" << endl;
    }
    cout << "];" << endl;

*/


    // Checking for crossing of cables
    for (int i=0;i<a.cols();i++)
    {
        for (int j=0;j<a.cols();j++)
        {
            if(a(0,i) < a(0,j) && B(0,i) > B(0,j))
            {
                val = -1.5;
                //cout << "Returned0 " << val << endl;
                return val;
            }
        }
        // Check that the tooltip is greater than the max value of a_x
        if (a(0,i) >= eaB(14))
        {
            val = -1.5;
            return val;
        }
    }

    if (abs(eaB(12) - eaB(13)) < 0.05 * length_scaffold)
    {
        val = -1.5;
        return val;
    }

    // Based on the tooltip, find the poses that the CG of the tool has to reach
    double dist_tooltip = eaB(14);
    Vector3d r_ee;
    r_ee << dist_tooltip, 0, 0;


    // Checking if the points in the taskspace are feasible across the given forces on the end effector
    vector<VectorXd>::iterator taskspace_iter;
    vector<Vector3d>::iterator f_ee_iter;


    for (taskspace_iter = taskspace.begin(); taskspace_iter!=taskspace.end(); ++taskspace_iter)
    {

        Matrix<double,5,1> r_ee_temp;
        // The taskspace orientation
        double alpha_y = (*taskspace_iter)(3,0);
        double alpha_z = (*taskspace_iter)(4,0);

        //cout << "alpha_y: " << alpha_y << "  alpha_z: " << alpha_z << endl;

        // Find the required position of the tool
        // Define rotation Matrix
        Matrix3d R_y, R_z, T_r;
        R_y << cos(alpha_y), 0, sin(alpha_y),
               0, 1, 0,
               -sin(alpha_y), 0, cos(alpha_y);
        R_z << cos(alpha_z), -sin(alpha_z), 0,
               sin(alpha_z), cos(alpha_z), 0,
               0, 0, 1;
        T_r = R_y * R_z;

        Vector3d r_ee_rotated = -T_r * r_ee;
        //cout << r_ee_rotated.transpose() << endl;

        //r_ee_temp << r_ee(0), r_ee(1), r_ee(2), 0.0, 0.0;
        r_ee_temp << r_ee_rotated(0,0), r_ee_rotated(1,0), r_ee_rotated(2,0), 0.0, 0.0;

        Matrix<double,5,1> taskspace_temp = (*taskspace_iter) + r_ee_temp/1000;
        //cout << taskspace_temp.transpose() * 1000 << ";" <<std::endl;

        for (f_ee_iter = f_ee_vec.begin(); f_ee_iter!=f_ee_vec.end(); ++f_ee_iter)
        {

            bool feasible_temp = false;
            Matrix<double,5,1> P;
            P << (taskspace_temp)(0,0), (taskspace_temp)(1,0), (taskspace_temp)(2,0), (taskspace_temp)(3,0), (taskspace_temp)(4,0);
            feasible_temp = feasible_pose(P, a/1000.0, B/1000.0, W, *f_ee_iter, r_ee/1000.0, t_min, t_max);
            if (!feasible_temp)
            {
                val = val - 1.0;
            }
        }
    }
    //cout << endl << dist_tooltip << endl << endl;

    if (val < 0.0)
    {
        //cout << val << endl;
        double val2 = double(val)/(double(taskspace.size()) * double(f_ee_vec.size()));
        //cout << "Returned1 " << val2 << endl;
        return val2;
    }


    // Calculate and return the zero wrench workspace if the taskspace condition is fulfilled.

    Vector3d zero_f_ee;
    zero_f_ee << 0,0,0;

    vector<Vector3d> zero_f_ee_vec;
    zero_f_ee_vec.push_back(zero_f_ee);

    dw_result dex_wp = dex_workspace(a/1000.0, B/1000.0, W, zero_f_ee_vec, r_ee/1000.0, phi_min, phi_max, t_min, t_max, length_scaffold/1000.0);

    val = dex_wp.size;
    //cout << "Returned2 " << val << endl;
    return val;
}