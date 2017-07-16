
#include <iostream>
#include <vector>
#include <cmath>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "Eigen/Eigen"

namespace cyclops
{
using namespace Eigen;

Matrix3d test_function(int num_vars, Vector3d position);

bool feasible_pose(Matrix<double, 5,1> P, Matrix<double,3,6> a, 
	               Matrix<double,3,6> B, Vector3d f_ee,
	               Vector3d r_ee, VectorXd t_min, VectorXd t_max);

}; //namespace cyclops