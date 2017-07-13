/*
This contains various functions required for calculating the workspace of the CYCLOPS.
*/

#include "cyclops.h"
#include "Eigen/Eigen"

using namespace Eigen;

Matrix3d cyclops::test_function(int num_vars, Vector3d position)
{
	Matrix3d Position_mat;
	Position_mat = position * position.transpose();
	return Position_mat;
}
