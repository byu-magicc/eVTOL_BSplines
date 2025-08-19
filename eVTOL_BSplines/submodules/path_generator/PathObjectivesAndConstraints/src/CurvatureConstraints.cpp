#include "CurvatureConstraints.hpp"

template <int D>
CurvatureConstraints<D>::CurvatureConstraints()
{

}

template <int D>
double CurvatureConstraints<D>::get_spline_curvature_constraint(double cont_pts[], 
    int num_cont_pts, double max_curvature)
{
    double constraint;
    constraint = curv_bound.get_spline_curvature_bound(cont_pts, num_cont_pts) - max_curvature;
    return constraint;
}

template <int D>
double* CurvatureConstraints<D>::get_interval_curvature_constraints(double cont_pts[], int num_cont_pts, 
    double max_curvature)
{
    int order = 3;
    int num_intervals = num_cont_pts - order;
    double* constraints = new double[num_intervals];
    double curvature_bound;
    for(unsigned int i = 0; i < num_intervals; i++)
    {
        Eigen::Matrix<double,D,4> interval_control_points = cbind_help.array_section_to_eigen(cont_pts, num_cont_pts, i);
        curvature_bound = curv_bound.evaluate_interval_curvature_bound(interval_control_points);
        constraints[i] = curvature_bound - max_curvature;
    }
    return constraints;
}

template <int D>
Eigen::VectorXd CrossTermBounds<D>::get_interval_curvature_bounds(double cont_pts[], int &num_control_points)
{
    int order = 3;
    int num_intervals = num_control_points - order;
    Eigen::VectorXd curvature_bounds(num_intervals);
    for (unsigned int i = 0; i < num_control_points-3; i++)
    {
        Eigen::Matrix<double,D,4> interval_control_points = cbind_help.array_section_to_eigen(cont_pts, num_control_points, i);
        curvature_bounds(i) = evaluate_interval_curvature_bound(interval_control_points);
    }
    return curvature_bounds;
}

// explicit instantiation
template class CurvatureConstraints<2>;
template class CurvatureConstraints<3>;