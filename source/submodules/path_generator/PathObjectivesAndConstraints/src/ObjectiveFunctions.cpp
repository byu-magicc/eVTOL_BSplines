#include "ObjectiveFunctions.hpp"
#include <iostream>

template <int D>
ObjectiveFunctions<D>::ObjectiveFunctions()
{

}

template <int D>
double ObjectiveFunctions<D>::minimize_jerk_magnitude(double cont_pts[], int num_cont_pts)
{
    int order = 3;
    double sum_of_integrals{0};
    int step = num_cont_pts-3;
    Eigen::MatrixXd control_points = cbind_help.array_to_eigen(cont_pts, num_cont_pts);
    Eigen::MatrixXd jerk_cont_pts = control_points.block(0,3,D,step) - 
        3*control_points.block(0,2,D,step) + 3*control_points.block(0,1,D,step) - 
        control_points.block(0,0,D,step);
    sum_of_integrals = (jerk_cont_pts.cwiseAbs2()).sum();
    return sum_of_integrals;
}


template <int D>
double ObjectiveFunctions<D>::minimize_acceleration_magnitude(double cont_pts[], int num_cont_pts)
{
    int order = 3;
    double sum_of_integrals{0};
    int step = num_cont_pts-2;
    Eigen::MatrixXd control_points = cbind_help.array_to_eigen(cont_pts, num_cont_pts);
    Eigen::MatrixXd acceleration_cont_pts = control_points.block(0,2,D,step)
         - 2*control_points.block(0,1,D,step) + control_points.block(0,0,D,step);
    sum_of_integrals = (acceleration_cont_pts.cwiseAbs2()).sum();
    return sum_of_integrals;
}

template <int D>
double ObjectiveFunctions<D>::minimize_velocity_magnitude(double cont_pts[], int num_cont_pts)
{
    int order = 3;
    double sum_of_integrals{0};
    int step = num_cont_pts-1;
    Eigen::MatrixXd control_points = cbind_help.array_to_eigen(cont_pts, num_cont_pts);
    Eigen::MatrixXd velocity_cont_pts = control_points.block(0,1,D,step)
        - control_points.block(0,0,D,step);
    sum_of_integrals = (velocity_cont_pts.cwiseAbs2()).sum();
    return sum_of_integrals;
}

//explicit instantiation

template class ObjectiveFunctions<2>;
template class ObjectiveFunctions<3>;