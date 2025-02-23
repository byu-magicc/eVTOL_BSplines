#include "WaypointConstraints.hpp"
#include <iostream>

template<int D>
WaypointConstraints<D>::WaypointConstraints()
{

}

template<int D>
double* WaypointConstraints<D>::velocity_constraint(double cont_pts[], int num_control_points, 
    double desired_velocity[], double inverse_scale_factor, bool isStartPoint)
{
    int order = 3;
    double* constraints = new double[D*2];
    Eigen::Matrix<double,D,1> scaled_vel_vector;
    if (isStartPoint)
    {
        scaled_vel_vector = get_initial_velocity(cont_pts, num_control_points, inverse_scale_factor);
    }
    else // is end point
    {
        scaled_vel_vector = get_final_velocity(cont_pts, num_control_points, inverse_scale_factor);
    }
    for(unsigned int i = 0; i < D; i++)
    {
            constraints[i] = scaled_vel_vector.coeff(i) - desired_velocity[i];
    }
    return constraints;
}

template<int D>
double* WaypointConstraints<D>::acceleration_constraint(double cont_pts[], int num_control_points, 
    double desired_acceleration[], double inverse_scale_factor, bool isStartPoint)
{
    int order = 3;
    double* constraints = new double[D*2];
    Eigen::Matrix<double,D,1> scaled_accel_vector;
    if (isStartPoint)
    {
        scaled_accel_vector = get_initial_acceleration(cont_pts, num_control_points, inverse_scale_factor);
    }
    else // is end point
    {
        scaled_accel_vector = get_final_acceleration(cont_pts, num_control_points, inverse_scale_factor);
    }
    for(unsigned int i = 0; i < D; i++)
    {
            constraints[i] = scaled_accel_vector.coeff(i) - desired_acceleration[i];
    }
    return constraints;
}

template<int D>
Eigen::Matrix<double,D,1> WaypointConstraints<D>::get_initial_velocity(double cont_pts[], int &num_cont_pts, double &inverse_scale_factor)
{
    unsigned int start_index = 0;
    Eigen::Matrix<double,D,4> initial_control_points = cbind_help.array_section_to_eigen(cont_pts, num_cont_pts, start_index);
    Eigen::Matrix<double,D,1> p1 = initial_control_points.block(0,0,D,1);
    Eigen::Matrix<double,D,1> p3 = initial_control_points.block(0,2,D,1);
    Eigen::Matrix<double,D,1> velocity_vector_0 = inverse_scale_factor*(p3 - p1)/2;
    return velocity_vector_0;
}

template<int D>
Eigen::Matrix<double,D,1> WaypointConstraints<D>::get_final_velocity(double cont_pts[], int &num_cont_pts, double &inverse_scale_factor)
{
    int order = 3;
    unsigned int start_index = num_cont_pts - order - 1;
    Eigen::Matrix<double,D,4> final_control_points = cbind_help.array_section_to_eigen(cont_pts, num_cont_pts, start_index);
    Eigen::Matrix<double,D,1> p2 = final_control_points.block(0,1,D,1);
    Eigen::Matrix<double,D,1> p4 = final_control_points.block(0,3,D,1);
    Eigen::Matrix<double,D,1> velocity_vector_f = inverse_scale_factor*(p4 - p2)/2;
    return velocity_vector_f;
}

template<int D>
Eigen::Matrix<double,D,1> WaypointConstraints<D>::get_initial_acceleration(double cont_pts[], int &num_cont_pts, double &inverse_scale_factor)
{
    unsigned int start_index = 0;
    Eigen::Matrix<double,D,4> initial_control_points = cbind_help.array_section_to_eigen(cont_pts, num_cont_pts, start_index);
    Eigen::Matrix<double,D,1> p1 = initial_control_points.block(0,0,D,1);
    Eigen::Matrix<double,D,1> p2 = initial_control_points.block(0,1,D,1);
    Eigen::Matrix<double,D,1> p3 = initial_control_points.block(0,2,D,1); 
    Eigen::Matrix<double,D,1> acceleration_vector_0 = inverse_scale_factor*inverse_scale_factor*(p1 - 2*p2 + p3);
    return acceleration_vector_0;
}

template<int D>
Eigen::Matrix<double,D,1> WaypointConstraints<D>::get_final_acceleration(double cont_pts[], int &num_cont_pts, double &inverse_scale_factor)
{
    int order = 3;
    unsigned int start_index = num_cont_pts - order - 1;
    Eigen::Matrix<double,D,4> final_control_points = cbind_help.array_section_to_eigen(cont_pts, num_cont_pts, start_index);
    Eigen::Matrix<double,D,1> p2 = final_control_points.block(0,1,D,1);
    Eigen::Matrix<double,D,1> p3 = final_control_points.block(0,2,D,1);
    Eigen::Matrix<double,D,1> p4 = final_control_points.block(0,3,D,1);
    Eigen::Matrix<double,D,1> acceleration_vector_f = inverse_scale_factor*inverse_scale_factor*(p2 - 2*p3 + p4);
    return acceleration_vector_f;
}


//explicit instantiation
template class WaypointConstraints<2>;
template class WaypointConstraints<3>;
