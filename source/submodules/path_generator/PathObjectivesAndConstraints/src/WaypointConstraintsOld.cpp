#include "WaypointConstraintsOld.hpp"
#include <iostream>

template<int D>
WaypointConstraintsOld<D>::WaypointConstraintsOld()
{

}

template<int D>
double* WaypointConstraintsOld<D>::velocity_at_waypoints_constraints(double cont_pts[], int num_cont_pts, 
            double scale_factor, double desired_velocities[], bool switches[])
{
    int order = 3;
    double* constraints = new double[D*2];
    double t0 = 0;
    double tf = scale_factor;
    unsigned int start_index = 0;
    unsigned int end_index = num_cont_pts - order - 1;
    if (switches[0] == true)
    {
        Eigen::Matrix<double,D,4> initial_control_points = cbind_help.array_section_to_eigen(cont_pts, num_cont_pts, start_index);
        // Eigen::Matrix<double,D,1> velocity_vector_0 = d_dt_eval.calculate_velocity_vector(t0, initial_control_points, scale_factor);
        Eigen::Matrix<double,D,1> p1 = initial_control_points.block(0,0,D,1);
        Eigen::Matrix<double,D,1> p3 = initial_control_points.block(0,2,D,1);
        Eigen::Matrix<double,D,1> velocity_vector_0 = (p3 - p1)/(2*scale_factor);
        for(unsigned int i = 0; i < D; i++)
        {
            constraints[2*i] = velocity_vector_0.coeff(i) - desired_velocities[2*i];
        }
    }
    else
    {
        for(unsigned int i = 0; i < D; i++)
        {
            constraints[2*i] = 0;
        }
    }
    if (switches[1] == true)
    {
        Eigen::Matrix<double,D,4> final_control_points = cbind_help.array_section_to_eigen(cont_pts, num_cont_pts, end_index);
        // Eigen::Matrix<double,D,1> velocity_vector_f = d_dt_eval.calculate_velocity_vector(tf, final_control_points, scale_factor);
        Eigen::Matrix<double,D,1> p2 = final_control_points.block(0,1,D,1);
        Eigen::Matrix<double,D,1> p4 = final_control_points.block(0,3,D,1);
        Eigen::Matrix<double,D,1> velocity_vector_f = (p4 - p2)/(2*scale_factor);
        for(unsigned int i = 0; i < D; i++)
        {
            constraints[2*i+1] = velocity_vector_f.coeff(i) - desired_velocities[2*i+1];
        }
    }
    else
    {
        for(unsigned int i = 0; i < D; i++)
        {
            constraints[2*i+1] = 0;
        }

    }
    return constraints;
}

template<int D>
double* WaypointConstraintsOld<D>::acceleration_at_waypoints_constraints(double cont_pts[], int num_cont_pts, 
            double scale_factor, double desired_accelerations[], bool switches[])
{
    int order = 3;
    double* constraints = new double[D*2];
    double t0 = 0;
    double tf = scale_factor;
    unsigned int start_index = 0;
    unsigned int end_index = num_cont_pts - order - 1;
    if (switches[0] == true)
    {
        Eigen::Matrix<double,D,4> initial_control_points = cbind_help.array_section_to_eigen(cont_pts, num_cont_pts, start_index);
        Eigen::Matrix<double,D,1> p1 = initial_control_points.block(0,0,D,1);
        Eigen::Matrix<double,D,1> p2 = initial_control_points.block(0,1,D,1);
        Eigen::Matrix<double,D,1> p3 = initial_control_points.block(0,2,D,1);
        Eigen::Matrix<double,D,1> acceleration_vector_0 = (p1 - 2*p2 + p3)/(scale_factor * scale_factor);
        for(unsigned int i = 0; i < D; i++)
        {
            constraints[2*i] = acceleration_vector_0.coeff(i) - desired_accelerations[2*i];
        }
    }
    else
    {
        for(unsigned int i = 0; i < D; i++)
        {
            constraints[2*i] = 0;
        }
    }
    if (switches[1] == true)
    {
        Eigen::Matrix<double,D,4> final_control_points = cbind_help.array_section_to_eigen(cont_pts, num_cont_pts, end_index);
        Eigen::Matrix<double,D,1> p2 = final_control_points.block(0,1,D,1);
        Eigen::Matrix<double,D,1> p3 = final_control_points.block(0,2,D,1);
        Eigen::Matrix<double,D,1> p4 = final_control_points.block(0,3,D,1);
        Eigen::Matrix<double,D,1> acceleration_vector_f = (p2 - 2*p3 + p4)/(scale_factor * scale_factor);
        for(unsigned int i = 0; i < D; i++)
        {
            constraints[2*i+1] = acceleration_vector_f.coeff(i) - desired_accelerations[2*i+1];
        }
    }
    else
    {
        for(unsigned int i = 0; i < D; i++)
        {
            constraints[2*i+1] = 0;
        }

    }
    return constraints;
}

//explicit instantiation
template class WaypointConstraintsOld<2>;
template class WaypointConstraintsOld<3>;