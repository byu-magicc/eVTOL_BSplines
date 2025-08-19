#ifndef WAYPOINTCONSTRAINTS_HPP
#define WAYPOINTCONSTRAINTS_HPP
#include "CBindHelperFunctions.hpp"
#include "RotationHelperFunctions.hpp"
template<int D>
class WaypointConstraints
{
    public:
        WaypointConstraints();
        double* velocity_constraint(double cont_pts[], int num_control_points, 
            double desired_velocity[], double inverse_scale_factor, bool isStartPoint);
        double* acceleration_constraint(double cont_pts[], int num_control_points, 
            double desired_acceleration[], double inverse_scale_factor, bool isStartPoint);

    private:
        CBindHelperFunctions<D> cbind_help{};
        RotationHelperFunctions<D> rot_help{};
        Eigen::Matrix<double,D,1> get_initial_velocity(double cont_pts[], int &num_control_points, double &scale_factor);
        Eigen::Matrix<double,D,1> get_final_velocity(double cont_pts[], int &num_control_points, double &scale_factor);
        Eigen::Matrix<double,D,1> get_initial_acceleration(double cont_pts[], int &num_control_points, double &scale_factor);
        Eigen::Matrix<double,D,1> get_final_acceleration(double cont_pts[], int &num_control_points, double &scale_factor);
};

extern "C"
{
    WaypointConstraints<2>* WaypointConstraints_2(){return new WaypointConstraints<2>();}
    double* velocity_constraint_2(WaypointConstraints<2>* obj, double cont_pts[], int num_control_points, 
            double desired_velocity[], double inverse_scale_factor, bool isStartPoint){return obj->velocity_constraint(
            cont_pts, num_control_points, desired_velocity, inverse_scale_factor, isStartPoint);}
    double* acceleration_constraint_2(WaypointConstraints<2>* obj, double cont_pts[], int num_control_points, 
            double desired_acceleration[], double inverse_scale_factor, bool isStartPoint){return obj->acceleration_constraint(
            cont_pts, num_control_points, desired_acceleration, inverse_scale_factor, isStartPoint);}

    WaypointConstraints<3>* WaypointConstraints_3(){return new WaypointConstraints<3>();}
    double* velocity_constraint_3(WaypointConstraints<3>* obj, double cont_pts[], int num_control_points, 
            double desired_velocity[], double inverse_scale_factor, bool isStartPoint){return obj->velocity_constraint(
            cont_pts, num_control_points, desired_velocity, inverse_scale_factor, isStartPoint);}
    double* acceleration_constraint_3(WaypointConstraints<3>* obj, double cont_pts[], int num_control_points, 
            double desired_acceleration[], double inverse_scale_factor, bool isStartPoint){return obj->acceleration_constraint(
            cont_pts, num_control_points, desired_acceleration, inverse_scale_factor, isStartPoint);}

}

#endif