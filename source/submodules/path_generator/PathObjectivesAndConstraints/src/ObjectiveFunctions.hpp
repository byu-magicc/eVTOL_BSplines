#ifndef OBJECTIVEFUNCTIONS_HPP
#define OBJECTIVEFUNCTIONS_HPP
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include "CBindHelperFunctions.hpp"

template <int D>
class ObjectiveFunctions
{
    public:
        ObjectiveFunctions();
        double minimize_jerk_magnitude(double cont_pts[], int num_cont_pts);
        double minimize_acceleration_magnitude(double cont_pts[], int num_cont_pts);
        double minimize_velocity_magnitude(double cont_pts[], int num_cont_pts);
    private:
        CBindHelperFunctions<D> cbind_help{};
};

extern "C"
{
    ObjectiveFunctions<2>* ObjectiveFunctions_2(){return new ObjectiveFunctions<2>();}
    double minimize_jerk_magnitude_2(ObjectiveFunctions<2>* obj, double cont_pts[], 
        int num_control_points){return obj->minimize_jerk_magnitude(
            cont_pts, num_control_points);}
    double minimize_acceleration_magnitude_2(ObjectiveFunctions<2>* obj, double cont_pts[], 
        int num_control_points){return obj->minimize_acceleration_magnitude(
            cont_pts, num_control_points);}
    double minimize_velocity_magnitude_2(ObjectiveFunctions<2>* obj, double cont_pts[], 
        int num_control_points){return obj->minimize_velocity_magnitude(
            cont_pts, num_control_points);}

    ObjectiveFunctions<3>* ObjectiveFunctions_3(){return new ObjectiveFunctions<3>();}
    double minimize_jerk_magnitude_3(ObjectiveFunctions<3>* obj, double cont_pts[], 
        int num_control_points){return obj->minimize_jerk_magnitude(
            cont_pts, num_control_points);}
    double minimize_acceleration_magnitude_3(ObjectiveFunctions<3>* obj, double cont_pts[], 
        int num_control_points){return obj->minimize_acceleration_magnitude(
            cont_pts, num_control_points);}
    double minimize_velocity_magnitude_3(ObjectiveFunctions<3>* obj, double cont_pts[], 
        int num_control_points){return obj->minimize_velocity_magnitude(
            cont_pts, num_control_points);}
}

#endif