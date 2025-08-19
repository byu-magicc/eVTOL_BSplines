#ifndef WAYPOINTCONSTRAINTSOLD_HPP
#define WAYPOINTCONSTRAINTSOLD_HPP
#include "DerivativeEvaluator.hpp"
#include "CBindHelperFunctions.hpp"
#include "CrossTermEvaluator.hpp"

template<int D>
class WaypointConstraintsOld
{
    public:
        WaypointConstraintsOld();
        double* velocity_at_waypoints_constraints(double cont_pts[], int num_control_points, 
            double scale_factor, double desired_velocities[], bool switches[]);

        double* acceleration_at_waypoints_constraints(double cont_pts[], int num_cont_pts, 
            double scale_factor, double desired_accelerations[], bool switches[]);

    private:
        CBindHelperFunctions<D> cbind_help{};
        CrossTermEvaluator<D> c_eval{};
        DerivativeEvaluator<D> d_dt_eval{};
};

extern "C"
{
    WaypointConstraintsOld<2>* WaypointConstraintsOld_2(){return new WaypointConstraintsOld<2>();}
    double* velocity_at_waypoints_constraints_2(WaypointConstraintsOld<2>* obj, double cont_pts[], int num_control_points,
            double scale_factor, double desired_velocities[], bool switches[]){return obj->velocity_at_waypoints_constraints(
            cont_pts, num_control_points, scale_factor, desired_velocities, switches);}
    double* acceleration_at_waypoints_constraints_2(WaypointConstraintsOld<2>* obj, double cont_pts[], int num_cont_pts, 
            double scale_factor, double desired_accelerations[], bool switches[]){return obj->acceleration_at_waypoints_constraints(
            cont_pts, num_cont_pts, scale_factor, desired_accelerations, switches);}

    WaypointConstraintsOld<3>* WaypointConstraintsOld_3(){return new WaypointConstraintsOld<3>();}
    double* velocity_at_waypoints_constraints_3(WaypointConstraintsOld<3>* obj, double cont_pts[], int num_control_points,
            double scale_factor, double desired_velocities[], bool switches[]){return obj->velocity_at_waypoints_constraints(
            cont_pts, num_control_points, scale_factor, desired_velocities, switches);}
    double* acceleration_at_waypoints_constraints_3(WaypointConstraintsOld<3>* obj, double cont_pts[], int num_cont_pts, 
            double scale_factor, double desired_accelerations[], bool switches[]){return obj->acceleration_at_waypoints_constraints(
            cont_pts, num_cont_pts, scale_factor, desired_accelerations, switches);}
}

#endif