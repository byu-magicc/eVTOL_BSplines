#include "InclinationConstraint.hpp"

template <int D>
InclinationConstraint<D>::InclinationConstraint()
{

}

template <int D>
double InclinationConstraint<D>::get_spline_incline_constraint(double cont_pts[], int num_cont_pts, 
            double scale_factor, double max_incline)
{
    int order = 3;
    int num_intervals = num_cont_pts - order;
    double incline_extrema = 0;
    for(unsigned int i = 0; i < num_intervals; i++)
    {
        Eigen::Matrix<double,D,4> interval_control_points = cbind_help.array_section_to_eigen(cont_pts, num_cont_pts, i);
        double inclination_bound = evaluate_bound(interval_control_points, scale_factor);
        if (incline_extrema < inclination_bound)
        {
            incline_extrema = inclination_bound;
        }
    }
    double constraint = incline_extrema - max_incline;
    return constraint;

}

template <int D>
double* InclinationConstraint<D>::get_interval_incline_constraints(double cont_pts[], int num_cont_pts, double scale_factor, double max_incline)
{
    int order = 3;
    int num_intervals = num_cont_pts - order;
    double* constraints = new double[num_intervals];
    for(unsigned int i = 0; i < num_intervals; i++)
    {
        Eigen::Matrix<double,D,4> interval_control_points = cbind_help.array_section_to_eigen(cont_pts, num_cont_pts, i);
        double inclination_bound = evaluate_bound(interval_control_points, scale_factor);
        constraints[i] = inclination_bound - max_incline;
    }
    return constraints;
}


template <int D>
double InclinationConstraint<D>::evaluate_bound(Eigen::Matrix<double,D,4> &control_points, double &scale_factor)
{
    unsigned int dimension = 3;
    double max_z_vel = d_bounds_3.find_max_velocity_magnitude_in_single_dimension(control_points, scale_factor, dimension);
    Eigen::Matrix<double,2,4> horizontal_control_pts = control_points.block(0,0,2,4);
    std::array<double,2> min_horizontal_vel_and_time = d_bounds_2.find_min_velocity_and_time(horizontal_control_pts, scale_factor);
    double min_horizontal_vel = min_horizontal_vel_and_time[0];
    double incline_bound = max_z_vel / min_horizontal_vel;
    return incline_bound;
}

//Explicit template instantiations
template class InclinationConstraint<3>;