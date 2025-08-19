#include "CrossTermBounds.hpp"
#include "CubicEquationSolver.hpp"
#include <iostream>
#include <stdexcept>

template <int D>
CrossTermBounds<D>::CrossTermBounds()
{

}

template <int D>
double CrossTermBounds<D>::get_spline_curvature_bound(double cont_pts[], int &num_control_points)
{
    double max_curvature{0};
    double curvature;
    for (unsigned int i = 0; i < num_control_points-3; i++)
    {
        Eigen::Matrix<double,D,4> interval_control_points = cbind_help.array_section_to_eigen(cont_pts, num_control_points, i);
        curvature = evaluate_interval_curvature_bound(interval_control_points);
        if (curvature > max_curvature)
        {
            max_curvature = curvature;
        }
    }
    return max_curvature;
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

template <int D>
double CrossTermBounds<D>::evaluate_interval_curvature_bound(Eigen::Matrix<double,D,4> &control_points)
{
    double scale_factor = 1;
    std::array<double,2> min_velocity_and_time = d_dt_bounds.find_min_velocity_and_time(control_points,scale_factor);
    double min_velocity = min_velocity_and_time[0];
    double time_at_min_velocity = min_velocity_and_time[1];
    double max_cross_term = find_maximum_cross_term(control_points,scale_factor);
    std::array<double,2> max_acceleration_and_time = d_dt_bounds.find_max_acceleration_and_time(control_points, scale_factor);
    double max_acceleration = max_acceleration_and_time[0];
    double curvature_bound;
    if (min_velocity <= 1.0e-6)
    {
        double acceleration_at_min_vel = 
            d_dt_eval.calculate_acceleration_magnitude(time_at_min_velocity,
                control_points, scale_factor);
        if(acceleration_at_min_vel == 0)
        {
            curvature_bound = 0;
        }
        else
        {
            curvature_bound = std::numeric_limits<double>::max();
        }
    }
    else
    {
        curvature_bound = max_acceleration/(min_velocity*min_velocity);
        double curvature = max_cross_term/(min_velocity*min_velocity*min_velocity);
        if (curvature < curvature_bound)
        {
            curvature_bound = curvature;
        }
    }
    return curvature_bound;
}

template <int D>
double CrossTermBounds<D>::find_maximum_cross_term(Eigen::Matrix<double,D,4> &control_points, double &scale_factor)
{
    Eigen::Vector4d coeficients;
    if (D == 2) {coeficients = c_prop.get_2D_cross_coefficients(control_points);}
    else {coeficients = c_prop.get_3D_cross_coefficients(control_points);}
    double a_term = coeficients(0);
    double b_term = coeficients(1);
    double c_term = coeficients(2);
    double d_term = coeficients(3);
    std::array<double,3> roots = CubicEquationSolver::solve_equation(a_term,
        b_term, c_term, d_term);
    double t0 = 0;
    double tf = 1.0;
    double max_cross_term = c_eval.calculate_cross_term_magnitude(t0,control_points,scale_factor);
    double cross_term_at_tf = c_eval.calculate_cross_term_magnitude(tf,control_points,scale_factor);
    if (cross_term_at_tf > max_cross_term)
    {
        max_cross_term = cross_term_at_tf;
    }
    for(int index = 0; index < 3; index++)
    {
        double root = roots[index];
        if(root > 0 && root < 1.0)
        {
            double cross_term =  c_eval.calculate_cross_term_magnitude(root, control_points,scale_factor);
            if (cross_term > max_cross_term)
            {
                max_cross_term = cross_term;
            }
        }
    }
    return max_cross_term;
}

//Explicit template instantiations
template class CrossTermBounds<2>;
template class CrossTermBounds<3>;