#ifndef INCLINATIONCONSTRAINT_HPP
#define INCLINATIONCONSTRAINT_HPP
#include "DerivativeBounds.hpp"
#include "CBindHelperFunctions.hpp"

template <int D> // D is the dimension of the spline
class InclinationConstraint
{
    public:
        InclinationConstraint();
        double  get_spline_incline_constraint(double cont_pts[], int num_cont_pts, 
            double scale_factor, double max_incline);
        double*  get_interval_incline_constraints(double cont_pts[], int num_cont_pts, 
            double scale_factor, double max_incline);

    private:
        DerivativeBounds<3> d_bounds_3{};
        DerivativeBounds<2> d_bounds_2{};
        CBindHelperFunctions<D> cbind_help{};
        
        double evaluate_bound(Eigen::Matrix<double,D,4> &control_points, double &scale_factor);
};

extern "C"
{
    InclinationConstraint<3>* InclinationConstraint_3(){return new InclinationConstraint<3>();}
    double get_spline_incline_constraint_3(InclinationConstraint<3>* obj, double cont_pts[], int num_cont_pts, 
        double scale_factor, double max_incline){return obj->get_spline_incline_constraint(
        cont_pts, num_cont_pts,scale_factor, max_incline);}
    double* get_interval_incline_constraints_3(InclinationConstraint<3>* obj, double cont_pts[], int num_cont_pts, 
        double scale_factor, double max_incline){return obj->get_interval_incline_constraints(
        cont_pts, num_cont_pts,scale_factor, max_incline);}
}

#endif