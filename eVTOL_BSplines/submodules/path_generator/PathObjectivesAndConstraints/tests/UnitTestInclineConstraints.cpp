#include "gtest/gtest.h"
#include "InclinationConstraint.hpp"

TEST(InclinationConstraintTests, InclineTest1)
{
    InclinationConstraint<3> inc_const{};
    double max_incline = 2;
    double true_constraint = 0.9938586931777061 - max_incline;
    int num_control_points = 4;
    double control_points[] =  {5,  8, 6, 7,
                               5, 11, 2, 7,
                               8,  7, 5, 9};
    double scale_factor = 1.5;
    double* constraints = inc_const.get_interval_incline_constraints(control_points, num_control_points, 
            scale_factor, max_incline);
    double tolerance = 0.000001;
    EXPECT_NEAR(true_constraint, constraints[0], tolerance);
}

TEST(InclinationConstraintTests, InclineTest2)
{
    InclinationConstraint<3> inc_const{};
    double max_incline = 2;
    double true_constraint_1 = 0.9938586931777061 - max_incline;
    double true_constraint_2 =  17.70564330959397 - max_incline;
    int num_control_points = 5;
    double control_points[] =  {5,  8, 6, 7, 3,
                               5, 11, 2, 7, 5,
                               8,  7, 5, 9, 2};
    double scale_factor = 1.5;
    double* constraints = inc_const.get_interval_incline_constraints(control_points, num_control_points, 
            scale_factor, max_incline);
    double tolerance = 0.000001;
    EXPECT_NEAR(true_constraint_1, constraints[0], tolerance);
    EXPECT_NEAR(true_constraint_2, constraints[1], tolerance);
}

TEST(InclinationConstraintTests, InclineTestWholeSpline)
{
    InclinationConstraint<3> inc_const{};
    double max_incline = 2;
    double true_constraint =  17.70564330959397 - max_incline;
    int num_control_points = 5;
    double control_points[] =  {5,  8, 6, 7, 3,
                               5, 11, 2, 7, 5,
                               8,  7, 5, 9, 2};
    double scale_factor = 1.5;
    double constraint = inc_const.get_spline_incline_constraint(control_points, num_control_points, 
            scale_factor, max_incline);
    double tolerance = 0.000001;
    EXPECT_NEAR(true_constraint, constraint, tolerance);
}