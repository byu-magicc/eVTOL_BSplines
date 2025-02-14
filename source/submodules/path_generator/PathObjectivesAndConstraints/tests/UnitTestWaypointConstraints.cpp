#include "gtest/gtest.h"
#include "WaypointConstraints.hpp"

TEST(WaypointConstraintTests, VelocitiesCase2DStart)
{
    WaypointConstraints<2> waypoint_const{};
    double control_points[] = {3, 0, 9, 5, 7, 8, 0, 1,
                               1, 6, 2, 2, 1, 4, 1, 4};
    int num_control_points = 8;
    double inverse_scale_factor = 1/2.0;
    double true_constraints[] = {0,0};
    double desired_velocity[] = {1.5, 0.25};
    bool isStartPoint = true;
    double* constraints = waypoint_const.velocity_constraint(control_points, 
        num_control_points, desired_velocity, inverse_scale_factor, isStartPoint);
    double tolerance = 0.0000001;
    EXPECT_NEAR(constraints[0], true_constraints[0], tolerance);
    EXPECT_NEAR(constraints[1], true_constraints[1], tolerance);
}

TEST(WaypointConstraintOldTests, VelocitiesCase3DEnd)
{
    WaypointConstraints<3> waypoint_const{};
    double control_points[] = {-3,  -4, -2, -.5, 1  ,   0,  2, 3.5,
                               .5, 3.5,  6, 5.5, 3.7,   2, -1,  2,
                                1, 3.2,  5,   0, 3.3, 1.5, -1, 2.5};
    int num_control_points = 8;
    double inverse_scale_factor = 1/2.0;
    double true_constraints[] = {0, 0, 0};
    double desired_velocity[] = {0.875, 0, 0.25};
    bool isStartPoint = false;
    double* constraints = waypoint_const.velocity_constraint(control_points, 
        num_control_points, desired_velocity, inverse_scale_factor, isStartPoint);
    double tolerance = 0.0000001;
    EXPECT_NEAR(constraints[0], true_constraints[0], tolerance);
    EXPECT_NEAR(constraints[1], true_constraints[1], tolerance);
    EXPECT_NEAR(constraints[2], true_constraints[2], tolerance);
}

TEST(WaypointConstraintOldTests, AccelerationCase2DEnd)
{
    WaypointConstraints<2> waypoint_const{};
    double control_points[] = {3, 0, 9, 5, 7, 8, 0, 1,
                               1, 6, 2, 2, 1, 4, 1, 4};
    int num_control_points = 8;
    double inverse_scale_factor = 1/2.0;
    double true_constraints[] = {0,0};
    double desired_acceleration[] = {2.25, 1.5};
    bool isStartPoint = false;
    double* constraints = waypoint_const.acceleration_constraint(control_points, 
        num_control_points, desired_acceleration, inverse_scale_factor, isStartPoint);
    double tolerance = 0.0000001;
    EXPECT_NEAR(constraints[0], true_constraints[0], tolerance);
    EXPECT_NEAR(constraints[1], true_constraints[1], tolerance);
}

TEST(WaypointConstraintOldTests, AccelerationCase3DStart)
{
    WaypointConstraints<3> waypoint_const{};
    double control_points[] = {-3,  -4, -2, -.5, 1  ,   0,  2, 3.5,
                               .5, 3.5,  6, 5.5, 3.7,   2, -1,  2,
                                1, 3.2,  5,   0, 3.3, 1.5, -1, 2.5};
    int num_control_points = 8;
    double inverse_scale_factor = 1/2.0;
    double true_constraints[] = {0,0,0};
    double desired_acceleration[] = {0.75, -0.125, -0.1};
    bool isStartPoint = true;
    double* constraints = waypoint_const.acceleration_constraint(control_points, 
        num_control_points, desired_acceleration, inverse_scale_factor, isStartPoint);
    double tolerance = 0.0000001;
    EXPECT_NEAR(constraints[0], true_constraints[0], tolerance);
    EXPECT_NEAR(constraints[1], true_constraints[1], tolerance);
    EXPECT_NEAR(constraints[2], true_constraints[2], tolerance);
}
