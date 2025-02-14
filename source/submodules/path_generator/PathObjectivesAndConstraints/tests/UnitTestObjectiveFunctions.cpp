#include "gtest/gtest.h"
#include "ObjectiveFunctions.hpp"


TEST(ObjectiveFunctionsTests, MinimizeJerkMagnitude)
{
    ObjectiveFunctions<3> objective_func{};
    double control_points[] = {2, 3, 2, 7, 8, 5, 6, 2,
                                6, 8, 3, 0, 3, 8, 2, 4,
                                9, 2, 0, 9, 8, 2, 9, 3};
    int num_control_points = 8;
    double scale_factor = 1;
    double true_objective = 2454;
    double objective = objective_func.minimize_jerk_magnitude(control_points, num_control_points);
    double tolerance = 0.0000001;
    EXPECT_NEAR(true_objective, objective, tolerance);
}

TEST(ObjectiveFunctionsTests, MinimizeAccelerationMagnitude)
{
    ObjectiveFunctions<2> objective_func{};
    double control_points[] = {0.46788998, -0.23394499,  0.46788998,  2.99996936,  5.53202098,  6.23398952, 5.53202096,
       -1.0601977,   0.02889528,  0.94461658,  1.25561647,  0.94470939,  0.02884888, -1.06010489};
    int num_control_points = 7;
    double scale_factor = 1;
    double true_objective = 11.818771273527938;
    double objective = objective_func.minimize_acceleration_magnitude(control_points, num_control_points);
    double tolerance = 0.0000001;
    EXPECT_NEAR(true_objective, objective, tolerance);
}

TEST(ObjectiveFunctionsTests, MinimizeVelocityMagnitude)
{
    ObjectiveFunctions<2> objective_func{};
    double control_points[] = {0.30056282, -0.14229062,  0.30944028,  2.99750765,  5.69037088,  6.14238504, 5.69924834,
                              -0.61897252, -0.06479341,  0.87814623,  0.90233775,  0.87285859, -0.06214961, -0.62426016};
    int num_control_points = 7;
    double scale_factor = 1;
    double true_objective = 17.665998040727704;
    double objective = objective_func.minimize_velocity_magnitude(control_points, num_control_points);
    double tolerance = 0.0000001;
    EXPECT_NEAR(true_objective, objective, tolerance);
}