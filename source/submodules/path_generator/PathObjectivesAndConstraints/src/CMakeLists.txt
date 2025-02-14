add_library(${PROJECT_NAME} SHARED 
    BsplineToMinvo.cpp
    BsplineToMinvo.hpp
    CubicEquationSolver.cpp
    CubicEquationSolver.hpp
    CBindHelperFunctions.cpp
    CBindHelperFunctions.hpp
    RotationHelperFunctions.cpp
    RotationHelperFunctions.hpp
    MDMAlgorithmClass.cpp
    MDMAlgorithmClass.hpp
    DerivativeEvaluator.cpp 
    DerivativeEvaluator.hpp
    WaypointConstraintsOld.cpp
    WaypointConstraintsOld.hpp
    WaypointConstraints.cpp
    WaypointConstraints.hpp
    ObjectiveFunctions.cpp
    ObjectiveFunctions.hpp
    ObstacleConstraints.cpp
    ObstacleConstraints.hpp
    CurvatureConstraints.cpp
    CurvatureConstraints.hpp
    DerivativeBounds.cpp
    DerivativeBounds.hpp
    CrossTermEvaluator.cpp
    CrossTermEvaluator.hpp
    CrossTermProperties.cpp
    CrossTermProperties.hpp
    CrossTermBounds.hpp
    CrossTermBounds.cpp
    InclinationConstraint.cpp
    InclinationConstraint.hpp
    SphereCollisionEvaluator.cpp
    SphereCollisionEvaluator.hpp)

target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")