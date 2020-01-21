# ADEKF
Automatic Differentiated Extended Kalman Filter (ADEKF) - This is a generic EKF Implementation that uses automatic differentiation to get rid of the need to define Jacobians. 

It is written as a C++ Header only Library. It runs on c++17

This implementation uses the ceres Jet implementation and their rotation header in a modified version.
http://ceres-solver.org/


Please cite xxx
if you use this implementation in a published work.


# Installation
## Dependencies
1.  Eigen 3 is required for the usage of this library
2.  Boost is an optional dependency. It is required if you want to use the later described ADEKF_MANIFOLD Macro

## Installation
Installation is not required since it is a header only library. You can simply #include the Headers.




# Usage 
To use the ADEKF you need:
1. a state description which defines, what you want to estimate

1. At least one of:
    1. a state transition model (dynamic model)  with optional input measurements (e.g. imu measurements)

    1. a measurement model (update model) with sensor measurements

    1. a pseudo measurement model (update model without a real sensor)
1. A start state
1. Covariance matrices for all of the above


## State Description
The ADEKF is based on Eigen. Thus as state you can use an Eigen Vector of any size. Additionally the ADEKF supports manifold states (e.q. quaternions). 
To use the ADEKF with a vector with DOF states simply create an Eigen Vector of the correct Size and pass it to the constructor. The ADEKF uses template argument deduction to deduce the correct internal types. 
Internally the ADEKF needs a representation of the state which works on different Scalar Types. Thus, you can choose your ScalarType yourself. 

The constructor takes the start state and the start covariance which needs to be an Eigen::Matrix with the squared dimension of the start state. 

```c++
#include <ADEKF.h>
#include <Eigen/Core>
int main(){
    // The size of the state
    constexpr int DOF=3; 
    // declare your start_state, optionally fill it with values
    Eigen::Matrix<double,DOF,1>state=state.Zero(); 
    // declare the start covariance
    Eigen::Matrix<double, DOF,DOF> covariance=covariance.Identity(); 
    // Set state's type as the state representation and its values as the start_state
    ADEKF::ADEKF adekf(state,covariance);  
}
```

Now we have declared an ADEKF with a vector state with 3 entries. This could be for example the velocity of a tracked object.


## State transition model
The ADEKF needs to evaluate models with different Scalar types. Thus, models need to be templated functors. An easy way to implement the models are the polymorphic lambdas of c++17.
```c++
... // declare adekf with your state (as in State Description)
auto dynamic model=[](auto & state, auto world_acceleration, double time_diff){
state+=acceleration*time_diff;
};
```
This code snippet declares 


