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
2.  Boost is an optional dependency. It is required if you want to use certain helper macros

## Installation
Installation is not required since it is a header only library. You can simply #include the headers.

## Prerequisites
To use this software you should know how to:
1. program in C++ 
1. design the state and model of an EKF
1. use Eigen library http://eigen.tuxfamily.org/index.php?title=Main_Page

We will not show the above explicitly in our tutorials. 


In the following we will describe the general usage of the ADKEF. For executable code with a sensefull meaning see our tutorials.

# Usage 
To use the ADEKF you need:
1. a state description which defines, what you want to estimate

1. At least one of:
    1. a state transition model (dynamic model)  with optional input measurements (e.g. imu measurements)

    1. a measurement model (update model) with sensor measurements

    1. a pseudo measurement model (update model without a real sensor)
1. A start state
1. Covariance matrices for all of the above


We will show you the complete way to define any of those and an easy way for each which is faster to code.

For Quickstarts who do not want to read this, this is how a code can look like. 

```c++
#include <ADEKF.h>
#include <Eigen/Core>
//Position estimate in m
int main(){
    // Set the start position to 0, 0 ,0
    adekf::ADEKF ekf(Eigen::Vector3d::Zero());
    //Create dynamic model which adds the velocity times time_diff to the position
    auto dynamic_model=[](auto & state, Eigen::Vector3d velocity, double time_diff){
        state=state+velocity*time_diff;
    };
    //create a measurement_model which gets a position measurement in mm
    auto measurement_model=[](auto state){
        return state*1000.;
    };
    //time step length
    double time_diff=0.01;
    // input measurement -> random velocity
    Eigen::Vector3d velocity=velocity.Random();
    //A position measurement -> random position;
    Eigen::Vector3d position=position.Random();
    //dynamic noise
    Eigen::Matrix3d process_noise=process_noise.Identity()*0.01;
    // measurement noise
    Eigen::Matrix3d measurement_noise=measurement_noise.Identity();
    //prediction:pass model, noise followed by the parameters you want in your dynamic _model (variadic acceptance)
    ekf.predict(dynamic_model,process_noise,velocity,time_diff);
    //measurement update: pass model , noise and a measurement followed by parameters you need for your model  (variadic acceptance)
    ekf.update(measurement_model,measurement_noise,position);
}
```



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
    // declare your start_state. Fill it with values
    Eigen::Matrix<double,DOF,1>state=state.Zero(); 
    // declare the start covariance
    Eigen::Matrix<double, DOF,DOF> covariance=covariance.Identity(); 
    // Set state's type as the state representation and its values as the start_state
    adekf::ADEKF ekf(state,covariance);  
}
```

Now we have declared an ADEKF with a vector state with 3 entries. This could be for example the position of a tracked object.
As a short version you can also write:

```c++
int main(){
    // Set state's type as the state representation and its values as the start_state
    adekf::ADEKF ekf(Eigen::Vector3d::Zero());  
}
```
Here we use the Eigen::Vector3d which is a vector with 3 entries. If no covariance is passed to the adekf, it will initialize it with the identity matrix.


## State transition model
The ADEKF needs to evaluate models with different Scalar types. Thus, models need to be templated functors.
The general structure of a functor looks like :

```c++

struct dynamic_model{
    template<typename T>
    void operator()(STATE_TYPE &state, ParameterPack ... params){
        state= ... //Implement dynamic model
    }
};
```
The functor requires the template T, which is the Scalar type. Type this anywhere where you would write double or float. This is required for the automatic differentiation. The ParameterPack is a variadic template. You can pass as any arguments to your dynamic model as you need. These do not have to be templated with type T. Take care that you write an & in front of state. The dynamic model has to overwrite the given state.
For example if you want a dynamic model which adds velocity to your position state then you define the functor as:


```c++
struct dynamic_model{
    template<typename T>
    void operator()(Eigen::Matrix<T,3,1> &state, Eigen::Vector3d velocity){
        state= state+velocity;
    }
};
```

To call the predict step of the adkef, we need the process noise which has the same dimension as the state.


```c++
int main(){
    // Set state's type as the state representation and its values as the start_state
    adekf::ADEKF ekf(Eigen::Vector3d::Zero());
    //Set processnoise
    Eigen::Matrix3d noise=noise.Identity();
    // set velocity vector in this case its just a vector of ones
    Eigen::Vector3d velocity=velocity.Ones();
    // call predict and pass velocity after the noise
    ekf.predict(dynamic_model(),noise,velocity);
}
```

At the call of predict you have to pass all the parameters you require in your dynamic model after the process noise. 

If you want to change your model to  multiply the velocity by the time difference before you add it to the position you could just add the time difference as another parameter:

```c++
struct dynamic_model_with_time_diff{
    template<typename T>
    void operator()(Eigen::Matrix<T,3,1> &state, Eigen::Vector3d velocity, double time_diff){
        state= state+velocity*time_diff;
    }
};
```
And pass the time_diff later to the predict function:

```c++
int main(){
    // Set state's type as the state representation and its values as the start_state
    adekf::ADEKF ekf(Eigen::Vector3d::Zero());
    //Set processnoise
    Eigen::Matrix3d noise=noise.Identity();
    // set velocity vector in this case its just a vector of ones
    Eigen::Vector3d velocity=velocity.Ones();
    //Set time diff to 1/100Hz
    double time_diff=0.01;
    // call predict and pass velocity after the noise
    ekf.predict(dynamic_model(),noise,velocity,time_diff);
}
```

Please read Pittfalls with local variables if you need them in your model.

The templated functors are currently the safest way to define models. 
An easier and faster way are c++17 polymorphic lambdas. The same dynamic model as before can be easily implemented by:

```c++
auto dynamic_model=[](auto & state, Eigen::Vector3d velocity, double time_diff){
    state=state+velocity*time_diff;
};
```
The auto keyword in front of state enables the automatic differentiation. Again take care to write & before state.
Please read Pitfalls with lambdas if you want to use them.

## Measurement models

Measurement models are similar to state transition models. The main difference is, that they return a measurement expection instead of changing the state
The general structure of a functor looks like :

```c++

struct measurement_model{
    template<typename T>
    MEASUREMENT_TYPE operator()(STATE_TYPE state, ParameterPack ... params){
        return= ... //Implement measurement model
    }
};
```

The template parameters have the same meaning as in the dynamic model. Again you can pass as many auxilary parameters as you need. Do not type an & in front of state, since your are not supposed to change the state
For example if you want a measurement model which where you measure the actual position in mm while your state contains m:


```c++
struct measurement_model{
    template<typename T>
    Eigen::Matrix<T,3,1> operator()(Eigen::Matrix<T,3,1>  state){
        return state*1000.;
    }
};
```

To call the update step of the adkef, we need the measurement noise which has the same dimension as the measurement and a measurement observation.


```c++
int main(){
    // Set state's type as the state representation and its values as the start_state
    adekf::ADEKF ekf(Eigen::Vector3d::Zero());
    //Set processnoise
    Eigen::Matrix3d noise=noise.Identity();
    // set position vector in this case its just a vector of ones
    Eigen::Vector3d position=position.Ones();
    // call predict and pass the model, the noise and the measurement observation. If your models requires auxilary parameters pass them after position
    ekf.update(measurement_model(),noise,position);
}
```

Please read Pittfalls with local variables if you need them in your model.


The measurement models can also be defined by  poylmorphic lambda expressions:

```c++
auto measurement_model=[](auto state){
    return state*1000.;
};
```

Please read Pitfalls with lambdas if you want to use them.

## Pitfalls with local variables 
 Be careful that you do not create variables with fixed scalar type inside the model:
Calls like:
```c++
struct dynamic_model_with_time_diff{
    template<typename T>
    void operator()(Eigen::Matrix<T,3,1> &state, Eigen::Vector3d velocity, double time_diff){
        Eigen::Vector3d pos_in_mm=state*1000.;
        state= (pos+velocity*time_diff)/1000.;
    }
};
```
will fail at compilation time. The operators are evaluated with your chosen Scalar type and the ceres::Jet type http://www.ceres-solver.org/automatic_derivatives.html . But the Jet type can not be assigned to classic scalar types.

If you need local variables use the template T  instead:

```c++
struct dynamic_model_with_time_diff{
    template<typename T>
    void operator()(Eigen::Matrix<T,3,1> &state, Eigen::Vector3d velocity, double time_diff){
        Eigen::Matrix<T,3,1> pos_in_mm=state*1000.;
        state= (pos+velocity*time_diff)/1000.;
    }
};
```


## Pitfalls with lambdas
Polymorphic lambdas allow much shorter code when defining the models. You dont have to care about templating since the auto keyword allows this automatically. But you have to be carefull, when you implemenent models. 


Unfortunately, Eigen and auto do not work well together. If you assign a calculation result  to a local variable make sure to call .eval() on it. The problem is, that Eigen uses lazy evaluation for  computation speed ups. Thus, the result of arithmetic operations is often an expression class instead of a real matrix or vector. The expression will be evaluated when you assign it to a matrix. But if you assign it to an auto variable, it may happen that the data the expression requires is already deleted.   See https://eigen.tuxfamily.org/dox/TopicPitfalls.html for details. As a workaround call .eval() on expressions before you assign them to an auto variable.

May cause state to not be updated:
```c++
auto dynamic_model=[](auto & state, Eigen::Vector3d velocity, double time_diff){
    auto state_copy=state+velocity*time_diff;
    state=state_copy;
};
```
But it works safely with eval:
```c++
auto dynamic_model=[](auto & state, Eigen::Vector3d velocity, double time_diff){
    auto state_copy=(state+velocity*time_diff).eval();
    state=state_copy;
};
```





