# ADEKF
Automatic Differentiated Extended Kalman Filter (ADEKF) - This is a generic EKF Implementation that uses automatic differentiation to get rid of the need to define Jacobians. 

It is written as a C++ Header only Library. It runs on c++17

This implementation uses the ceres Jet implementation and their rotation header in a modified version.
http://ceres-solver.org/

 
Please cite us as: Lukas Post, Tom L. Koller "An Automatically Differentiating Extended Kalman Filter on Boxplus Manifolds", 2020 https://github.com/TomLKoller/ADEKF
if you use this implementation in a published work.


# Installation
## Dependencies
1.  Required: GCC >= 7 
1.  Required: Eigen 3 is used for matrix computations 
1.  Optional: Boost >= 1.37 It is required if you want to use certain helper macros
1.  Optional: CMake  >=3.0 It is used to compile the examples and to install the header files.

The ADEKF was originally tested on a ubuntu 18.04.4 LTS with :
1. GCC 7.4
1. Eigen 3
1. Boost 1.65.1
1. CMake 3.10.2



## Installation
Installation is not required since it is a header only library. You can simply #include the headers.

However, if you want convenient access to the library you can use the provided cmake file to install the headers via:
```
mkdir build
cd build
cmake ..
sudo make install 
```

You can then use the package ADEKF as it is shown in the [example cmake file](examples/velocity_without_orientation_tracking/CMakeLists.txt)


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
    adekf::ADEKF ekf(Eigen::Vector3d::Zero(),Eigen::Matrix3d::Identity());
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
    adekf::ADEKF ekf(Eigen::Vector3d::Zero(), Eigen::Matrix3d::Identity());  
}
```
Here we use the Eigen::Vector3d which is a vector with 3 entries. 


## State transition model
The ADEKF needs to evaluate models with different Scalar types. Thus, models need to be templated functors.
The general structure of a functor looks like :

```c++

struct dynamic_model{
    template<typename T>
    void operator()(STATE_TYPE<T> &state, ParameterPack ... params){
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
     adekf::ADEKF ekf(Eigen::Vector3d::Zero(), Eigen::Matrix3d::Identity());  
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
    adekf::ADEKF ekf(Eigen::Vector3d::Zero(), Eigen::Matrix3d::Identity());
    //Set processnoise
    Eigen::Matrix3d noise=noise.Identity();
    // set velocity vector in this case its just a vector of ones
    Eigen::Vector3d velocity=velocity.Ones();
    //Set time diff to 1/100Hz
    double time_diff=0.01;
    // call predict and pass velocity after the noise
    ekf.predict(dynamic_model_with_time_diff(),noise,velocity,time_diff);
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

If you use lambdas as models, take care to omit the () since you can not instantiate the lambda.

```c++
int main(){
     //Declare the dynamic model
    auto dynamic_model=[](auto & state, Eigen::Vector3d velocity, double time_diff){
        state=state+velocity*time_diff;
    };

    // Set state's type as the state representation and its values as the start_state
    adekf::ADEKF ekf(Eigen::Vector3d::Zero(), Eigen::Matrix3d::Identity());
    //Set processnoise
    Eigen::Matrix3d noise=noise.Identity();
    // set velocity vector in this case its just a vector of ones
    Eigen::Vector3d velocity=velocity.Ones();
    //Set time diff to 1/100Hz
    double time_diff=0.01;
    // call predict and pass velocity after the noise
    ekf.predict(dynamic_model,noise,velocity,time_diff);

}
```

## Measurement models

Measurement models are similar to state transition models. The main difference is, that they return a measurement expection instead of changing the state
The general structure of a functor looks like :

```c++

struct measurement_model{
    template<typename T>
    MEASUREMENT_TYPE operator()(STATE_TYPE<T> state, ParameterPack ... params){
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
     adekf::ADEKF ekf(Eigen::Vector3d::Zero(), Eigen::Matrix3d::Identity());
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
        state= (pos_in_mm+velocity*time_diff)/1000.;
    }
};
```
will fail at compilation time. The operators are evaluated with your chosen Scalar type and the ceres::Jet type http://www.ceres-solver.org/automatic_derivatives.html . But the Jet type can not be assigned to classic scalar types.
The corresponding error is the "Jet can not be assigned to double"

If you need local variables use the template T  instead:

```c++
struct dynamic_model_with_time_diff{
    template<typename T>
    void operator()(Eigen::Matrix<T,3,1> &state, Eigen::Vector3d velocity, double time_diff){
        Eigen::Matrix<T,3,1> pos_in_mm=state*1000.;
        state= (pos_in_mm+velocity*time_diff)/1000.;
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
This problem definetly appears at the use of SO3: 

So instead of
```c++
auto dynamic_model=[](auto & state, Eigen::Vector3d velocity, double time_diff){
    auto state_copy=state+state.orientation*velocity*time_diff;
    state=state_copy;
};
```
use .eval():
```c++
auto dynamic_model=[](auto & state, Eigen::Vector3d velocity, double time_diff){
    auto state_copy=(state+state.orientation*velocity*time_diff).eval();
    state=state_copy;
};
```


Since lambdas automatically deduce their return type, the same problem occurs at the return statement.


This call will fail :
```c++
auto measurement_model=[](auto & state, Eigen::Vector3d velocity, double time_diff){
    return state+state.orientation*velocity*time_diff;
};
```

Whereas this works:

```c++
auto measurement_model=[](auto & state, Eigen::Vector3d velocity, double time_diff){
    return (state+state.orientation*velocity*time_diff).eval();
};
```










## Manifolds (e.g. Quaternions)
In state estimation you often have the problem that you require something like a rotation quaternion in the state space.
An important property of the rotation quaternion is, that it always have to be normed. In the standard EKF this constraint will be destroyed which will end in a failed estimation of the orientation.

The ADEKF implements a general solution for this problem.  A quaternion is a so called boxplus manifold. For now, just think of a boxplus manifold of some mathematical structure that is represented by N numbers. The manifold now defines two operations, boxplus and boxminus, which define how to use these numbers in a Kalman Filter. 
The first operation boxplus, defines how a small vectorial change is applied on the manifold. We use quaternions as an example:

A quaternion consists of 4 numbers so N=4. Now you can rotate a quaternion with the measurement of a gyrometer which is an Euler-Angle Axis ( a vector with 3 entries). The boxplus of the quaternion now handles that the quaternion is rotated by the axis. 

In contrast, boxminus defines what the difference of two manifolds is. So for quaternions, we receive from  w=q2 boxminus q1 an Euler-Angle Axis w that defines how you have to rotate q1 to get to q2. So, q1 boxplus w =q2.


This principle can be used in an Kalman Filter to teach the filter how to interact with the manifold without destroying its properties. The ADEKF generically supports manifolds. For now, the quaternion (SO3) and the 2d rotation (SO2) are implemented. So to estimate the orientation you can just use the adekf::SO3 class as the state:

```c++
#include <ADEKF.h>
#include <Eigen/Core>
#include <types/SO3.h>
//Position estimate in m
int main(){
   // Set the start orientation to identity
    adekf::ADEKF ekf(adekf::SO3d(),Eigen::Matrix3d::Identity());
    //Create dynamic model which adds the angular velocity times time_diff to the position
    auto dynamic_model=[](auto & state, Eigen::Vector3d angular_rate, double time_diff){
        state= state + angular_rate * time_diff;
    };
    //create a measurement_model which gets a orientation measurement in mm
    auto measurement_model=[](auto state){
        return state;
    };
    //time step length
    double time_diff=0.01;
    // input measurement -> random velocity
    Eigen::Vector3d angular_rate=angular_rate.Random();
    //A orientation measurement -> random orientation;
    adekf::SO3d orientation(Eigen::Quaterniond::UnitRandom());
    //dynamic noise
    Eigen::Matrix3d process_noise=process_noise.Identity()*0.01;
    // measurement noise
    Eigen::Matrix3d measurement_noise=measurement_noise.Identity();
    //prediction:pass model, noise followed by the parameters you want in your dynamic _model (variadic acceptance)
    ekf.predict(dynamic_model, process_noise*time_diff, angular_rate, time_diff);
    //measurement update: pass model , noise and a measurement followed by parameters you need for your model  (variadic acceptance)
    ekf.update(measurement_model, measurement_noise, orientation);
}
```




The Interesting part here is the line:


```c++
    adekf::ADEKF ekf(adekf::SO3d(),Eigen::Matrix3d::Identity());    
```
While the quaternion requires 4 numbers, the passed covariance is only of dimension 3. The ADEKF requires only the real degrees of freedom of the quaternion (which are 3) in the covariance. Thus, the covariance has to be size 3.

## Compound States
States often consist of multiple sub states. For example, estimating the position of an object with an IMU requires not only to track the position, but also the orientation and the velocity.
The ADEKF ships with Macro support to setup those compound states.

Usually, a state can be divided into manifold substates m1 ... mn  and vector substates v1 ... vk.
The ADEKF provides the Macro ADEKF_MANIFOLD() to setup up the compound states. It takes 3 arguments:

```c++
ADEKF_MANIFOLD(CompoundName,Manifolds,Vectors)
```

The CompoundName is the name of the resulting Compound State Class. 
The Argument Manifolds declares all sub manifolds of the CompoundState It takes tuples of (ManifoldType, ManifoldName)
The ManifoldType is the class of the Manifold, e.g. adekf::SO3, whereas the name is the unique identifier to access the substate. Each tuple requires an additional bracket pair around it.  So the call:

```c++
ADEKF_MANIFOLD(DoubleQuaternion,((adekf::SO3,q1))((adekf::SO3,q2)))
```
Creates a compound state DoubleQuaternion which consists of 2SO3 parameterizations. 
The substates can be accessed via their names:
```c++
DoubleQuaternion<double> dq;
adekf::SO3d q1=dq.q1;
adekf::SO3d q2=dq.q2;
```


The Vectors parameter takes tuples aswell but this time of the form (VectorSize,VectorName). They do not require additional brackets and are concatenated by a "," (the parameter is variadic).


```c++
ADEKF_MANIFOLD(Pose3D,((adekf::SO3,orientation)),(3,position),(3,velocity) )
```

While the Vectors parameter is optional the Manifolds parameter is not. Thus, only compound states with Manifolds are supported.

The created CompoundState can be used as any other manifold as the state of the ADEKF, since the macro creates required operators and attributes.

To initialize a compound manifold just pass a value for each of the substates or use the default constructor of each:

```c++
Pose3D pose(adekf::SO3d(),Eigen::Vector3d::Zero(),Eigen::Vector3d()::Zero())
Pose3D pose_default();
```









## Non-additive noise
Sometimes it is usefull to have noise which applies on certain parameters. For example, if you have an gyrometer  the noise model is:

measured_angular_rate=angular_rate+noise

To use this, the ADEKF supports non additive noise with the functions predictWithNonAdditiveNoise() and
updateWithNonAdditiveNoise.
Those functions change two things in contrast to predict() and update.
1. They accept a slightly different model
    instead of the usual dynamic model predictWithNonAdditiveNoise accepts
    
```c++
struct dynamic_model{
    template<typename T>
    void operator()(STATE_TYPE<T> &state,Eigen::Vector<T,NOISE_DOF> noise,   ParameterPack ... params){
        state= ... //Implement dynamic model
    }
};
```
In your dynamic model simply add the noise vector where it applies.
The change for the update model is similar.
1. when calling predictWithNonAdditiveNoise and updateWithNonAdditiveNoise you set the dimension of the noise vector with the dimension of your passed covariance matrix.


So when we try to write the orientation estimation of the previous section we will use this code instead:


```c++
#include "ADEKF.h"
#include <Eigen/Core>
#include "types/SO3.h"
//Position estimate in m
void test(){
    // Set the start orientation to identity
    adekf::ADEKF ekf(adekf::SO3d(),Eigen::Matrix3d::Identity());
    //Create dynamic model which adds the angular velocity times time_diff to the position
    // With non-additive noise the ADEKF gives you a noise vector as second argument which you can use to tell it where the noise applies
    auto dynamic_model=[](auto & state, auto noise, Eigen::Vector3d angular_rate, double time_diff){
        state= state + (angular_rate+noise )* time_diff;
    };
    //create a measurement_model which gets a orientation measurement in mm
    auto measurement_model=[](auto state){
        return state;
    };
    //time step length
    double time_diff=0.01;
    // input measurement -> random velocity
    Eigen::Vector3d angular_rate=angular_rate.Random();
    //A orientation measurement -> random orientation;
    adekf::SO3d orientation(Eigen::Quaterniond::UnitRandom());
    //dynamic noise
    Eigen::Matrix3d process_noise=process_noise.Identity()*0.01;
    // measurement noise
    Eigen::Matrix3d measurement_noise=measurement_noise.Identity();
    //prediction:pass model, noise followed by the parameters you want in your dynamic _model (variadic acceptance)
    ekf.predictWithNonAdditiveNoise(dynamic_model, process_noise, angular_rate, time_diff);
    //measurement update: pass model , noise and a measurement followed by parameters you need for your model  (variadic acceptance)
    ekf.update(measurement_model, measurement_noise, orientation);
}
```


Now we dont have to multiply our process noise with time_diff anymore since it is multiplied with it in the dynamic model.



