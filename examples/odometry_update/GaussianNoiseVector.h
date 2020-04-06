//
// Created by tomlucas on 02.04.20.
//

#ifndef ADEKF_GAUSSIANNOISEVECTOR_H
#define ADEKF_GAUSSIANNOISEVECTOR_H

#include <Eigen/Core>
#include <random>

namespace adekf {

    /**
     * Class to create noise vectors with a gaussian random distribution.
     *
     * Use it like a usual eigen vector in computations. Call poll() to sample a new random vector.
     * @tparam Scalar type of the vector
     * @tparam size size of the vector
     */
    template<typename Scalar, int size>
    class GaussianNoiseVector : public Eigen::Matrix<Scalar, size, 1> {
        //Type of the vector
        typedef Eigen::Matrix<Scalar, size, 1> BASE;
        //generator for random numbers.
        static inline std::default_random_engine generator;
        //normal distribution to sample from
        std::normal_distribution<Scalar> distribution;
    public:
        //Mean and standard deviation of the gaussian distribution
        const Scalar mu, sigma;

        /**
         * Initialises the vectors random distribution and assigns a random vector
         * @param mu mean of the gaussian
         * @param sigma standard deviation of the gaussian
         */
        GaussianNoiseVector(Scalar mu, Scalar sigma) : distribution(mu, sigma), mu(mu), sigma(sigma) {
            this->poll();
        }
        /**
         * Samples a new gaussian random vector into this and also returns this
         */
        BASE poll() {
            //Call move assign of the base class
            BASE::operator=(this->unaryExpr([this](int) { return this->distribution(generator); }));
            return *this;
        }


    };

}

#endif //ADEKF_GAUSSIANNOISEVECTOR_H
