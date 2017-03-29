/*!
 * \file utils/random.cpp
 * \brief Random utilities implementation
 * \author Quentin Jouet
 */

#ifndef FLOE_UTILS_RANDOM_CPP
#define FLOE_UTILS_RANDOM_CPP

#include "floe/utils/random.hpp"
#include <chrono>


namespace floe { namespace random
{


unsigned long get_unique_seed(){
    // random seed for not having same distrib (microseconds since epoch)
    // (even for quasi-simultaneous execs, where time in seconds is not sufficient)
    return std::chrono::system_clock::now().time_since_epoch() / std::chrono::microseconds(1);
}


void initiate_generator(std::default_random_engine& generator){
    generator.seed(get_unique_seed());
}


std::default_random_engine get_uniquely_seeded_generator(){
    std::default_random_engine generator;
    initiate_generator(generator);
    return generator;
}


}} // namespace floe::random


#endif // FLOE_UTILS_RANDOM_CPP