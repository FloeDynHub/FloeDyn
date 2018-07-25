/*!
 * \file utils/random.hpp
 * \brief Random utilities declaration
 * \author Quentin Jouet
 */

#ifndef FLOE_UTILS_RANDOM_HPP
#define FLOE_UTILS_RANDOM_HPP

#include <random>

namespace floe { namespace random
{

unsigned long get_unique_seed();

void initiate_generator(std::default_random_engine& generator);

std::default_random_engine get_uniquely_seeded_generator();

std::string gen_random(const int len);

}} // namespace floe::random


#endif // FLOE_UTILS_RANDOM_HPP