/*!
 * \file ope/proximity_detector.hpp
 * \brief Proximity Detector
 * \author Quentin Jouet
 */

#ifndef PROBLEM_PROXIMITY_DETECTOR_HPP
#define PROBLEM_PROXIMITY_DETECTOR_HPP



//TODO
// #include "floe/ope/proxy-detector "
//TODO

namespace floe { namespace ope
{

/*! ProximityDetector
 *
 * It represents the whole problem of moving N floes in interval time [0, T].
 *
 * \tparam TFloe  
 * \tparam TDetector 
 * \tparam TProxymityDetector   
 *
 */

template <
    typename TDetector_h
>
class ProximityDetector
{

public:

    TDetector_h m_detector_h;

private:

};

}} // namespace floe::ope


#endif // PROBLEM_PROXIMITY_DETECTOR_HPP