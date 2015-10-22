/*!
 * \file ope/proximity_detector.hpp
 * \brief Proximity Detector
 * \author Quentin Jouet
 */

#ifndef PROBLEM_PROXIMITY_DETECTOR_HPP
#define PROBLEM_PROXIMITY_DETECTOR_HPP
 

namespace floe { namespace ope
{

/*! ProximityDetector
 *
 * \tparam TDetector_h 	Discrete level detector type
 *
 */

template <
    typename TDetector_h
>
class ProximityDetector
{

public:
    using detector_h_type = TDetector_h;

    TDetector_h m_detector_h;

private:

};

}} // namespace floe::ope


#endif // PROBLEM_PROXIMITY_DETECTOR_HPP