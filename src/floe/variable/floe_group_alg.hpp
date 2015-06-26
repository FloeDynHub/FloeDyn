/*!
 * \file variable/floe_group_alg.hpp
 * \brief Algebraic Floe Configuration class
 * \author Quentin Jouet
 */

#ifndef VARIABLE_FLOES_ALG_HPP
#define VARIABLE_FLOES_ALG_HPP


#include "floe/variable/floe_alg.hpp"
#include <vector>
#include "floe/lcp/lcp.hpp"



namespace floe { namespace variable
{

/*! FloeGroup_alg
 *
 * It represents the floe configuration at the algebraic level.
 *
 */
template <
    typename TFloe_alg,
    typename TLCP = floe::lcp::LCP<double> // TODO find another way to get this type
>
class FloeGroup_alg
{

public:

    using LCP_type = TLCP;

    //! Default constructor.
    // FloeGroup_alg()

    std::vector<TFloe_alg*> m_list_floe_alg;

    inline void add_floe( TFloe_alg& floe )
        {
        }

    inline void add_LCP( TLCP& lcp){
        m_list_LCP.push_back( lcp );
    }

    inline std::vector<TLCP>& get_list_LCP(){ return m_list_LCP; }
    inline void empty_list_LCP() { m_list_LCP.clear(); }

private:

    std::vector<TLCP> m_list_LCP;

};

}} // namespace floe::variable


#endif // VARIABLE_FLOES_ALG_HPP
