/**
 * @file fem_problem.hpp
 * @author Silouane 
 * @brief fem problem, assembles and solves a linear elasticity problem 
 * @version 0.1
 * @date 2023-10-13
 * 
 * 
 */

#ifndef FEM_PROBLEM
#define FEM_PROBLEM

#include <type_traits>
#include "floe/integration/quadrature.hpp"
#include "floe/floes/static_floe.hpp"
#include "floe/floes/floe_exception.hpp"

#include "floe/state/space_time_state.hpp"
#include "floe/geometry/frame/frame_transformers.hpp"

#include "floe/floes/floe_h.hpp"
#include "floe/geometry/arithmetic/dot_product.hpp"
#include "floe/geometry/arithmetic/arithmetic.hpp"

#include "floe/floes/floe_interface.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>


namespace floe { namespace fem {

/**
 * @brief 
 * 
 * @tparam TFloeGroup 
 * @tparam TProxymityDetector 
 * @tparam TCollisionManager 
 * @tparam TDynamicsManager 
 * @tparam TDomain 
 */
template <
    typename TFloeGroup
    // typename TProxymityDetector,
    // typename TCollisionManager,
    // typename TDynamicsManager,
    // typename TDomain
>
class FemProblem
{
public :
    
    using real_type = typename TFloeGroup::real_type;
    using point_type = typename TFloeGroup::point_type;
    // using floe_group_type = TFloeGroup; 
    // using time_scale_manager_type = TDomain::TimeScaleManager<typename TProxymityDetector::proximity_data_type>;
    using floe_type = typename TFloeGroup::floe_type;
    typedef typename floe_type::mesh_type     mesh_type;

    //! Default constructor.
    FemProblem();
    FemProblem(floe_type * floe): m_floe{floe}{}
    
    bool prepare();
    bool assembleK();
    bool assembleF()
    {
        return true; 
    };
    real_type computeKE()
    {
        return true; 
    };
    bool solve();

private : 
    floe_type * m_floe;
    size_t m_nE;
    size_t m_nN;
    size_t m_nDof; // number of degrees of freedom. 2 translations for instance 
    Eigen::SparseMatrix<real_type> m_K;
    Eigen::SparseMatrix<real_type> m_F;
    Eigen::SparseMatrix<real_type> m_Sol;
    
};


template < typename TFloeGroup>
bool
FemProblem<TFloeGroup>::assembleK()
{
    return true; 
};


template < typename TFloeGroup>
bool
FemProblem<TFloeGroup>::prepare()
{
    m_nE = m_floe->mesh().get_n_cells();
    m_nN = m_floe->mesh().get_n_nodes();
    m_nDof = 2; // 2 displacements on each node. 
    
    if (!assembleK())
        return false;
    if (!assembleF())
        return false; 
    return true; 
};


}} // namespace floe::fem

#endif
