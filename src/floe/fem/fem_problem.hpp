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


#define WHEREAMI std::cout << std::endl << "no crash until line " << __LINE__ << " in the file " __FILE__ << std::endl;



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
 * @tparam TFloe
 * @tparam TProxymityDetector 
 * @tparam TCollisionManager 
 * @tparam TDynamicsManager 
 * @tparam TDomain 
 */
template <
    typename TFloe
    // typename TProxymityDetector,
    // typename TCollisionManager,
    // typename TDynamicsManager,
    // typename TDomain
>
class FemProblem
{
public :
    
    // using point_type = typename TFloeGroup::point_type;
    // // using floe_group_type = TFloeGroup; 
    // // using time_scale_manager_type = TDomain::TimeScaleManager<typename TProxymityDetector::proximity_data_type>;
    // using floe_type = typename TFloeGroup::floe_type;
    // typedef typename floe_type::mesh_type     mesh_type;

    using real_type = typename TFloe::real_type;
    using point_type = typename TFloe::point_type;
    // using floe_type = typename TFloe;
    typedef TFloe floe_type;
    typedef typename floe_type::mesh_type     mesh_type;
    using multi_point_type = typename mesh_type::multi_point_type;
    using connectivity_type = typename mesh_type::connectivity_type;


    //! Default constructor.
    FemProblem(): m_floe{nullptr}, m_is_prepared{false}, m_largest_value{0} {}
    FemProblem(floe_type * floe): 
        m_floe{floe}, 
        m_E{9000000000.0}, // ice, from https://tc.copernicus.org/articles/17/3883/2023/tc-17-3883-2023.pdf 
        // m_E{10000000000.0}, // for verification purposes, to compare with analytical solution from mecagora 
        m_nu{0.3}, // from https://tc.copernicus.org/articles/17/3883/2023/tc-17-3883-2023.pdf  
        m_nE{m_floe->mesh().get_n_cells()},
        m_nN{m_floe->mesh().get_n_nodes()},
        m_nDof{2}, // 2 displacements at each node 
        m_nNodesPerElt{3}, // linear triangle elements only, for now... 
        m_Solution{Eigen::SparseMatrix<real_type> (1, 1)},
        m_Stress{Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic> (m_nE, 3)},
        m_elasticEnergies{Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic> (m_nE, 1)},
        m_is_prepared{false}, 
        m_stress_is_computed{false},
        m_largest_value{0},
        m_last_total_impulse{0},
        m_fracture_points {0,0}
    {}
    
    
    bool prepare();
    bool assembleK();
    bool addDirichlet(std::vector<size_t> gamma_d, std::vector<point_type> values);
    bool addContactDirichlet(std::vector<size_t> gamma_d, std::vector<point_type> values);

    bool assembleF()
    {
        m_F = Eigen::SparseMatrix<real_type> (m_nN*m_nDof, 1); // this way, F is initialized at zero 
        return true; 
    };
   
    bool solve();
    bool performComputation(std::vector<size_t> gamma_d, std::vector<point_type> values);
    inline Eigen::SparseMatrix<real_type> get_solution() const {return m_Solution;};

    /**
     * @brief Get the solution vector object. Same content as m_Solution, but in an std::vector instead of Eigen::SparseMatrix
     * 
     * @return std::vector<real_type>. The first m_nN components correspond to u (displacement along e_x), the last m_nN ones to v (displacement along e_y). Nope, not any more. 
     */
    std::vector<real_type> get_solution_vector() const 
    {
        std::vector<real_type> v;
        v.resize(m_nN*2);
        if (m_Solution.rows() < 2*m_nN)
        {
            std::cout << "problem with matrix size" << std::endl;
            WHEREAMI
            return v;
        }
        for (size_t iPoint = 0; iPoint < m_nN*2 ; ++iPoint)
        {
            v[iPoint] = m_Solution.coeff(iPoint,0);
        }
        return v;
    };

    std::vector<real_type> get_stress_vector() const 
    {
        std::vector<real_type> v;
        if ((m_Solution.rows() != 2*m_nN) || (!m_stress_is_computed))
            return v;
        v.resize(m_nE*3);
        for (size_t iElem = 0 ; iElem < m_nE ; ++iElem)
        {
            v[iElem] = m_Stress(iElem, 0);
            v[iElem+m_nE] = m_Stress(iElem, 1);
            v[iElem+2*m_nE] = m_Stress(iElem, 2);
        }
        return v;
    };

    bool compute_stress_vector() 
    {
        if (m_Solution.rows() != 2*m_nN)
            return false;
        Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic> stress(m_nE, 3);
        Eigen::Matrix<real_type, 3,1> sigma;
        Eigen::Matrix<real_type, 6,1> sol_elem;
        connectivity_type connect; 
        connect = m_floe->mesh().connectivity();
        Eigen::Matrix<real_type, 3, 3> H = computeH(); 
        Eigen::Matrix<real_type, 3, 6> B;
        for (size_t iElem = 0 ; iElem < m_nE ; ++iElem)
        {
            B = computeB(iElem); 
            sol_elem(0,0) = m_Solution.coeff(connect[iElem][0]*2, 0);
            sol_elem(1,0) = m_Solution.coeff(connect[iElem][0]*2+1, 0);
            sol_elem(2,0) = m_Solution.coeff(connect[iElem][1]*2, 0);
            sol_elem(3,0) = m_Solution.coeff(connect[iElem][1]*2+1, 0);
            sol_elem(4,0) = m_Solution.coeff(connect[iElem][2]*2, 0);
            sol_elem(5,0) = m_Solution.coeff(connect[iElem][2]*2+1, 0);
            sigma = H*B*sol_elem;
            stress.block(iElem, 0, 1, 3) = sigma.transpose();
        } 
        m_Stress = stress;
        return true;
    };

    bool compute_elementary_energies()
    {
        if (!m_stress_is_computed)
            return false;
        m_elasticEnergies = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>(m_nE, 1);

        connectivity_type connect = m_floe->mesh().connectivity();
        floe::integration::RefGaussLegendre<real_type,2,2> i;
        Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic> trianglePoints = i.pointsAndWeights(); 
        size_t nPoints = trianglePoints.rows();

        Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic> coord = trianglePoints.block(0,0,nPoints, 2);
        Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic> weights = trianglePoints.block(0,2,nPoints, 1);

        for (size_t iElem = 0; iElem < m_nE; iElem++)
        {
            Eigen::Matrix<real_type, 1,1> e ;
            e(0,0) = 0.0;
            real_type detJac = m_floe->mesh().get_jacobian(iElem);
            Eigen::Matrix<real_type, 6, 1> u;
            u(0,0) = m_Solution.coeff(connect[iElem][0]*2,   0);
            u(1,0) = m_Solution.coeff(connect[iElem][0]*2+1, 0);
            u(2,0) = m_Solution.coeff(connect[iElem][1]*2,   0);
            u(3,0) = m_Solution.coeff(connect[iElem][1]*2+1, 0);
            u(4,0) = m_Solution.coeff(connect[iElem][2]*2,   0);
            u(5,0) = m_Solution.coeff(connect[iElem][2]*2+1, 0);
            Eigen::Matrix<real_type, 3, 6> B = computeB(iElem); 
            Eigen::Matrix<real_type, 3, 1> epsilon;
            epsilon = B*u;
            Eigen::Matrix<real_type, 1, 3> stress = m_Stress.block(iElem, 0, 1, 3);
            e = 0.5*detJac*0.5*stress*epsilon; 
            // Epe = 1/2*\int_Omega^e sigma*epsilon dOmega^e. Ici detJac = surface de l'élément*2. Stress est constant sur l'élément, pas besoin de plus de points. 
            // Plein d'autre manières de faire, stress*H.inverse()*stress.transpose(), epsilon.transpose()*H*epsilon, etc. 
            m_elasticEnergies(iElem,0) = e.coeff(0,0);
        }
        return true ;
    }

    real_type compute_elastic_energy(std::vector<size_t> elems) 
    { 
        real_type e(0.0);
        if (!m_stress_is_computed)
            return e;
        for (size_t iElem = 0; iElem < elems.size() ; ++iElem)
        {
            e+=m_elasticEnergies.coeff(elems[iElem],0);
        }
        return e;
    }; 

    inline bool unset_prepared() {m_is_prepared = false;};
    Eigen::Matrix<double, 3, 6> computeB(size_t iElem) const;
    Eigen::Matrix<double, 3, 3> computeH() const;

    bool look_for_fracture()
    {
        // au menu ici : 
        // * calculer les énergies élémentaires
        // * tester plein de droites
        // * faire des listes d'éléments 
        // * calculer les énergies 
        // * si on tombe sur une énergie plus faible avec fracture, on ajoute les deux noeuds dans m_fracture_points et on renvoit true
        // * sinon, on met deux fois zero dans m_fracture points et on renvoit false
        // * question : on casse sur des noeuds ? Des endroits arbitraires ?  
        m_fracture_points.clear();
        m_fracture_points.push_back(0);
        m_fracture_points.push_back(0);
        return false;
    }

private : 
    floe_type * m_floe;
    real_type m_E; // Young modulus 
    real_type m_nu; // poisson coefficient 
    size_t m_nE; // number of elements in the mesh 
    size_t m_nN; // number of nodes in the mesh 
    size_t m_nDof; // number of degrees of freedom. 2 translations for instance 
    size_t m_nNodesPerElt; // number of nodes per element. A priori only triangles will be used but on ne sait jamais. 
    std::vector<Eigen::Triplet<real_type>> m_KTriplet;
    std::vector<Eigen::Triplet<real_type>> m_DirichletTriplet;
    std::vector<Eigen::Triplet<real_type>> m_FTriplet;
    Eigen::SparseMatrix<real_type> m_K;
    Eigen::SparseMatrix<real_type> m_F;
    Eigen::SparseMatrix<real_type> m_Solution;
    Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic> m_Stress;
    Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic> m_elasticEnergies;
    bool m_is_prepared;
    bool m_stress_is_computed;
    bool m_energies_are_computed;
    real_type m_largest_value; 
    real_type m_last_total_impulse; 
    std::vector<size_t> m_fracture_points;
};


template < typename TFloe>
bool
FemProblem<TFloe>::assembleK()
{
    Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic> Ke;
    real_type detJac(0);		
    floe::integration::RefGaussLegendre<double,2,2> i;
    Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic> gp = i.pointsAndWeights();// getGauss(22, 5);// gauss points and weights for current element
    Eigen::Matrix<real_type, 3, 6> B; // idem for gradient matrix
    Eigen::Matrix<real_type, 3,3> H;

    size_t ngp = gp.rows(); // number of integration points 
    size_t gw(2); // indicates the column of gp associated with the weight. (ex: for a 2D case, the columns 0 and 1 contain x_i and y_i coordinates, column 2 contains w_i )
    multi_point_type coordinates;
    connectivity_type connect; 
    coordinates = m_floe->mesh().points();
    connect = m_floe->mesh().connectivity();

    m_KTriplet.resize(m_nE*m_nNodesPerElt*m_nDof*m_nNodesPerElt*m_nDof);
    size_t k(0);
    H = computeH(); 
    // loop over the elements, construction of elementary matrices 
    for (size_t iElem = 0 ; iElem < m_nE ; ++iElem)
    { 
        detJac = m_floe->mesh().get_jacobian(iElem);
        if (detJac == 0) {std::cerr << "in FemProblem::prepare, could not get jacobian of element" << iElem << std::endl; }

        B = computeB(iElem); 
        Ke = detJac/2*B.transpose()*H*B; 
        for (unsigned int i = 0 ; i < m_nNodesPerElt ; i++)
        {
            size_t globalI = m_nDof*connect[iElem][i];
            for (unsigned int j = 0; j < m_nNodesPerElt ; j++ )
            {
                size_t globalJ = m_nDof*connect[iElem][j];
                m_KTriplet[k++] = Eigen::Triplet<real_type>(globalI, globalJ, Ke.coeff(2*i,2*j)); 
                m_KTriplet[k++] = Eigen::Triplet<real_type>(globalI+1, globalJ, Ke.coeff(2*i+1,2*j)); 
                m_KTriplet[k++] = Eigen::Triplet<real_type>(globalI+1, globalJ+1, Ke.coeff(2*i+1,2*j+1)); 
                m_KTriplet[k++] = Eigen::Triplet<real_type>(globalI, globalJ+1, Ke.coeff(2*i,2*j+1)); 
                if((isnan((Ke)(i,j))) || (isnan((Ke)(i+1,j))) || (isnan((Ke)(i+1,j+1))) || (isnan((Ke)(i+1,j+1))))
                {
                    std::cerr << "Found a nan in Ke, element " << iElem << ". Stopping right now. " <<std:: endl << Ke << std::endl;
                    return false; 	
                }
                // preparing Dirichlet condition 
                if (m_largest_value < abs(Ke.maxCoeff()))
                    m_largest_value = abs(Ke.maxCoeff());
            }
        }
    }
    
    return true; 
};

template < typename TFloe>
bool
FemProblem<TFloe>::addDirichlet(std::vector<size_t> gamma_d, std::vector<point_type> values)
{
    if (m_largest_value == 0)
    {
        std::cerr << "largest_value has not been initialized" << std::endl;
        return false; // this should have been set to a larger value in the 
    }
    size_t n_dirichlet(gamma_d.size());
    if (values.size() != n_dirichlet)
    {
        std::cerr << "Incoherent input" << std::endl;
        return false; 
    }

    real_type very_big_stuff=10000*m_largest_value;
    real_type theta(m_floe->get_frame().theta());
    m_FTriplet.clear();
    m_FTriplet.resize(n_dirichlet*2);
    m_DirichletTriplet.clear();
    m_DirichletTriplet.resize(n_dirichlet*2);

    for (size_t iDir = 0; iDir < n_dirichlet ; ++iDir)
    {
        // inverse rotation to transform the CL vector in the reference frame of the floe 
        m_FTriplet[iDir*2] = Eigen::Triplet<real_type>(gamma_d[iDir]*2, 0, very_big_stuff*(values[iDir].x*cos(theta) + values[iDir].y*sin(theta)));
        m_FTriplet[iDir*2+1] = Eigen::Triplet<real_type>(gamma_d[iDir]*2+1, 0, very_big_stuff*(values[iDir].x*sin(theta)*(-1) + values[iDir].y*cos(theta)));

        m_DirichletTriplet[iDir*2] = Eigen::Triplet<real_type>(gamma_d[iDir]*2, gamma_d[iDir]*2, very_big_stuff);
        m_DirichletTriplet[iDir*2+1] = Eigen::Triplet<real_type>(gamma_d[iDir]*2+1, gamma_d[iDir]*2+1, very_big_stuff);
    }
    // pour la solution analytique, on met une force sur le noeud en haut à gauche de la poutre 
    real_type P(300000);
    m_FTriplet.push_back(Eigen::Triplet<real_type>(1*2+1, 0, P));
    return true;
};



/**
 * @brief constructs a Dirichlet BC that 'looks like' an impact : 3 consecutive nodes with amplitudes {values/2, values, values/2} respectively.  
 * 
 * @tparam TFloe
 * @param std::vector<size_t> gamma_d : list of nodes on which there is an impact 
 * @param std::vector<point_type> values : 'amplitudes' of impact
 * @return bool : true if everything went fine 
 */
template < typename TFloe>
bool
FemProblem<TFloe>::addContactDirichlet(std::vector<size_t> gamma_d, std::vector<point_type> values)
{
    if (m_largest_value == 0)
    {
        std::cerr << "largest_value has not been initialized" << std::endl;
        return false; // this should have been set to a larger value in the constructor 
    }
    size_t n_dirichlet(gamma_d.size());
    if (values.size() != n_dirichlet)
    {
        std::cerr << "Incoherent input" << std::endl;
        return false; 
    }

    real_type very_big_stuff=10000*m_largest_value;
    real_type theta(m_floe->get_frame().theta());
    m_FTriplet.clear();
    m_DirichletTriplet.clear();

    connectivity_type connect; 
    connect = m_floe->mesh().connectivity();

    for (size_t iDir = 0; iDir < n_dirichlet ; ++iDir)
    {
        // looking for the two surrounding points around the impact point : 
        // building the list of all elements containing the current node 
        std::vector<size_t> surr_elements;
        std::vector<size_t> surr_nodes_temp;
        std::vector<size_t> surr_nodes;
        for (size_t iElem = 0; iElem < connect.size(); iElem++)
        {
            if ((connect[iElem][0] == gamma_d[iDir]) || (connect[iElem][1] == gamma_d[iDir]) || (connect[iElem][2] == gamma_d[iDir]))
                surr_elements.push_back(iElem);
        }
        // * if the list contains only one element : both other nodes are the surrounding ones 
        if (surr_elements.size() == 1)
        {
            if (connect[surr_elements[0]][0] == iDir)
            {
                surr_nodes.push_back(connect[surr_elements[0]][1]); 
                surr_nodes.push_back(connect[surr_elements[0]][2]);
            }
            else if (connect[surr_elements[0]][1] == iDir)
            {
                surr_nodes.push_back(connect[surr_elements[0]][0]); 
                surr_nodes.push_back(connect[surr_elements[0]][2]);
            }
            else
            {
                surr_nodes.push_back(connect[surr_elements[0]][0]); 
                surr_nodes.push_back(connect[surr_elements[0]][1]);
            }
        }
        // * otherwise, we build the list of all nodes contained in these elements. All nodes appearing twice in the list are removed. There should remain only two nodes. There you go. 
        else 
        {

            for (size_t iElem = 0; iElem < surr_elements.size(); iElem++)
            {
                surr_nodes_temp.push_back(connect[surr_elements[iElem]][0]);
                surr_nodes_temp.push_back(connect[surr_elements[iElem]][1]);
                surr_nodes_temp.push_back(connect[surr_elements[iElem]][2]);
            }
            
            for (size_t iNode = 0; iNode < surr_nodes_temp.size(); iNode++)
            {
                surr_nodes.push_back(surr_nodes_temp[iNode]);
                size_t count(0);
                for (size_t iNode2 = 0; iNode2 < surr_nodes_temp.size(); iNode2++)
                {
                    if (surr_nodes_temp[iNode2] == surr_nodes_temp[iNode])
                    {
                        count++;
                    }
                }
                if (count > 1)
                {
                    surr_nodes.pop_back(); 
                }
            }
        }
    
        // checking that the surrounding nodes are not already in the triplet m_FTriplet s
        for (size_t iDir2 = 0; iDir2 < n_dirichlet ; ++iDir2)
            if (iDir2 != iDir)
            {
                if (iDir == surr_nodes[0])
                    surr_nodes.erase(surr_nodes.begin()); 
                if (iDir == surr_nodes[1])
                    surr_nodes.pop_back(); 
            } 
        // pour l'instant : on y met une dirichlet avec values/2, à remplacer lorsque la condition dirichlet aura été déterminée 
        m_FTriplet.push_back(Eigen::Triplet<real_type>(gamma_d[iDir]*2, 0, very_big_stuff*(values[iDir].x*cos(theta) + values[iDir].y*sin(theta))));
        m_FTriplet.push_back(Eigen::Triplet<real_type>(gamma_d[iDir]*2+1, 0, very_big_stuff*(values[iDir].x*sin(theta)*(-1) + values[iDir].y*cos(theta))));
        m_DirichletTriplet.push_back(Eigen::Triplet<real_type>(gamma_d[iDir]*2, gamma_d[iDir]*2, very_big_stuff));
        m_DirichletTriplet.push_back(Eigen::Triplet<real_type>(gamma_d[iDir]*2+1, gamma_d[iDir]*2+1, very_big_stuff));

        for (size_t iPoint = 0 ; iPoint < surr_nodes.size() ; ++iPoint)
        {
            size_t point = surr_nodes[iPoint];
            m_FTriplet.push_back(Eigen::Triplet<real_type>(point*2, 0, 0.5*very_big_stuff*(values[iDir].x*cos(theta) + values[iDir].y*sin(theta))));
            m_FTriplet.push_back(Eigen::Triplet<real_type>(point*2+1, 0, 0.5*very_big_stuff*(values[iDir].x*sin(theta)*(-1) + values[iDir].y*cos(theta))));

            m_DirichletTriplet.push_back(Eigen::Triplet<real_type>(point*2, point*2, very_big_stuff));
            m_DirichletTriplet.push_back(Eigen::Triplet<real_type>(point*2+1, point*2+1, very_big_stuff));
        }
        // Idée d'accélération à faire lorsque la fonction de contact dirichlet aura été définie :
        // * enregistrer une fois pour toutes pour chaque noeud une liste des plus proches voisins
        // * combien de voisin pour un contact ? faut voir sur combien de points on l'étale. 
        // * en effet, la recherche de voisins doit aller vite car faite à chaque nouveau contact. 
    }
    return true;
};

template < typename TFloe>
bool
FemProblem<TFloe>::solve()
{
    if (!m_is_prepared)
        if (!prepare())
            return false; 
    m_stress_is_computed = false; 
    std::vector<Eigen::Triplet<real_type>> stiffness_and_bc;
    stiffness_and_bc = m_KTriplet;
    stiffness_and_bc.insert(stiffness_and_bc.end(), m_DirichletTriplet.begin(), m_DirichletTriplet.end() );
    // ATTENTION 
    // at this point, we should check that the max indices in triplets do not go beyond the matrix size. Otherwise ça crashe de toute façon. 
    for (size_t iPoint = 0; iPoint < stiffness_and_bc.size(); iPoint++)
    {
        if (stiffness_and_bc[iPoint].row() >= m_nDof*m_nN)
            std::cerr << "not supposed to find in K.row() " << stiffness_and_bc[iPoint].row()  << " " << stiffness_and_bc[iPoint].col() << " " << stiffness_and_bc[iPoint].value()<< std::endl;
        if (stiffness_and_bc[iPoint].col() >= m_nDof*m_nN)
            std::cerr << "not supposed to find in K.col() " << stiffness_and_bc[iPoint].row()  << " " << stiffness_and_bc[iPoint].col() << " " << stiffness_and_bc[iPoint].value()<< std::endl;
    }
    for (size_t iPoint = 0; iPoint < m_FTriplet.size(); iPoint++)
    {
        if (m_FTriplet[iPoint].row() >= m_nDof*m_nN)
            std::cerr << "not supposed to find in F.row() " << m_FTriplet[iPoint].row() << " " << m_FTriplet[iPoint].col() << " " << m_FTriplet[iPoint].value()<< std::endl;
        if (m_FTriplet[iPoint].col() >= 1)
            std::cerr << "not supposed to find in F.col() " << m_FTriplet[iPoint].row() << " " << m_FTriplet[iPoint].col() << " " << m_FTriplet[iPoint].value()<< std::endl;
    }
    m_K.setFromTriplets(stiffness_and_bc.begin(), stiffness_and_bc.end());
    m_F.setFromTriplets(m_FTriplet.begin(), m_FTriplet.end());

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<real_type>> solver;
    solver.compute(m_K);
    if(solver.info()!=Eigen::Success) 
    {// decomposition
        std::cerr << "LDLt Decomposition failed" << std::endl;
        return false;
    }
    m_Solution = solver.solve(m_F);
    if(solver.info()!=Eigen::Success) 
    {// descente - remontee
        std::cerr << "LDLt Solve failed" << std::endl;
        return false;
    }
    return true;
};


/**
 * @brief matrix of the shape function gradients. doc http://ressources.unit.eu/cours/Mecagora4/accueil.html, (calcul de poutre 2D en flexion, phase d'assemblage )
 * 
 * @tparam TFloe 
 * @param iElem 
 * @return Eigen::Matrix<double, 3, 6> 
 */
template < typename TFloe>
Eigen::Matrix<double, 3, 6>
FemProblem<TFloe>::computeB(size_t iElem) const
{
    Eigen::Matrix<double, 3, 6> B;
    
    multi_point_type coordinates;
    connectivity_type connect; 
    coordinates = m_floe->mesh().points();
    connect = m_floe->mesh().connectivity();


    double x0 = coordinates[connect[iElem][0]][0];
    double x1 = coordinates[connect[iElem][1]][0];
    double x2 = coordinates[connect[iElem][2]][0];
    double y0 = coordinates[connect[iElem][0]][1];
    double y1 = coordinates[connect[iElem][1]][1];
    double y2 = coordinates[connect[iElem][2]][1];
    double detJac = m_floe->mesh().get_jacobian(iElem);

    B << y1-y2, 0.0, y2-y0, 0.0, y0 - y1, 0.0,
        0.0, x2-x1, 0.0, x0-x2, 0.0, x1-x0,
        x2-x1, y1-y2, x0-x2, y2-y0, x1-x0, y0-y1;
    B *=1/(detJac);
    return B;
    // maybe store this function in triangle_mesh ? There is also a computeN there... 
};



/**
 * @brief elastic constants matrix, see for instance http://ressources.unit.eu/cours/Mecagora4/EXEMPLE/FICHES/F4/CTPAF4.htm 
 * 2D, contraintes planes 
 * 
 * @tparam TFloe 
 * @return Eigen::Matrix<double, 3, 3> 
 */
template < typename TFloe>
Eigen::Matrix<double, 3, 3>
FemProblem<TFloe>::computeH() const
{
    Eigen::Matrix<double, 3, 3> H;
    H << 1.0, m_nu, 0.0,
        m_nu, 1.0, 0.0,
        0.0, 0.0, (1.0-m_nu)/2;
    H *= m_E/(1.0-pow(m_nu,2));
    return H;
};
        

/**
 * @brief prepares the FEM computation : initializes all necessary variables and assembles the stiffness matrix and sollicitation vector 
 * 
 * @tparam TFloe 
 * @return bool, true if everything went fine 
 */
template < typename TFloe>
bool
FemProblem<TFloe>::prepare()
{
    if ( m_floe == nullptr ) 
    {
        std::cerr << "In FemProblem::prepare, m_floe has not been set to an initialized floe" << std::endl;
        return false; 
    }
    if (m_is_prepared)
        return true;
    m_nE = m_floe->mesh().get_n_cells();
    m_nN = m_floe->mesh().get_n_nodes();
    m_nDof = 2; // 2 displacements on each node. 
    m_nNodesPerElt = 3; // triangle elements 
    m_Solution = Eigen::SparseMatrix<real_type> (m_nN*m_nDof, 1);
    m_largest_value = 0;
    m_K = Eigen::SparseMatrix<real_type> (m_nN*m_nDof, m_nN*m_nDof);

    if (!assembleK())
    {
        std::cerr << "In FemProblem::prepare, could not assemble K" << std::endl;
        return false;
    }
    if (!assembleF())
    {
        std::cerr << "In FemProblem::prepare, could not assemble F" << std::endl;
        return false; 
    }
    m_is_prepared = true; 
    return true; 
};


/**
 * @brief performs all computation steps if required : setting new dirichlet condition, solve, compute stress, compute elementary energy, looks for fracture.   
 * 
 * @tparam TFloe 
 * @return bool, true if everything went fine 
 */
template < typename TFloe>
bool
FemProblem<TFloe>::performComputation(std::vector<size_t> gamma_d, std::vector<point_type> values)
{
    if (m_floe->total_received_impulse() == m_last_total_impulse)
    {
        // No need to solve, no evolution since last time
        return false; 
    }
    real_type amplitude = m_floe->total_received_impulse() - m_last_total_impulse;
    m_last_total_impulse = m_floe->total_received_impulse();
    for (size_t iDir = 0; iDir < values.size() ; ++iDir)
    {
        amplitude = 1; // pour l'instant, pour s'assurer que la condition est bien respectée. 
        point_type a(amplitude*values[iDir].x, amplitude*values[iDir].y);
        values[iDir] = a;
    }
    // if (!addDirichlet(gamma_d, values))
    // {
    //     std::cout << "Could not build dirichlet condition" << std::endl;
    //     std::cerr << "In FemProblem::performComputation, could not build dirichlet condition" << std::endl;
    //     return false; 
    // }

    if (!addContactDirichlet(gamma_d, values))
    {
        std::cout << "Could not build dirichlet condition" << std::endl;
        std::cerr << "In FemProblem::performComputation, could not build dirichlet condition" << std::endl;
        return false; 
    }

    if (!solve())
    {
        std::cout << "Could not solve " << std::endl;
        std::cerr << "In FemProblem::performComputation, could not solve " << std::endl;
        return false; 
    }

    m_stress_is_computed = compute_stress_vector();
    if (!m_stress_is_computed)
    {
        std::cout << "Could not compute stress " << std::endl;
        std::cerr << "In FemProblem::performComputation, could not compute stress " << std::endl;
        return false; 
    }

    m_energies_are_computed = compute_elementary_energies();
    if (!m_energies_are_computed)
    {
        std::cout << "Could not compute elastic energies " << std::endl;
        std::cerr << "In FemProblem::performComputation, could not compute elastic energies " << std::endl;
        return false; 
    }

    // recherche de fractures 
    std::vector<size_t> elems;
    for (size_t i = 0; i < m_nE; i++)
        elems.push_back(i);
    std::cout << "énergie totale : " << compute_elastic_energy(elems) << std::endl; 
    
    return true; 
};

}} // namespace floe::fem

#endif
