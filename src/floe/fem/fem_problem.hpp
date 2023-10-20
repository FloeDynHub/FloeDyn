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
        m_E{9300000000.0}, // dixit wikipedia 
        m_nu{0.4}, // no idea what it should be for ice, I took the one from the glass comme ça, à pouf. 
        m_nE{m_floe->mesh().get_n_cells()},
        m_nN{m_floe->mesh().get_n_nodes()},
        m_nDof{2}, // 2 displacements at each node 
        m_nNodesPerElt{3}, // triangle elements only... 
        m_Sol{Eigen::SparseMatrix<real_type> (1, 1)},
        m_is_prepared{false}, 
        m_largest_value{0},
        m_last_total_impulse{0}
    {}
    
    
    bool prepare();
    bool assembleK();
    bool addDirichlet(std::vector<size_t> gamma_d, std::vector<point_type> values);
    bool assembleF()
    {
        m_F = Eigen::SparseMatrix<real_type> (m_nN*m_nDof, 1); // F is initialized at zero 
        return true; 
    };
    real_type computeKE()
    {
        return 0.0; 
    };
    bool solve();
    bool performComputation(std::vector<size_t> gamma_d, std::vector<point_type> values);
    inline Eigen::SparseMatrix<real_type> get_solution() const {return m_Sol;};
    std::vector<real_type> get_solution_vector() const 
    {
        std::vector<real_type> v;
        v.resize(m_nN);
        if (m_Sol.rows() < 2*m_nN)
        {
            std::cout << "problem with matrix size" << std::endl;
            WHEREAMI
        }
        for (size_t iPoint = 0; iPoint < m_nN ; ++iPoint)
            v[iPoint] = m_Sol.coeff(2*iPoint,0);
        return v;
    };

    // std::vector<std::vector<real_type>> get_stress_vector() const 
    std::vector<real_type> get_stress_vector() const 
    // nouvelle tentative : on met bout à bout les trois composantes, donc le vecteur est de taille 3*nE par 1. 
    {
        // returns a vector of stress components. 
        // std::vector<std::vector<real_type>> v;
        std::vector<real_type> v;
        if (m_Sol.rows() != 2*m_nN)
            return v;
        v.resize(m_nE*3);
        Eigen::Matrix<real_type, 3,1> sigma;
        Eigen::Matrix<real_type, 6,1> sol_elem;
        connectivity_type connect; 
        connect = m_floe->mesh().connectivity();
        Eigen::Matrix<real_type, 3, 3> H = computeH(); 
        Eigen::Matrix<real_type, 3, 6> B;
        for (size_t iElem = 0 ; iElem < m_nE ; ++iElem)
        {
            WHEREAMI
            // v[iElem].resize(3); // bad malloc, chais pas pourquoi... replacing this by calls to push_back
            B = computeB(iElem); 
            // WHEREAMI
            // v[iElem].clear();
            sol_elem(0,0) = m_Sol.coeff(connect[iElem][0]*2, 0);
            sol_elem(1,0) = m_Sol.coeff(connect[iElem][0]*2+1, 0);
            sol_elem(2,0) = m_Sol.coeff(connect[iElem][1]*2, 0);
            sol_elem(3,0) = m_Sol.coeff(connect[iElem][1]*2+1, 0);
            sol_elem(4,0) = m_Sol.coeff(connect[iElem][2]*2, 0);
            sol_elem(5,0) = m_Sol.coeff(connect[iElem][2]*2+1, 0);
            sigma = H*B*sol_elem;
            // WHEREAMI
            // v[iElem].push_back(sigma(0,0));
            // WHEREAMI
            // v[iElem].push_back(sigma(1,0));
            // WHEREAMI
            // std::cout << "là j'essaie d'écrire " << sigma(2,0) << " à l'emplacement " << iElem << " de taille " << v[iElem].size()<< std::endl;
            // WHEREAMI
            // v[iElem].push_back(sigma(2,0));
            // WHEREAMI
            v[iElem] = sigma(0,0);
            // WHEREAMI
            v[iElem+m_nE] = sigma(1,0);
            // WHEREAMI 
            v[iElem+2*m_nE] = sigma(2,0);
            // WHEREAMI
            // v[iElem][0] = sigma(0,0);
            // WHEREAMI
            // v[iElem][1] = sigma(1,0);
            // WHEREAMI 
            // v[iElem][2] = sigma(2,0);
            // std::cout << "At elem " << iElem << " sigma = " << sigma << std::endl;
        }    
        return v;
    };
    inline bool unset_prepared() {m_is_prepared = false;};
    Eigen::Matrix<double, 3, 6> computeB(size_t iElem) const;
    Eigen::Matrix<double, 3, 3> computeH() const;

private : 
    floe_type * m_floe;
    real_type m_E; // dixit wikipedia 
    real_type m_nu; // no idea what it should be for ice, I took the one from the glass comme ça, à pouf. 
    size_t m_nE;
    size_t m_nN;
    size_t m_nDof; // number of degrees of freedom. 2 translations for instance 
    size_t m_nNodesPerElt; // number of nodes per element. A priori only triangles will be used but on ne sait jamais. 
    std::vector<Eigen::Triplet<real_type>> m_KTriplet;
    std::vector<Eigen::Triplet<real_type>> m_DirichletTriplet;
    std::vector<Eigen::Triplet<real_type>> m_FTriplet;
    Eigen::SparseMatrix<real_type> m_K;
    Eigen::SparseMatrix<real_type> m_F;
    Eigen::SparseMatrix<real_type> m_Sol;
    bool m_is_prepared;
    real_type m_largest_value; 
    real_type m_last_total_impulse; 

};


template < typename TFloe>
bool
FemProblem<TFloe>::assembleK()
// FemProblem::assembleK()
{
    // initialization 
    // std::vector<Eigen::Triplet<real_type>> m_KTriplet;
    

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
    // loop over the elements, construction of elementary matrices 
    for (size_t iElem = 0 ; iElem < m_nE ; ++iElem)
    { 
        // real_type x0 = coordinates[connect[iElem][0]][0];
        // real_type x1 = coordinates[connect[iElem][1]][0];
        // real_type x2 = coordinates[connect[iElem][2]][0];
        // real_type y0 = coordinates[connect[iElem][0]][1];
        // real_type y1 = coordinates[connect[iElem][1]][1];
        // real_type y2 = coordinates[connect[iElem][2]][1];

        // std::cout << "element " << iElem << " : " << connect[iElem][0] << " " << connect[iElem][1] << " " << connect[iElem][2] << std::endl;
        size_t gotcha(3);
        if (connect[iElem][0] == gotcha || connect[iElem][1] == gotcha || connect[iElem][2] == gotcha )
            std::cout << "Point " << gotcha << " is found in element " << iElem << std::endl; 
        detJac = m_floe->mesh().get_jacobian(iElem);
        if (detJac == 0) {std::cerr << "in FemProblem::prepare, could not get jacobian of element" << iElem << std::endl; }

        

        B = computeB(iElem); 
        H = computeH(); 

        // std::cout << "B : " << B << std::endl << "H : " << H << std::endl; 
        Ke = Eigen::MatrixXd::Constant(m_nNodesPerElt*m_nDof, m_nNodesPerElt*m_nDof , 0);
        for (int iG = 0; iG < ngp; iG++) // for linear elements, B is constant, no need to integrate over several points you stupid son of a beach 
        {
            // std::cout << " detJac : " << detJac << "; weight : " << gp.coeff(iG, gw) << "; Bt H B : " << B.transpose()*H*B << std::endl; 
            Ke = Ke + detJac*gp.coeff(iG, gw)*B.transpose()*H*B; 
        }
        // std::cout << "Ke : " << Ke << std::endl << std::endl; 
        for (unsigned int i = 0 ; i < m_nNodesPerElt ; i++)
        {
            size_t globalI = m_nDof*connect[iElem][i];
            for (unsigned int j = 0; j < m_nNodesPerElt ; j++ )
            {
                size_t globalJ = m_nDof*connect[iElem][j];
                m_KTriplet[k++] = Eigen::Triplet<real_type>(globalI, globalJ, Ke.coeff(i,j)); 
                // std::cout << "adding " << Ke.coeff(i,j) << " at location " << globalI << " by " << globalJ << std::endl;
                m_KTriplet[k++] = Eigen::Triplet<real_type>(globalI+1, globalJ, Ke.coeff(i+1,j)); 
                m_KTriplet[k++] = Eigen::Triplet<real_type>(globalI+1, globalJ+1, Ke.coeff(i+1,j+1)); 
                m_KTriplet[k++] = Eigen::Triplet<real_type>(globalI, globalJ+1, Ke.coeff(i,j+1)); 
                if((isnan((Ke)(i,j))) || (isnan((Ke)(i+1,j))) || (isnan((Ke)(i+1,j+1))) || (isnan((Ke)(i+1,j+1))))
                {
                    std::cout << "Found a nan in Ke, element " << iElem << ". Stopping right now. " <<std:: endl << Ke << std::endl;
                    return false; 	
                }
                // preparing Dirichlet condition 
                if (m_largest_value < abs(Ke.coeff(i,j)))
                    m_largest_value = abs(Ke.coeff(i,j));
                if (m_largest_value < abs(Ke.coeff(i+1,j)))
                    m_largest_value = abs(Ke.coeff(i+1,j));
                if (m_largest_value < abs(Ke.coeff(i+1,j+1)))
                    m_largest_value = abs(Ke.coeff(i+1,j+1));
                if (m_largest_value < abs(Ke.coeff(i,j+1)))
                    m_largest_value = abs(Ke.coeff(i,j+1));
            }
        }
    }
    // std::cout << "during assembly, m_largest_value = " << m_largest_value << std::endl;
    // for (size_t iPoint = 0; iPoint < m_KTriplet.size() ; ++iPoint)
    //     if (abs(m_KTriplet[iPoint].value()) > m_largest_value)
    //         m_largest_value = m_KTriplet[iPoint].value();
    //     // else 
    //     //     std::cout << abs(m_KTriplet[iPoint].value()) << "is smaller than " << m_largest_value << std::endl;
    // std::cout << "another attempt, m_largest_value = " << m_largest_value << std::endl;
    
    // size_t maxi(0); 
    // size_t maxj(0); 
    // size_t mini(0); 
    // size_t minj(0); 
    // for (size_t iNode = 0; iNode < m_KTriplet.size() ; ++iNode)
    // {
    //     if (m_KTriplet[iNode].row() > maxi)
    //         maxi = m_KTriplet[iNode].row();
    //     if (m_KTriplet[iNode].col() > maxj)
    //         maxj = m_KTriplet[iNode].col();
    //     if (m_KTriplet[iNode].row() < mini)
    //         mini = m_KTriplet[iNode].row();
    //     if (m_KTriplet[iNode].col() < minj)
    //         minj = m_KTriplet[iNode].col();
    // }
    // std::cout << "maxi " << maxi << "; maxj " << maxj << "mini " << mini << "; minj " << minj << std::endl;
    // std::cout << "size of KTriplet : " << m_KTriplet.size() << ". m_nE*m_nDof*m_nNodesPerElt*m_nDof*m_nNodesPerElt : " << m_nE*m_nDof*m_nNodesPerElt*m_nDof*m_nNodesPerElt<< std::endl;
    // std::cout << m_KTriplet.end()- m_KTriplet.begin() << std::endl; 
    
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

    m_FTriplet.clear();
    m_FTriplet.resize(n_dirichlet*2);
    m_DirichletTriplet.clear();
    m_DirichletTriplet.resize(n_dirichlet*2);

    for (size_t iDir = 0; iDir < n_dirichlet ; ++iDir)
    {
        m_FTriplet[iDir*2] = Eigen::Triplet<real_type>(gamma_d[iDir]*2, 0, very_big_stuff*values[iDir].x);
        m_FTriplet[iDir*2+1] = Eigen::Triplet<real_type>(gamma_d[iDir]*2+1, 0, very_big_stuff*values[iDir].y);
        // std::cout << "adding " << very_big_stuff*values[iDir] << " at location " << gamma_d[iDir]*2 << " * " << 0 << std::endl;
        m_DirichletTriplet[iDir*2] = Eigen::Triplet<real_type>(gamma_d[iDir]*2, gamma_d[iDir]*2, very_big_stuff);
        m_DirichletTriplet[iDir*2+1] = Eigen::Triplet<real_type>(gamma_d[iDir]*2+1, gamma_d[iDir]*2+1, very_big_stuff);
        // std::cout << "adding " << very_big_stuff << " at location " << gamma_d[iDir]*2 << " * " << gamma_d[iDir]*2 << std::endl;
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
    std::vector<Eigen::Triplet<real_type>> stiffness_and_bc;
    stiffness_and_bc = m_KTriplet;
    stiffness_and_bc.insert(stiffness_and_bc.end(), m_DirichletTriplet.begin(), m_DirichletTriplet.end() );
    m_K.setFromTriplets(stiffness_and_bc.begin(), stiffness_and_bc.end());
    m_F.setFromTriplets(m_FTriplet.begin(), m_FTriplet.end());

    // std::cout << "size of stiffness and bc : " << stiffness_and_bc.size() << std::endl;

    // size_t N(80);
    // std::cout << "----------- KTriplets ------------ " << std::endl;
    // for (size_t i = 0; i < N ; ++i)
    //     std::cout << m_KTriplet[i].value() << " at location " << m_KTriplet[i].row() << ", " << m_KTriplet[i].col() << std::endl; 
    // std::cout << std::endl;

    // std::cout << "----------- stiffness_and_bc ------------ " << std::endl;
    // // for (size_t i = 0; i < N ; ++i)
    // //     std::cout << stiffness_and_bc[i].value() << " at location " << stiffness_and_bc[i].row() << ", " << stiffness_and_bc[i].col() << std::endl; 
    // for (size_t i = 0; i < stiffness_and_bc.size() ; ++i)
    //     std::cout << stiffness_and_bc[i].row() << ", " << stiffness_and_bc[i].col() << " : " << stiffness_and_bc[i].value() <<  std::endl; 
    // std::cout << std::endl;

    // std::cout << "----------- K ------------ " << std::endl;
    // for (size_t i = 0; i < N ; ++i)
    // {
    //     for (size_t j = 0; j < N ; ++j)
    //         std::cout << std::setw(3) << std::setprecision(3) << m_K.coeff(i,j) << " ";
    //     std::cout << std::endl;    
    // }

    // std::cout << "----------- F ------------ " << std::endl;
    // std::cout << std::endl;    
    // for (size_t i = 0; i < N ; ++i)
    //     std::cout << m_F.coeff(N,1) << " "; 

    // std::cout << "m_K contains " << m_K.nonZeros() << " nonzeros, and F " << m_F.nonZeros() << " ones." << std::endl;

    // on devrait pouvoir utiliser LLT, c'est pas complexe et c'est symmétrique 
    std::cout << "Performing FEM solve" << std::endl; 
    Eigen::SparseLU<Eigen::SparseMatrix<real_type>, Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern(m_K);
    solver.factorize(m_K);
    m_Sol = solver.solve(m_F); 
    if(solver.info()!=Eigen::Success) 
    {
        return false;
    }	
    // std::cout << "-> solve successful" << std::endl << m_Sol << std::endl; 
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
};

/**
 * @brief elastic constants matrix, see for instance http://ressources.unit.eu/cours/Mecagora4/accueil.html, (calcul de poutre 2D en flexion, phase d'assemblage)
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
        


template < typename TFloe>
bool
FemProblem<TFloe>::prepare()
// FemProblem::prepare()
{
    // WHEREAMI
    if ( m_floe == nullptr ) 
    {
        std::cerr << "In FemProblem::prepare, m_floe has not been set to an initialized floe" << std::endl;
        return false; 
    }
    // WHEREAMI
    if (m_is_prepared)
        return true;
    // WHEREAMI
    m_nE = m_floe->mesh().get_n_cells();
    m_nN = m_floe->mesh().get_n_nodes();
    m_nDof = 2; // 2 displacements on each node. 
    m_nNodesPerElt = 3; // triangle elements 
    m_Sol = Eigen::SparseMatrix<real_type> (m_nN*m_nDof, 1);
    m_largest_value = 0;
    m_K = Eigen::SparseMatrix<real_type> (m_nN*m_nDof, m_nN*m_nDof);
    // std::cout << "nE : " << m_nE <<  "; nN : " << m_nN <<  "; nDofs : " << m_nDof <<  "; nNpe : " << m_nNodesPerElt << std::endl; 

    if (!assembleK())
    {
        std::cerr << "In FemProblem::prepare, could not assemble K" << std::endl;
        return false;
    }
    std::cout << "-> Assembly successful" << std::endl;
    if (!assembleF())
    {
        std::cerr << "In FemProblem::prepare, could not assemble F" << std::endl;
        return false; 
    }
    std::cout << "-> RHS successful" << std::endl;
    m_is_prepared = true; 
    // WHEREAMI
    return true; 
};


template < typename TFloe>
bool
FemProblem<TFloe>::performComputation(std::vector<size_t> gamma_d, std::vector<point_type> values)
{
    // WHEREAMI
    // std::vector<size_t> merdouilles = gamma_d;
    // std::vector<point_type> merdouillesValues = values;

    if (m_floe->total_received_impulse() == m_last_total_impulse)
    {
        std::cout << "No need to solve, nothing new since last time" << std::endl; 
        return false; 
    }
    real_type amplitude = m_floe->total_received_impulse() - m_last_total_impulse;
    m_last_total_impulse = m_floe->total_received_impulse();
    std::cout << "Performing computation with amplitude = " << amplitude << std::endl;

    for (size_t iDir = 0; iDir < values.size() ; ++iDir)
    {
        amplitude = 1; // pour l'instant, pour s'assurer que la condition est bien respectée. 
        point_type a(amplitude*values[iDir].x, amplitude*values[iDir].y);
        values[iDir] = a;
    }

    if (!addDirichlet(gamma_d, values))
    {
        std::cerr << "In FemProblem::performComputation, could not build dirichlet condition" << std::endl;
        return false; 
    }
    std::cout << "-> Dirichlet successful" << std::endl;
    if (!solve())
    {
        std::cerr << "In FemProblem::performComputation, could not solve " << std::endl;
        return false; 
    }
    std::cout << "-> Solve successful" << std::endl;
    // WHEREAMI
    return true; 
};


}} // namespace floe::fem

#endif
