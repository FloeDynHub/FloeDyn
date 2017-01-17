/*!
 * \file floe/collision/matlab/mpi_detector.hpp
 * \brief Collision detector (MPI version)
 * \see MatlabDetector for more explanations.
 * \author Quentin Jouet
 */

#ifndef FLOE_COLLISION_MPI_DETECTOR_HPP
#define FLOE_COLLISION_MPI_DETECTOR_HPP

#include <mpi.h>
#include "floe/collision/matlab/detector.h"
#include <math.h>

namespace floe { namespace collision { namespace matlab
{


template <
    typename TFloeGroup,
    typename TProximityData = ProximityData<TFloeGroup, OptimizedFloe<typename TFloeGroup::floe_type>>,
    typename TContact = ContactPoint<typename TFloeGroup::floe_type> 
>
class MPIMatlabDetector : public MatlabDetector<TFloeGroup>
{

public:
    using base_class = MatlabDetector<TFloeGroup>;
    using floe_group_type = TFloeGroup;
    using value_type = typename base_class::value_type;
    using point_type = typename base_class::point_type;
    using optim_type = typename base_class::optim_type;
    using floe_interface_type = typename base_class::floe_interface_type;
    using optim_interface_type = typename base_class::optim_interface_type;
    using floe_distrib_type = std::map<int, std::vector<std::size_t>>;

    MPIMatlabDetector() : base_class() {}

    //! Deleted copy constructor
    MPIMatlabDetector( MPIMatlabDetector<TFloeGroup, TContact> const& ) = delete;

    //! Deleted copy operator
    MPIMatlabDetector<TFloeGroup, TContact>& operator= (MPIMatlabDetector const&) = delete;

    inline floe_distrib_type const& floe_process_distribution() const { return m_floe_process_distribution; };
    inline floe_distrib_type const& border_floe_process_distribution() const { return m_border_floe_process_distribution; };
    inline floe_distrib_type const& crossing_floe_process_distribution() const { return m_crossing_floe_process_distribution; };
    
    /* MASTER PART */

    void distribute_floes(){
        resize_window();
        // calc grid dimension
        int nb_process;
        MPI_Comm_size(MPI_COMM_WORLD, &nb_process);
        int square_grid_dim = 0;
        while (4 * (square_grid_dim + 1) * ((square_grid_dim + 1) - 1) + 2 <= nb_process) square_grid_dim++;
        m_dim_grid = {{square_grid_dim, square_grid_dim}};
        // Clear previous distribution
        m_floe_process_distribution.clear();
        m_border_floe_process_distribution.clear();
        m_crossing_floe_process_distribution.clear();
        // Ocean division
        for (std::size_t i = 0; i < this->data().get_optims().size(); ++i){
            auto const& optim = this->data().get_optim(i);
            distribute(optim, i);
        }
    }

    void display_floe_distrib(){
        std::cout << "Floe distrib " << "(" << this->grid_dim_x() << "x" << this->grid_dim_y() << ") : ";
        for (auto const& iter : m_floe_process_distribution)
            std::cout << "s" << iter.first << "(" << iter.second.size() << ") ";
        std::cout << "| Border : ";
        for (auto const& iter : m_border_floe_process_distribution)
            std::cout << "b" << iter.first << "(" << iter.second.size() << ") ";
        std::cout << " | Cross : ";
        for (auto const& iter : m_crossing_floe_process_distribution)
            std::cout << "c" << iter.first << "(" << iter.second.size() << ") ";
        std::cout << std::endl;
    }

    virtual void set_floe_group(floe_group_type const& floe_group) override {
        base_class::set_floe_group(floe_group);
        set_max_floe_radius();
    }

    void set_max_floe_radius(){
        value_type max_radius = 0.;
        for (auto const* optim : this->data().get_optims()){
            max_radius = std::max(max_radius, optim->global_disk().radius);
        }
        m_max_floe_radius = max_radius;
    }

    /* WORKER PART */

    void prepare_optims() override {
        for (std::size_t i=0; i< this->data().nb_floes(); ++i)
            this->get_optim(i).update();
    }

private:
    floe_distrib_type m_floe_process_distribution;
    floe_distrib_type m_border_floe_process_distribution;
    floe_distrib_type m_crossing_floe_process_distribution;
    std::array<value_type, 4> m_window;
    std::array<int, 2> m_dim_grid;
    value_type m_max_floe_radius;

    inline int grid_dim_x() const { return m_dim_grid[0]; }
    inline int grid_dim_y() const { return m_dim_grid[1]; }

    void resize_window(){
        // Update floe optimizers (but not their local disks)
        for ( auto optim_ptr : this->m_prox_data.get_optims() ) optim_ptr->update(false);
        // compute minimal rectangle ocean area
        value_type mg = 1e-8; // margin
        value_type min_x, min_y, max_x, max_y;
        min_x = min_y = std::numeric_limits<value_type>::max();
        max_x = max_y = - std::numeric_limits<value_type>::max();
        for (auto const* optim : this->data().get_optims()){
            max_x = std::max(optim->global_disk().center.x, max_x);
            min_x = std::min(optim->global_disk().center.x, min_x);
            max_y = std::max(optim->global_disk().center.y, max_y);
            min_y = std::min(optim->global_disk().center.y, min_y);
        }
        m_window = {{ min_x - mg, max_x + mg, min_y - mg, max_y + mg }};
    }

void distribute(optim_type const& optim, std::size_t floe_id){
        // Calc ocean 2D parcel dimension
        value_type parcel_width = (m_window[1] - m_window[0]) / this->grid_dim_x();
        value_type parcel_height = (m_window[3] - m_window[2]) / this->grid_dim_y();
        // Coordinates relative to grid
        auto X_grid = (optim.global_disk().center.x - m_window[0]) / parcel_width;     
        auto Y_grid = (optim.global_disk().center.y - m_window[2]) / parcel_height;
        // indices of floe parcel
        int X_id{(int)floor(X_grid)};
        int Y_id{(int)floor(Y_grid)};
        // indices of floe closest border
        int X_border_id = (int)round(X_grid);
        int Y_border_id = (int)round(Y_grid);
        // Is floe in some border region ?
        bool x_border_floe{
            X_border_id != 0 && X_border_id != this->grid_dim_x() && std::abs(X_grid - X_border_id) * parcel_width - optim.global_disk().radius < m_max_floe_radius
        };
        bool y_border_floe{
            Y_border_id != 0 && Y_border_id != this->grid_dim_y() && std::abs(Y_grid - Y_border_id) * parcel_height - optim.global_disk().radius < m_max_floe_radius
        };
        bool xy_border_floe{x_border_floe && y_border_floe};
        // assigning floe to related processes
        int NUM_PROC = 1; // Start worker counting at 1 (0 is the master)
        int worker_id = Y_id * this->grid_dim_x() + X_id + NUM_PROC ;
        m_floe_process_distribution[worker_id].push_back(floe_id); // order : left->right, down->up
        NUM_PROC += this->grid_dim_x() * this->grid_dim_y();
        if (x_border_floe){
            worker_id = Y_id * (this->grid_dim_x() - 1) + (X_border_id -1) + NUM_PROC;
            m_border_floe_process_distribution[worker_id].push_back(floe_id);
        }
        NUM_PROC += this->grid_dim_y() * (this->grid_dim_x() - 1);
        if (y_border_floe){
            worker_id = X_id * (this->grid_dim_y() - 1) + (Y_border_id - 1) + NUM_PROC;
            m_border_floe_process_distribution[worker_id].push_back(floe_id);
        }
        NUM_PROC += this->grid_dim_x() * (this->grid_dim_y() - 1);
        if (xy_border_floe){
            worker_id = (Y_border_id - 1) * (this->grid_dim_x() - 1) + (X_border_id -1) + NUM_PROC;
            m_crossing_floe_process_distribution[worker_id].push_back(floe_id);
        }
    }
};

}}} // namespace floe::collision::matlab

#endif // FLOE_COLLISION_MPI_DETECTOR_HPP
