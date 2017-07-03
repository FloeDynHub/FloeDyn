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
    using real_type = typename base_class::real_type;
    using point_type = typename base_class::point_type;
    using optim_type = typename base_class::optim_type;
    using floe_interface_type = typename base_class::floe_interface_type;
    using optim_interface_type = typename base_class::optim_interface_type;
    using floe_distrib_type = std::map<int, std::vector<std::size_t>>;
    using process_list_type = std::vector<int>;
    using process_partition_type = std::map<std::string, process_list_type>;
    using floe_multi_partition_type = std::map<std::string, floe_distrib_type>;

    MPIMatlabDetector() : base_class() {
        this->init_grid_dimension();
        this->partition_processes();
        for (auto key : this->collision_process_partition_keys()){
            // init keys
            this->m_process_partition[key];
        }
    }

    //! Deleted copy constructor
    MPIMatlabDetector( MPIMatlabDetector<TFloeGroup, TContact> const& ) = delete;

    //! Deleted copy operator
    MPIMatlabDetector<TFloeGroup, TContact>& operator= (MPIMatlabDetector const&) = delete;

    inline floe_distrib_type const& floe_process_distribution() const {
        // return m_floe_process_distribution;
        // return m_floe_process_distributions["grid"];
        return m_floe_process_distrib;
    };
    inline process_partition_type const& process_partition() const {
        return m_process_partition;
    };
    inline process_list_type all_worker_processes() const {
        process_list_type resp(m_nb_workers);
        std::iota(std::begin(resp), std::end(resp), 1);
        return resp;
    };
    inline process_list_type const& base_grid_processes() const {
        return m_process_partition.at("grid");
    };
    inline std::vector<std::string> collision_process_partition_keys() const {
        // list partition keys
        return {{ "grid", "x_border", "y_border", "crossing" }};
    };
    inline process_list_type borders_processes() const {
        process_list_type resp;
        for (auto& key : { "x_border", "y_border", "crossing" }){
            resp.insert(
                resp.end(),
                this->process_partition().at(key).begin(),
                this->process_partition().at(key).end()
            );
        }
        return resp;
    }
    // inline floe_distrib_type const& border_floe_process_distribution() const {
    //     // return m_border_floe_process_distribution;
    // };
    // inline floe_distrib_type const& x_border_floe_process_distribution() const {
    //     // return m_border_floe_process_distribution;
    //     return m_floe_process_distributions["x_border"];
    // };
    // inline floe_distrib_type const& y_border_floe_process_distribution() const {
    //     // return m_border_floe_process_distribution;
    //     return m_floe_process_distributions["x_border"];
    // };
    // inline floe_distrib_type const& crossing_floe_process_distribution() const {
    //     return m_crossing_floe_process_distribution;
    // };
    
    /* MASTER PART */

    void distribute_floes(){
        this->resize_window();
        // Clear previous distribution
        // m_floe_process_distribution.clear();
        // m_border_floe_process_distribution.clear();
        // m_crossing_floe_process_distribution.clear();
        m_floe_process_distrib.clear();
        // Ocean division
        for (std::size_t i = 0; i < this->data().get_optims().size(); ++i){
            auto const& optim = this->data().get_optim(i);
            this->distribute(optim, i);
        }
    }

    void display_floe_distrib(){
        std::cout << "Floe distrib " << "(" << this->grid_dim_x() << "x" << this->grid_dim_y() << ") : ";
        for (auto const& p_id : m_process_partition["grid"])
        {
            std::cout << "s" << p_id << "(" << this->m_floe_process_distrib[p_id].size() << ") ";
        }
        std::cout << "| Border : ";
        for (auto const& p_id : m_process_partition["x_border"])
        {
            std::cout << "b" << p_id << "(" << this->m_floe_process_distrib[p_id].size() << ") ";
        }
        for (auto const& p_id : m_process_partition["y_border"])
        {
            std::cout << "b" << p_id << "(" << this->m_floe_process_distrib[p_id].size() << ") ";
        }
        std::cout << " | Cross : ";
        for (auto const& p_id : m_process_partition["crossing"])
        {
            std::cout << "c" << p_id << "(" << this->m_floe_process_distrib[p_id].size() << ") ";
        }
        std::cout << std::endl;
    }

    virtual void set_floe_group(floe_group_type const& floe_group) override {
        base_class::set_floe_group(floe_group);
        this->set_max_floe_radius();
    }

    void set_max_floe_radius(){
        real_type max_radius = 0.;
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
    // floe_distrib_type m_floe_process_distribution;
    // floe_distrib_type m_border_floe_process_distribution;
    // floe_distrib_type m_crossing_floe_process_distribution;
    floe_distrib_type m_floe_process_distrib;
    process_partition_type m_process_partition;
    // floe_multi_partition_type m_floe_process_distributions;
    std::array<real_type, 4> m_window;
    std::array<int, 2> m_dim_grid;
    int m_nb_workers;
    real_type m_max_floe_radius;

    inline int grid_dim_x() const { return m_dim_grid[0]; }
    inline int grid_dim_y() const { return m_dim_grid[1]; }

    void init_grid_dimension(){
        // calc grid dimension
        int nb_process;
        MPI_Comm_size(MPI_COMM_WORLD, &nb_process);
        int square_grid_dim = 0;
        while (4 * (square_grid_dim + 1) * ((square_grid_dim + 1) - 1) + 2 <= nb_process) square_grid_dim++;
        m_dim_grid = {{square_grid_dim, square_grid_dim}};
        m_nb_workers = 4 * square_grid_dim * (square_grid_dim - 1) + 1;
    }

    void partition_processes(){
        int BASE_NUM_PROC = 1; // Start worker counting at 1 (0 is the master)
        int size_grid = this->grid_dim_x() * this->grid_dim_y();
        for (int i = 0; i < size_grid; i++){
            m_process_partition["grid"].push_back(BASE_NUM_PROC + i);
        }
        BASE_NUM_PROC += size_grid;
        int size_x_border = this->grid_dim_y() * (this->grid_dim_x() - 1);
        for (int i = 0; i < size_x_border; i++){
            m_process_partition["x_border"].push_back(BASE_NUM_PROC + i);
        }
        BASE_NUM_PROC += size_x_border;
        int size_y_border = this->grid_dim_x() * (this->grid_dim_y() - 1);
        for (int i = 0; i < size_y_border; i++){
            m_process_partition["y_border"].push_back(BASE_NUM_PROC + i);
        }
        BASE_NUM_PROC += size_x_border;
        int size_crossing = (this->grid_dim_x() - 1) * (this->grid_dim_y() - 1);
        for (int i = 0; i < size_crossing; i++){
            m_process_partition["crossing"].push_back(BASE_NUM_PROC + i);
        }
    }

    void resize_window(){
        // Update floe optimizers (but not their local disks)
        for ( auto optim_ptr : this->m_prox_data.get_optims() ) optim_ptr->update(false);
        // compute minimal rectangle ocean area
        real_type mg = 1e-8; // margin
        real_type min_x, min_y, max_x, max_y;
        min_x = min_y = std::numeric_limits<real_type>::max();
        max_x = max_y = - std::numeric_limits<real_type>::max();
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
        real_type parcel_width = (m_window[1] - m_window[0]) / this->grid_dim_x();
        real_type parcel_height = (m_window[3] - m_window[2]) / this->grid_dim_y();
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
        // int NUM_PROC = 1; // Start worker counting at 1 (0 is the master)
        // int worker_id = Y_id * this->grid_dim_x() + X_id + NUM_PROC ;
        int worker_id = Y_id * this->grid_dim_x() + X_id;
        m_floe_process_distrib[m_process_partition["grid"][worker_id]].push_back(floe_id); // order : left->right, down->up
        // m_floe_process_distributions["grid"][worker_id].push_back(floe_id);
        // m_floe_process_distribution[worker_id].push_back(floe_id); // order : left->right, down->up
        // NUM_PROC += this->grid_dim_x() * this->grid_dim_y();
        if (x_border_floe){
            // worker_id = Y_id * (this->grid_dim_x() - 1) + (X_border_id -1) + NUM_PROC;
            worker_id = Y_id * (this->grid_dim_x() - 1) + (X_border_id -1);
            m_floe_process_distrib[m_process_partition["x_border"][worker_id]].push_back(floe_id);
            // m_border_floe_process_distribution[worker_id].push_back(floe_id);
            // m_floe_process_distributions["x_border"][worker_id].push_back(floe_id);
        }
        // NUM_PROC += this->grid_dim_y() * (this->grid_dim_x() - 1);
        if (y_border_floe){
            // worker_id = X_id * (this->grid_dim_y() - 1) + (Y_border_id - 1) + NUM_PROC;
            worker_id = X_id * (this->grid_dim_y() - 1) + (Y_border_id - 1);
            m_floe_process_distrib[m_process_partition["y_border"][worker_id]].push_back(floe_id);
            // m_border_floe_process_distribution[worker_id].push_back(floe_id);
            // m_floe_process_distributions["y_border"][worker_id].push_back(floe_id);
        }
        // NUM_PROC += this->grid_dim_x() * (this->grid_dim_y() - 1);
        if (xy_border_floe){
            // worker_id = (Y_border_id - 1) * (this->grid_dim_x() - 1) + (X_border_id -1) + NUM_PROC;
            worker_id = (Y_border_id - 1) * (this->grid_dim_x() - 1) + (X_border_id -1);
            m_floe_process_distrib[m_process_partition["crossing"][worker_id]].push_back(floe_id);
            // m_crossing_floe_process_distribution[worker_id].push_back(floe_id);
            // m_floe_process_distributions["crossing"][worker_id].push_back(floe_id);
        }
    }
};

}}} // namespace floe::collision::matlab

#endif // FLOE_COLLISION_MPI_DETECTOR_HPP
