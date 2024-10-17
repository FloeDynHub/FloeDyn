/*!
 * \file floe/collision/matlab/mpi_detector.hpp
 * \brief Collision detector (MPI Periodic version)
 * \see MatlabDetector for more explanations.
 * \author Quentin Jouet
 */

#ifndef FLOE_COLLISION_MPI_PERIODIC_DETECTOR_HPP
#define FLOE_COLLISION_MPI_PERIODIC_DETECTOR_HPP

#include <mpi.h>
#include "floe/collision/matlab/mpi_detector.hpp"
#include <math.h>

namespace floe { namespace collision { namespace matlab
{


template <
    typename TFloeGroup,
    typename TProximityData = ProximityData<TFloeGroup, OptimizedFloe<typename TFloeGroup::floe_type>>,
    typename TContact = ContactPoint<typename TFloeGroup::floe_type>,
    typename TDetector = MatlabDetector<TFloeGroup>
>
class MPIPeriodicDetector : public MPIMatlabDetector<TFloeGroup, TProximityData, TContact, TDetector>
{

public:
    using base_class = MPIMatlabDetector<TFloeGroup, TProximityData, TContact, TDetector>;
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

    MPIPeriodicDetector() : base_class() {
        std::cout << "MPIPeriodicDetector init" << std::endl;
        this->init_grid_dimension();
        this->m_process_partition.clear(); // TODO because base_class already inits it with no PBC
        this->partition_processes();
        for (auto key : this->collision_process_partition_keys()){
            // init keys
            this->m_process_partition[key];
        }
    }

    inline std::vector<std::string> collision_process_partition_keys() const {
        // list partition keys
        return {{ "grid", "x_border", "y_border", "crossing" }};
    };

    void distribute_floes(){
        this->resize_window();
        // Clear previous distribution
        this->m_floe_process_distrib.clear();
        // Ocean division
        for (std::size_t i = 0; i < this->data().get_optims().size(); ++i){
            auto const& optim = this->data().get_optim(i);
            this->distribute(optim, i);
        }
        std::cout << "END distribute_floes " << std::endl;
    }

private:

    void init_grid_dimension(){
        // calc grid dimension
        int nb_process;
        MPI_Comm_size(MPI_COMM_WORLD, &nb_process);
        int square_grid_dim = 0;
        while (4 * square_grid_dim * square_grid_dim + 1 < nb_process) square_grid_dim++;
        this->m_dim_grid = {{square_grid_dim, square_grid_dim}};
        this->m_nb_workers = 4 * square_grid_dim * square_grid_dim;
    }

    void partition_processes(){
        // base class
        int BASE_NUM_PROC = 1; // Start worker counting at 1 (0 is the master)
        int size_grid = this->grid_dim_x() * this->grid_dim_y();
        for (int i = 0; i < size_grid; i++){
            this->m_process_partition["grid"].push_back(BASE_NUM_PROC + i);
        }
        BASE_NUM_PROC += size_grid;
        int size_x_border = this->grid_dim_y() * (this->grid_dim_x() - 1);
        for (int i = 0; i < size_x_border; i++){
            this->m_process_partition["x_border"].push_back(BASE_NUM_PROC + i);
        }
        BASE_NUM_PROC += size_x_border;
        int size_y_border = this->grid_dim_x() * (this->grid_dim_y() - 1);
        for (int i = 0; i < size_y_border; i++){
            this->m_process_partition["y_border"].push_back(BASE_NUM_PROC + i);
        }
        BASE_NUM_PROC += size_x_border;
        int size_crossing = (this->grid_dim_x() - 1) * (this->grid_dim_y() - 1);
        for (int i = 0; i < size_crossing; i++){
            this->m_process_partition["crossing"].push_back(BASE_NUM_PROC + i);
        }
        // end base class
        // int BASE_NUM_PROC = 1; // Start worker counting at no-PBC version limit
        // Add periodic border processes
        BASE_NUM_PROC += size_crossing;
        int size_x_pbc_border = this->grid_dim_x();
        for (int i = 0; i < size_x_pbc_border; i++){
            this->m_process_partition["x_border"].push_back(BASE_NUM_PROC + i);
        }
        BASE_NUM_PROC += size_x_pbc_border;
        int size_y_pbc_border = this->grid_dim_y();
        for (int i = 0; i < size_y_pbc_border; i++){
            this->m_process_partition["y_border"].push_back(BASE_NUM_PROC + i);
        }
        BASE_NUM_PROC += size_y_pbc_border;
        int size_pbc_crossing = this->grid_dim_x() + this->grid_dim_y() -1;
        for (int i = 0; i < size_pbc_crossing; i++){
            this->m_process_partition["crossing"].push_back(BASE_NUM_PROC + i);
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
        this->m_window = {{ min_x - mg, max_x + mg, min_y - mg, max_y + mg }};
    }

void distribute(optim_type const& optim, std::size_t floe_id){
        // base class
        // Calc ocean 2D parcel dimension
        real_type parcel_width = (this->m_window[1] - this->m_window[0]) / this->grid_dim_x();
        real_type parcel_height = (this->m_window[3] - this->m_window[2]) / this->grid_dim_y();
        // Coordinates relative to grid
        auto X_grid = (optim.global_disk().center.x - this->m_window[0]) / parcel_width; // TODO /!\ global disk not updated !
        auto Y_grid = (optim.global_disk().center.y - this->m_window[2]) / parcel_height;
        // indices of floe parcel
        int X_id{(int)floor(X_grid)};
        int Y_id{(int)floor(Y_grid)};
        // indices of floe closest border
        int X_border_id = (int)round(X_grid);
        int Y_border_id = (int)round(Y_grid);
        // Is floe in some border region ?
        bool x_border_floe{
            X_border_id != 0 && X_border_id != this->grid_dim_x() && std::abs(X_grid - X_border_id) * parcel_width - optim.global_disk().radius < this->m_max_floe_radius
        };
        bool y_border_floe{
            Y_border_id != 0 && Y_border_id != this->grid_dim_y() && std::abs(Y_grid - Y_border_id) * parcel_height - optim.global_disk().radius < this->m_max_floe_radius
        };
        bool xy_border_floe{x_border_floe && y_border_floe};
        // assigning floe to related processes
        int worker_id = Y_id * this->grid_dim_x() + X_id;
        if (floe_id == 0) std::cout << "POS ZERO " << optim.global_disk().center.x << " " << optim.global_disk().center.y << " -> worker_id " << worker_id << std::endl;
        this->m_floe_process_distrib[this->m_process_partition["grid"][worker_id]].push_back(floe_id); // order : left->right, down->up
        if (x_border_floe){
            worker_id = Y_id * (this->grid_dim_x() - 1) + (X_border_id -1);
            this->m_floe_process_distrib[this->m_process_partition["x_border"][worker_id]].push_back(floe_id);
        }
        if (y_border_floe){
            worker_id = X_id * (this->grid_dim_y() - 1) + (Y_border_id - 1);
            this->m_floe_process_distrib[this->m_process_partition["y_border"][worker_id]].push_back(floe_id);
        }
        if (xy_border_floe){
            worker_id = (Y_border_id - 1) * (this->grid_dim_x() - 1) + (X_border_id -1);
            this->m_floe_process_distrib[this->m_process_partition["crossing"][worker_id]].push_back(floe_id);
        }
        // end base class
        int size_x_border = this->grid_dim_y() * (this->grid_dim_x() - 1);
        int size_y_border = this->grid_dim_x() * (this->grid_dim_y() - 1);
        int size_crossing = (this->grid_dim_x() - 1) * (this->grid_dim_y() - 1);
        bool X_pbc_border_floe {
          (X_border_id == 0 || X_border_id == this->grid_dim_x()) && std::abs(X_grid - X_border_id) * parcel_width - optim.global_disk().radius < this->m_max_floe_radius
        };
        if (X_pbc_border_floe) {
          worker_id = size_x_border + Y_id; // TODO VERIFIER (formule 1er jet)
        //   std::cout << "X_pbc_border_floe " << worker_id << " " << size_x_border << " " << Y_id << std::endl;
        //   if (worker_id > 3) std::cout << "PB distrib X_pbc_border_floe" << worker_id << std::endl;
          this->m_floe_process_distrib[this->m_process_partition["x_border"][worker_id]].push_back(floe_id);
        }
        bool Y_pbc_border_floe {
          (Y_border_id == 0 || Y_border_id == this->grid_dim_y()) && std::abs(Y_grid - Y_border_id) * parcel_height - optim.global_disk().radius < this->m_max_floe_radius
        };
        if (Y_pbc_border_floe) {
          worker_id = size_y_border + X_id; // TODO VERIFIER (formule 1er jet)
        //   std::cout << "Y_pbc_border_floe " << worker_id << std::endl;
        //   if (worker_id > 3) std::cout << "PB distrib Y_pbc_border_floe" << worker_id << std::endl;
          this->m_floe_process_distrib[this->m_process_partition["y_border"][worker_id]].push_back(floe_id);
        }
        bool XY_pbc_border_floe = X_pbc_border_floe && Y_pbc_border_floe;
        if (XY_pbc_border_floe) {
          worker_id = size_crossing + X_id + Y_id; // TODO VERIFIER (formule 1er jet)
        //   std::cout << "XY_pbc_border_floe " << worker_id << std::endl;
        //   if (worker_id > 3) std::cout << "PB distrib XY_pbc_border_floe" << worker_id << std::endl;
          this->m_floe_process_distrib[this->m_process_partition["crossing"][worker_id]].push_back(floe_id);
        }
    }
};

}}} // namespace floe::collision::matlab

#endif // FLOE_COLLISION_MPI_PERIODIC_DETECTOR_HPP
