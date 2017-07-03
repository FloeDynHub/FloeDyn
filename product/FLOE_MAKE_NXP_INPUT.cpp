#include <iostream>
#include <cassert>
#include "../product/config/config.hpp"

/*
Create an hdf5 unit input file (domain = 1m * 1m) from an input file
by loading it N * P times and translating floes.
*/


int main( int argc, char* argv[] )
{   
    using namespace std;
    using namespace types;

    if ( argc < 4 )
    {
        cout << "Usage: " << argv[0] << " <input_file_name> <nb_rows> <nb_cols>" << endl;
        return 1;
    }

    std::string input_list_floe_filename = argv[1];
    int nb_rows = atoi(argv[2]);
    int nb_cols = atoi(argv[3]);

    problem_type P;
    auto& floes = P.get_floe_group().get_floes();
    P.load_config(input_list_floe_filename);
    int nb_floes = floes.size();
    auto w = P.get_floe_group().get_initial_window();
    auto win_width = (w[1] - w[0]);
    auto win_height = (w[3] - w[2]);
    for (int k=0; k<nb_cols*nb_rows-1; k++){
        P.load_config(input_list_floe_filename);
    }
    decltype(w) new_window{
        w[0] - ((nb_cols - 1) * win_width / 2),
        w[1] + ((nb_cols - 1) * win_width / 2),
        w[2] - ((nb_rows - 1) * win_height / 2),
        w[3] + ((nb_rows - 1) * win_height / 2)
    };
    int id_floe = 0;
    for (int i=0; i<nb_rows; i++){
        for (int j=0; j<nb_cols; j++){
            point_type trans{
                new_window[0]+ (2*j + 1) * win_width / 2,
                new_window[2]+ (2*i + 1) * win_height / 2
            };
            for (int l=0; l < nb_floes; l++){
                auto& floe = floes[id_floe];
                floe.state().pos += trans;
                id_floe++;
            }
        }
    }
    P.get_floe_group().set_initial_window(new_window);
    P.make_input_file();

    return 0;
}