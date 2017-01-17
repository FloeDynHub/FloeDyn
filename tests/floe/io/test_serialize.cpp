#include "../tests/catch.hpp"
#include <iostream>
#include <vector>

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <sstream>

struct MyData{

    MyData() : m_x{0}, m_y{0.3}, m_v{}, m_z{0} {}
    MyData(int x, double y, std::vector<double> v, long z) : m_x{x}, m_y{y}, m_v{v}, m_z{z} {}

    // This method lets cereal know which data members to serialize
    template<class Archive>
    void serialize(Archive & archive)
    {
    archive( m_x, m_y, m_v, m_z ); // serialize things by passing them to the archive
    }

    long z(){ return m_z; }

    int m_x;
    double m_y;
    std::vector<double> m_v;
private:
    long m_z;
};


TEST_CASE( "Test serialization with cereal", "[io]" ) {

    std::stringstream ss; // any stream can be used

    {
        cereal::BinaryOutputArchive oarchive(ss); // Create an output archive
        int m1{1}, m2{2}, m3{3};
        MyData md{4, 1.23, {{5, 6, 7}}, 9};
        oarchive(m1, m2, m3, md); // Write the data to the archive
    } // archive goes out of scope, ensuring all contents are flushed

    {
        cereal::BinaryInputArchive iarchive(ss); // Create an input archive

        int m1, m2, m3;
        MyData md;
        iarchive(m1, m2, m3, md); // Read the data from the archive

        std::cout << m1 << m2 << m3 << " " << md.m_x << " " << md.m_y << " " << md.z() << std::endl;
        for (auto d : md.m_v) std::cout << d;
        std::cout << std::endl;
    }

}
