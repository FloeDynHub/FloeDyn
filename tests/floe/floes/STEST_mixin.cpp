/*!
 * \file floe/floes/TEST_mixin.cpp
 * \brief Test file for mixin technics.
 * \author Roland Denis
 */

#include <iostream>
#include <cstddef>

template <typename T>
class Base
{
    public:
        T test;
        Base(T test) : test(test) {}
        T get_test() { return test; };
};

template < typename T>
struct Identity
{
    T toto;
};

template <
    typename T,
    typename Base
>
class Mixin : public Base, public Identity<T>
{
    public:
        using Base::Base;
};

int main()
{
    Mixin<std::size_t, Base<double>> mixin_base(3.);
    mixin_base.toto = 1;
    mixin_base.test = 0.;
    auto t = mixin_base.get_test();

    return 0;
}
