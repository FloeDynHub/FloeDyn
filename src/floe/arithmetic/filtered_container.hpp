/*!
 * \file floe/floes/partial_floe_group.hpp
 * \brief Partial Floe Configuration class (subset of floes)
 * \author Quentin Jouet
 */

#ifndef FILTERED_CONTAINER_HPP
#define FILTERED_CONTAINER_HPP

#include <vector>


template<typename TContainer>
struct SubIterator {
    SubIterator(TContainer& obj, std::size_t n) : m_obj{&obj}, m_id{n} {};
    TContainer* m_obj;
    std::size_t m_id;
    SubIterator& operator++() { ++m_id; return *this; } //prefix increment
    SubIterator operator++(int) { auto old_this = *this; ++m_id; return old_this; }
    typename TContainer::value_type& operator*() { return (*m_obj)[m_id]; }
    bool operator==(SubIterator const& it) const { return it.m_id == this->m_id; }
    bool operator!=(SubIterator const& it) const { return !(it==*this); }
};

template<typename TContainer>
struct ConstSubIterator {
    ConstSubIterator(TContainer const& obj, std::size_t n) : m_obj{&obj}, m_id{n} {};
    TContainer const* m_obj;
    std::size_t m_id;
    ConstSubIterator& operator++() { ++m_id; return *this; } //prefix increment
    ConstSubIterator operator++(int) { auto old_this = *this; ++m_id; return old_this; }
    typename TContainer::value_type const& operator*() const { return (*m_obj)[m_id]; }
    bool operator==(ConstSubIterator const& it) const { return it.m_id == this->m_id; }
    bool operator!=(ConstSubIterator const& it) const { return !(it==*this); }
};

template<typename T>
class FilteredVector : public std::vector<T>{
public:
    using base_class = std::vector<T>;
    using value_type = T;
    using iterator_type = SubIterator<FilteredVector<T>>;
    using const_iterator_type = ConstSubIterator<FilteredVector<T>>;
    FilteredVector() : base_class(), m_filter{false} {}
    FilteredVector(std::initializer_list<value_type> il) : base_class(il), m_filter{false} {}
    T const& operator[](std::size_t n) const { std::size_t id = m_filter ? m_ids[n] : n; return base_class::operator[](id); }
    T& operator[](std::size_t n) { std::size_t id = m_filter ? m_ids[n] : n; return base_class::operator[](id); }
    T& operator()(std::size_t n) { return base_class::operator[](n); }
    T const& operator()(std::size_t n) const { return base_class::operator[](n); }
    void update_ids(std::vector<std::size_t> ids){ m_filter = true; m_ids = ids; }
    std::size_t absolute_id(std::size_t n) const { return m_filter ? m_ids[n] : n; }
    const_iterator_type begin() const { return const_iterator_type(*this, 0); }
    const_iterator_type end() const { return const_iterator_type(*this, size()); }
    iterator_type begin() { return iterator_type(*this, 0); }
    iterator_type end() { return iterator_type(*this, size()); }
    std::size_t size() const {
        return m_filter ? m_ids.size() : base_class::size();
    }
private:
    bool m_filter;
    std::vector<std::size_t> m_ids;
};


#endif // FILTERED_CONTAINER_HPP
