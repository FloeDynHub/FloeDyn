#ifndef FLOE_STATE_STATE_OBSERVER_HPP
#define FLOE_STATE_STATE_OBSERVER_HPP

#include <vector>
#include <cstdef>

namespace flo { namespace state
{

/*! Observable pattern for state
 */
class StateObservable
{
    void attach( StateObserver 
private:
    std::vector<StateObserver*> m_observers;
};

}} // namespace floe::state
#endif // FLOE_STATE_STATE_OBSERVER_HPP

