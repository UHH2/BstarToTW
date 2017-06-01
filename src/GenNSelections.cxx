#include "UHH2/BstarToTW/include/GenNSelections.h"
#include "UHH2/common/include/NSelections.h"

using namespace uhh2;
using namespace std;

namespace {
    
template<typename T>
bool passes_minmax(const vector<T> & objects, int nmin, int nmax, const Event & event, const boost::optional<std::function<bool (const T &, const Event & )>> & object_id){
    int n_objects = objects.size();
    if(object_id){
        n_objects = 0;
        for(const auto & obj : objects){
            if((*object_id)(obj, event)) ++n_objects;
        }
    }
    return n_objects >= nmin && (nmax < 0 || n_objects <= nmax);
}
    
}

NGenTopJetSelection::NGenTopJetSelection(int nmin, int nmax,
    const boost::optional<GenTopJetId> & gentopjetid,
    const boost::optional<Event::Handle<std::vector<GenTopJet> > > & gentopjetcollection) :
    m_nmin(nmin), m_nmax(nmax), m_gentopjetid(gentopjetid), m_gentopjetcollection(gentopjetcollection){}

bool NGenTopJetSelection::passes(const Event & event){
    const auto & jets = m_gentopjetcollection ? event.get(*m_gentopjetcollection) : *event.gentopjets;
    return passes_minmax(jets, m_nmin, m_nmax, event, m_gentopjetid);
}
