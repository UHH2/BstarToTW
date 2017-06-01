#pragma once

#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/ObjectIdUtils.h"

#include "UHH2/BstarToTW/include/HOTVRGenIds.h"

class NGenTopJetSelection: public uhh2::Selection {
public:
    explicit NGenTopJetSelection(int nmin, int nmax = -1,
        const boost::optional<GenTopJetId> & gentopjetid = boost::none,
				 const boost::optional<uhh2::Event::Handle<std::vector<GenTopJet> > > & gentopjetcollection = boost::none);
    virtual bool passes(const uhh2::Event & event);
private:
    int m_nmin, m_nmax;
    boost::optional<GenTopJetId> m_gentopjetid;
    boost::optional<uhh2::Event::Handle<std::vector<GenTopJet> > > m_gentopjetcollection;
};
