#include "UHH2/BstarToTW/include/BstarToTWGen.h"

using namespace std;
using namespace uhh2;

/* ------------------------------------------------------------------
 * BstarToTWGen - Interprets an b* -> tW event by going through the
 * GenParticle list and identifying all GenParticles from the hard
 * process. It also gets a list of all stable GenParticles
 * (status() == 1).
 */

// possible to write this more efficient..
BstarToTWGen::BstarToTWGen(const Event & event)
{
  vector<GenParticle> genparticles = *event.genparticles;
  vector<GenJet> genjets = *event.genjets;
  // Get first 11 genparticles (this should be all particles from the hard process) and match them
  vector<GenParticle> hardGenparticles(&genparticles[0], &genparticles[10]);

  is_singlemuo = false;
  is_singleele = false;
  bool has_bstar = false;
  for(unsigned int i = 0; i < genparticles.size(); ++i) 
    {
      const GenParticle & genp = genparticles[i];
      // is b*
      if (abs(genp.pdgId()) == 1005)
	{
	  m_bstar = genp.v4();
	  has_bstar = true;
	}
      // is top
      else if (abs(genp.pdgId()) == 6)
	{
	  m_tbstar = genp.v4();
	  auto tw = genp.daughter(&genparticles, 1);
	  auto tb = genp.daughter(&genparticles, 2);
	  if (abs(tw->pdgId()) != 24) std::swap(tw, tb);
	  // if W from t decay was skipped find daughters and reconstruct it
	  if (abs(tw->pdgId()) != 24)
	    {
	      if (abs(tw->pdgId()) == 5) std::swap(tw, tb);
	      auto twd1 = tw;
	      auto twd2 = genp.daughter(&genparticles, 3);
	      m_tWdecay1 = twd1->v4();
	      m_tWdecay2 = twd2->v4();
	      m_tb = tb->v4();
	      m_tW = m_tWdecay1 + m_tWdecay2;
	      if (abs(twd1->pdgId()) >= 11 && abs(twd2->pdgId()) >= 11 && abs(twd1->pdgId()) < 15 && abs(twd2->pdgId()) < 15)
		{
		  is_tophad = false;
		  if (twd1->charge() != 0)
		    {
		      m_cl = m_tWdecay1;
		      m_nl = m_tWdecay2;
		    }
		  else
		    {
		      m_cl = m_tWdecay1;
		      m_nl = m_tWdecay2;
		    }
		  if (abs(twd1->pdgId()) == 11 || abs(twd2->pdgId()) == 11) is_singleele = true;
		  else if (abs(twd1->pdgId()) == 13 || abs(twd2->pdgId()) == 13) is_singlemuo = true;
		}
	      else
		{
		  is_tophad = true;
		}
	    }
	  else
	    {
	      auto twd1 = tw->daughter(&genparticles, 1);
	      auto twd2 = tw->daughter(&genparticles, 2);
	      m_tWdecay1 = twd1->v4();
	      m_tWdecay2 = twd2->v4();
	      m_tb = tb->v4();
	      m_tW = tw->v4();
	      if (abs(twd1->pdgId()) >= 11 && abs(twd2->pdgId()) >= 11 && abs(twd1->pdgId()) < 15 && abs(twd2->pdgId()) < 15)
		{
		  is_tophad = false;
		  if (twd1->charge() != 0)
		    {
		      m_cl = m_tWdecay1;
		      m_nl = m_tWdecay2;
		    }
		  else
		    {
		      m_cl = m_tWdecay1;
		      m_nl = m_tWdecay2;
		    }
		  if (abs(twd1->pdgId()) == 11 || abs(twd2->pdgId()) == 11) is_singleele = true;
		  else if (abs(twd1->pdgId()) == 13 || abs(twd2->pdgId()) == 13) is_singlemuo = true;
		}
	      else is_tophad = true;
	    }  
	}
      // is W from b*
      else if (abs(genp.pdgId()) == 24 && abs(genp.mother(&genparticles, 1)->pdgId()) != 6)
	{
	  m_Wbstar = genp.v4();
	  auto wd1 = genp.daughter(&genparticles, 1);
	  auto wd2 = genp.daughter(&genparticles, 2);
	  m_Wdecay1 = wd1->v4();
	  m_Wdecay2 = wd2->v4();
	  if (abs(wd1->pdgId()) >= 11 && abs(wd2->pdgId()) >= 11 && abs(wd1->pdgId()) < 15 && abs(wd2->pdgId()) < 15)
	    {
	      is_whad = false;
	      if (wd1->charge() != 0)
		{
		  m_cl = m_Wdecay1;
		  m_nl = m_Wdecay2;
		}
	      else
		{
		  m_cl = m_Wdecay1;
		  m_nl = m_Wdecay2;
		}
	      if (abs(wd1->pdgId()) == 11 || abs(wd2->pdgId()) == 11) is_singleele = true;
	      else if (abs(wd1->pdgId()) == 13 || abs(wd2->pdgId()) == 13) is_singlemuo = true;
	    }
	  else is_whad = true;
	}
    }
  // if b* was skipped reconstruct it from daughters
  if (!has_bstar)
    {
      m_bstar = m_tbstar + m_Wbstar;
    }

  // Get stable Particles (except neutrinos)
  for (unsigned int i = 0; i < genparticles.size(); ++i){
    const GenParticle & genp = genparticles[i];
    if ( genp.status() == 1 ){
      if (abs(genp.pdgId() != 12) || abs(genp.pdgId() != 14) || abs(genp.pdgId() != 16) ) m_stable.push_back(genp.v4());
      if (abs(genp.pdgId()) == 11) m_ele.push_back(genp.v4());
      else if (abs(genp.pdgId()) == 13) m_muo.push_back(genp.v4());
    }
  } 
}


LorentzVector BstarToTWGen::bstar() const{
  return m_bstar;
}

LorentzVector BstarToTWGen::Wbstar() const{
  return m_Wbstar;
}

LorentzVector BstarToTWGen::tbstar() const{
  return m_tbstar;
}

LorentzVector BstarToTWGen::Wdecay1() const{
  return m_Wdecay1;
} 

LorentzVector BstarToTWGen::Wdecay2() const{
  return m_Wdecay2;
} 

LorentzVector BstarToTWGen::tW() const{
  return m_tW;
} 

LorentzVector BstarToTWGen::tb() const{
  return m_tb;
}

LorentzVector BstarToTWGen::tWdecay1() const{
  return m_tWdecay1;
} 

LorentzVector BstarToTWGen::tWdecay2() const{
  return m_tWdecay2;
} 

LorentzVector BstarToTWGen::ChargedLepton() const{
  return m_cl;
}

LorentzVector BstarToTWGen::Neutrino() const{
  return m_nl;
}


vector<LorentzVector> BstarToTWGen::stable() const{
  return m_stable;
}

vector<LorentzVector> BstarToTWGen::electrons() const{
  return m_ele;
}

vector<LorentzVector> BstarToTWGen::muons() const{
  return m_muo;
}


bool BstarToTWGen::IsSemiLeptonicDecay() const{
  return (is_tophad != is_whad);
}

bool BstarToTWGen::IsTopHadronicDecay() const{
  return is_tophad;
}

bool BstarToTWGen::IsWHadronicDecay() const{
  return is_whad;
}

bool BstarToTWGen::IsMuonDecay() const{
  return is_singlemuo;
}

bool BstarToTWGen::IsElectronDecay() const{
  return is_singleele;
}

BstarToTWGenProducer::BstarToTWGenProducer(uhh2::Context & ctx, const std::string & name)
{
  h_BstarToTWgen = ctx.get_handle<BstarToTWGen>(name);
}

bool BstarToTWGenProducer::process(Event & event){
  event.set(h_BstarToTWgen, BstarToTWGen(event));
  return true;
}
