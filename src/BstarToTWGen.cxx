#include "UHH2/BstarToTW/include/BstarToTWGen.h"

using namespace std;
using namespace uhh2;

BstarToTWGen::BstarToTWGen(const vector<GenParticle> & genparticles, bool throw_on_failure): m_type(e_notfound) {    
  int n_bstar = 0;
  for(unsigned int i=0; i<genparticles.size(); ++i) {
    const GenParticle & genp = genparticles[i];
    if (abs(genp.pdgId()) == 1005){
      auto w = genp.daughter(&genparticles, 1);
      auto t = genp.daughter(&genparticles, 2);
      if(!w || !t){
	if(throw_on_failure) throw runtime_error("BstarToTWGen: bstar has not ==2 daughters");
	return;
      }
      if(abs(w->pdgId()) != 24){
	std::swap(w, t);
	if(abs(w->pdgId()) != 24){
	  if(throw_on_failure) throw runtime_error("BstarToTWGen: bstar has no W daughter");
	  return;
	}
      }
      
      // NOTE: here, we could skip over intermediate W bosons. However,
      // this Pythia8-related problem is now fixed when creating ntuples already,
      // so this should not be necessary.
      
      if(abs(t->pdgId()) != 6){
	if(throw_on_failure) throw runtime_error("BstarToTWGen: bstar has no t daughter");
	return;
      }

      // get W daughters:
      auto wd1 = w->daughter(&genparticles, 1);
      auto wd2 = w->daughter(&genparticles, 2);
      if(!wd1 || !wd2){
	if(throw_on_failure) throw runtime_error("BstarToTWGen: W from Bstar decay has not ==2 daughters");
	return;
      }

      // get t daughters:
      auto tw = t->daughter(&genparticles, 1);
      auto tb = t->daughter(&genparticles, 2);      
      if(!tw || !tb){
	if(throw_on_failure) throw runtime_error("BstarToTWGen: t has not ==2 daughters");
	return;
      }
      if(abs(tw->pdgId()) != 24){
	std::swap(tw, tb);
	if(abs(tw->pdgId()) != 24){
	  if(throw_on_failure) throw runtime_error("BstarToTWGen: t has no W daughter");
	  return;
	}
      }
      // get W daughters of tW:
      auto twd1 = tw->daughter(&genparticles, 1);
      auto twd2 = tw->daughter(&genparticles, 2);
      if(!twd1 || !twd2){
	if(throw_on_failure) throw runtime_error("BstarToTWGen: W from t decay has not ==2 daughters");
	return;
      }

      // now that we collected everything, fill the member variables. 
      m_bstar = genp;
      m_Wbstar = *w;
      m_tbstar = *t;
      m_Wdecay1 = *wd1;
      m_Wdecay2 = *wd2;
      m_tW = *tw;
      m_tb = *tb;
      m_tWdecay1 = *twd1;
      m_tWdecay2 = *twd2;
      ++n_bstar;
    }
  }
  // check if bstar was in the event
  if (n_bstar != 1) {
    if(throw_on_failure)  throw runtime_error("BstarToTWGen: did not find exactly one bstar in the event");
    return;
  }

  // calculate decay channel by counting the number of charged leptons
  // in the W daughters:
  int n_e = 0, n_m = 0, n_t = 0;
  for(const auto & wd : {m_Wdecay1, m_Wdecay2, m_tWdecay1, m_tWdecay2}){
    int id = abs(wd.pdgId());
    if(id == 11) ++n_e;
    else if(id == 13) ++n_m;
    else if(id == 15) ++n_t;
  }

  // dilepton channels:
  if(n_e == 2){
    m_type = e_ee;
  }
  else if(n_e == 1 && n_m == 1){
    m_type = e_emu;
  }
  else if(n_e == 1 && n_t == 1){
    m_type = e_etau;
  }
  else if(n_m == 2){
    m_type = e_mumu;
  }
  else if(n_m == 1 && n_t == 1){
    m_type = e_mutau;
  }
  else if(n_t == 2){
    m_type = e_tautau;
  }

  // lepton+jet channels:
  else if(n_e == 1){
    m_type = e_ehad;
  }
  else if(n_m == 1){
    m_type = e_muhad;
  }
  else if(n_t == 1){
    m_type = e_tauhad;
  }

  // hadronic:
  else{
    m_type = e_had;
  }
}   
 

GenParticle BstarToTWGen::bstar() const{
  return m_bstar;
}

GenParticle BstarToTWGen::Wbstar() const{
  return m_Wbstar;
}

GenParticle BstarToTWGen::tbstar() const{
  return m_tbstar;
}

GenParticle BstarToTWGen::Wdecay1() const{
  return m_Wdecay1;
} 

GenParticle BstarToTWGen::Wdecay2() const{
  return m_Wdecay2;
} 

GenParticle BstarToTWGen::tW() const{
  return m_tW;
} 

GenParticle BstarToTWGen::tb() const{
  return m_tb;
}

GenParticle BstarToTWGen::tWdecay1() const{
  return m_tWdecay1;
} 

GenParticle BstarToTWGen::tWdecay2() const{
  return m_tWdecay2;
} 

BstarToTWGen::E_DecayChannel BstarToTWGen::DecayChannel()  const{  
  return m_type;
}

bool BstarToTWGen::IsPrimaryHadronicDecay()  const{
  return abs(m_Wdecay1.pdgId()) <= 5;
}

bool BstarToTWGen::IsSecundaryHadronicDecay()  const{
  return abs(m_tWdecay1.pdgId()) <= 5;
}

bool BstarToTWGen::IsSemiLeptonicDecay() const{
  return m_type == e_ehad ||  m_type == e_muhad;  
}


namespace {
    
  bool is_charged_lepton(const GenParticle & gp){
    int id = abs(gp.pdgId());
    return id == 11 || id == 13 || id == 15;
  }

  bool is_neutrino(const GenParticle & gp){
    int id = abs(gp.pdgId());
    return id == 12 || id == 14 || id == 16;
  }

}

GenParticle BstarToTWGen::ChargedLepton() const{
  if (m_type != e_ehad &&  m_type != e_muhad && m_type!= e_tauhad){
    throw runtime_error("BstarToTWGen::ChargedLepton called, but this is no l+jets bstar event!");
  }
  for(const auto & wd : {m_Wdecay1, m_Wdecay2, m_tWdecay1, m_tWdecay2}){
    if(is_charged_lepton(wd)) return wd;
  }
  throw logic_error("logic error in BstarToTWGen::ChargedLepton");
}

GenParticle BstarToTWGen::Neutrino() const{
  if (m_type != e_ehad &&  m_type != e_muhad && m_type!= e_tauhad){
    throw runtime_error("BstarToTWGen::ChargedLepton called, but this is no l+jets ttbar event!");
  }
  for(const auto & wd : {m_Wdecay1, m_Wdecay2, m_tWdecay1, m_tWdecay2}){
    if(is_neutrino(wd)) return wd;
  }
  throw logic_error("logic error in BstarToTWGen::Neutrino");
}

// GenParticle BstarToTWGen::TopLep() const{
//   if(ChargedLepton().charge()>0) return Top();
//   else return Antitop();
// }

// GenParticle BstarToTWGen::TopHad() const{
//   if(ChargedLepton().charge()<0) return Top();
//   else return Antitop();
// }

// GenParticle BstarToTWGen::BLep() const{
//   if(ChargedLepton().charge()>0) return bTop();
//   else return bAntitop();
// }

// GenParticle BstarToTWGen::BHad() const{
//   if(ChargedLepton().charge()<0) return bTop();
//   else return bAntitop();
// }

// GenParticle BstarToTWGen::WLep() const{
//   if(ChargedLepton().charge()>0) return WTop();
//   else return WAntitop();
// }

// GenParticle BstarToTWGen::WHad() const{
//   if(ChargedLepton().charge()<0) return WTop();
//   else return WAntitop();
// }

// GenParticle BstarToTWGen::Q1() const{
//   if(ChargedLepton().charge()>0) return WMinusdecay1();
//   else return Wdecay1();
// }

// GenParticle BstarToTWGen::Q2() const{
//   if(ChargedLepton().charge()>0) return WMinusdecay2();
//   else return Wdecay2();
// }


BstarToTWGenProducer::BstarToTWGenProducer(uhh2::Context & ctx, const std::string & name, bool throw_on_failure_): throw_on_failure(throw_on_failure_){
  h_BstarToTWgen = ctx.get_handle<BstarToTWGen>(name);
}

bool BstarToTWGenProducer::process(Event & event){
  event.set(h_BstarToTWgen, BstarToTWGen(*event.genparticles, throw_on_failure));
  return true;
}
