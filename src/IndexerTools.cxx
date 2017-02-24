#include "UHH2/BstarToTW/include/IndexerTools.h"

#include <vector>
#include <iostream>

using namespace std;
using namespace uhh2;



/** 
 * Remove indices with pt < pt_min and eta > eta_max
 */
PtEtaTopIndexCleaner::PtEtaTopIndexCleaner(Context &ctx, double pt_min_, double eta_max_) {
  h_TopTagIndexer = ctx.get_handle<TopTagIndexer>("TopTagIndexer");
  pt_min = pt_min_;
  eta_max = eta_max_;
}

bool PtEtaTopIndexCleaner::process(Event &event) {
  vector<TopJet> *topjets = event.topjets;
  TopTagIndexer indexer = event.get(h_TopTagIndexer);
  vector<int> ind = indexer.GetIndex(); 
  vector<int> selInd;
  for (int i : ind)
    { 
      TopJet topjet = topjets->at(i);
      double pt = topjet.v4().pt();
      double eta = topjet.v4().eta();
      if ( (pt > pt_min) && (abs(eta) < eta_max) ) selInd.push_back(i);
    }
  indexer.SetIndex(selInd);
  event.set(h_TopTagIndexer, indexer);
  return (selInd.size() > 0);
}

/** 
 * Remove indices that are not within m_min < m < m_max
 */
MTopIndexCleaner::MTopIndexCleaner(Context &ctx, double m_min_, double m_max_) {
  h_TopTagIndexer = ctx.get_handle<TopTagIndexer>("TopTagIndexer");
  m_min = m_min_;
  m_max = m_max_;
}

bool MTopIndexCleaner::process(Event &event) {
  vector<TopJet> *topjets = event.topjets;
  TopTagIndexer indexer = event.get(h_TopTagIndexer);
  vector<int> ind = indexer.GetIndex();
  vector<int> selInd;
  for (int i : ind)
    { 
      TopJet topjet = topjets->at(i);
      double m = topjet.v4().M();
      if (m_min < m && m < m_max) selInd.push_back(i);
    }
  indexer.SetIndex(selInd);
  event.set(h_TopTagIndexer, indexer);
  return (selInd.size() > 0);
}

/*
 * Remove indices with N_subjets < nsub_min
 */
NSubTopIndexCleaner::NSubTopIndexCleaner(Context &ctx, unsigned int nsub_min_) {
  h_TopTagIndexer = ctx.get_handle<TopTagIndexer>("TopTagIndexer");
  nsub_min = nsub_min_;
}

bool NSubTopIndexCleaner::process(Event & event) {
  vector<TopJet> *topjets = event.topjets;
  TopTagIndexer indexer = event.get(h_TopTagIndexer);
  vector<int> ind = indexer.GetIndex();
  vector<int> selInd;
  for (int i : ind)
    { 
      vector<Jet> subjets = topjets->at(i).subjets();
      unsigned int nsub = subjets.size();
      if (nsub_min <= nsub) selInd.push_back(i);
    }
  indexer.SetIndex(selInd);
  event.set(h_TopTagIndexer, indexer);
  return (selInd.size() > 0);
}

/*
 * Remove indices with fpt > fpt_max
 */
FptTopIndexCleaner::FptTopIndexCleaner(Context &ctx, double fpt_max_) {
  h_TopTagIndexer = ctx.get_handle<TopTagIndexer>("TopTagIndexer");
  fpt_max = fpt_max_;
}

bool FptTopIndexCleaner::process(Event &event) {
  vector<TopJet> *topjets = event.topjets;
  TopTagIndexer indexer = event.get(h_TopTagIndexer);
  vector<int> ind = indexer.GetIndex();
  vector<int> selInd;
  for (int i : ind)
    { 
      vector<Jet> subjets = topjets->at(i).subjets();
      double pt = topjets->at(i).v4().pt();
      double ptsub = subjets.at(0).v4().pt();
      if (subjets.size() < 1) continue;
      double fpt = ptsub/pt;
      if (fpt < fpt_max) selInd.push_back(i);
    }
  indexer.SetIndex(selInd);
  event.set(h_TopTagIndexer, indexer);
  return (selInd.size() > 0);
}

/*
 * Remove indices with mpair < mpair_min
 */
MpairTopIndexCleaner::MpairTopIndexCleaner(Context &ctx, double mpair_min_) {
  h_TopTagIndexer = ctx.get_handle<TopTagIndexer>("TopTagIndexer");
  mpair_min = mpair_min_;
}

bool MpairTopIndexCleaner::process(Event &event) {
  vector<TopJet> *topjets = event.topjets;
  TopTagIndexer indexer = event.get(h_TopTagIndexer);
  vector<int> ind = indexer.GetIndex();
  vector<int> selInd;
  for (int i : ind)
    { 
      vector<Jet> subjets = topjets->at(i).subjets();
      if (subjets.size() < 3) continue;
      double m12 = (subjets.at(0).v4() + subjets.at(1).v4()).M();
      double m13 = (subjets.at(0).v4() + subjets.at(2).v4()).M();
      double m23 = (subjets.at(1).v4() + subjets.at(2).v4()).M();
      double mpair = min(min(m12, m13), m23);
      if (mpair_min < mpair) selInd.push_back(i);
    }
  indexer.SetIndex(selInd);
  event.set(h_TopTagIndexer, indexer);
  return (selInd.size() > 0);
}

/*
 * Remove indices with tau3/tau2 > t32_max
 */
Tau32TopIndexCleaner::Tau32TopIndexCleaner(Context &ctx, double t32_max_) {
  h_TopTagIndexer = ctx.get_handle<TopTagIndexer>("TopTagIndexer");
  t32_max = t32_max_;
}

bool Tau32TopIndexCleaner::process(Event &event) {
  vector<TopJet> *topjets = event.topjets;
  TopTagIndexer indexer = event.get(h_TopTagIndexer);
  vector<int> ind = indexer.GetIndex();
  vector<int> selInd;
  for (int i : ind)
    { 
      TopJet topjet = topjets->at(i);
      if (topjet.tau3_groomed()/topjet.tau2_groomed() < t32_max) selInd.push_back(i);
    }
  indexer.SetIndex(selInd);
  event.set(h_TopTagIndexer, indexer);
  return (selInd.size() > 0);
}

/*
 * Remove indices with delta Phi between topjet and muon < dphi_min
 */
DeltaPhiTopIndexCleaner::DeltaPhiTopIndexCleaner(Context &ctx, double dphi_min_) {
  h_TopTagIndexer = ctx.get_handle<TopTagIndexer>("TopTagIndexer");
  dphi_min = dphi_min_;
}

bool DeltaPhiTopIndexCleaner::process(Event &event) {
  vector<TopJet> *topjets = event.topjets;
  Muon muon = event.muons->at(0);
  TopTagIndexer indexer = event.get(h_TopTagIndexer);
  vector<int> ind = indexer.GetIndex();
  vector<int> selInd;
  for (int i : ind)
    { 
      TopJet topjet = topjets->at(i);
      if (dphi_min < deltaPhi(topjet.v4(), muon.v4())) selInd.push_back(i);
    }
  indexer.SetIndex(selInd);
  event.set(h_TopTagIndexer, indexer);
  return (selInd.size() > 0);
}




/*
 * Selection for events with == 1 topjet
 */
NTopIndexSelection::NTopIndexSelection(Context &ctx, unsigned int n_min_, unsigned int n_max_) {
  h_TopTagIndexer = ctx.get_handle<TopTagIndexer>("TopTagIndexer");
  n_min = n_min_;
  n_max = n_max_;
}

bool NTopIndexSelection::passes(const Event & event) {
  TopTagIndexer indexer = event.get(h_TopTagIndexer);
  unsigned int n = indexer.GetIndex().size();
  return (n_min <= n && n <= n_max);
}

// PtEta
PtEtaTopIndexSelection::PtEtaTopIndexSelection(Context &ctx, double pt_min_, double eta_max_) {
  h_TopTagIndexer = ctx.get_handle<TopTagIndexer>("TopTagIndexer");
  pt_min = pt_min_;
  eta_max = eta_max_;
}

bool PtEtaTopIndexSelection::passes(const Event &event) {
  vector<TopJet> *topjets = event.topjets;
  TopTagIndexer indexer = event.get(h_TopTagIndexer);
  vector<int> ind = indexer.GetIndex();
  unsigned int n = 0;
  for (int i : ind)
    { 
      TopJet topjet = topjets->at(i);
      double pt = topjet.v4().pt();
      double eta = topjet.v4().eta();
      if ( (pt > pt_min) && (abs(eta) < eta_max) ) ++n;
    }
  return (n > 0);
}

// MTop
MTopIndexSelection::MTopIndexSelection(Context &ctx, double m_min_, double m_max_) {
  h_TopTagIndexer = ctx.get_handle<TopTagIndexer>("TopTagIndexer");
  m_min = m_min_;
  m_max = m_max_;
}

bool MTopIndexSelection::passes(const Event &event) {
  vector<TopJet> *topjets = event.topjets;
  TopTagIndexer indexer = event.get(h_TopTagIndexer);
  vector<int> ind = indexer.GetIndex();
  unsigned int n = 0;
  for (int i : ind)
    { 
      TopJet topjet = topjets->at(i);
      double m = topjet.v4().M();
      if (m_min < m && m < m_max) ++n;
    }
  return (n > 0);
}

// Nsub
NSubTopIndexSelection::NSubTopIndexSelection(Context &ctx, unsigned int nsub_min_) {
  h_TopTagIndexer = ctx.get_handle<TopTagIndexer>("TopTagIndexer");
  nsub_min = nsub_min_;
}

bool NSubTopIndexSelection::passes(const Event & event) {
  vector<TopJet> *topjets = event.topjets;
  TopTagIndexer indexer = event.get(h_TopTagIndexer);
  vector<int> ind = indexer.GetIndex();
  unsigned int n = 0;
  for (int i : ind)
    { 
      vector<Jet> subjets = topjets->at(i).subjets();
      unsigned int nsub = subjets.size();
      if (nsub_min <= nsub) ++n;
    }
  return (n > 0);
}

// Fpt
FptTopIndexSelection::FptTopIndexSelection(Context &ctx, double fpt_max_) {
  h_TopTagIndexer = ctx.get_handle<TopTagIndexer>("TopTagIndexer");
  fpt_max = fpt_max_;
}

bool FptTopIndexSelection::passes(const Event &event) {
  vector<TopJet> *topjets = event.topjets;
  TopTagIndexer indexer = event.get(h_TopTagIndexer);
  vector<int> ind = indexer.GetIndex();
  unsigned int n = 0;
  for (int i : ind)
    { 
      vector<Jet> subjets = topjets->at(i).subjets();
      double pt = topjets->at(i).v4().pt();
      double ptsub = subjets.at(0).v4().pt();
      if (subjets.size() < 1) continue;
      double fpt = ptsub/pt;
      if (fpt < fpt_max) ++n;
    }
  return (n > 0);
}

//Mpair
MpairTopIndexSelection::MpairTopIndexSelection(Context &ctx, double mpair_min_) {
  h_TopTagIndexer = ctx.get_handle<TopTagIndexer>("TopTagIndexer");
  mpair_min = mpair_min_;
}

bool MpairTopIndexSelection::passes(const Event &event) {
  vector<TopJet> *topjets = event.topjets;
  TopTagIndexer indexer = event.get(h_TopTagIndexer);
  vector<int> ind = indexer.GetIndex();
  unsigned int n = 0;
  for (int i : ind)
    { 
      vector<Jet> subjets = topjets->at(i).subjets();
      if (subjets.size() < 3) continue;
      double m12 = (subjets.at(0).v4() + subjets.at(1).v4()).M();
      double m13 = (subjets.at(0).v4() + subjets.at(2).v4()).M();
      double m23 = (subjets.at(1).v4() + subjets.at(2).v4()).M();
      double mpair = min(min(m12, m13), m23);
      if (mpair_min < mpair) ++n;
    }
  return (n > 0);
}

// tau32
Tau32TopIndexSelection::Tau32TopIndexSelection(Context &ctx, double t32_max_) {
  h_TopTagIndexer = ctx.get_handle<TopTagIndexer>("TopTagIndexer");
  t32_max = t32_max_;
}

bool Tau32TopIndexSelection::passes(const Event &event) {
  vector<TopJet> *topjets = event.topjets;
  TopTagIndexer indexer = event.get(h_TopTagIndexer);
  vector<int> ind = indexer.GetIndex();
  unsigned int n = 0;
  for (int i : ind)
    { 
      TopJet topjet = topjets->at(i);
      if (topjet.tau3_groomed()/topjet.tau2_groomed() < t32_max) ++n;
    }
  return (n > 0);
}

//delta phi
DeltaPhiTopIndexSelection::DeltaPhiTopIndexSelection(Context &ctx, double dphi_min_) {
  h_TopTagIndexer = ctx.get_handle<TopTagIndexer>("TopTagIndexer");
  dphi_min = dphi_min_;
}

bool DeltaPhiTopIndexSelection::passes(const Event &event) {
  vector<TopJet> *topjets = event.topjets;
  Muon muon = event.muons->at(0);
  TopTagIndexer indexer = event.get(h_TopTagIndexer);
  vector<int> ind = indexer.GetIndex();
  unsigned int n = 0;
  for (int i : ind)
    { 
      TopJet topjet = topjets->at(i);
      if (dphi_min < deltaPhi(topjet.v4(), muon.v4())) ++n;
    }
  return (n > 0);
}
