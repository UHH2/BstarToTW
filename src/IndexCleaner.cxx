#include "UHH2/BstarToTW/include/IndexCleaner.h"

#include <vector>
#include <iostream>

using namespace std;
using namespace uhh2;

/* 
 * Index cleaner for pt_top
 */
PtTopIndexCleaner::PtTopIndexCleaner(Context &ctx, double pt_min_) {
  h_TopTagIndexer = ctx.get_handle<TopTagIndexer>("TopTagIndexer");
  pt_min = pt_min_;
}

bool PtTopIndexCleaner::process(Event &event) {
  vector<TopJet> *topjets = event.topjets;
  TopTagIndexer indexer = event.get(h_TopTagIndexer);
  vector<int> ind = indexer.GetIndex(); 
  vector<int> selInd;
  for (int i : ind)
    { 
      TopJet topjet = topjets->at(i);
      double pt = topjet.v4().pt();
      if (pt > pt_min) selInd.push_back(i);
    }
  indexer.SetIndex(selInd);
  event.set(h_TopTagIndexer, indexer);
  return (selInd.size() > 0);
}

/* 
 * Index cleaner for M_top
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
 * Index cleaner for eta_top
 */
EtaTopIndexCleaner::EtaTopIndexCleaner(Context &ctx, double eta_max_) {
  h_TopTagIndexer = ctx.get_handle<TopTagIndexer>("TopTagIndexer");
  eta_max = eta_max_;
}

bool EtaTopIndexCleaner::process(Event & event) {
  vector<TopJet> *topjets = event.topjets;
  TopTagIndexer indexer = event.get(h_TopTagIndexer);
  vector<int> ind = indexer.GetIndex();
  vector<int> selInd;
  for (int i : ind)
    { 
      TopJet topjet = topjets->at(i);
      double eta = abs(topjet.v4().eta());
      if (eta < eta_max) selInd.push_back(i);
    }
  indexer.SetIndex(selInd);
  event.set(h_TopTagIndexer, indexer);
  return (selInd.size() > 0);
}

/*
 * Index cleaner for number of subjets
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
 * Index cleaner for pt fraction of first subjet
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
 * Index cleaner for pairwise mass of first three subjets
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
 * Index cleaner for tau32 of topjet
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
 * Index cleaner for HOTVR criteria
 */
HotvrTopIndexCleaner::HotvrTopIndexCleaner(Context &ctx, double pt_min_, double eta_max_, unsigned int nsub_min_, double fpt_max_, double mpair_min_, double tau32_max_) {
  h_TopTagIndexer = ctx.get_handle<TopTagIndexer>("TopTagIndexer");
  pt_min = pt_min_;
  eta_max = eta_max_;
  nsub_min = nsub_min_;
  fpt_max = fpt_max_;
  mpair_min = mpair_min_;
  tau32_max = tau32_max_;
}

bool HotvrTopIndexCleaner::process(Event &event) {
  vector<TopJet> *topjets = event.topjets;
  TopTagIndexer indexer = event.get(h_TopTagIndexer);
  vector<int> ind = indexer.GetIndex();
  vector<int> selInd;
  for (int i : ind)
    { 
      TopJet topjet = topjets->at(i);
      vector<Jet> subjets = topjet.subjets();
      double pt = topjet.v4().pt();
      double eta = abs(topjet.v4().eta());
      unsigned int nsub = subjets.size();
      if (!( (pt_min < pt) && (eta < eta_max) && (nsub_min <= nsub) )) continue;
      double ptsub = subjets.at(0).v4().pt();
      double fpt   = ptsub/pt;
      double m12   = (subjets.at(0).v4() + subjets.at(1).v4()).M();
      double m13   = (subjets.at(0).v4() + subjets.at(2).v4()).M();
      double m23   = (subjets.at(1).v4() + subjets.at(2).v4()).M();
      double mpair = min(min(m12, m13), m23);
      double tau32 = topjet.tau3_groomed()/topjet.tau2_groomed();
      if ((fpt < fpt_max) && (mpair_min < mpair) && (tau32 < tau32_max)) selInd.push_back(i);
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
