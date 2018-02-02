#include "UHH2/BstarToTW/include/TopTagPerformanceHists.h"

#include "TH1F.h"
#include "TH2F.h"

using namespace std;
using namespace uhh2;


TopTagPerformanceHists::TopTagPerformanceHists(Context & ctx, const string & dirname, bool is_qcd, TopJetId id , const boost::optional<Event::Handle<std::vector<TopJet> > > &topjetcollection): 
  Hists(ctx, dirname),
  b_is_qcd(is_qcd),
  m_id(id),
  h_topjetcollection(topjetcollection) 
{
  double intervall[] = {0, 200, 400, 600, 800, 1000, 1500, 2000};
  N_evt = book<TH1F>("N", "N", 1, 0.5, 1.5);
  N_mat = book<TH1F>("N_mat", "N", 1, 0.5, 1.5);
  pt_top = book<TH1F>("pt_gentop", "p_{T, topjet}", 7, intervall);
  pt_matched = book<TH1F>("pt_top_matched", "p_{T, topjet}", 7, intervall);
  pt_mismatched = book<TH1F>("pt_top_mismatched", "p_{T, topjet}", 7, intervall);
  pt_gen = book<TH1F>("pt_gen", "p_{T, gen}", 7, intervall);
  dR_max = book<TH2F>("dR_max", "#Delta R_{max}", 7, intervall, 60, 0, 3);
}

void TopTagPerformanceHists::fill(const Event & event){  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // double weight = event.weight;
  // N_evt->Fill(1., weight);

  // const vector<TopJet> & topjets = h_topjetcollection ? event.get(*h_topjetcollection) : *event.topjets;
  // const vector<GenParticle> genparticles = *event.genparticles;
  // if (m_id)
  //   {
  //     for (TopJet topjet : topjets)
  // 	{
  // 	  if ((*m_id)(topjet, event)) 
  // 	    {
  // 	      double pt = 0;
  // 	      GenParticle matched_gen;
  // 	      for (GenParticle gen : genparticles)
  // 		{
  // 		  if (deltaR(gen.v4(), topjet.v4()) < 0.2 && gen.pt() > pt)
  // 		    {
  // 		      pt = gen.pt();
  // 		      matched_gen = gen;
  // 		    }
  // 		}
  // 	      if (abs(matched_gen.pdgId()) == 6 )	  
  // 		{
  // 		  bool top_had = false;
  // 		  GenParticle gend1 = *matched_gen.daughter(&genparticles,1);
  // 		  GenParticle gend2 = *matched_gen.daughter(&genparticles,2);
		  
  // 		  if (abs(gend1.pdgId()) == 5)
  // 		    {
  // 		      const GenParticle *gendd = gend2.daughter(&genparticles,1);
  // 		      if (gendd  && abs(gendd->pdgId()) < 6) top_had=true;
  // 		    }
  // 		  else if (abs(gend2.pdgId()) == 5)
  // 		    {
  // 		      const GenParticle *gendd = gend1.daughter(&genparticles,1);
  // 		      if (gendd && abs(gendd->pdgId()) < 6) top_had=true;
  // 		    }
  // 		  if (top_had)
  // 		    {
  // 		      if (pt < 2000) pt_matched->Fill(pt, weight);
  // 		      else pt_matched->Fill(1001., weight);
  // 		      N_mat->Fill(1., weight);
  // 		    }
  // 		}
  // 	      else 
  // 		{
  // 		  if (pt < 2000) pt_mismatched->Fill(pt, weight);
  // 		  else pt_mismatched->Fill(1001., weight);
  // 		}
  // 	    }
  // 	}
  //   }
  // else
  //   {
  //     for (GenParticle gen : genparticles)
  // 	{
  // 	  if (gen.status() == 23) pt_gen->Fill(gen.pt(), weight);		  
  // 	  if (abs(gen.pdgId()) == 6)
  // 	    {
  // 	      bool top_had = false;
  // 	      GenParticle gend1 = *gen.daughter(&genparticles,1);
  // 	      GenParticle gend2 = *gen.daughter(&genparticles,2);
		  
  // 	      if (abs(gend1.pdgId()) == 5)
  // 		{
  // 		  const GenParticle *gendd = gend2.daughter(&genparticles,1);
  // 		  if (gendd  && abs(gendd->pdgId()) < 6) top_had=true;
  // 		}
  // 	      else if (abs(gend2.pdgId()) == 5)
  // 		{
  // 		  const GenParticle *gendd = gend1.daughter(&genparticles,1);
  // 		  if (gendd && abs(gendd->pdgId()) < 6) top_had=true;
  // 		}
  // 	      if (top_had) pt_top->Fill(gen.pt(), weight);
  // 	    }
  // 	}
  //   }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  double weight = event.weight;
  const vector<TopJet> & topjets = h_topjetcollection ? event.get(*h_topjetcollection) : *event.topjets;

  double matching_dist = 0.4;
  
  if (b_is_qcd)
    {
      const vector<Particle> genjets = *event.genjets;
      for (Particle gen : genjets)
	{
	  pt_gen->Fill(gen.pt(), weight);
	}
      for (TopJet topjet : topjets)
	{
	  double pt = -1;
	  for (Particle gen : genjets)
	    {
	      bool is_tagged = m_id(topjet, event);
	      if (is_tagged) 
		{	  
		  if (deltaR(gen.v4(), topjet.v4()) < matching_dist && gen.pt() > pt)
		    {
		      pt = gen.pt();
		    }
		}
	    }
	  if (pt < 2000) pt_mismatched->Fill(pt, weight);
	  else pt_mismatched->Fill(1999., weight);
	}
    }
  else
    {
      const vector<GenParticle> genparticles = *event.genparticles;

      GenParticle top;
      bool top_had = false;
      for (GenParticle gen : genparticles)
	{		  
	  if (abs(gen.pdgId()) == 6)
	    {
	      GenParticle gend1 = *gen.daughter(&genparticles,1);
	      GenParticle gend2 = *gen.daughter(&genparticles,2);
	      double deltaR_max = -1;
		  
	      if (abs(gend1.pdgId()) == 5)
		{
		  const GenParticle *gendd = gend2.daughter(&genparticles,1);
		  if (gendd  && abs(gendd->pdgId()) <= 6)
		    {
		      top_had=true;
		      deltaR_max =  max( max(deltaR(gend1.v4(), gendd->v4()), deltaR(gend1.v4(), gend2.daughter(&genparticles,2)->v4())), deltaR(gendd->v4(), gend2.daughter(&genparticles,2)->v4()));
		    }
		}
	      else if (abs(gend2.pdgId()) == 5)
		{
		  const GenParticle *gendd = gend1.daughter(&genparticles,1);
		  if (gendd && abs(gendd->pdgId()) <= 6) 
		    {
		      top_had=true;
		      deltaR_max =  max( max(deltaR(gend2.v4(), gendd->v4()), deltaR(gend2.v4(), gend1.daughter(&genparticles,2)->v4())), deltaR(gendd->v4(), gend1.daughter(&genparticles,2)->v4()));
		    }
		}
	      if (top_had) 
		{
		  top = gen;
		  pt_gen->Fill(gen.pt(), weight);
		  dR_max->Fill(gen.pt(), deltaR_max, weight);
		}
	    }
	}
      if (!top_had) return;
      for (TopJet topjet : topjets)
	{
	  if (m_id(topjet, event)) 
	    {
	      
	      if (deltaR(top.v4(), topjet.v4()) < matching_dist)
		{
		  if (top.pt() < 2000) pt_matched->Fill(top.pt(), weight);
		  else pt_matched->Fill(1999., weight);
		}
	    }
	}
    }
}


TopTagPerformanceHists::~TopTagPerformanceHists() {}

