#include "UHH2/BstarToTW/include/CutflowHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

CutflowHists::CutflowHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  cutflow = book<TH1F>("Cutflow", "Cutflow", 10, 0, 10);

}


void CutflowHists::fill(const Event & event){

}

void CutflowHists::fill(const Event & event, int i){
  double weight = event.weight;
  cutflow->Fill(i, weight);
}


CutflowHists::~CutflowHists(){}
