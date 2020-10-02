import ROOT

if __name__ == "__main__":

    base_dir = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/PreSelection/Electron/"
    # base_dir = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/PreSelection/Muon/"

    years = {"2016":"_2016v3", "2017":"_2017v2", "2018": "_2018"}

    # base_dir = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/Electron/"
    # base_dir = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/Muon/"

    # years = {"all":""}

    file_prefix = "uhh2.AnalysisModuleRunner."

    cuts = ["Cleaning", "1LepCut", "1JetCut", "100HTCut", "50METCut", "400STCut", "deltaR_jetlep", "deltaPhiMETLep", "1TopJetCut"]
    # cuts = ["PreSel", "Trigger", "BadHCAL", "1btag", "1btagdr", "1btag1toptag","1btag1toptag20chi2", "2btag", "2btag1toptag"]
    n_cuts = len(cuts)

    processes = ["TT", "ST", "DATA", "BstarToTW1000LH", "BstarToTW2000LH", "BstarToTW3000LH"]

    hist_name = "_Counter/NEvt"

    hists = {}
    
    for process in processes:
        hists[process] = ROOT.TH1F("cutflow_"+process, "cutflow_"+process, n_cuts, 0, n_cuts)
        for year, year_postfix in years.iteritems():
            process_prefix = "MC."
            if process == "DATA": process_prefix = "DATA."
            file_path = base_dir + year + "/NOMINAL/" + file_prefix + process_prefix + process + year_postfix + ".root"
            fin = ROOT.TFile(file_path)
            ind = 1 
            for cut in cuts:
                hists[process].Fill(ind - 0.5,fin.Get(cut+hist_name).Integral())
                hists[process].GetXaxis().SetBinLabel(ind, cut)
                ind += 1
            fin.Close()
        hists[process].Scale(1/hists[process].GetBinContent(1))
        


    c = ROOT.TCanvas("c","c", 1200, 600)
    hists["DATA"].Draw("HIST")
    hists["TT"].SetLineColor(632)
    hists["TT"].Draw("HIST SAME")    
    hists["ST"].SetLineColor(800)
    hists["ST"].Draw("HIST SAME")
    hists["BstarToTW1000LH"].SetLineColor(921)
    hists["BstarToTW1000LH"].Draw("HIST SAME")
    hists["BstarToTW2000LH"].SetLineColor(922)
    hists["BstarToTW2000LH"].Draw("HIST SAME")
    hists["BstarToTW3000LH"].SetLineColor(923)
    hists["BstarToTW3000LH"].Draw("HIST SAME")
    leg = ROOT.TLegend(0.6,0.7,0.9,0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(hists["DATA"], "data","l")
    leg.AddEntry(hists["TT"], "tt", "l")
    leg.AddEntry(hists["ST"], "single t", "l")
    leg.AddEntry(hists["BstarToTW1000LH"], "1TeV signal", "l")
    leg.AddEntry(hists["BstarToTW2000LH"], "2TeV signal", "l")
    leg.AddEntry(hists["BstarToTW3000LH"], "3TeV signal", "l")
    leg.Draw()
    c.Print("cutflow_selection_lin.eps")
    ROOT.gPad.SetLogy(True)
    c.Print("cutflow_selection_log.eps")

