import os
import ROOT
import numpy as np
from collections import defaultdict

def write_histogram(pname, hin, hout, fout):
    """
    get histogram from input file and write as to output file
    """
    print "write {} {} {}".format(pname, hin, hout)
    fin = ROOT.TFile(pname)
    h = fin.Get(hin)
    # rate = h.Integral()
    rate = -1
    fout.cd()
    h.Write(hout)
    fin.Close()
    return rate

def draw_channel(fout, channels, signal, backgrounds):
    """
    Load histograms from output file and draw histogram. 
    This can be used as a cross check and to easily create pre-fit distribution plots.
    """
    colors = {"TT":810, 
               "ST":797,
               "WJets":600, 
               "DYJets":834,
               "Diboson":435,
               "QCD":867,
               "Other":858}
    
    for channel in channels:
        h_signal = fout.Get(channel+"_"+signal)
        h_stack = ROOT.THStack("h_stack","")
        h_list = range(len(backgrounds))
        for i in reversed(range(0,len(backgrounds))):
            h_tmp = fout.Get(channel+"_"+backgrounds[i])
            h_list[i] = (h_tmp.Clone())
            h_list[i].SetFillColor(colors[backgrounds[i]])
            h_list[i].SetFillStyle(1001)
            h_stack.Add(h_list[i])
            print i
            print backgrounds[i]
        c = ROOT.TCanvas("c","c",600,600)
        ROOT.gPad.SetLogy(True)
        h_stack.Draw("hist")
        h_signal.Draw("SAME")
        c.Print(signal+"_"+channel+".eps")

if __name__ == "__main__":
    """
    create Combine datacards
    """

    seperator = "-"*64
    ##########################################
    #                Settings                #
    ##########################################

    file_prefix = "uhh2.AnalysisModuleRunner."
    # channels = {"muon_region_1b1t":"1btag1toptag20chi2_reco/Bstar_reco_M_rebin", "muon_region_2b1t":"2btag1toptag_reco/Bstar_reco_M_rebin"}
    channels = {"muon_region_2b1t":"2btag1toptag_reco/Bstar_reco_M_rebin",
                "electron_region_2b1t":"2btag1toptag_reco/Bstar_reco_M_rebin"}
    signals = ["BstarToTW0700LH", "BstarToTW0800LH", "BstarToTW0900LH", "BstarToTW1000LH", "BstarToTW1100LH",
               "BstarToTW1200LH", "BstarToTW1400LH", "BstarToTW1600LH", "BstarToTW1800LH", "BstarToTW2000LH",
               "BstarToTW2200LH", "BstarToTW2400LH", "BstarToTW2600LH", "BstarToTW2800LH", "BstarToTW3000LH",
               "BstarToTW3200LH", "BstarToTW3400LH", "BstarToTW3600LH", "BstarToTW3800LH", "BstarToTW4000LH", "BstarToTW4200LH"]
    # signals = ["BstarToTW0700RH", "BstarToTW0800RH", "BstarToTW0900RH", "BstarToTW1000RH", "BstarToTW1100RH",
    #            "BstarToTW1200RH", "BstarToTW1400RH", "BstarToTW1600RH", "BstarToTW1800RH", "BstarToTW2000RH",
    #            "BstarToTW2200RH", "BstarToTW2400RH", "BstarToTW2600RH", "BstarToTW2800RH", "BstarToTW3000RH",
    #            "BstarToTW3200RH", "BstarToTW3400RH", "BstarToTW3600RH", "BstarToTW3800RH", "BstarToTW4000RH", "BstarToTW4200RH"]
    # signals = ["BstarToTW0700VL", "BstarToTW0800VL", "BstarToTW0900VL", "BstarToTW1000VL", "BstarToTW1100VL",
    #            "BstarToTW1200VL", "BstarToTW1400VL", "BstarToTW1600VL", "BstarToTW1800VL", "BstarToTW2000VL",
    #            "BstarToTW2200VL", "BstarToTW2400VL", "BstarToTW2600VL", "BstarToTW2800VL", "BstarToTW3000VL",
    #            "BstarToTW3200VL", "BstarToTW3400VL", "BstarToTW3600VL", "BstarToTW3800VL", "BstarToTW4000VL", "BstarToTW4200VL"]

    # signals = ["BprimeTToTW0700RH", "BprimeTToTW0800RH", "BprimeTToTW0900RH", "BprimeTToTW1000RH", "BprimeTToTW1100RH","BprimeTToTW1200RH",
    #            "BprimeTToTW1300RH", "BprimeTToTW1400RH", "BprimeTToTW1500RH", "BprimeTToTW1600RH", "BprimeTToTW1700RH","BprimeTToTW1800RH"]
    # signals = ["BprimeTToTW0700LH", "BprimeTToTW0800LH", "BprimeTToTW0900LH", "BprimeTToTW1000LH", "BprimeTToTW1100LH","BprimeTToTW1200LH",
    #            "BprimeTToTW1300LH", "BprimeTToTW1400LH", "BprimeTToTW1500LH", "BprimeTToTW1600LH", "BprimeTToTW1700LH","BprimeTToTW1800LH"]
    backgrounds = ["TT", "ST", "Other"]
    # backgrounds = ["TT", "ST", "WJets", "DYJets", "Diboson", "QCD"]
    # uncertainties
    lumierr = 1.025
    mc_norm = {"TT":1.056,
               "ST":1.1,
               "WJets":1.1, 
               "DYJets":1.1,
               "Diboson":1.2,
               "QCD":2.0}

    # shorthand for systematics shared by regions
    all_syst = []
    shape_systematics = {}
    common_systematics = [
        'PDF',
        'PU',
        # 'BTAG_bc',
        # 'BTAG_udsg',
        'BTAG_bc_2016',
        'BTAG_bc_2017',
        'BTAG_bc_2018',
        'BTAG_udsg_2016',
        'BTAG_udsg_2017',
        'BTAG_udsg_2018',
        # 'TOPTAG',
        'TOPTAG_2016',
        'TOPTAG_2017',
        'TOPTAG_2018',
        'Prefiring',
        # 'JECs',
        # 'JERs'
        'JECs_2016',
        'JECs_2017',
        'JECs_2018',
        'JERs_2016',
        'JERs_2017',
        'JERs_2018'
        ]
    muon_systematics = [
        # 'MUOID',
        # 'MUOTR',
        # 'MUOISO',
        'MUOID_2016',
        'MUOID_2017',
        'MUOID_2018',
        'MUOISO_2016',
        'MUOISO_2017',
        'MUOISO_2018',
        'MUOTR_2016',
        'MUOTR_2017',
        'MUOTR_2018'
        ]

    electron_systematics = [
        # 'ELEID',
        # 'ELEREC',
        # 'ELETR',
        'ELEID_2016',
        'ELEID_2017',
        'ELEID_2018',
        'ELEREC_2016',
        'ELEREC_2017',
        'ELEREC_2018',
        'ELETR_2016',
        'ELETR_2017',
        'ELETR_2018',
        ]

    all_syst = common_systematics + muon_systematics + electron_systematics + ["BGEST"] + ["TOPPT_alpha", "TOPPT_beta"] + ["SCALE_TT", "SCALE_ST", "SCALE_SIGNAL"]

    shape_systematics = {
        "muon_region_1b1t":{
            "TT": common_systematics + muon_systematics + ["TOPPT_alpha", "TOPPT_beta", "SCALE_TT"],
            "ST": common_systematics + muon_systematics + ["SCALE_ST"],
            "signal": common_systematics + muon_systematics + ["SCALE_SIGNAL"],
            "Other": ["BGEST"]
            },
        "muon_region_2b1t":{
            "TT": common_systematics + muon_systematics + ["TOPPT_alpha", "TOPPT_beta", "SCALE_TT"],
            "ST": common_systematics + muon_systematics + ["SCALE_ST"],
            "signal": common_systematics + muon_systematics + ["SCALE_SIGNAL"],
            "Other": ["BGEST"]
            },
        "electron_region_1b1t":{
            "TT": common_systematics + electron_systematics + ["TOPPT_alpha", "TOPPT_beta", "SCALE_TT"],
            "ST": common_systematics + electron_systematics + ["SCALE_ST"],
            "signal": common_systematics + electron_systematics + ["SCALE_SIGNAL"],
            "Other": ["BGEST"]
            },
        "electron_region_2b1t":{
            "TT": common_systematics + electron_systematics + ["TOPPT_alpha", "TOPPT_beta", "SCALE_TT"],
            "ST": common_systematics + electron_systematics + ["SCALE_ST"],
            "signal": common_systematics + electron_systematics + ["SCALE_SIGNAL"],
            "Other": ["BGEST"]
            }
        }
    
    ##########################################

    # loop over all signal processes and create separate root files
    # open output file
    fout_name = "BstarToTW_combined.root"
    fout = ROOT.TFile(fout_name, "RECREATE")
    observed = {}
    rates = defaultdict(dict)
    
    for lepton in ['Muon','Electron']:
        analysis_dir = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/%s/all/"%lepton
        for channel, hname in channels.iteritems():    
            if lepton.lower() not in channel: continue
            # get nominal histograms and rates
            # signal
            for signal in signals:
                pname = analysis_dir + "NOMINAL/" + file_prefix + "MC." + signal + ".root"
                rates[channel][signal] = write_histogram(pname, hname, channel+"__"+signal, fout)
                for syst in shape_systematics[channel]["signal"]:
                    pname = analysis_dir + syst+"_up/" + file_prefix + "MC." + signal + ".root"
                    write_histogram(pname, hname, channel+"__"+signal+"_"+syst+"Up", fout)
                    pname = analysis_dir + syst+"_down/" + file_prefix + "MC." + signal + ".root"
                    write_histogram(pname, hname, channel+"__"+signal+"_"+syst+"Down", fout)
            # background
            for background in backgrounds: 
                pname = analysis_dir + "NOMINAL/" + file_prefix + "MC." + background + ".root"
                rates[channel][background] = write_histogram(pname, hname, channel+"__"+background, fout)
                for syst in shape_systematics[channel][background]:
                    pname = analysis_dir + syst+"_up/" + file_prefix + "MC." + background + ".root"
                    write_histogram(pname, hname, channel+"__"+background+"_"+syst+"Up", fout)
                    pname = analysis_dir + syst+"_down/" + file_prefix + "MC." + background + ".root"
                    write_histogram(pname, hname, channel+"__"+background+"_"+syst+"Down", fout)

            #  data
            pname = analysis_dir + "NOMINAL/" + file_prefix + "DATA.DATA.root"
            observed[channel] = write_histogram(pname, hname, channel+"__data_obs", fout)

        # draw_channel(fout,channels,signal,backgrounds);
    fout.Close()
    print rates
        # create datacard
    for signal in signals:
        fcard_name = "datacard_" + signal + ".txt"
        with open(fcard_name, "w") as fcard:
            fcard.write("# automatically generated datacard for {}\n".format(signal))
            fcard.write(seperator+"\n")
            fcard.write("imax {}\n".format(len(channels)))
            fcard.write("jmax {}\n".format(len(backgrounds)))
            fcard.write("kmax *\n")
            fcard.write(seperator+"\n")
            fcard.write("shapes * * {} $CHANNEL__$PROCESS $CHANNEL__$PROCESS_$SYSTEMATIC \n".format(fout_name))
            fcard.write(seperator+"\n")
            fcard.write("#{:^62}#\n".format('CHANNELS'))
            fcard.write(seperator+"\n")
            fcard.write("bin {}\n".format(" ".join([channel for channel in channels])))
            fcard.write("observation {}\n".format(" ".join([str(observed[channel]) for channel in channels])))
            fcard.write(seperator+"\n")
            fcard.write("#{:^62}#\n".format('PROCESSES'))
            fcard.write(seperator+"\n")
            fcard.write("bin {}\n".format(" ".join([(" "+channel)*(len(backgrounds) +1) for channel in channels])))
            fcard.write("process {}\n".format((" "+signal+" "+" ".join([background for background in backgrounds]))*len(channels)))
            fcard.write("process {}\n".format(" ".join([str(x) for x in range(0,len(backgrounds) + 1) * len(channels)])))
            fcard.write("rate {}\n".format(" ".join([str(rates[channel][signal]) +" " + " ".join([str(rates[channel][background]) for background in backgrounds]) for channel in channels])))
            fcard.write(seperator+"\n")
            fcard.write("#{:^62}#\n".format('SYSTEMATICS'))
            fcard.write(seperator+"\n")
            fcard.write("lumi lnN {}\n".format(" ".join([" - " if process == "Other" else " "+str(lumierr)+" " for channel in channels for process in (["signal"] + backgrounds) ])))
            for background in backgrounds:
                if not background in mc_norm.keys(): continue
                fcard.write("{}mc lnN {}\n".format(background, (" - "+" -"*backgrounds.index(background) + " "+str(mc_norm[background]) + " -"*(len(backgrounds) - backgrounds.index(background) - 1))*len(channels)))

            for syst in all_syst:
                fcard.write("{} shape {}\n".format(syst, " ".join([" 1 " if syst in shape_systematics[channel][process] else " - " for channel in channels for process in (["signal"] + backgrounds) ])))
            fcard.write("* autoMCStats 0 0 1")
