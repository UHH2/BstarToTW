import os.path
import ROOT

b_debug = False
b_warnings = True
def get_histogram(path, hin_name, hout_name):
    if b_debug:
        print "get {} from {} as {}".format(hin_name, path, hout_name)
    fin = ROOT.TFile(path, "read")
    h = fin.Get(hin_name).Clone(hout_name)
    h.SetDirectory(0)
    return h

def write_histogram(fout, hout):
    """
    write histogram to TFile
    """
    if b_debug:
        print "write {} to {}".format(hout.GetName(), fout.GetName())
    fout.cd()
    hout.Write()

def create_template_file(analysis_path, channel, year, processes, variations, regions, workspace_name):
    """
    collect histograms for one analysis channel and write them to a file
    """    

    # helper variables
    file_prefix = "uhh2.AnalysisModuleRunner."
    variation_postfixes = []
    name_postfix = ""
    process_prefix = ""
    process_postfix = ""
    year_postfixes = {"2016":"_2016v3","2017":"_2017v2","2018":"_2018"}

    # Create a new TFile
    file_name = workspace_name + "_" + channel + "_" + year + ".root"
    if b_debug:
        print "create template file {}".format(file_name)
    template_file = ROOT.TFile(file_name, "recreate")

    for process in processes: 
        # set file prefix for data or MC
        if process == "DATA": process_prefix = "DATA."
        else: process_prefix = "MC." 
        # add year specifit file postfix to signal samples
        if "BstarToTW" in process: process_postfix = year_postfixes[year]
        else: process_postfix = ""
        for variation, b_correlated in variations.iteritems():
            for region, hist_name in regions.iteritems():
                work_path = analysis_path + "/" + channel + "/" + year + "/"
                if process == "DATA" and variation != "NOMINAL": continue
                if variation == "NOMINAL": variation_postfixes = [""]
                else: variation_postfixes = ["_up","_down"]
                if variation == "BGEST": name_postfix = "_"+channel+"_"+region
                elif variation == "SCALE": name_postfix = "_"+process
                else: name_postfix = ""
                if not b_correlated : name_postfix += "_" + year

                for variation_postfix in variation_postfixes:
                    fin_path = work_path + variation + variation_postfix + "/" + file_prefix + process_prefix + process + process_postfix + ".root" 
                    if (os.path.exists(fin_path)) : 
                        if variation == "NOMINAL":
                            hout_name = channel + "_" + region + "__" + process
                        else:
                            hout_name = channel + "_" + region + "__" + process + "_" + variation + name_postfix + variation_postfix.replace("_up","Up").replace("_down","Down")
                        if process == "DATA":
                            hout_name = channel + "_" + region + "__" + "data_obs"
                        hist = get_histogram(fin_path, hist_name, hout_name)
                        write_histogram(template_file,hist)
                    elif b_debug or b_warnings:
                        print "[WARNING] {} does not exist".format(fin_path)

if __name__ == "__main__":
    """
    get histograms from root files and create combine compatible root file
    """
   
    # Setup
    workspace_name = "BstarToTW"
    analysis_path = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis"
    channels = ["Electron", "Muon"]
    years = ["2016", "2017", "2018"]
    processes = ["DATA", "TT","ST","Other",
                 "BstarToTW0700LH", "BstarToTW0800LH", "BstarToTW0900LH", "BstarToTW1000LH", "BstarToTW1100LH",
                 "BstarToTW1200LH", "BstarToTW1400LH", "BstarToTW1600LH", "BstarToTW1800LH", "BstarToTW2000LH",
                 "BstarToTW2200LH", "BstarToTW2400LH", "BstarToTW2600LH", "BstarToTW2800LH", "BstarToTW3000LH",
                 "BstarToTW3200LH", "BstarToTW3400LH", "BstarToTW3600LH", "BstarToTW3800LH", "BstarToTW4000LH", 
                 "BstarToTW4200LH"]
    variations = {"NOMINAL"     : True,
                  "PDF"         : True,
                  "PU"          : True,        
                  "Prefiring"   : True,
                  "SCALE"       : True,
                  "MUOID"       : False,
                  "MUOISO"      : False,
                  "MUOTR"       : False,
                  "ELEID"       : False,
                  "ELEREC"      : False,
                  "ELETR"       : False,
                  "JECs"        : False,
                  "JERs"        : False,
                  "BTAG_bc"     : False,
                  "BTAG_udsg"   : False,
                  "TOPTAG"      : False,
                  "TOPPT_alpha" : True,
                  "TOPPT_beta"  : True,
                  "BGEST"       : False,
                  }
    
    hist_name = "_reco/Bstar_reco_M_rebin" # path to hist in root file is region + hist_name_base
    regions = {"region_1b1t":"1btag1toptag20chi2"+hist_name, "region_2b1t":"2btag1toptag"+hist_name}

    # collect all the templates and create .root files
    for channel in channels:
        for year in years:
            create_template_file(analysis_path, channel, year, processes, variations, regions, workspace_name)
            
