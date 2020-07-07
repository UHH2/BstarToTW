from collections import defaultdict
import ROOT
import numpy as np

def get_limits(signal_dir, signals, b_blind):
    """
    shitty method to read limits from combine root output. super not optimized
    """
    limits = defaultdict(list)
    keys = ["low_95", "low_68", "central", "high_68", "high_95", "observed"]
    for signal in signals:
        f = ROOT.TFile(signal_dir+"/higgsCombine"+signal+".AsymptoticLimits.mH120.root")
        t = f.Get("limit")
        values = defaultdict(list)
        n = 0
        mod = 5 if b_blind else 6
        for entry in t:
            ind = n % mod
            values[keys[ind]].append(entry.limit)
            n += 1
        for i in range(0,len(values["central"])):
            values["low_95"][i] = abs(values["low_95"][i] - values["central"][i])
            values["low_68"][i] = abs(values["low_68"][i] - values["central"][i])
            values["high_68"][i] = abs(values["high_68"][i] - values["central"][i])
            values["high_95"][i] = abs(values["high_95"][i] - values["central"][i])
        for key in keys:
            if values[key]:
                limits[key].append(sum(values[key]) / len(values[key]))
    print signal_dir
    print limits
    return limits

def draw_limits(masspoints,limits,dirs):
    """
    make nice limit plots
    """

    ROOT.gStyle.SetPalette(ROOT.kBlackBody)
    masspoints = np.array(masspoints, dtype=float)
    zero = np.zeros(len(masspoints), dtype=float)
    g_expected = ROOT.TMultiGraph()
    for i in range(0,len(limits)):
        central = np.array(limits[i]["central"], dtype=float)
        g_tmp = ROOT.TGraph(len(masspoints), masspoints, central)
        g_tmp.SetLineWidth(2)
        g_tmp.SetLineStyle(1)
        g_tmp.SetTitle(dirs[i])
        g_expected.Add( g_tmp, "PL" )
    
    maximum = 10
    minimum = 0.01
    
    c = ROOT.TCanvas("c", "combine limits", 800,600)
    ROOT.gPad.SetLogy()

    g_expected.Draw("A pmc plc")

    leg = c.BuildLegend(0.35,0.65,0.85,0.92,"")
    leg.SetBorderSize(0)
    leg.SetTextSize(0.035)
    leg.SetFillStyle(0)
    leg.SetLineColor(1)
    leg.SetTextFont(42)
    leg.SetHeader("expected limit at 95% C.L.")
    leg.Draw()
    
    h = g_expected.GetHistogram()
    h.GetXaxis().SetRangeUser(masspoints[0], masspoints[-1])
    h.SetXTitle("M_{b*} [GeV]")
    h.SetYTitle("#sigma(bg#rightarrow b* #rightarrow tW) [pb]")

    h.Draw("AXIS SAME")

    c.SaveAs("Combine_Limits.eps")

if __name__ == "__main__":

    dirs = ['combined_extrapolation_vjetsnlo','combined_mconly_vjetsnlo']
    limits = []

    signals = ["BstarToTW0700LH", "BstarToTW0800LH", "BstarToTW0900LH", "BstarToTW1000LH", "BstarToTW1100LH",
               "BstarToTW1200LH", "BstarToTW1400LH", "BstarToTW1600LH", "BstarToTW1800LH", "BstarToTW2000LH",
               "BstarToTW2200LH", "BstarToTW2400LH", "BstarToTW2600LH", "BstarToTW2800LH", "BstarToTW3000LH"]

    masspoints = [700, 800, 900, 1000, 1100, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000]


    for d in dirs:
        limits.append(get_limits(d,signals,True))

    draw_limits(masspoints, limits, dirs)
