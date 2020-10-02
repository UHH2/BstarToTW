import argparse
import sys, os
import subprocess
import signal
import time
from collections import defaultdict
import ROOT
import numpy as np

def run_combine(signal_dir, signals, max_jobs):
    """
    spawn parallel combine jobs from datacards
    """
    
    n_jobs = len(signals)
    n_running = 0
    n_completed = 0
    processes = []
    for signal in signals:
        
        b_wait = (n_running >= max_jobs)
        while b_wait:
            n_running = 0
            n_completed = 0
            for proc in processes:
                if proc[0].poll() == None: n_running += 1
                else: 
                    n_completed += 1
                    if not proc[1].closed:
                        proc[1].close()
                        print 'Job "{}" has finished.'.format(proc[1].name)
            percent = round(float(n_completed)/float(n_jobs)*100, 1)
            sys.stdout.write( '{0:d} of {1:d} ({2:4.2f}%) jobs done.\r'.format(n_completed, n_jobs, percent))
            sys.stdout.flush()
            time.sleep(10)
            b_wait = (n_running >= max_jobs)

        print 'Spawning job: {}'.format(signal)        
        n_running += 1
        f = open("log/{}.log".format(signal),'w')
        command = "nice -n {} combine {}/datacard_{}.txt -M AsymptoticLimits -n {}".format(10, signal_dir, signal, signal)
        # command = "nice -n {} combine {}/datacard_{}.txt -M HybridNew --LHCmode LHC-limits -n {}".format(10, signal_dir, signal, signal)
        processes.append((subprocess.Popen(command, stdout=f, shell=True),f))
    b_wait = (n_completed < n_jobs) 
    while b_wait:
        n_running = 0
        n_completed = 0
        for proc in processes:
            if proc[0].poll() == None: n_running += 1
            else: 
                n_completed += 1
                if not proc[1].closed:
                    proc[1].close()
                    print 'Job "%s" has finished.' % proc[1].name
        percent = float(n_completed)/float(n_jobs)*100
        sys.stdout.write( '{0:d} of {1:d} ({2:4.2f} %) jobs done.\r'.format(n_completed, n_jobs, percent))
        sys.stdout.flush()
        time.sleep(10)
        b_wait = (n_completed < n_jobs)
    print ''
    print 'Done'
    return

def get_limits(signal_dir, signals, b_blind):
    """
    shitty method to read limits from combine root output. super not optimized
    """
    limits = defaultdict(list)
    keys = ["low_95", "low_68", "central", "high_68", "high_95", "observed"]
    for signal in signals:
        # f = ROOT.TFile("higgsCombine"+signal+".AsymptoticLimits.mH120.123456.root")
        f = ROOT.TFile("higgsCombine"+signal+".AsymptoticLimits.mH120.root")
        # f = ROOT.TFile("higgsCombine"+signal+".HybridNew.mH120.root")
        t = f.Get("limit")
        values = defaultdict(list)
        n = 0
        # mod = 6 if b_blind else 5
        mod = 6
        for entry in t:
            # print entry.limit
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
    return limits

def draw_limits(masspoints,limits, predicted=None):
    """
    make nice limit plots
    """
    masspoints = np.array(masspoints, dtype=float)
    central = np.array(limits["central"], dtype=float)


    low_95 = np.array(limits["low_95"], dtype=float)
    low_68 = np.array(limits["low_68"], dtype=float)
    high_68 = np.array(limits["high_68"], dtype=float)
    high_95 = np.array(limits["high_95"], dtype=float)
    zero = np.zeros(len(masspoints), dtype=float)

    g_expected = ROOT.TGraph(len(masspoints), masspoints, central)
    g_expected_68 = ROOT.TGraphAsymmErrors(len(masspoints), masspoints, central, zero, zero, low_68, high_68)
    g_expected_95 = ROOT.TGraphAsymmErrors(len(masspoints), masspoints, central, zero, zero, low_95, high_95)

    
    g_expected.SetLineWidth(2)
    g_expected.SetLineStyle(1)
    g_expected.SetLineColor(ROOT.kBlack)

    g_expected_68.SetFillStyle(1001)
    g_expected_68.SetFillColor(ROOT.kGreen + 1) # recommended color

    g_expected_95.SetFillStyle(1001)
    g_expected_95.SetFillColor(ROOT.kOrange) # recommended color

    maximum = 10
    minimum = 0.001
    
    g_expected_95.SetMaximum(maximum*1.1)
    g_expected_95.SetMinimum(minimum*0.33)
    g_expected_95.SetTitle("")

    c = ROOT.TCanvas("c", "combine limits", 600,600)
    ROOT.gPad.SetLogy()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    
    g_expected_95.Draw("A3")
    g_expected_68.Draw("SAME3")
    g_expected.Draw("SAME") 
    
    observed = np.array(limits["observed"], dtype=float)
    g_observed = ROOT.TGraph(len(masspoints), masspoints, observed)
    g_observed.Draw("SAME")

    leg = ROOT.TLegend(0.45,0.6,0.85,0.85,"")
    leg.SetBorderSize(0)
    leg.SetTextSize(0.035)
    leg.SetFillStyle(0)
    leg.SetLineColor(1)
    leg.SetTextFont(42)
    leg.SetHeader("expected limit at 95% C.L.")
    leg.AddEntry(g_expected, "expected","L")
    leg.AddEntry(g_expected_68, "68% expected","F")
    leg.AddEntry(g_expected_95, "95% expected","F")
    # if predicted is not None:
    #     predicted = np.array(predicted, dtype=float)
    #     g_predicted = ROOT.TGraph(len(masspoints), masspoints, predicted)
    #     g_predicted.SetLineWidth(2)
    #     g_predicted.SetLineStyle(9)
    #     g_predicted.SetLineColor(ROOT.kBlack)
    #     g_predicted.Draw("SAME")
    #     leg.AddEntry(g_predicted, "predicted","L")
    leg.Draw()
    
    h = g_expected_95.GetHistogram()
    h.GetXaxis().SetRangeUser(masspoints[0], masspoints[-1])
    h.SetXTitle("M_{b*} [GeV]")
    h.SetYTitle("#sigma(bg#rightarrow b* #rightarrow tW) [pb]")
    h.GetYaxis().CenterTitle()
    h.GetYaxis().SetTitleOffset(1.25)
    h.Draw("AXIS SAME")

    c.SaveAs("Combine_Limits.eps")

if __name__ == "__main__":
    """
    run datacards for all masspoints with Combine
    """
    parser = argparse.ArgumentParser(description="run combine for datacards in given directory")
    parser.add_argument("signal_dir", help="directory of the datacards")
    parser.add_argument("-j", "--max_jobs", default=16, type=int,
                        help="maximum number of parallel jobs")
    parser.add_argument("--unblind", dest='blind', action='store_true',
                        help="unblind analysis")


    args = parser.parse_args()
    max_jobs = args.max_jobs
    signal_dir = args.signal_dir
    b_blind = args.blind

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


    masspoints = [700, 800, 900, 1000, 1100, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200]
    # masspoints = [700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800]
    predicted = [32.8, 17.06, 9.46, 5.93, 3.23, 1.98, 0.81, 0.35, 0.17, 0.081, 0.041, 0.022, 0.012, 0.0067, 0.0038, 0.001861, 0.001064, 0.0006169, 0.0003632, 0.0002179, 0.0001324]
    # predicted =  [2 *32.8, 2 *17.06, 2 *9.46, 2 *5.93, 2 *3.23, 2 *1.98, 2 *0.81, 2 *0.35, 2 *0.17, 2 *0.081, 2 *0.041, 2 *0.022, 2 *0.012, 2 *0.0067, 2 *0.0038, 2*0.001861, 2*0.001064, 2*0.0006169, 2*0.0003632, 2*0.0002179, 2*0.0001324]
    predicted = []
    run_combine(signal_dir, signals, max_jobs)
    # print b_blind
    limits = get_limits(signal_dir, signals, b_blind)
    draw_limits(masspoints, limits, predicted)
