import ROOT
import AnalysisPlotter.Plotter as Plotter

class LimitPlotter(Plotter):
    """
    Plotter to create nice limit plots
    """
    def __init__(self):
        super().__init__()    

    def draw_limits(masspoints,limits, predicted=None):
        """
        make nice limit plots
        """
        print "start drawing..."
        canvas_text = []
        maximum = 10
        minimum = 0.001
        masspoints = np.array(masspoints, dtype=float)
        central = np.array(limits["central"], dtype=float)

        low_95 = np.array(limits["low_95"], dtype=float)
        low_68 = np.array(limits["low_68"], dtype=float)
        high_68 = np.array(limits["high_68"], dtype=float)
        high_95 = np.array(limits["high_95"], dtype=float)
        # observed = np.array(limits["observed"], dtype=float)
        zero = np.zeros(len(masspoints), dtype=float)

        g_expected = ROOT.TGraph(len(masspoints), masspoints, central)
        g_expected_68 = ROOT.TGraphAsymmErrors(len(masspoints), masspoints, central, zero, zero, low_68, high_68)
        g_expected_95 = ROOT.TGraphAsymmErrors(len(masspoints), masspoints, central, zero, zero, low_95, high_95)
        # g_observed = ROOT.TGraph(len(masspoints), masspoints, observed)
        # g_observed.Draw("SAME")
        
        g_expected.SetLineWidth(2)
        g_expected.SetLineStyle(1)
        g_expected.SetLineColor(ROOT.kBlack)
        g_expected_68.SetFillStyle(1001)
        g_expected_68.SetFillColor(ROOT.kGreen + 1) # recommended color
        g_expected_95.SetFillStyle(1001)
        g_expected_95.SetFillColor(ROOT.kOrange) # recommended color
        g_expected_95.SetMaximum(maximum*1.1)
        g_expected_95.SetMinimum(minimum*0.33)
        g_expected_95.SetTitle("")

        c = self.get_canvas()
        ROOT.gPad.SetLogy()
        ROOT.gPad.SetTickx()
        ROOT.gPad.SetTicky()
    
        g_expected_95.Draw("A3")
        g_expected_68.Draw("SAME3")
        g_expected.Draw("SAME") 
        
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
        if predicted is not None:
            predicted = np.array(predicted, dtype=float)
            g_predicted = ROOT.TGraph(len(masspoints), masspoints, predicted)
            g_predicted.SetLineWidth(2)
            g_predicted.SetLineStyle(9)
            g_predicted.SetLineColor(ROOT.kBlack)
            g_predicted.Draw("SAME")
            leg.AddEntry(g_predicted, "predicted","L")
        leg.Draw()
    
        h = g_expected_95.GetHistogram()
        h.GetXaxis().SetRangeUser(masspoints[0], masspoints[-1])
        h.SetXTitle("M_{b*} [GeV]")
        h.SetYTitle("#sigma(bg#rightarrow b* #rightarrow tW) [pb]")
        h.GetYaxis().CenterTitle()
        h.GetYaxis().SetTitleOffset(1.25)
        h.Draw("AXIS SAME")

        canvas_text+=self.draw_cmstext(extratext="Preliminary")
        canvas+text+=self.draw_infotext("137.2 fb^{-1} (13 TeV)")
        self.save_plot(c)
