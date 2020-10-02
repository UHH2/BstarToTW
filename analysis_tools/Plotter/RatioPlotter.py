import ROOT
import AnalysisPlotter.Plotter as Plotter

class RatioPlotter(Plotter):
    """
    Plotter to create nice ratio plots
    """
    def __init__(self):
        super().__init__()    
        

    def get_ratio_pads():
        yplot = 0.65
        yratio = 0.34

        y3 = 0.99 
        y2 = y3 - yplot
        y1 = y2 - yratio
        x1 = 0.01
        x2 = 0.99

        plot_pad = ROOT.TPad("pad1", "Control Plots 2", x1, y2, x2, y3)
        ratio_pad = ROOT.TPad("rp1", "Ratio2", x1, y1, x2, y2)
        return [plot_pad, ratio_pad]

    def draw_ratio(h_num, h_den):
        """
        Draws ratio plot of hists h_num and h_den.
        h_num is the numerator and h_den is the denominator
        """
        canvas_items = [] # keep everything drawn on the canvas here to avoid things going out of scope
        c = super().get_canvas("ratio_plot")
        pads = get_ratio_pads()

        # begin drawing the histograms
        pad[0].cd()
        

        # begin drawing the ratio
        pad[1].cd()
