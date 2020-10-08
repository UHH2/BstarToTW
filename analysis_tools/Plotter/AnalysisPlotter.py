import CMSPlotStyle
import ROOT
import yaml

class Plotter(object):
    """
    Generic Plotter to create a variety of Plots to be used in AN / Paper / Thesis
    """
    def __init__(self, config_file):
        """
        initialize
        """
        self.style = CMSPlotStyle.get_style()
        self.config = {}
        with open(config_file, 'r') as stream:
            try:
                self.config = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
        
    def save_plot(self,c):
        if "b_save_eps" in config and config["b_save_eps"]:
            c.Print(c.GetName()+".eps")
        if "b_save_pdf" in config and config["b_save_pdf"]:
            c.Print(c.GetName()+".pdf")
                
    def get_canvas(self,name, title= ""):
         """
         create a new TCavas
         """
         canvas = ROOT.TCanvas(name, title)
         return canvas

    def get_pad(self):
        y1 = 0.1 
        y2 = 0.99
        x1 = 0.01
        x2 = 0.99
        pad = ROOT.TPad("pad", "", x1, y1, x2, y3)
        return pad

    def draw_cmstext(self, position="default",extratext=""):
        """
        Draws "CMS" on current active canvas.
        Additionally extra text (e.g. "Preliminary") can be drawn by 
        Use position to choose from a number of predefined positions for the text.
        """
        
        x = left_margin + 0.04
        y = 1 - top_margin - top_margin * 0.5
        align =13

        if "default" in position:
            x = left_margin + 0.04
            y = 1 - top_margin - top_margin * 0.5
            align =13
        elif "right" in position:
            x = 1 - right_margin - 0.04
            y = 1 - top_margin - top_margin * 0.5
            align =33
    
        cms_text = ROOT.TLatex(3.5, 24, "CMS")
        cms_text.SetNDC()
        cms_text.SetTextAlign(align)
        cms_text.SetX(x)
        cms_text.SetY(y)
        cms_text.SetTextFont(62)
        cms_text.SetTextSize(0.75*top_margin)
        cms_text.Draw()

        prelim_text = ROOT.TLatex(3.5, 24, extratext);
        prelim_text.SetNDC()
        prelim_text.SetTextAlign(align)
        prelim_text.SetX(x)
        prelim_text.SetY(y-top_margin*0.7)
        prelim_text.SetTextFont(52)
        prelim_text.SetTextSize(0.76*text1.GetTextSize())
        prelim_text.Draw()
        
        return [cms_text,prelim_text]

    def draw_info(self, infotext):
        """
        Draws infotext (e.g. luminosity / sqrt(s) ) on the top right corner of the plot.
        """
        
        x = 1 - right_margin - 0.01
        y = 1 - top_margin + 0.01
        
        text = ROOT.TLatex(3.5,24, infotext)
        text.SetNDC()
        text.SetTextAlign(31)
        text.SetX(x)
        text.SetY(y)
        text.SetTextFont(42)
        text.SetTextSize(0.6*top_margin)
        text.Draw()
        
        return [text]

        
