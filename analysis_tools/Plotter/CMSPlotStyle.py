import ROOT

""" Global Variables for Styling """
top_margin = 0.07
bottom_margin = 0.12
left_margin = 0.17
right_margin = 0.06

width = 600
height = 600

extratext = "Preliminary"

def get_style():
    """
    Get plot style used in CMS
    uses style guidline defined in 
    https://twiki.cern.ch/twiki/bin/viewauth/CMS/Internal/FigGuidelines
    especially definitions from 
    https://twiki.cern.ch/twiki/pub/CMS/Internal/FigGuidelines/tdrstyle.C
    """
    
    cmsstyle = ROOT.TStyle('CMS_Style', 'Style for CMS publishing')

    # canvas
    cmsstyle.SetCanvasBorderMode(0)
    cmsstyle.SetCanvasColor(0)
    cmsstyle.SetCanvasDefH(height)
    cmsstyle.SetCanvasDefW(width)
    cmsstyle.SetCanvasDefX(0)
    cmsstyle.SetCanvasDefY(0)
    
    # pad
    cmsstyle.SetPadBorderMode(0)
    cmsstyle.SetPadColor(0)
    cmsstyle.SetPadGridX(False)
    cmsstyle.SetPadGridY(False)
    cmsstyle.SetGridColor(0)
    cmsstyle.SetGridStyle(3)
    cmsstyle.SetGridWidth(1)
    
    # frame
    cmsstyle.SetFrameBorderMode(0)
    cmsstyle.SetFrameBorderSize(1)
    cmsstyle.SetFrameFillColor(0)
    cmsstyle.SetFrameFillStyle(0)
    cmsstyle.SetFrameLineColor(1)
    cmsstyle.SetFrameLineStyle(1)
    cmsstyle.SetFrameLineWidth(1)

    # stat box
    cmsstyle.SetOptFile(0)
    cmsstyle.SetOptStat(0)
    
    # margins
    cmsstyle.SetPadTopMargin(top_margin)
    cmsstyle.SetPadBottomMargin(bottom_margin)
    cmsstyle.SetPadLeftMargin(left_margin)
    cmsstyle.SetPadRightMargin(right_margin)

    # global title
    cmsstyle.SetOptTitle(0)
    cmsstyle.SetTitleFont(42)
    cmsstyle.SetTitleColor(1)
    cmsstyle.SetTitleTextColor(1)
    cmsstyle.SetTitleFillColor(10)
    cmsstyle.SetTitleFontSize(0.05)
    
    # axis
    cmsstyle.SetAxisColor(1, "XYZ")
    cmsstyle.SetStripDecimals(True)
    cmsstyle.SetTickLength(0.03, "XYZ")
    cmsstyle.SetNdivisions(510, "XYZ")
    cmsstyle.SetPadTickX(1)
    cmsstyle.SetPadTickY(1)

    # axis title
    cmsstyle.SetTitleColor(1, "XYZ")
    cmsstyle.SetTitleFont(42, "XYZ") # Set relative font size 
    cmsstyle.SetTitleSize(0.06, "XYZ")
    cmsstyle.SetTitleXOffset(0.95)
    cmsstyle.SetTitleYOffset(1.4)

    # axis labels
    cmsstyle.SetLabelColor(1, "XYZ")
    cmsstyle.SetLabelFont(42, "XYZ") # Set relative font size 
    cmsstyle.SetLabelOffset(0.007, "XYZ")
    cmsstyle.SetLabelSize(0.05, "XYZ")

    # postscript
    cmsstyle.SetPaperSize(20.,20.)
    cmsstyle.SetHatchesLineWidth(5)
    cmsstyle.SetHatchesSpacing(0.05)
    return cmsstyle
 
