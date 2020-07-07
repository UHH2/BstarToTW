#include "TStyle.h"

// tdrGrid: Turns the grid lines on (true) or off (false)

// void cmsGrid(bool gridOn) {
//   cmsStyle->SetPadGridX(gridOn);
//   cmsStyle->SetPadGridY(gridOn);
// }

// fixOverlay: Redraws the axis

void fixOverlay() {
  gPad->RedrawAxis();
}

void setCMSStyle() {
  TStyle *cmsStyle = new TStyle("cmsStyle","Style for CMS publications");

  // For the canvas:
  cmsStyle->SetCanvasBorderMode(0);
  cmsStyle->SetCanvasColor(kWhite);
  cmsStyle->SetCanvasDefH(600); //Height of canvas
  cmsStyle->SetCanvasDefW(600); //Width of canvas
  cmsStyle->SetCanvasDefX(0);   //POsition on screen
  cmsStyle->SetCanvasDefY(0);

  // For the Pad:
  cmsStyle->SetPadBorderMode(0);
  // cmsStyle->SetPadBorderSize(Width_t size = 1);
  cmsStyle->SetPadColor(kWhite);
  cmsStyle->SetPadGridX(false);
  cmsStyle->SetPadGridY(false);
  cmsStyle->SetGridColor(0);
  cmsStyle->SetGridStyle(3);
  cmsStyle->SetGridWidth(1);

  // For the frame:
  cmsStyle->SetFrameBorderMode(0);
  cmsStyle->SetFrameBorderSize(1);
  cmsStyle->SetFrameFillColor(0);
  cmsStyle->SetFrameFillStyle(0);
  cmsStyle->SetFrameLineColor(1);
  cmsStyle->SetFrameLineStyle(1);
  cmsStyle->SetFrameLineWidth(1);
  
  // For the histo:
  // cmsStyle->SetHistFillColor(1);
  // cmsStyle->SetHistFillStyle(0);
  cmsStyle->SetHistLineColor(1);
  cmsStyle->SetHistLineStyle(0);
  cmsStyle->SetHistLineWidth(1);
  // cmsStyle->SetLegoInnerR(Float_t rad = 0.5);
  // cmsStyle->SetNumberContours(Int_t number = 20);

  cmsStyle->SetEndErrorSize(2);
  // cmsStyle->SetErrorMarker(20);
  //cmsStyle->SetErrorX(0.);
  
  cmsStyle->SetMarkerStyle(20);
  
  //For the fit/function:
  cmsStyle->SetOptFit(1);
  cmsStyle->SetFitFormat("5.4g");
  cmsStyle->SetFuncColor(2);
  cmsStyle->SetFuncStyle(1);
  cmsStyle->SetFuncWidth(1);

  //For the date:
  cmsStyle->SetOptDate(0);
  // cmsStyle->SetDateX(Float_t x = 0.01);
  // cmsStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  cmsStyle->SetOptFile(0);
  cmsStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  cmsStyle->SetStatColor(kWhite);
  cmsStyle->SetStatFont(42);
  cmsStyle->SetStatFontSize(0.025);
  cmsStyle->SetStatTextColor(1);
  cmsStyle->SetStatFormat("6.4g");
  cmsStyle->SetStatBorderSize(1);
  cmsStyle->SetStatH(0.1);
  cmsStyle->SetStatW(0.15);
  // cmsStyle->SetStatStyle(Style_t style = 1001);
  // cmsStyle->SetStatX(Float_t x = 0);
  // cmsStyle->SetStatY(Float_t y = 0);

  // Margins:
  cmsStyle->SetPadTopMargin(0.07);
  cmsStyle->SetPadBottomMargin(0.12);
  cmsStyle->SetPadLeftMargin(0.13);
  cmsStyle->SetPadRightMargin(0.06);

  // For the Global title:

  cmsStyle->SetOptTitle(0);
  cmsStyle->SetTitleFont(42);
  cmsStyle->SetTitleColor(1);
  cmsStyle->SetTitleTextColor(1);
  cmsStyle->SetTitleFillColor(10);
  cmsStyle->SetTitleFontSize(0.05);
  // cmsStyle->SetTitleH(0); // Set the height of the title box
  // cmsStyle->SetTitleW(0); // Set the width of the title box
  // cmsStyle->SetTitleX(0); // Set the position of the title box
  // cmsStyle->SetTitleY(0.985); // Set the position of the title box
  // cmsStyle->SetTitleStyle(Style_t style = 1001);
  // cmsStyle->SetTitleBorderSize(2);

  // For the axis titles:

  cmsStyle->SetTitleColor(1, "XYZ");
  cmsStyle->SetTitleFont(42, "XYZ");
  cmsStyle->SetTitleSize(0.05, "XYZ");
  // cmsStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // cmsStyle->SetTitleYSize(Float_t size = 0.02);
  cmsStyle->SetTitleXOffset(1.1);
  cmsStyle->SetTitleYOffset(1.1);
  // cmsStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:

  cmsStyle->SetLabelColor(1, "XYZ");
  cmsStyle->SetLabelFont(42, "XYZ");
  cmsStyle->SetLabelOffset(0.007, "XYZ");
  cmsStyle->SetLabelSize(0.04, "XYZ");

  // For the axis:

  cmsStyle->SetAxisColor(1, "XYZ");
  cmsStyle->SetStripDecimals(kTRUE);
  cmsStyle->SetTickLength(0.03, "XYZ");
  cmsStyle->SetNdivisions(510, "XYZ");
  cmsStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  cmsStyle->SetPadTickY(1);

  // Change for log plots:
  cmsStyle->SetOptLogx(0);
  cmsStyle->SetOptLogy(0);
  cmsStyle->SetOptLogz(0);

  // Postscript options:
  cmsStyle->SetPaperSize(20.,20.);
  // cmsStyle->SetLineScalePS(Float_t scale = 3);
  // cmsStyle->SetLineStyleString(Int_t i, const char* text);
  // cmsStyle->SetHeaderPS(const char* header);
  // cmsStyle->SetTitlePS(const char* pstitle);

  // cmsStyle->SetBarOffset(Float_t baroff = 0.5);
  // cmsStyle->SetBarWidth(Float_t barwidth = 0.5);
  // cmsStyle->SetPaintTextFormat(const char* format = "g");
  // cmsStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // cmsStyle->SetTimeOffset(Double_t toffset);
  // cmsStyle->SetHistMinimumZero(kTRUE);

  cmsStyle->SetHatchesLineWidth(5);
  cmsStyle->SetHatchesSpacing(0.05);

  cmsStyle->cd();

}
