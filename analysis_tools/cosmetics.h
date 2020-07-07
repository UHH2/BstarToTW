

TPad* SetupPad() {
  Float_t yplot = 0.65;
  Float_t yratio = 0.34;

  //  coordinates:
  // set up the coordinates of the two pads:    //  			 
  Float_t y1, y2, y3;                           //  y3 +-------------+	
  y3 = 0.99;                                    //     |             |	
  y2 = y3-yplot;                                //     |     pad1    |	
  y1 = y2-yratio;                               //  y2 |-------------|	
  Float_t x1, x2;                               //     |     rp1     |	
  x1 = 0.01;                                    //  y1 +-------------+	
  x2 = 0.99;                                    //     x1            x2	
                                                // 			
                                                // No Pad 2!            
                                              
  TPad* m_pad = new TPad("pad1", "Control Plots 2", x1, y1, x2, y3);             

  m_pad->SetTopMargin(0.05);
  m_pad->SetBottomMargin(0.12);
  m_pad->SetLeftMargin(0.14);
  m_pad->SetRightMargin(0.05);

  return m_pad;
}

TPad* SetupPad2D() {
  Float_t yplot = 0.65;
  Float_t yratio = 0.34;

  //  coordinates:
  // set up the coordinates of the two pads:    //  			 
  Float_t y1, y2, y3;                           //  y3 +-------------+	
  y3 = 0.99;                                    //     |             |	
  y2 = y3-yplot;                                //     |     pad1    |	
  y1 = y2-yratio;                               //  y2 |-------------|	
  Float_t x1, x2;                               //     |     rp1     |	
  x1 = 0.01;                                    //  y1 +-------------+	
  x2 = 0.99;                                    //     x1            x2	
                                                // 			
                                                // No Pad 2!            
                                              
  TPad* m_pad = new TPad("pad1", "Control Plots 2", x1, y1, x2, y3);             

  m_pad->SetTopMargin(0.05);
  m_pad->SetBottomMargin(0.12);
  m_pad->SetLeftMargin(0.1);
  m_pad->SetRightMargin(0.14);

  return m_pad;
}

TPad* SetupRatioPad() {
  Float_t yplot = 0.65;
  Float_t yratio = 0.34;

  //  coordinates:
  // set up the coordinates of the two pads:    //  			 
  Float_t y1, y2, y3;                           //  y3 +-------------+	
  y3 = 0.99;                                    //     |             |	
  y2 = y3-yplot;                                //     |     pad1    |	
  y1 = y2-yratio;                               //  y2 |-------------|	
  Float_t x1, x2;                               //     |     rp1     |	
  x1 = 0.01;                                    //  y1 +-------------+	
  x2 = 0.99;                                    //     x1            x2	
                                                // 			
                                                // No Pad 2!            
                      
  TPad* m_rp = new TPad("rp1", "Ratio2", x1, y1, x2, y2);
                              
  m_rp->SetTopMargin(0.0);
  m_rp->SetBottomMargin(0.35);
  m_rp->SetLeftMargin(0.19);
  m_rp->SetRightMargin(0.05);

  return m_rp;
}

TPad* SetupRatioPadTop() {
  Float_t yplot = 0.65;
  Float_t yratio = 0.34;

  //  coordinates:
  // set up the coordinates of the two pads:    //  			 
  Float_t y1, y2, y3;                           //  y3 +-------------+	
  y3 = 0.99;                                    //     |             |	
  y2 = y3-yplot;                                //     |     pad1    |	
  y1 = y2-yratio;                               //  y2 |-------------|	
  Float_t x1, x2;                               //     |     rp1     |	
  x1 = 0.01;                                    //  y1 +-------------+	
  x2 = 0.99;                                    //     x1            x2	
                                                // 			
                                                // No Pad 2!            
                      
  TPad* m_rp_top = new TPad("pad1", "Control Plots 2", x1, y2, x2, y3);
 
                              
  m_rp_top->SetTopMargin(0.065);
  m_rp_top->SetBottomMargin(0.0);
  m_rp_top->SetLeftMargin(0.19);
  m_rp_top->SetRightMargin(0.05);
 
  return m_rp_top;
}

void Hist_Cosmetics(TH1* hist, bool ratio = false) {
 /* hist->SetLineWidth(2); */
  /* hist->SetMarkerStyle(8); */
  /* hist->SetMarkerSize(0.7); */
  
  // X label
  hist->GetXaxis()->SetLabelFont(43);
  hist->GetXaxis()->SetLabelSize(18);
  // X title
  hist->GetXaxis()->SetTitleFont(43);
  hist->GetXaxis()->SetTitleSize(18);
    
  // Y label
  hist->GetYaxis()->SetLabelFont(43);
  hist->GetYaxis()->SetLabelSize(18);
  hist->GetYaxis()->SetNdivisions(505);
  // Y title
  hist->GetYaxis()->SetTitleFont(43);
  hist->GetYaxis()->SetTitleSize(18);

  hist->GetXaxis()->SetNdivisions(505);
  hist->SetMinimum(0.01); 

  // offset
  hist->GetYaxis()->SetTitleOffset(1.4);
  if (ratio)
    {
      hist->GetXaxis()->SetTitleOffset(3);
    } 
  else
    {
      hist->GetXaxis()->SetTitleOffset(1.3);
    }
}

void Hist_Cosmetics(TGraph* hist, bool ratio = false) {
  hist->SetLineWidth(2);
  hist->SetMarkerStyle(8);
  hist->SetMarkerSize(0.7);
  
  // X label
  hist->GetXaxis()->SetLabelFont(43);
  hist->GetXaxis()->SetLabelSize(18);
  // X title
  hist->GetXaxis()->SetTitleFont(43);
  hist->GetXaxis()->SetTitleSize(18);
    
  // Y label
  hist->GetYaxis()->SetLabelFont(43);
  hist->GetYaxis()->SetLabelSize(18);
  // Y title
  hist->GetYaxis()->SetTitleFont(43);
  hist->GetYaxis()->SetTitleSize(18);

  hist->GetXaxis()->SetNdivisions(505);
  hist->SetMinimum(0.01); 

  // offset
  hist->GetYaxis()->SetTitleOffset(1.4);
  if (ratio)
    {
      hist->GetXaxis()->SetTitleOffset(3);
    } 
}


void Hist_Cosmetics(TMultiGraph* hist, bool ratio = false) {
  
  // X label
  hist->GetXaxis()->SetLabelSize(.05);
  hist->GetXaxis()->SetNdivisions(505);
  // X title
  hist->GetXaxis()->SetTitleSize(.05);
  hist->GetXaxis()->SetTitleOffset(1.1);
  // Y label
  hist->GetYaxis()->SetLabelSize(.05);
  hist->GetYaxis()->SetNdivisions(505);
  // Y title
  hist->GetYaxis()->SetTitleSize(.05);
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->SetMinimum(0.01); 
  if (ratio)
    {
      hist->GetXaxis()->SetTitleOffset(1.6);
    } 
}

void Hist_Cosmetics_2D(TH1* hist) {
  
  hist->SetTitle("");
  // X label
  hist->GetXaxis()->SetLabelSize(.05);
  hist->GetXaxis()->SetNdivisions(505);
  // X title
  hist->GetXaxis()->SetTitleSize(.05);
  hist->GetXaxis()->SetTitleOffset(1.1);
  // Y label
  hist->GetYaxis()->SetLabelSize(.05);
  hist->GetYaxis()->SetNdivisions(505);
  // Y title
  hist->GetYaxis()->SetTitleSize(.05);
  hist->GetYaxis()->SetTitleOffset(0.9);
  hist->SetMinimum(0.01); 
  // Z label
  hist->GetZaxis()->SetLabelSize(.05);
  hist->GetZaxis()->SetNdivisions(505);

}

