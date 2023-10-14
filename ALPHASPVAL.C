//   Author: Thabo Pilusa

#include <TGraphErrors.h>
void ALPHASPVAL() 
{

     // creating a Canvas 
   TCanvas* canvas = new TCanvas("canvas", "Multiple Line Graphs", 1800, 1800);
    //Creating multiple TGRaph objects and filling them with data
   Int_t n = 9; 

   //Create graph and open the txt files
   TGraphErrors* graph1 = new TGraphErrors( "atlaspval25.txt");
  TGraphErrors* graph2 = new TGraphErrors("atlaspval27.txt");
   TGraphErrors* graph3 = new TGraphErrors("atlaspval29.txt");
  //  TGraphErrors* graph4 = new TGraphErrors("atlaspval33.txt");
     TGraphErrors* graph5 = new TGraphErrors("atlascombpval.txt");
   // TGraph* graph1 = new TGraph("cmspval27.txt");
   
  
   
   // Customize the plots
  graph1->SetLineColor(kRed);
   graph2->SetLineColor(kBlue);
   graph3->SetLineColor(kOrange);
//   graph4->SetLineColor(kGreen);
   graph5->SetLineColor(kBlack);
   
  graph1->SetLineWidth(3);
  graph2->SetLineWidth(3);
  graph3->SetLineWidth(3);
//  graph4->SetLineWidth(3);
  graph5->SetLineWidth(4);

   //name ze axis name
    graph1->GetXaxis()->SetTitle("m_{Y} [GeV]");
   graph1->GetYaxis()->SetTitle("p-value");
   graph1->SetTitle(" ");
   graph1->GetXaxis()->SetTitleOffset(1.);
   graph1->GetXaxis()->SetTitleSize(0.04);
   graph1->GetYaxis()->SetTitleOffset(1.);
   graph1->GetYaxis()->SetTitleSize(0.04);
  
   //Draw first graph
   graph1->Draw("AC");
   
    
   //PLot the other graphs
   graph2->Draw("CSAME");
   graph3->Draw("CSAME");
 //  graph4->Draw("CSAME");
   graph5->Draw("CSAME");
   
   
   
   //graph1->GetXaxis()->SetRangeUser(65,110);
   //graph1->GetXaxis()->SetLimits(68,112);
   graph1->GetYaxis()->SetRangeUser(0.00001,100);//remove 2
    //graph1->GetYaxis()->SetLimits(0.0000001,10000);
  
   
   
   
   
   // Add a legend
   TLegend* legend = new TLegend(0.9, 0.9, 0.7, 0.7); // set the legend position
   legend->AddEntry(graph1, " 0.24 < #alpha < 0.26", "l");
   legend->AddEntry(graph2, " 0.26 < #alpha < 0.28", "l");
   legend->AddEntry(graph3, " 0.28 < #alpha < 0.30 ", "l");
//  legend->AddEntry(graph4, "  0.32 < #alpha < 0.34 ", "l");
   legend->AddEntry(graph5, " Combination", "l");
  
   legend->Draw(); //Draw ze legend
   
   
   
  //Adding ze straightlins (Beginner level style  🥲️)
  TLine *line = new TLine(2899.191, 1.008, 4100.809, 1.008);
  line->SetLineWidth(2);
  line->Draw();
  
  TLine *line1 = new TLine(2899.191, 0.331, 4100.809, 0.321);
  line1->SetLineWidth(1);
  line1->SetLineStyle(9);
  line1->Draw();
  
  TLine *line2 = new TLine(2899.191, 0.042,4100.809, 0.042);
  line2->SetLineWidth(1);
  line2->SetLineStyle(9);
  line2->Draw();
  
  TLine *line3 = new TLine(2899.191, 0.0028,4100.809, 0.0028);
  line3->SetLineWidth(1);
  line3->SetLineStyle(9);
  line3->Draw();
  
  TLine *line4 = new TLine(2899.191, 0.000061, 4100.809, 0.000061);
  line4->SetLineWidth(1);
  line4->SetLineStyle(9);
  line4->Draw();
  
  
  TLine *line5 = new TLine(2.899, 0.00000058, 4.099, 0.00000058);
  line5->SetLineWidth(1);
  line5->SetLineStyle(9);
 // line5->Draw();
  //LOG SCALE (Seems to work on TPad not TGraph, so set log y on the pad on your own)
   //graph1->SetLogy();
   
  TLine *line6 = new TLine(66.404, 0.0000000020, 113.960, 0.0000000030);
  line6->SetLineWidth(1);
  line6->SetLineStyle(9);
 // line6->Draw();
  //LOG SCALE (Seems to work on TPad not TGraph, so set log y on the pad on your own)
   //graph1->SetLogy();
  
  
  
  
  
  
  
  //This is for adding any text to the plot
  
    TLatex * text = new TLatex();
    text->SetNDC();
    text->SetTextFont(42);
    text->SetTextSize(gStyle->GetTextSize() * 0.75);
    text->SetTextColor(1);
  
    text->DrawLatex(0.91, 0.60, "1#sigma");
    text->DrawLatex(0.91, 0.50 , "2#sigma");
    text->DrawLatex(0.91, 0.37 , "3#sigma");
    text->DrawLatex(0.9, 0.18 , "4#sigma");
    //text->DrawLatex(0.9, 0.15, "5#sigma");
    //text->DrawLatex(0.9, 0.12, "6#sigma");
 //    THE END!!
    
}
