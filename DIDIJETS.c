/*            #################################################
             #              🏃️🏃‍♂️️🏃‍♀️️💨️                         #
            #                   By T. Pilusa                #
           #                                                #
          #  Di-dijets macro : 2023 - 2024 WITS ICPP       #
         #  @copyrights reserved                          #
        #                                                #
       ##################################################


*/

#include <iostream>
#include <cmath>
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TCanvas.h" 

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

void DIDIJETS()
{
  
//##################Load events For SIGNAL#################  
  
  // Load the Delphes library
  gSystem->Load("libDelphes");

  // Create a chain of root trees
  TChain chain("Delphes");
  //load the events from the root file
  chain.Add("/home/electronicslab/Desktop/MG5_aMC_v3_5_3/DIDIJETS_NLO_v3/Events/run_02/tag_1_delphes_events.root");

  // Create an object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");



//#######################Do the same For BACKGROUND#########
// Create a chain for background events
TChain backgroundChain("Delphes");
backgroundChain.Add("/home/electronicslab/Desktop/MG5_aMC_v3_5_3/DIDIJETSBACKGROUNDpp2wpm/Events/run_02/tag_2_delphes_events.root");

// Create an object of class ExRootTreeReader for background events
ExRootTreeReader *backgroundTreeReader = new ExRootTreeReader(&backgroundChain);
Long64_t numberOfBackgroundEntries = backgroundTreeReader->GetEntries();





  // Book Histograms
     //..........Eta.............
  
  TH1F *histetaj1 = new TH1F("histEtaj1", "#eta_{j1};#eta_{j1} [rad];Event Count", 100, -3.15, 3.15); // Adjust the range as needed
  TH1F *histetaj2 = new TH1F("histEtaj2", "#eta_{j2};#eta_{j2} [rad];Event Count", 100, -3.15, 3.15); // Adjust the range as needed
  TH1F *histetaj3 = new TH1F("histEtaj3", "#eta_{j3};#eta_{j3} [rad];Event Count", 100, -3.15, 3.15); // Adjust the range as needed
  TH1F *histetaj4 = new TH1F("histEtaj4", "#eta_{j4};#eta_{j4} [rad];Event Count", 100, -3.15, 3.15); // Adjust the range as needed
   
          //.......Phi......... 
  
  TH1F *histpj1 = new TH1F("histPhij1", "#phi_{j1};#phi_{j1} [rad];Event Count", 100, -4, 4); // Adjust the range as needed
  TH1F *histpj2 = new TH1F("histPhij2", "#phi_{j2};#phi_{j2} [rad];Event Count", 100, -4, 4); // Adjust the range as needed
  TH1F *histpj3 = new TH1F("histPhij3", "#phi_{j3};#phi_{j3} [rad];Event Count", 100, -4, 4); // Adjust the range as needed
  TH1F *histpj4 = new TH1F("histPhij4", "#phi_{j4};#phi_{j4} [rad];Event Count", 100, -4, 4); // Adjust the range as needed
  
       //-------Delta-eta
           
           
 TH1F *histetaj1j2 = new TH1F("histEtaj1", "#Delta#eta_{j1j2};#Delta#eta_{j1j2} ;Event Count", 100, -3.15, 3.15); // Adjust the range as needed
  TH1F *histetaj3j4 = new TH1F("histEtaj2", "#Delta#eta_{j3j4};#Delta#eta_{j3j4} ;Event Count", 100, -3.15, 3.15); // Adjust the range as needed
  TH1F *histetaj1j3 = new TH1F("histEtaj3", "#Delta#eta_{j1j3};#Delta#eta_{j1j3} ;Event Count", 100, -3.15, 3.15); // Adjust the range as needed
  TH1F *histetaj2j4 = new TH1F("histEtaj4", "#Delta#eta_{j2j4};#Delta#eta_{j2j4} ;Event Count", 100, -3.13, 3.15); // Adjust the range as needed
  TH1F *histetaj1j4 = new TH1F("histEtaj3", "#Delta#eta_{j1j4};#Delta#eta_{j1j4} ;Event Count", 100, -3.15, 3.15); // Adjust the range as needed
  TH1F *histetaj2j3 = new TH1F("histEtaj3", "#Delta#eta_{j2j3};#Delta#eta_{j2j3};Event Count", 100, -3.15, 3.15); // Adjust the range as needed
  //-----------Delta-phi------------
  
  TH1F *histphij1j2 = new TH1F("histEtaj1", "#Delta#phi_{j1j2};#Delta#phi_{j1j2} ;Event Count", 100, -4, 4); // Adjust the range as needed
  TH1F *histphij3j4 = new TH1F("histEtaj2", "#Delta#phi_{j3j4};#Delta#phi_{j3j4} ;Event Count", 100, -4, 4); // Adjust the range as needed
  TH1F *histphij1j3 = new TH1F("histEtaj3", "#Delta#phi_{j1j3};#Delta#phi_{j1j3} ;Event Count", 100, -4, 4); // Adjust the range as needed
  TH1F *histphij2j4 = new TH1F("histEtaj4", "#Delta#phi_{j2j4};#Delta#phi_{j2j4} ;Event Count", 100,-4, 4); // Adjust the range as needed
  TH1F *histphij1j4 = new TH1F("histEtaj3", "#Delta#phi_{j1j4};#Delta#phi_{j1j4} ;Event Count", 100, -4, 4); // Adjust the range as needed
  TH1F *histphij2j3 = new TH1F("histEtaj3", "#Delta#phi_{j2j3};#Delta#phi_{j2j3};Event Count", 100, -4, 4); // Adjust the range as needed
  
  
       //..........ΔR1 and ΔR2.............
  
  TH1F *histRj1j2 = new TH1F("histR1", "#DeltaR_{j1j2};#DeltaR_{j1j2};Event Count", 100, -4, 4); // Adjust the range as needed
  TH1F *histRj3j4 = new TH1F("histR2", "#DeltaR_{j3j4};#DeltaR_{j3j4}  ;Event Count", 100, -4, 4); // Adjust the range as needed
  TH1F *histRj1j3 = new TH1F("histR1", "#DeltaR_{j1j3};#DeltaR_{j1j3};Event Count", 100, -4, 4); // Adjust the range as needed
  TH1F *histRj2j4 = new TH1F("histR2", "#DeltaR_{j2j4};#DeltaR_{j2j4} ;Event Count", 100, -4, 4); // Adjust the range as needed
  TH1F *histRj1j4 = new TH1F("histR1", "#DeltaR_{j1j4};#DeltaR_{j1j4} ;Event Count", 100, -4, 4); // Adjust the range as needed
  TH1F *histRj2j3 = new TH1F("histR2", "#DeltaR_{j2j3};#DeltaR_{j2j3} ;Event Count", 100, -4, 4); // Adjust the range as needed
  
       //........................ΔR...........
          
  TH1F *histRalpha = new TH1F("histRalpha", "#DeltaR_{Alpha};#DeltaR ;Event Count", 100, -6, 6); // Adjust the range as needed
  TH1F *histRbeta = new TH1F("histbeta", "#DeltaR_{Beta} ;#DeltaR ;Event Count", 100,-6, 6); // Adjust the range as needed
  TH1F *histRgamma = new TH1F("histRgamma", "#DeltaR_{Gamma};#DeltaR ;Event Count", 100,-6, 6); // Adjust the range as needed
                 // also histogram for the minimum ΔR

  TH1F *histminR= new TH1F("histRgamma", "Minimum #DeltaR;#DeltaR [rad];Event Count", 100, -6, 6); // Adjust the range as needed      
  
         //....invariant mass ......
       
              
  TH1F *histMj1j2 = new TH1F("histMj1j2", "M_{j1j2}; M_{j1j2} [GeV] ;Event Count", 100, 0, 3000);
  TH1F *histMj1j3 = new TH1F("histMj1j3", "M_{j1j3}; M_{j1j3} [GeV] ;Event Count", 100, 0, 3000);  
  TH1F *histMj1j4 = new TH1F("histMj1j4", "M_{j1j4}; M_{j1j4} [GeV] ;Event Count", 100, 0, 3000);    
  TH1F *histMj2j3 = new TH1F("histMj2j3", "M_{j2j3}; M_{j2j3} [GeV] ;Event Count", 100, 0, 3000);  
  TH1F *histMj2j4 = new TH1F("histMj2j4", "M_{j2j4}; M_{j2j4} [GeV] ;Event Count", 100, 0, 3000);  
  TH1F *histMj3j4 = new TH1F("histMj3j4", "M_{j3j4}; M_{j3j4} [GeV] ;Event Count", 100, 0, 3000); 
  TH1F *histM4j = new TH1F("histM4j", "M_{4j}; M_{4j} [GeV] ;Event Count", 100, 0, 5000); // Adjust the range as needed
  
  
                     //......asymetry........
                     
    TH1F *histAsy = new TH1F("histAsy", "Asymetry;Assymetry ;Event Count", 100, -1, 1); // Adjust the range as needed
    
                   //..............PLot ∆η = |η1 − η2 | < 1.1..............
                   

    TH1F *histSETA = new TH1F("histSETA", "#Delta#eta ;#Delta#eta ;Event Count", 100, -5, 5); // Adjust the range as needed
  
  
                   
//--------------------------ALPHA--------------------

    TH1F *histAlpha = new TH1F("histAlpha", "#alpha ;#alpha ;Event Count", 10, 0, 1); // Adjust the range as needed
 
  
  
               //.........acceptance........
  
  TH1F *histAcc = new TH1F("histAcc", "M_{Z} = 3.5 TeV, M_{W} = 1 TeV ; #acceptance ;Event Count", 10, 0.07796559502133332, 0.3611505733363813);
  
  //------------Grand Finale[ Acceptance nd Alpha]
  
  
  TH1F *histAccAlpha = new TH1F("histAccAlpha", "ACCEPTANCE ; #alpha ;Acceptance", 12, 0.07796559502133332, 0.3611505733363813); 
  
  
  
    
  // [[[[[[[[[[[[[[[[[[[[[[[[[[-----------------END OF BOOKING HISTOGRAMS----------------]]]]]]]]]]]]]]]]]]]]]]]]]]]
  
 
 
 
//-----------------------------------------------------------------
 
  // Signal Event counters
  Int_t totalEvents = 0;
  Int_t passedBothCuts = 0;
  Int_t passedPtJetCount = 0;
  Int_t passedEtaJetCount = 0;
  Int_t passedRCount = 0;
  Int_t passedSETACount = 0;
  Int_t passedAsyCount = 0;
  Int_t fourjets = 0;
  // Background Event counters
  Int_t totalBackgroundEvents = 0 ;
  Int_t passedPtBackgroundJetCount = 0 ;
  Int_t passedEtaBackgroundJetCount = 0 ;
  Int_t passedBothBackgroundCuts = 0;
  Int_t passedRCount_bg = 0;
  Int_t passedSETACount_bg = 0;
  Int_t passedAsyCount_bg = 0;
  Int_t fourbackgroundjets = 0; 
//-----------------------------------------------------------
  
  
/*
---------------------------------------------

          Insert ALPHA HERE 

  -----------------------     --------------
                           ||
                           \/       
 
 
 */   
 
   Double_t desiredAlpha25 = 0.25;   
   Double_t desiredAlpha27 = 0.27;
   Double_t desiredAlpha29 = 0.29;
   Double_t desiredAlpha31 = 0.31; 
   Double_t desiredAlpha33 = 0.33;       



 
  

   
//------------Event selection for SIGNAL starts here -------------------------
          
  // Loop over all signal events
  for (Long64_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from the specified event
    treeReader->ReadEntry(entry);
    Event *event = (Event *)branchEvent->At(0); 

    // Access jets from the Jet branch
    TClonesArray *jets = branchJet;
    Int_t numberofjets = jets->GetEntries();

    // Increment the total event counter
    totalEvents++;

    // Check the number of jets in each event
    if (jets->GetEntries() >= 4)
    {  fourjets++;
        
        
        /* for (Int_t i = 0; i < numberofjets; ++i) //iz what's causing the events to misbehave
         {
           Jet *jet = (Jet*)jets->At(i); 
            // Apply pT cut for each jet
            if (jet->PT > 60)
            {
            passedPtJetCount++;
            }
            
            // Apply eta cut for each jet
            if (std::abs(jet->Eta) < 2.4)
            {
            passedEtaJetCount++;
            }
            
         }

*/
        // Access the first four leading jets
        Jet *jet1 = (Jet *)jets->At(0);
        Jet *jet2 = (Jet *)jets->At(1);
        Jet *jet3 = (Jet *)jets->At(2);
        Jet *jet4 = (Jet *)jets->At(3);

      // Check if the event passes both pT and eta cuts FOR SIGNAL
      if (jet1->PT > 60 && jet2->PT > 60 && jet3->PT > 60 && jet4->PT > 60)
      {
        passedPtJetCount++;
       if (std::abs(jet1->Eta) < 2.4 && std::abs(jet2->Eta) < 2.4 && std::abs(jet3->Eta) < 2.4 && std::abs(jet4->Eta) < 2.4)
        {
         passedEtaJetCount++;

        // Declaring Lorentz Vectors for jet 1,2,3,4
        TLorentzVector vecJet1;
        TLorentzVector vecJet2;
        TLorentzVector vecJet3;
        TLorentzVector vecJet4;



        // Assign Lorentz Vectors to 4-momenta
        vecJet1 = jet1->P4();
        vecJet2 = jet2->P4();
        vecJet3 = jet3->P4();
        vecJet4 = jet4->P4();

        // Calculate eta for each jet
        Double_t etaj1 = vecJet1.Eta();
        Double_t etaj2 = vecJet2.Eta();
        Double_t etaj3 = vecJet3.Eta();
        Double_t etaj4 = vecJet4.Eta();

        histetaj1->Fill(etaj1);
        histetaj2->Fill(etaj2);
        histetaj3->Fill(etaj3);
        histetaj4->Fill(etaj4);

        // Calculate phi for each jet
        Double_t phij1 = vecJet1.Phi();
        Double_t phij2 = vecJet2.Phi();
        Double_t phij3 = vecJet3.Phi();
        Double_t phij4 = vecJet4.Phi();

        histpj1->Fill(phij1);
        histpj2->Fill(phij2);
        histpj3->Fill(phij3);
        histpj4->Fill(phij4);
        //.......LET'S CALCULATE DELTA-ETA AND DELTA-PHI
               //--Eta
        Double_t etaj1j2 = std::abs(std::abs(etaj1) -  std::abs(etaj2));
        Double_t etaj3j4 = std::abs(std::abs(etaj3) -  std::abs(etaj4));
        Double_t etaj1j3 = std::abs(std::abs(etaj1) -  std::abs(etaj3));
        Double_t etaj2j4 = std::abs(std::abs(etaj2) -  std::abs(etaj4));
        Double_t etaj1j4 = std::abs(std::abs(etaj1) -  std::abs(etaj4));
        Double_t etaj2j3 = std::abs(std::abs(etaj2) - std::abs(etaj3));
        

        
        
           histetaj1j2->Fill(etaj1j2);
           histetaj3j4->Fill(etaj3j4);
           histetaj1j3->Fill(etaj1j3);
           histetaj2j4->Fill(etaj2j4);
           histetaj1j4->Fill(etaj1j4);
           histetaj2j3->Fill(etaj2j3);        
       
          //---phi----
          
          
        Double_t phij1j2 = std::abs(std::abs(phij1) - std::abs(phij2));
        Double_t phij3j4 = std::abs(std::abs(phij3) - std::abs(phij4));
        Double_t phij1j3 = std::abs(std::abs(phij1) - std::abs(phij3));
        Double_t phij2j4 = std::abs(std::abs(phij2) - std::abs(phij4));
        Double_t phij1j4 = std::abs(std::abs(phij1) - std::abs(phij4));
        Double_t phij2j3 = std::abs(std::abs(phij2) - std::abs(phij3));
        
           histphij1j2->Fill(phij1j2);
           histphij3j4->Fill(phij3j4);
           histphij1j3->Fill(phij1j3);
           histphij2j4->Fill(phij2j4);
           histphij1j4->Fill(phij1j4);
           histphij2j3->Fill(phij2j3);
        
//......THIS IS WHERE we START TO MAKE COMB1, COMB2, COMB3 event-by-event instead of doing one combination for all events..........
  //......We basically try to see which combination is the best for a certain event, instead of taking the best combination for the overall events........
     //......This means each event will have it's own comination.........

              
        //calculate ΔR1 and ΔR2 for different jet combinations
        
        Double_t Rj1j2 =  std::sqrt(std::pow((etaj1 - etaj2), 2) + std::pow((phij1 - phij2), 2));

     
        Double_t Rj3j4 = std::sqrt(std::pow((etaj3 - etaj4), 2) + std::pow((phij3 - phij4), 2));

        
        Double_t Rj1j3 = std::sqrt(std::pow((etaj1 - etaj3), 2) + std::pow((phij1 - phij3), 2));
        
      
        Double_t Rj2j4 = std::sqrt(std::pow((etaj2 - etaj4), 2) + std::pow((phij2 - phij4), 2));
        
        
        Double_t Rj1j4 = std::sqrt(std::pow((etaj1 - etaj4), 2) + std::pow((phij1 - phij4), 2));
               
        
        Double_t Rj2j3 = std::sqrt(std::pow((etaj2 - etaj3), 2) + std::pow((phij2 - phij3), 2));
        
            //fill histograms
            
            histRj1j2->Fill(Rj1j2);
            histRj3j4->Fill(Rj3j4);
            histRj1j3->Fill(Rj1j3);
            histRj2j4->Fill(Rj2j4);
            histRj1j4->Fill(Rj1j4);
            histRj2j3->Fill(Rj2j3);
            
           //calculate R for the 3 combinations 
         Double_t Ralpha = std::abs((Rj1j2 - 0.8)) + std::abs((Rj3j4 - 0.8));   
         Double_t Rbeta = std::abs((Rj1j3 - 0.8)) + std::abs((Rj2j4 - 0.8));  
         Double_t Rgamma = std::abs((Rj1j4 - 0.8)) + std::abs((Rj2j3 - 0.8)); 

    
  
            //Fill histograms
                
            histRalpha->Fill(Ralpha);
            histRbeta->Fill(Rbeta);
            histRgamma->Fill(Rgamma);

         
     /*       
   // Loop over all possible combinations of changeRi
          Double_t minChangeR = std::numeric_limits<Double_t>::max();

            // Loop over all possible combinations of changeRi for the current event
            for (int i = 1; i <= 6; ++i)
            {
                Double_t changeRi1 = 0.0;
                Double_t changeRi2 = 0.0;

                // Calculate changeRi based on the combination i
                switch (i)
                {
                case 1:
                    changeRi1 = std::sqrt(std::pow(etaj1 - etaj2, 2) + std::pow(phij1 - phij2, 2));
                    changeRi2 = std::sqrt(std::pow(etaj3 - etaj4, 2) + std::pow(phij3 - phij4, 2));
                    break;
                case 2:
                    changeRi1 = std::sqrt(std::pow(etaj3 - etaj4, 2) + std::pow(phij3 - phij4, 2));
                    changeRi2 = std::sqrt(std::pow(etaj1 - etaj2, 2) + std::pow(phij1 - phij2, 2));
                    break;
                case 3:
                    changeRi1 = std::sqrt(std::pow(etaj1 - etaj3, 2) + std::pow(phij1 - phij3, 2));
                    changeRi2 = std::sqrt(std::pow(etaj2 - etaj4, 2) + std::pow(phij2 - phij4, 2));
                    break;
                case 4:
                    changeRi1 = std::sqrt(std::pow(etaj2 - etaj4, 2) + std::pow(phij2 - phij4, 2));
                    changeRi2 = std::sqrt(std::pow(etaj1 - etaj3, 2) + std::pow(phij1 - phij3, 2));
                    break;
                case 5:
                    changeRi1 = std::sqrt(std::pow(etaj1 - etaj4, 2) + std::pow(phij1 - phij4, 2));
                    changeRi2 = std::sqrt(std::pow(etaj2 - etaj3, 2) + std::pow(phij2 - phij3, 2));
                    break;
                case 6:
                    changeRi1 = std::sqrt(std::pow(etaj2 - etaj3, 2) + std::pow(phij2 - phij3, 2));
                    changeRi2 = std::sqrt(std::pow(etaj1 - etaj4, 2) + std::pow(phij1 - phij4, 2));
                    break;
                }

                // Calculate changeR based on the formula
                Double_t changeR = std::abs((changeRi1 - 0.8)) + std::abs((changeRi2 - 0.8));
              
                // Update minChangeR if the current combination has a smaller value
                if (changeR < minChangeR)
                {
                    minChangeR = changeR;
                      histR->Fill(minChangeR);
                }
            }

// Now minChangeR contains the minimum changeR value



        // Apply the R cut
        
        */
        
//      >>>>-------------If Rgamma pass the R cut combination---------------
        if (Rj1j4 < 2 && Rj2j3 < 2)
        { passedRCount++;

         
         
          //------ apply ∆η = |η1 − η2 | < 1.1 cut----------
            // η1 & η2 are the deltaetas from above
           
        Double_t ETA1 = std::abs(etaj1 - etaj4);
        Double_t ETA2 = std::abs(etaj2 - etaj3);
            
        Double_t DELTAETA = std::abs(ETA1 - ETA2);
         histSETA->Fill(DELTAETA);
              
          // Apply SETA cut
          if (DELTAETA < 1.1 )
          {
            passedSETACount++;

          // Calculate the Invariant Mass
          Double_t Mj1j2 = (vecJet1 + vecJet2).M();
          Double_t Mj1j3 = (vecJet1 + vecJet3).M(); 
          Double_t Mj1j4 = (vecJet1 + vecJet4).M();         
          Double_t Mj2j3 = (vecJet2 + vecJet3).M();
          Double_t Mj2j4 = (vecJet2 + vecJet4).M();
          Double_t Mj3j4 = (vecJet3 + vecJet4).M();
          
          Double_t M4j = (vecJet1 + vecJet4 + vecJet2 + vecJet3).M();
          
          Double_t asy = std::abs(Mj1j4 - Mj2j3) / (Mj1j4 + Mj2j3);
          
           
            histMj1j2->Fill(Mj1j2); 
            histMj1j3->Fill(Mj1j3); 
            histMj1j4->Fill(Mj1j4); 
            histMj2j3->Fill(Mj2j3); 
            histMj2j4->Fill(Mj2j4); 
            histMj3j4->Fill(Mj3j4); 
            histM4j->Fill(M4j); 
            histAsy->Fill(asy);
          //Apply the asymmetry cut
           if (asy < 0.1)
           {
          passedAsyCount++;
	  
           //Calculate Alpha
          Double_t Alpha = (Mj1j4 + Mj2j3)/(2*M4j);
          histAlpha->Fill(Alpha);
          
          
          //-histogram for alpha vs acceptance
        
        histAccAlpha->Fill(Alpha);
        
        

          
          }
        }
      }
     }
    }
  }
}

//-----------End of Events selection for SIGNAL---------------




// Get pointers to branches used in this analysis for background events
TClonesArray *backgroundBranchJet = backgroundTreeReader->UseBranch("Jet");
TClonesArray *backgroundBranchEvent = backgroundTreeReader->UseBranch("Event");



  // Now, you can analyze the background events similar to signal events
for (Long64_t backgroundEntry = 0; backgroundEntry < numberOfBackgroundEntries; ++backgroundEntry) {
    
    
 // Load selected branches with data from the specified background event
    backgroundTreeReader->ReadEntry(backgroundEntry);
    Event *backgroundEvent = (Event *)backgroundBranchEvent->At(0);

    // Access jets from the Jet branch
    TClonesArray *backgroundJets = backgroundBranchJet;

    // Increment the total background event counter
    totalBackgroundEvents++;
    
// Check the number of jets in each background event
    if (backgroundJets->GetEntries() >= 4)
    { fourbackgroundjets++;
        // Count jets that pass the pT cut
       /* for (Int_t iJet = 0; iJet < backgroundJets->GetEntries(); ++iJet)
        {
            Jet *backgroundJet = (Jet *)backgroundJets->At(iJet);

            // Apply pT cut
            if (backgroundJet->PT > 60.0)
            {
                passedPtBackgroundJetCount++;
             }
                
                 // Apply eta cut
            if (std::abs(backgroundJet->Eta) < 2.4)
                {
                passedEtaBackgroundJetCount++;
                }
        }
        */
        
    // Access the first four leading jets for background
    Jet *jet1_bg = (Jet *)backgroundJets->At(0);
    Jet *jet2_bg = (Jet *)backgroundJets->At(1);
    Jet *jet3_bg = (Jet *)backgroundJets->At(2);
    Jet *jet4_bg = (Jet *)backgroundJets->At(3);
        
        
         // Check if the event passes both pT and eta cuts
     if (jet1_bg->PT > 60 && jet2_bg->PT > 60 && jet3_bg->PT > 60 && jet4_bg->PT > 60)
      {passedPtBackgroundJetCount++;
        
          if (std::abs(jet1_bg->Eta) < 2.4 && std::abs(jet2_bg->Eta) < 2.4 && std::abs(jet3_bg->Eta) < 2.4 &&std::abs(jet4_bg->Eta) < 2.4)
           {passedEtaBackgroundJetCount++;
             // Declaring Lorentz Vectors for jet 1,2,3,4
             TLorentzVector vecJet1_bg;
             TLorentzVector vecJet2_bg;
             TLorentzVector vecJet3_bg;
             TLorentzVector vecJet4_bg;


    // Assign Lorentz Vectors to 4-momenta for background
    vecJet1_bg = jet1_bg->P4();
    vecJet2_bg = jet2_bg->P4();
    vecJet3_bg = jet3_bg->P4();
    vecJet4_bg = jet4_bg->P4();

    // Calculate eta for each jet for background
    Double_t etaj1_bg = vecJet1_bg.Eta();
    Double_t etaj2_bg = vecJet2_bg.Eta();
    Double_t etaj3_bg = vecJet3_bg.Eta();
    Double_t etaj4_bg = vecJet4_bg.Eta();


// Calculate phi for each jet
Double_t phij1_bg = vecJet1_bg.Phi();
Double_t phij2_bg = vecJet2_bg.Phi();
Double_t phij3_bg = vecJet3_bg.Phi();
Double_t phij4_bg = vecJet4_bg.Phi();


// Calculate delta-eta and delta-phi
//--Eta
Double_t etaj1j2_bg = std::abs(std::abs(etaj1_bg) -  std::abs(etaj2_bg));
Double_t etaj3j4_bg = std::abs(std::abs(etaj3_bg) -  std::abs(etaj4_bg));
Double_t etaj1j3_bg = std::abs(std::abs(etaj1_bg) -  std::abs(etaj3_bg));
Double_t etaj2j4_bg = std::abs(std::abs(etaj2_bg) -  std::abs(etaj4_bg));
Double_t etaj1j4_bg = std::abs(std::abs(etaj1_bg) -  std::abs(etaj4_bg));
Double_t etaj2j3_bg = std::abs(std::abs(etaj2_bg) - std::abs(etaj3_bg));
 
 // --PHI

Double_t phij1j2_bg = std::abs(std::abs(phij1_bg) - std::abs(phij2_bg));
Double_t phij3j4_bg = std::abs(std::abs(phij3_bg) - std::abs(phij4_bg));
Double_t phij1j3_bg = std::abs(std::abs(phij1_bg) - std::abs(phij3_bg));
Double_t phij2j4_bg = std::abs(std::abs(phij2_bg) - std::abs(phij4_bg));
Double_t phij1j4_bg = std::abs(std::abs(phij1_bg) - std::abs(phij4_bg));
Double_t phij2j3_bg = std::abs(std::abs(phij2_bg) - std::abs(phij3_bg));


// Calculate ΔR for different jet combinations
Double_t Rj1j2_bg = std::sqrt(std::pow((etaj1_bg - etaj2_bg), 2) + std::pow((phij1_bg - phij2_bg), 2));
Double_t Rj3j4_bg = std::sqrt(std::pow((etaj3_bg - etaj4_bg), 2) + std::pow((phij3_bg - phij4_bg), 2));
Double_t Rj1j3_bg = std::sqrt(std::pow((etaj1_bg - etaj3_bg), 2) + std::pow((phij1_bg - phij3_bg), 2));
Double_t Rj2j4_bg = std::sqrt(std::pow((etaj2_bg - etaj4_bg), 2) + std::pow((phij2_bg - phij4_bg), 2));
Double_t Rj1j4_bg = std::sqrt(std::pow((etaj1_bg - etaj4_bg), 2) + std::pow((phij1_bg - phij4_bg), 2));
Double_t Rj2j3_bg = std::sqrt(std::pow((etaj2_bg - etaj3_bg), 2) + std::pow((phij2_bg - phij3_bg), 2));

// Calculate R for the three combinations
Double_t Ralpha_bg = std::abs((Rj1j2_bg - 0.8)) + std::abs((Rj3j4_bg - 0.8));   
Double_t Rbeta_bg = std::abs((Rj1j3_bg - 0.8)) + std::abs((Rj2j4_bg - 0.8));  
Double_t Rgamma_bg = std::abs((Rj1j4_bg - 0.8)) + std::abs((Rj2j3_bg - 0.8)); 


        // Apply the R cut
        
        
if (Rj1j4_bg < 2 && Rj2j3_bg < 2)
{
    passedRCount_bg++;

    // Apply ∆η = |η1 − η2| < 1.1 cut
    // η1 & η2 are the deltaetas from above
    Double_t ETA1_bg = std::abs(etaj1_bg - etaj4_bg);
    Double_t ETA2_bg = std::abs(etaj2_bg - etaj3_bg);
            
    Double_t DELTAETA_bg = std::abs(ETA1_bg - ETA2_bg);
    
      
    
    
  if (DELTAETA_bg < 1.1)
{
    passedSETACount_bg++;
    
     // Calculate the Invariant Mass
    Double_t Mj1j2_bg = (vecJet1_bg + vecJet2_bg).M();
    Double_t Mj1j3_bg = (vecJet1_bg + vecJet3_bg).M(); 
    Double_t Mj1j4_bg = (vecJet1_bg + vecJet4_bg).M();         
    Double_t Mj2j3_bg = (vecJet2_bg + vecJet3_bg).M();
    Double_t Mj2j4_bg = (vecJet2_bg + vecJet4_bg).M();
    Double_t Mj3j4_bg = (vecJet3_bg + vecJet4_bg).M();
          
    Double_t M4j_bg = (vecJet1_bg + vecJet4_bg + vecJet2_bg + vecJet3_bg).M();
          
    Double_t asy_bg = std::abs(Mj1j4_bg - Mj2j3_bg) / (Mj1j4_bg + Mj2j3_bg);
     if (asy_bg < 0.1)
    {
        passedAsyCount_bg++;
        // Calculate Alpha
        Double_t Alpha_bg = (Mj1j4_bg + Mj2j3_bg) / (2 * M4j_bg);
        
             
      }
     }
    }
   }
  }
 }    
}
 
//__________Done doing analysis for background___________











//Now you can show resulting histograms, just uncomment the one you wish to plot :)

  //.........Plot Eta...............
/*
  TCanvas *canvas = new TCanvas("Eta", "Eta");
  canvas->Divide(2, 2);
  canvas->cd(1);
  histetaj1->Draw();
    canvas->cd(2);
  histetaj2->Draw();
    canvas->cd(3);
  histetaj3->Draw();
    canvas->cd(4);
  histetaj4->Draw();
  */
  
 /*
  
    //.........Plot Phi...............
  TCanvas *canvas = new TCanvas("Phi", "Phi");
  canvas->Divide(2, 2);
  canvas->cd(1);
  histpj1->Draw();
    canvas->cd(2);
  histpj2->Draw();
    canvas->cd(3);
  histpj3->Draw();
    canvas->cd(4);
  histpj4->Draw();

  */
  
  //--------Plot Delta-Eta----------

  
 /*
  TCanvas *canvas = new TCanvas("delta-eta", "delta-eta",1000, 1000);
  canvas->Divide(3, 2, 0.001, 0.01);
  canvas->cd(1);
   histetaj1j2->Draw();
  canvas->cd(2);
    histetaj3j4->Draw();
  canvas->cd(3);
    histetaj1j3->Draw();
  canvas->cd(4);
    histetaj2j4->Draw();
  canvas->cd(5);
     histetaj1j4->Draw();
  canvas->cd(6);
  
     histetaj2j3->Draw();
  */
  
 
 /*
 
   //---------plot Delta-Phi---------
 TCanvas *canvas = new TCanvas("delta-phi", "dtlta-phi",1000, 1000);
  canvas->Divide(3, 2, 0.001, 0.01);
  canvas->cd(1);
   histphij1j2->Draw();
  canvas->cd(2);
    histphij3j4->Draw();
  canvas->cd(3);
    histphij1j3->Draw();
  canvas->cd(4);
    histphij2j4->Draw();
  canvas->cd(5);
     histphij1j4->Draw();
  canvas->cd(6);
     histphij2j3->Draw();
  */   
    
  
  
 //.........Plot ΔR1 and ΔR2 .........
/*
  TCanvas *canvas = new TCanvas("R1", "R2",1000, 1000);
  canvas->Divide(3, 2, 0.001, 0.01);
  canvas->cd(1);
   histRj1j2->Draw();
  canvas->cd(2);
    histRj3j4->Draw();
  canvas->cd(3);
    histRj1j3->Draw();
  canvas->cd(4);
    histRj2j4->Draw();
  canvas->cd(5);
     histRj1j4->Draw();
  canvas->cd(6);
     histRj2j3->Draw();
*/
/* 
  //.........Plot ΔR for different dijet pairs

    TCanvas *canvas = new TCanvas("R1", "R2",1000, 1000);
  canvas->Divide(3, 2, 0.001, 0.01);
  canvas->cd(1);
   histRalpha->Draw();
  canvas->cd(2);
    histRbeta->Draw();
  canvas->cd(3);
    histRgamma->Draw();
 */

    //.........Plot Sum of Eta1&Eta2, Eta3&Eta4, and Sum of all.........
  /*
   TCanvas *canvas = new TCanvas("Sum of Eta", "Sum of Eta");
   canvas->Divide(1, 1);
   canvas->cd(1);
   histSETA->Draw();
 
  */
 //.......plot invariant mass

     //------plot first 3, then last 3, finally 4j------

//first 3

/*
TCanvas *canvas = new TCanvas("M", "M");
canvas->Divide(3, 3);

canvas->cd(1);
histMj1j2->Draw();

canvas->cd(2);
histMj1j3->Draw();

canvas->cd(3);
histMj1j4->Draw();
*/
// last 3

/*
TCanvas *canvas = new TCanvas("M", "M");
canvas->Divide(3, 3);


canvas->cd(1);
histMj2j3->Draw();

canvas->cd(2);
histMj2j4->Draw();

canvas->cd(3);
histMj3j4->Draw();
*/
/*
TCanvas *canvas = new TCanvas("M", "M");
histM4j->Draw();
*/
    

   /*
      //.........Plot asymmetry
  
    TCanvas *canvas = new TCanvas("asymetry", "asymetry");
   canvas->Divide(1, 1);
   canvas->cd(1);
   histAsy->Draw(); 
   
 */
 
 
 /*  
//--------ALPHA------------
  
    TCanvas *canvas = new TCanvas("Alpha", "Alpha");
   canvas->Divide(1, 1);
   canvas->cd(1);
   histAlpha->Draw(); 
   
*/

  
//        >>>>>>>>>>>> This part is for the acceptance >>>>>>>>>

  /*  TCanvas *canvas = new TCanvas("Alpha", "Alpha");
   canvas->Divide(1, 1);
   canvas->cd(1);
   histAcc->Draw(); */


//>>>>>>>>>>>ACCEPTANCE ND ALPHA>>>>>>>>>>>
  TCanvas *canvas = new TCanvas("AlphavsAcc", "AlphavsACC");
   canvas->Divide(1, 1);
   canvas->cd(1);
   histAccAlpha->Scale(1/static_cast<Double_t>(totalEvents));
   //histAccAlpha->Draw(); 
  
  
  // ...............................Done with Histograms..................

// ..histAccAlpha is the histogram for Alpha vs acceptance


// Access the bin corresponding to Alpha 

Int_t binIndex25 = histAccAlpha->GetXaxis()->FindBin(desiredAlpha25);
Int_t binIndex27 = histAccAlpha->GetXaxis()->FindBin(desiredAlpha27);
Int_t binIndex29 = histAccAlpha->GetXaxis()->FindBin(desiredAlpha29);
Int_t binIndex31 = histAccAlpha->GetXaxis()->FindBin(desiredAlpha31);
Int_t binIndex33 = histAccAlpha->GetXaxis()->FindBin(desiredAlpha33);

// Get the acceptance value for Alpha
Double_t acceptance25 = histAccAlpha->GetBinContent(binIndex25);
Double_t acceptance27 = histAccAlpha->GetBinContent(binIndex27);
Double_t acceptance29 = histAccAlpha->GetBinContent(binIndex29);
Double_t acceptance31 = histAccAlpha->GetBinContent(binIndex31);
Double_t acceptance33 = histAccAlpha->GetBinContent(binIndex33);


//-------Calculate the significance: significance = signal / sqrt( signal + bacground events)
//Double_t significance = passedAsyCount/std::sqrt(passedAsyCount_bg + passedAsyCount);
Double_t S =  passedAsyCount;
Double_t B =  passedAsyCount_bg;
const double ε = 0.033; //efficiency


Double_t significance = sqrt(2*((B+(ε*S))*log(1 + (ε*S)/B) -(ε*S)));
//Double_t significance = S/std::sqrt(S+B);


  std::cout <<"RESULTS FOR SIGNAL " << std::endl;
  std::cout << "Total SIGNAL events: " << totalEvents << std::endl;
  std::cout << "SignalJets: " << fourjets << std::endl;
  std::cout << "Events passing pT cut: " << passedPtJetCount << std::endl;
  std::cout << "Events passing eta cut: " << passedEtaJetCount << std::endl;
  std::cout << "Events passing R < 2: " << passedRCount << std::endl;
  std::cout << "Events passing Sum of eta < 1.1: " << passedSETACount << std::endl;
  std::cout << "Events passing Asymetry < 0.1: " << passedAsyCount << std::endl;
  
  std::cout <<"RESULTS FOR BACKGROUND:" << std::endl;
    
  std::cout << "Total BG events: " << totalBackgroundEvents << std::endl;
  std::cout << "backgroundjets:" << fourbackgroundjets << std::endl;
  std::cout << "BG Events passing pT cut: " << passedPtBackgroundJetCount << std::endl;
  std::cout << "BG Events passing eta cut: " << passedEtaBackgroundJetCount << std::endl;
  //std::cout << "BG Events passing both pT and eta cuts: " << passedBothBackgroundCuts<< std::endl;
  std::cout << "BG Events passing R < 2: " <<passedRCount_bg << std::endl;
  std::cout << "BG Events passing Sum of eta < 1.1: " << passedSETACount_bg << std::endl;
  std::cout << "BG Events passing Asymetry < 0.1: " << passedAsyCount_bg << std::endl;
  

  
    std::cout <<"ACCEPTANCE FOR SIGNAL:" << std::endl;
    
  //std::cout << "if Alpha = " << desiredAlpha25 <<" Acceptance = " << acceptance25 << std::endl;
  std::cout << "if Alpha = "<< desiredAlpha27 <<" Acceptance = " << acceptance27 << std::endl;
  std::cout << "if Alpha = "<< desiredAlpha29 <<" Acceptance = " << acceptance29 << std::endl;
  std::cout << "if Alpha = "<< desiredAlpha31 <<" Acceptance = " << acceptance31 << std::endl;
  std::cout << "if Alpha = "<< desiredAlpha33 <<" Acceptance = " << acceptance33 << std::endl;

  std::cout <<"THE SIGNIFICANCE:" << std::endl;
  std::cout <<"significance = " << significance << std::endl;
 
  
  
    
 std::cout <<"WELL DONE :)" << std::endl;

  
 



}

//The End Game

// 2023 - 2024 @copyrights reserved
