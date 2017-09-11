#include "classes/DelphesClasses.h"

#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootTreeWriter.h"
#include "external/ExRootAnalysis/ExRootTreeBranch.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "external/ExRootAnalysis/ExRootUtilities.h"


void DrawHists(const char *inputFileList)
{
    
  void MyMessage(const char * msg, double num, bool show);
  TChain *chain= new TChain("Delphes");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);  
  if(!FillChain(chain, inputFileList)) return;
  // switch off reading unuseful branches
  chain->SetBranchStatus("*",0); //disable all branches
  chain->SetBranchStatus("Particle*",1);
  chain->SetBranchStatus("Electron*",1);  
  chain->SetBranchStatus("Muon*",1);  
  chain->SetBranchStatus("Jet*",1);  

  ExRootResult result;
  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");

  Double_t delta,mean,nen;
  TCanvas*c1=new TCanvas("c1","c1",400*3,300*2);
  c1->Divide(3,2);
  TH1::SetDefaultSumw2();  
  TH1D *hist[7];
  for (int i=0; i<7; i++) hist[i]=NULL;


              // 0) Truth level positive charged lepton pt
               // From here on only reco level quantities:
               // 1) positive charged lepton pt
               // 2) pt of the l+l- system
               // 3) M(l+l-)
               // 4) E(l+) + E(l-)
               // 5) pt(p+) + pt(l-)
               // From here on new observables:
               // # skip 6) | y(l+) -  y(l-) | 
               // 7) pt of difference l+ - l-
  
    hist[0] = new TH1D("hist_0", "Truth level positive charged lepton pt",  50, 15, 200.0);
    hist[1] = new TH1D("hist_1", "positive charged lepton pt", 50, 15.0, 200.0);
    hist[2] = new TH1D("hist_2", "pt of the l+l- system", 50, 0.0, 200.0);
    hist[3] = new TH1D("hist_3", "M(l+l-)", 50, 0.0, 250.0);
    hist[4] = new TH1D("hist_4", "E(l+) + E(l-)", 50, 40, 500.0);
    hist[5] = new TH1D("hist_5", "pt(p+) + pt(l-)", 50, 40, 400.0);
    hist[6] = new TH1D("hist_6", "pt of difference l+ - l-", 50, 0, 250.0);

    string xTitles[7]={"p_{T}(l^{+}) (GeV)",
                       "p_{T}(l^{+}) (GeV)",
                       "p_{T}(l^{+}+l^{-}) (GeV)",
                       "m(l^{+}l^{-}) (GeV)",
                       "E(l^{+}) + E(l^{-}) (GeV)",
                       "p_{T}(l^{+}) + pt_{T}(l^{-}) (GeV)",
                       "p_{T}(l^{+}-l^{-}) (GeV)"};

  // new way, taken from Example3
  
  Long64_t allEntries = treeReader->GetEntries(), goodEntries=0;

  
  cout << "** Chain contains " << allEntries << " events" << endl;

  TLorentzVector ele_part,muon_part;
  Int_t muon_charge,ele_charge;
  Electron *electron;
  Muon *muon;
  Jet *jet;

  TLorentzVector momentum,elecand,muoncand;

  Bool_t goodEle,goodMuon,debug=false;
 

  Long64_t entry;

  Int_t i, j,jmom,jobs, pdgCode,goodJets,goodBJets;
  
  Int_t nobs;
  std::vector<double> observables; // check always that "hist" has a coherent size!

  // Loop over all events. 
  // Not yet looking for btag...FIXME
  for(entry = 0; entry < ( (debug) ? 10 : allEntries); ++entry){
      // Load selected branches with data from specified event
      if (((entry+1)%5000)==0) cout << "Analysing entry: "<<entry<<endl;      
      treeReader->ReadEntry(entry);
      
      goodJets=0; 
      goodBJets=0; 
      goodEle=false;
      goodMuon=false;

      // Loop over all jets in the event
      for(i = 0; i < branchJet->GetEntriesFast(); ++i) // take electron with highest pt This doesn't necessarily give the correct lepton FIXME
      {
          jet = (Jet*) branchJet->At(i);
          if (jet->PT > 30 && TMath::Abs(jet->Eta) < 2.4 ) {
              goodJets++;
              if (jet->BTag>0) goodBJets++;
          }
      }
    
      if (goodJets<2 || goodBJets<1) continue;

      // Loop over all electrons in event
      elecand.SetPtEtaPhiM(0,0,0,0);
      for(i = 0; i < branchElectron->GetEntriesFast(); ++i) // take electron with highest pt This doesn't necessarily give the correct lepton FIXME
      {
          electron = (Electron*) branchElectron->At(i);
          if (electron->PT > elecand.Pt()) {
              elecand=electron->P4();
              ele_charge=electron->Charge;
              ele_part = ((GenParticle*) electron->Particle.GetObject())->P4();
              goodEle=( elecand.Pt() > 20.0 && TMath::Abs(elecand.Eta() ) < 2.4);
              MyMessage("## Found electron candidate with Pt",elecand.Pt(),debug );
              MyMessage("#### Eta",elecand.Eta(),debug );
              MyMessage("#### Charge",ele_charge,debug );
          }
      }
      
      // Loop over all muons in event, taking the hardest one. This doesn't necessarily give the correct lepton FIXME
      muoncand.SetPtEtaPhiM(0,0,0,0);
      for(i = 0; i < branchMuon->GetEntriesFast(); ++i) {
          muon = (Muon*) branchMuon->At(i);
          if (muon->PT > muoncand.Pt()) {
              muoncand=muon->P4();
              muon_charge=muon->Charge;
              muon_part = ( (GenParticle*) muon->Particle.GetObject() )->P4();
              goodMuon=( muoncand.Pt() > 20.0 && TMath::Abs(muoncand.Eta() ) < 2.4);
              MyMessage("## Found muon candidate with Pt",muoncand.Pt(),debug );
              MyMessage("#### Eta",muoncand.Eta(),debug );
              MyMessage("#### Charge",muon_charge,debug );
          }
       }
       
       if (!(goodEle && goodMuon)) continue;
       else { // we have two kinematically allowed leptons. 
           if ( (ele_charge * muon_charge) > 0) continue;
           else { // two good, opposite sign leptons: go on and save all relevant observables
               MyMessage("Found two good leptons!",1,debug);
               goodEntries++;
               //save positive lepton PT
                    // save other observables ...
               
               // 0) Truth level positive charged lepton pt
               // From here on only reco level quantities:
               // 1) positive charged lepton pt
               // 2) pt of the l+l- system
               // 3) M(l+l-)
               // 4) E(l+) + E(l-)
               // 5) pt(p+) + pt(l-)
               // From here on new observables:
               // # skip 6) | y(l+) -  y(l-) | 
               // 7) pt of difference l+ - l-
               // # skip 8) M of difference l+ - l-
               observables.clear();
               
               if (ele_charge > 0) {
                   observables.push_back(ele_part.Pt());
                   observables.push_back(elecand.Pt());
               }
               else {
                   observables.push_back(muon_part.Pt());
                   observables.push_back(muoncand.Pt());
               }
               momentum=elecand+muoncand;
               observables.push_back(momentum.Pt()); // 2
               observables.push_back(momentum.M() ); // 3
               observables.push_back(elecand.E()+muoncand.E()); // 4 
               observables.push_back(elecand.Pt()+muoncand.Pt()); // 5
               // observables.push_back(elecand.Eta()-muoncand.Eta()); // 6
               momentum=elecand-muoncand;
               observables.push_back(momentum.Pt()); // 7
               // observables.push_back(momentum.M()); // 8
               
               // compute moments of pt 
                    // and other observables
               nobs=(Int_t) observables.size();
               for (int iobs=0; iobs<nobs; iobs++){
                       hist[iobs]->Fill(observables[iobs]) ;
                       MyMessage("#  Filling observable ",iobs,debug);
               }
           }
       } // end good events relevant observables
  } // end loop over entries
  
  double acceptance=((double) goodEntries)/allEntries;
  // Now normalise histograms and save moments 
  MyMessage("good/all ",acceptance,true);
  for (int iobs=1; iobs<7; iobs++){ // nobs has the last value set, hopefully not buggy
      if (!hist[iobs]) {
        MyMessage("Failing with histo",iobs,true);
        continue;
      }
      nen =  hist[iobs]->Integral(1,hist[iobs]->GetNbinsX()  ) ; 
      hist[iobs]->Scale(1.0/nen); // this should scale both contents and errors 
      c1->cd(iobs);
      hist[iobs]->GetXaxis()->SetTitle(xTitles[iobs].c_str());
      hist[iobs]->GetXaxis()->SetTitleSize(0.06);
      hist[iobs]->GetYaxis()->SetTitleSize(0.06);
      hist[iobs]->GetYaxis()->SetTitle("Norm. entries");
      hist[iobs]->Draw("hist");
      hist[iobs]->SetDirectory(gROOT);
  }
} 


void MyMessage(const char * msg, double num, bool show)
{
    if (show) cout << msg << num << endl;
}

