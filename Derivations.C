#include "classes/DelphesClasses.h"

#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootTreeWriter.h"
#include "external/ExRootAnalysis/ExRootTreeBranch.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "external/ExRootAnalysis/ExRootUtilities.h"



void MyMessage(const char * msg, double num, bool show)
{
    if (show) cout << msg << num << endl;
}

void Derivations(const char *inputFile) 
{
    
  void MyMessage(const char * msg, double num, bool show);
  
  TFile *oldfile = new TFile(inputFile);
  TTree *chain = (TTree*)oldfile->Get("Delphes");

  TClonesArray *branchParticle=0,*branchElectron=0,*branchMuon=0,*branchJet=0;

  chain->SetBranchAddress("Particle",&branchParticle);
  chain->SetBranchAddress("Electron",&branchElectron);
  chain->SetBranchAddress("Muon",&branchMuon);
  chain->SetBranchAddress("Jet",&branchJet);
  
  
  Long64_t allEntries = chain->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;
  
  
  std::string fileout(inputFile);
  fileout.erase(fileout.find_last_of("."), string::npos);
  fileout+=".der.root";
  TFile *newfile = new TFile(fileout.c_str(),"recreate");
  TTree *chaincopy = chain->CloneTree(0);
  TH1D*cutflow=new TH1D("cutflow","Efficiency of cuts/detector effects: ",2,0,2);
 
  TLorentzVector ele_part,muon_part;
  Int_t muon_charge,ele_charge;
  Electron *electron;
  Muon *muon;

  TLorentzVector momentum,elecand,muoncand;

  Bool_t goodEle,goodMuon,debug=false;

  Long64_t entry;

  Int_t i, j, pdgCode;
  
  Int_t nobs;
  std::vector<double> observables; // check always that "hist" has a coherent size!

  // Loop over all events. 
  // Not yet looking for btag...FIXME
  for(entry = 0; entry < allEntries; ++entry){
  //for(entry = 0; entry < 10; ++entry){
      // Load selected branches with data from specified event
      if ((entry%5000)==0) cout << "Analysing entry: "<<entry<<endl;
      chain->GetEntry(entry);
      
      goodEle=false;
      goodMuon=false;
      
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
               chaincopy->Fill();
               branchParticle->Clear();
               branchElectron->Clear();
               branchMuon->Clear();
               branchJet->Clear();
           }
       } // end good events relevant observables
  } // end loop over entries
  
  
  cutflow->Fill(1.0, ((double) chaincopy->GetEntries() ) / allEntries );
  cout << "Writing "<<fileout  <<endl;
  //chaincopy->Print();
  chaincopy->AutoSave();
  cutflow->Write();
  delete oldfile;
  delete newfile;
}
