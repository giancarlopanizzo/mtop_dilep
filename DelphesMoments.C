#include "classes/DelphesClasses.h"

#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootTreeWriter.h"
#include "external/ExRootAnalysis/ExRootTreeBranch.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "external/ExRootAnalysis/ExRootUtilities.h"

void DelphesMoments(Int_t combine,Int_t lhaid, Int_t fixedmasspoint){ // combine: 0-> all in one, 1->only one mass point and exit, 2-> combine previous
  
  double moments[9][4][20],errors[9][4][20];
  double masses[]={165.0,167.0,169.0,171.0,173.0,175.0,177.0,179.0,181.0};

  void GetMoments(const char *inputFileList,double moments[9][4][20],double errors[9][4][20],int masspoint);
  void GetMomentsFromFile(const char *inputFile,double moments[9][4][20],double errors[9][4][20],int masspoint);
  void PrintCorrMatrix(const char *inputFile);

  stringstream str_i;
  int startat=0; int stopat=0;
  if (combine==1) {startat=fixedmasspoint; stopat=fixedmasspoint;}
  else {startat=0; stopat=7;}

  int numberofpoints=stopat-startat+1;
  for (int i=startat; i<=stopat; i++){
	str_i.str("");
	str_i << lhaid <<".list.";
	if (i<10) str_i<<"0";
        str_i<< i;
	  if (combine!=2){
		GetMoments( ("mt.moments."+str_i.str()).c_str(),moments,errors,i-startat);
	  }
 	 else {  // assume combine 2, so just collect results:open relevant file and pick what needed
		GetMomentsFromFile( ("mt.moments."+str_i.str()+".result").c_str(),moments,errors,i-startat);
  	 }
  }

  if (combine==1) {// in this case str_i still contains correct value 
      PrintCorrMatrix(  ("mt.moments."+str_i.str()+".result").c_str() );
      return;
  }
  TH1::SetDefaultSumw2();
  TCanvas *c1=new TCanvas("c1","Moments");
  
  string functionfit;
  TGraphErrors grMom[9][4];
  for (int lev=0; lev<9; lev++){
    for (int i=0; i<4; i++){
            grMom[lev][i]=TGraphErrors(numberofpoints,masses,moments[lev][i],0,errors[lev][i]);
            grMom[lev][i].GetXaxis()->SetTitle("m_{t} (GeV)");
            grMom[lev][i].GetYaxis()->SetTitle("#mu^{1} (GeV)");
            grMom[lev][i].GetYaxis()->SetTitleOffset(1.2);
            grMom[lev][i].SetTitle("");

            str_i.str("");
            str_i<< i+1;
            functionfit="[0] * pow(173,"+str_i.str()+") + [1] * pow(x,"+str_i.str()+")";
            TF1 *fFrix = new TF1 ("fFrix", functionfit.c_str(), 165.0, 181.0);
            grMom[lev][i].Fit(fFrix);
            grMom[lev][i].Draw("AP");
            str_i<<"."<<lhaid<<".lev"<<lev;
            c1->SaveAs((str_i.str()+".png").c_str());
    }
  }
}


void GetMoments(const char *inputFileList,double moments[9][4][20],double errors[9][4][20],int masspoint)
{
    
  void MyMessage(const char * msg, double num, bool show);
  TChain *chain= new TChain("Delphes");

  if(!FillChain(chain, inputFileList)) return;
  // switch off reading unuseful branches
  chain->SetBranchStatus("*",0); //disable all branches
  chain->SetBranchStatus("Particle*",1);
  chain->SetBranchStatus("Electron*",1);  
  chain->SetBranchStatus("Muon*",1);  

  ExRootResult result;
  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");

  Double_t delta,mean,nen;

  TH1::SetDefaultSumw2();  
  TH1D *hist[9];
  TH2D *corrmatr=(TH2D*)result.AddHist2D("corrmatr","Correlation matrix entries","i","j",9*4,0.0,9.0*4, 9*4,0.0,9.0*4);
  
  stringstream str_i;
  for (int iobs=0; iobs<9; iobs++){
    str_i.str("hist_");
    str_i << iobs;
    hist[iobs] = (TH1D*) result.AddHist1D(str_i.str().c_str(), ( "negative leptons moments, "+str_i.str() ).c_str(), "negative leptons moments, relevant units", "number of entries", 4, 1.0, 5.0);
  }

  // new way, taken from Example3
  
  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

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
  for(entry = 0; entry < ( (debug) ? 10 : allEntries); ++entry){
      // Load selected branches with data from specified event
      if (((entry+1)%5000)==0) cout << "Analysing entry: "<<entry<<endl;      
      treeReader->ReadEntry(entry);
      
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
               // 6) | y(l+) -  y(l-) | 
               // 7) pt of difference l+ - l-
               // 8) M of difference l+ - l-
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
               observables.push_back(elecand.Eta()-muoncand.Eta()); // 6
               momentum=elecand-muoncand;
               observables.push_back(momentum.Pt()); // 7
               observables.push_back(momentum.M()); // 8
               
               // compute moments of pt 
                    // and other observables
               nobs=(Int_t) observables.size();
               for (int iobs=0; iobs<nobs; iobs++){
                   for (int imom=1; imom<=4; imom++){ // here UNnormalised moments! remember somewhere to normalise them!
                       hist[iobs]->Fill(imom,TMath::Power(observables[iobs],imom)) ;
                       MyMessage("#  Filling observable ",iobs,debug);
                       MyMessage("#  With content ",TMath::Power(observables[iobs],imom),debug);
                       for (int jobs=iobs; jobs<nobs; jobs++) for (int jmom=imom; jmom<=4; jmom++) corrmatr->Fill(iobs+imom-1,jobs+jmom-1,TMath::Power(observables[iobs],imom)*TMath::Power(observables[jobs],jmom));
                   }
               }
           }
       } // end good events relevant observables
  } // end loop over entries
  
  // Now normalise histograms and save moments 
  nen = (double) ( hist[0]->GetEntries() / 4); // 4 moments, entries are 4* goodevents. Assuming all histograms have same nen ...
  MyMessage("Number of entries: ",nen,debug);
  corrmatr->Scale(1.0/nen);
  for (int iobs=0; iobs<nobs; iobs++){ // nobs has the last value set, hopefully not buggy
      
      hist[iobs]->Scale(1.0/nen); // this should scale both contents and errors
      
      MyMessage("** Computing moments for observable ",iobs,true);
      for (int imom=1; imom<=4; imom++){// useful to leave last one as first moments for correlation matrix computation 
          mean=hist[iobs]->GetBinContent(imom);
          moments[iobs][imom-1][masspoint]=mean;
          cout <<imom<<":"<< mean << " +- ";

	  // Now correct the errors: they now are sqrt(  <x^2>  /N )
	  // They should be sqrt(  ( <x^2>  - <x>^2 ) /N )
          delta=hist[iobs]->GetBinError(imom);
	  delta=TMath::Sqrt( delta*delta - mean*mean/nen );
	  errors[iobs][imom-1][masspoint]=delta;
          hist[iobs]->SetBinError(imom,delta); // necessary to allow errors to be read by following copies
          cout << delta << endl;
      
          // now correlation matrix. "delta" contains sigma/sqrt(nen) of iobs imom, hopefully
          for (int jobs=0; jobs<=iobs; jobs++) for (int jmom=0; jmom<=imom; jmom++) { // refine corrmatr step by step
              mean=corrmatr->GetBinContent(jobs+jmom,iobs+imom);
              // MyMessage("Corr matrix was: ",mean,debug);
              mean-=hist[iobs]->GetBinContent(imom)*hist[jobs]->GetBinContent(jmom); // now mean is <iobs jobs> - <iobs><jobs>
              MyMessage("Corr matrix is: ",mean,debug);
              mean/=nen*(delta*errors[jobs][jmom-1][masspoint]); // we are computing the corrmatrix of -this- mass point . . .
              corrmatr->SetBinContent(jobs+jmom,iobs+imom,mean); // done! 
          }
      }
  }
  
  
  // old way, through TTreeDraw. 
/*
  stringstream str_i;
  string momdraw,cutdraw;
  
  for (int i=1; i<=4; i++){
        str_i.str("");
        str_i<< i;
        momdraw = str_i.str() + ">> histtemp";
        // gen level:
        cutdraw = "TMath::Power(Particle.PT,"+str_i.str()+")*( ( Particle.PID==11 || Particle.PID==13) && Particle.Status==1&& (TMath::Abs(Particle.Eta)<2.4 && Particle.PT>20))";

  	nen=chain->Draw(momdraw.c_str(),cutdraw.c_str());
	temp=histtemp->GetBinContent(i)/nen;
	cout <<i<<":"<< temp << " +- ";
	moments[0][i-1][masspoint]=temp;
	histtemp->SetBinContent(i,temp);
	temp=histtemp->GetBinError(i)/nen;
	errors[0][i-1][masspoint]=temp;
        cout <<	temp <<	endl;
        histtemp->SetBinError(i,temp);
        hist[0]->Add(histtemp);
        
        // reco level:electrons
        cutdraw = "TMath::Power(Electron.PT,"+str_i.str()+")*( Muon_size==1 && Electron.Charge>0 && (TMath::Abs(Electron.Eta)<2.4 && Electron.PT>20))";
        nen=(double) chain->Draw(momdraw.c_str(),cutdraw.c_str());
        
        cutdraw = "TMath::Power(Muon.PT,"+str_i.str()+")*( Electron_size==1 && Muon.Charge>0 && (TMath::Abs(Muon.Eta)<2.4 && Muon.PT>20))";
        momdraw = str_i.str() + ">>+ histtemp"; // add to electrons 
        nen+=(double) chain->Draw(momdraw.c_str(),cutdraw.c_str());
        
        
	temp=histtemp->GetBinContent(i)/nen;
	cout <<i<<":"<< temp << " +- ";
	moments[1][i-1][masspoint]=temp;
	histtemp->SetBinContent(i,temp);
	temp=histtemp->GetBinError(i)/nen;
	errors[1][i-1][masspoint]=temp;
        cout <<	temp <<	endl;
        histtemp->SetBinError(i,temp);
        hist[1]->Add(histtemp);
  
  }
*/
  std::string fileout(inputFileList);
  fileout+=".result";
  cout << "Writing "<<fileout  <<endl;
  // result.Write(fileout.str().c_str());
  result.Write(fileout.c_str());
  delete treeReader;
  delete chain;
}

void GetMomentsFromFile(const char *inputFile,double moments[9][4][20],double errors[9][4][20], int masspoint)
{
  TFile*f1=new TFile(inputFile,"READ");
  TH1D *hist[9];
  
  stringstream str_i;
  
  for (int j=0; j<9; j++){
      str_i.str("hist_");
      str_i << j;
      hist[j] = (TH1D*) f1->Get( str_i.str().c_str() );
	  for (int i=1; i<=4; i++){
	       moments[j][i-1][masspoint]=hist[j]->GetBinContent(i);
	       errors[j][i-1][masspoint]=hist[j]->GetBinError(i);
	  }
  }
  delete f1;
}

void PrintCorrMatrix(const char *inputFile)
{
  TFile*f1=new TFile(inputFile,"READ");
  TH2D *corrmatr=(TH2D*)f1->Get("corrmatr");
  
  cout << "Correlation matrix for file "<<inputFile<<endl;
  for (int j=0; j<9*4; j++){
	  for (int i=0; i<9*4; i++){
		printf("%5.3f ",corrmatr->GetBinContent(j+1,i+1) );
	  }
	  cout << endl;
  }
  delete f1;
}

void MyMessage(const char * msg, double num, bool show)
{
    if (show) cout << msg << num << endl;
}
