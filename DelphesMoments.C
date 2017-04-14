#include "classes/DelphesClasses.h"

#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootTreeWriter.h"
#include "external/ExRootAnalysis/ExRootTreeBranch.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "external/ExRootAnalysis/ExRootUtilities.h"

void DelphesMoments(){

  double moments[2][4][20],errors[2][4][20];
  double masses[]={165.0,167.0,169.0,171.0,173.0,175.0,177.0,179.0,181.0};

  void GetMoments(const char *inputFileList,double moments[2][4][20],double errors[2][4][20],int masspoint);

  stringstream str_i;
  int startat=0; int stopat=7;
  int numberofpoints=stopat-startat+1;
  for (int i=startat; i<=stopat; i++){
	str_i.str("");
	if (i<10) str_i<<"0";
        str_i<< i;
	GetMoments( ("mt.moments.list."+str_i.str()).c_str(),moments,errors,i-startat);
  }


  TH1::SetDefaultSumw2();
  TCanvas *c1=new TCanvas("c1","Moments");
  
  string functionfit;
  TGraphErrors grMom[2][4];
  for (int lev=0; lev<2; lev++){
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
            str_i<<".lev"<<lev;
            c1->SaveAs((str_i.str()+".pdf").c_str());
    }
  }
}


void GetMoments(const char *inputFileList,double moments[2][4][20],double errors[2][4][20],int masspoint)
{
  TChain *chain= new TChain("Delphes");

  if(!FillChain(chain, inputFileList)) return;

  ExRootResult result;

  Double_t temp,nen;

  TH1::SetDefaultSumw2();
  TH1 *histtemp = result.AddHist1D("histtemp", "negative leptons moments", "negative leptons moments, relevant units", "number of entries", 4, 1.0, 5.0);
  TH1 *hist[2];
  hist[0] = result.AddHist1D("hist", "negative leptons moments, gen", "negative leptons moments, gen, relevant units", "number of entries", 4, 1.0, 5.0);
  hist[1] = result.AddHist1D("hist", "negative leptons moments, reco", "negative leptons moments, reco, relevant units", "number of entries", 4, 1.0, 5.0);

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
        nen=chain->Draw(momdraw.c_str(),cutdraw.c_str());
        
        cutdraw = "TMath::Power(Muon.PT,"+str_i.str()+")*( Electron_size==1 && Muon.Charge>0 && (TMath::Abs(Muon.Eta)<2.4 && Muon.PT>20))";
        momdraw = str_i.str() + ">>+ histtemp"; // add to electrons 
        nen+=chain->Draw(momdraw.c_str(),cutdraw.c_str());
        
        
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

  result.Print(); 
}
