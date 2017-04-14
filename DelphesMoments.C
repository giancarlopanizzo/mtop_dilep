#include "classes/DelphesClasses.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"

void DelphesMoments(){

  double moments[4][20],errors[4][20];
  double masses[]={165.0,167.0,169.0,171.0,173.0,175.0,177.0,179.0,181.0};

  void GetMoments(const char *inputFileList,double moments[4][20],double errors[4][20],int masspoint);

  stringstream str_i;
  int startat=0; int stopat=8;
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
  TGraphErrors grMom[4];
  for (int i=0; i<4; i++){
  	grMom[i]=TGraphErrors(numberofpoints,masses,moments[i],0,errors[i]);
        grMom[i].GetXaxis()->SetTitle("m_{t} (GeV)");
        grMom[i].GetYaxis()->SetTitle("#mu^{1} (GeV)");
        grMom[i].GetYaxis()->SetTitleOffset(1.2);
        grMom[i].SetTitle("");
       


	str_i.str("");
        str_i<< i+1;
	functionfit="[0] * pow(173,"+str_i.str()+") + [1] * pow(x,"+str_i.str()+")";
	TF1 *fFrix = new TF1 ("fFrix", functionfit.c_str(), 165.0, 181.0);
	grMom[i].Fit(fFrix);
        grMom[i].Draw("AP");
         c1->SaveAs((str_i.str()+".pdf").c_str());
  }
}


void GetMoments(const char *inputFileList,double moments[4][20],double errors[4][20],int masspoint)
{
  TChain *chain= new TChain("Delphes");

  if(!FillChain(chain, inputFileList)) return;

  ExRootResult result;

  Double_t temp,nen;

  TH1::SetDefaultSumw2();
  TH1 *histtemp = result.AddHist1D("histtemp", "negative leptons moments", "negative leptons moments, relevant units", "number of entries", 4, 1.0, 5.0);
  TH1 *hist = result.AddHist1D("hist", "negative leptons moments", "negative leptons moments, relevant units", "number of entries", 4, 1.0, 5.0);

  stringstream str_i;
  string momdraw,cutdraw;

  for (int i=1; i<=4; i++){
        str_i.str("");
        str_i<< i;
        momdraw = str_i.str() + ">> histtemp";
        cutdraw = "TMath::Power(Particle.PT,"+str_i.str()+")*( (TMath::Abs(Particle.PID)==11||TMath::Abs(Particle.PID)==13) && (TMath::Abs(Particle.Eta)<2.4 && Particle.PT>20))";

  	nen=chain->Draw(momdraw.c_str(),cutdraw.c_str());
	temp=histtemp->GetBinContent(i)/nen;
	cout <<i<<":"<< temp << " +- ";
	moments[i-1][masspoint]=temp;
	histtemp->SetBinContent(i,temp);
	temp=histtemp->GetBinError(i)/nen;
	errors[i-1][masspoint]=temp;
        cout <<	temp <<	endl;
        histtemp->SetBinError(i,temp);
        hist->Add(histtemp);
  }

  result.Print(); 
}
