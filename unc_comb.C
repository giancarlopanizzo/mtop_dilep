
void MyMessage(const char * msg, double num, bool show)
{
    if (show) cout << msg << num << endl;
}

void unc_comb(const char *inputFileList) 
{
    
  void MyMessage(const char * msg, double num, bool show);
  
  //adapted from FillChain, ExRootAnalysis
  ifstream infile(inputFileList);
  string buffer;
  TFile*fitres;
  TF1 *fFrix;
  TGraphErrors*gr;
  Int_t nobs=7,nmoms=4;
  Double_t temp[2];
  std::vector<std::vector<double>> parcomb(nobs*nmoms, std::vector<double>(2)); // will contain alpha and beta combined
  std::vector<std::vector<double>> parerrcomb(nobs*nmoms, std::vector<double>(2)); // will contain their error - associated with this particular variation -
  
  for (int i=0; i<2; i++) for (metaobs=0; metaobs<nobs*nmoms; metaobs++) {
      parcomb[metaobs][i]=0;
      parerrcomb[metaobs][i]=0;
  }

  if(!infile.is_open())
  {
    cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << endl;
    return;
  }
   
  Int_t n_members=0;
  while(1)
  {
    infile >> buffer;
    if(!infile.good()) break;
    n_members++;
    fitres=new TFile(buffer.c_str(),"READ"); // new file means new member in the combination
                                             // for each member get parameters of each observable and add to the total
    
    for (metaobs=0; metaobs<nobs*nmoms; metaobs++){
        str_i.str("");
        str_i<<"gr"<<metaobs/nobs<<"_"<<metaobs%nobs+1;
        gr=(TGraphErrors*) fitres->Get(str_i));
        
        str_i.str("");
        str_i<< metaobs%nobs+1;
        fFrix =   (TF1*)   ( gr->GetFunction( ("fFrix"+str_i.str() ).c_str() ) );
        fFrixGetParameters(temp);
        for (int i=0; i<2; i++) {
            parcomb[metaobs][i]+=temp[i];
            parerrcomb[metaobs][i]+=temp[i]*temp[i];
        }
    }
    delete fitres;
  }// end of reading input files
 infile.close();
 
 for (int i=0; i<2; i++) for (metaobs=0; metaobs<nobs*nmoms; metaobs++) {
      parcomb[metaobs][i]/=n_members;
      parerrcomb[metaobs][i]/=n_members;
      parerrcomb[metaobs][i]-=parcomb[metaobs][i]*parcomb[metaobs][i]; // variance= <x^2>-<x>^2
      if (parerrcomb[metaobs][i]>0) parerrcomb[metaobs][i]=TMath::Sqrt(parerrcomb[metaobs][i]);
      else cout << "Error in variance computation, please check your algorithm."<< endl;
  } 
   
  ofstream outfile("combresults.log"); 
  // now print, for convenience, Delta( alpha ) / beta
  outfile << "Delta(alpha) / beta :"<< endl;
  for (metaobs=0; metaobs<nobs*nmoms; metaobs++) cout << << parerrcomb[metaobs][0]/parcomb[metaobs][1]<< endl;
  outfile.close();
}
