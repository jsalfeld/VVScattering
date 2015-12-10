double effSF_m_25_medium[10] = {0.893181,0.903896,0.918628,0.970606,0.961617,0.988883,0.981609,0.988133,0.997004,0.980483};
double effSF_e_25_medium[10] = {0.891279,0.897481,0.915062,0.973597,0.925746,1.038128,0.942930,1.031350,0.965295,0.961879};
double effSF_m_25_tight[10]  = {0.890531,0.912268,0.916243,0.974802,0.959355,0.994539,0.978811,0.993689,0.993997,0.981905};
double effSF_e_25_tight[10]  = {0.888321,0.855963,0.905287,0.940666,0.915090,1.028653,0.938059,1.028201,0.957405,0.955777};

double fake_rate_m_25_medium[5][5] = {
0.289,0.199,0.183,0.168,0.164,
0.311,0.218,0.199,0.188,0.191,
0.345,0.260,0.241,0.227,0.225,
0.373,0.299,0.284,0.266,0.261,
0.376,0.304,0.294,0.293,0.289
};
double fake_rate_e_25_medium[5][5] = {
0.233,0.219,0.187,0.202,0.188,
0.236,0.221,0.199,0.165,0.185,
0.293,0.248,0.215,0.191,0.195,
0.281,0.242,0.236,0.256,0.279,
0.284,0.262,0.268,0.291,0.317
};
double fake_rate_m_25_tight[5][5] = {
0.297,0.204,0.188,0.174,0.175,
0.315,0.220,0.201,0.190,0.196,
0.348,0.260,0.241,0.228,0.229,
0.377,0.300,0.284,0.267,0.265,
0.380,0.304,0.295,0.295,0.294
};
double fake_rate_e_25_tight[5][5] = {
0.136,0.111,0.101,0.111,0.097,
0.131,0.118,0.113,0.087,0.101,
0.170,0.144,0.115,0.107,0.102,
0.205,0.177,0.154,0.151,0.154,
0.206,0.165,0.154,0.168,0.184
};

double weightEWKCorr(float pt, int type){
  double parWZ08[2] = { 2.85714,-0.05714};
  double parZZ08[2] = {-4.57143,-0.06857};
  double parWZ14[3] = {3.69800,-0.0726117,0.0000318044};
  double parZZ14[3] = {-0.586985000,-0.099845900,0.0000445083};
  double corrA = 0.0;
  double corrB = 0.0;
  if     (type == 0){ // WZ13
    corrA = (parWZ08[0]+parWZ08[1]*pt)/100.;
    corrB = (parWZ14[0]+parWZ14[1]*pt+parWZ14[2]*pt*pt)/100.;
  }
  else if(type == 1){ // ZZ13
    corrA = (parZZ08[0]+parZZ08[1]*pt)/100.;
    corrB = (parZZ14[0]+parZZ14[1]*pt+parZZ14[2]*pt*pt)/100.;
  }
  double corr = corrB - (corrB-corrA)/6.;

  if(corr >= 0.0) return 1.0;
  return (1.0+corr);
}

char **strsplit(const char* str, const char* delim, size_t* numtokens) {

    // copy the original string so that we don't overwrite parts of it

    // (don't do this if you don't need to keep the old line,

    // as this is less efficient)

    char *s = strdup(str);

    // these three variables are part of a very common idiom to
    // implement a dynamically-growing array

    size_t tokens_alloc = 1;

    size_t tokens_used = 0;

    char **tokens = (char**)calloc(tokens_alloc, sizeof(char*));

    char *token, *strtok_ctx;

    for (token = strtok_r(s, delim, &strtok_ctx);
            token != NULL;
            token = strtok_r(NULL, delim, &strtok_ctx)) {
        // check if we need to allocate more space for tokens
        if (tokens_used == tokens_alloc) {
            tokens_alloc *= 2;
            tokens = (char**)realloc(tokens, tokens_alloc * sizeof(char*));
        }
        tokens[tokens_used++] = strdup(token);
    }

    // cleanup
    if (tokens_used == 0) {
        free(tokens);
        tokens = NULL;
    } else {
        tokens = (char**)realloc(tokens, tokens_used * sizeof(char*));
    }
    *numtokens = tokens_used;
    free(s);
    return tokens;

}

double nPUScaleFactor(TH1D *fhDPU, float npu){
  double mynpu = TMath::Min(npu,(float)39.999);
  Int_t npuxbin = fhDPU->GetXaxis()->FindBin(mynpu);
  return fhDPU->GetBinContent(npuxbin);
}

double ratioFactor(TH1D *fhDVar, float var){
  Int_t nbin = fhDVar->GetXaxis()->FindBin(var);
  return fhDVar->GetBinContent(nbin);
}

double selectIdIsoCut(TString type, int pdgId, double pt, double eta, double iso, int selBits){
  bool isEB = TMath::Abs(eta) < 1.479;
  double isoCut = 0.;
  bool idCut = false;
  if     (TMath::Abs(pdgId) == 13) {
    isoCut = 0.12;
    if     (type == "medium")  idCut = (selBits & BareLeptons::LepMediumIP) == BareLeptons::LepMediumIP;
    else if(type == "tight")   idCut = (selBits & BareLeptons::LepTightIP)  == BareLeptons::LepTightIP;
    else if(type == "default") idCut = (selBits & BareLeptons::LepMediumIP) == BareLeptons::LepMediumIP;
  }
  else if(TMath::Abs(pdgId) == 11) {
    if     (type == "medium")  isoCut = (isEB ? 0.0766 : 0.0678);
    else if(type == "tight")   isoCut = (isEB ? 0.0354 : 0.0646);
    else if(type == "default") isoCut = (isEB ? 0.0354 : 0.0646);
    if     (type == "medium")  idCut = (selBits & BareLeptons::LepMedium) == BareLeptons::LepMedium;
    else if(type == "tight")   idCut = (selBits & BareLeptons::LepTight)  == BareLeptons::LepTight;
    else if(type == "default") idCut = (selBits & BareLeptons::LepTight)  == BareLeptons::LepTight;
  }
  else {
    printf("Problem with selectIsoCut!\n");
    assert(0);
  }
  return (idCut && iso/pt < isoCut);
}

void InitializeJetIdCuts(Float_t fMVACut[4][4])
{
  float cutValues[4][4] = {
    -0.95, -0.96 ,-0.94, -0.95,
    -0.95, -0.96 ,-0.94, -0.95,
    -0.15, -0.26 ,-0.16, -0.16,
    -0.15, -0.26 ,-0.16, -0.16
  };
  
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      fMVACut[i][j] = cutValues[i][j];
    }
  }

}

bool passJetId(Float_t fMVACut[4][4], double mva, double pt, double eta){

  int lPtId = 3;
  if     (pt < 10.)
    lPtId = 0;
  else if(pt < 20.)
    lPtId = 1;
  else if(pt < 30.)
    lPtId = 2;

  int lEtaId = 3;
  if     (eta < 2.50)
    lEtaId = 0;
  else if(eta < 2.75)
    lEtaId = 1;
  else if(eta < 3.00)
    lEtaId = 2;

  if (mva > fMVACut[lPtId][lEtaId])
    return true;
  
  return false;

}

double effScaleFactor(double pt, double eta, int nsel, int period, TString type){
  int iPt = -1;
  if	 (pt < 15) iPt = 0;
  else if(pt < 20) iPt = 1;
  else if(pt < 25) iPt = 2;
  else if(pt < 30) iPt = 3;
  else  	   iPt = 4;

  int iEta = -1;
  if	 (TMath::Abs(eta) < 1.5) iEta = 0;
  else  			 iEta = 1;

  int iPoint = -1;

  if     (iPt==0&&iEta==0) iPoint = 0;
  else if(iPt==0&&iEta==1) iPoint = 1;
  else if(iPt==1&&iEta==0) iPoint = 2;
  else if(iPt==1&&iEta==1) iPoint = 3;
  else if(iPt==2&&iEta==0) iPoint = 4;
  else if(iPt==2&&iEta==1) iPoint = 5;
  else if(iPt==3&&iEta==0) iPoint = 6;
  else if(iPt==3&&iEta==1) iPoint = 7;
  else if(iPt==4&&iEta==0) iPoint = 8;
  else if(iPt==4&&iEta==1) iPoint = 9;
  else assert(0);

  if     (TMath::Abs(nsel) == 13 && period == 1 && (type== "medium" || type== "default"))  return effSF_m_25_medium[iPoint];
  else if(TMath::Abs(nsel) == 11 && period == 1 &&  type== "medium")                       return effSF_e_25_medium[iPoint];
  else if(TMath::Abs(nsel) == 13 && period == 1 &&  type== "tight")                        return effSF_m_25_tight [iPoint];
  else if(TMath::Abs(nsel) == 11 && period == 1 && (type== "tight" || type== "default"))   return effSF_e_25_tight [iPoint];

  assert(0);

  return 0.0;
}

double fakeRateFactor(double pt, double eta, int nsel, int period, TString type){
  int iPt = -1;
  if	 (pt < 15) iPt = 0;
  else if(pt < 20) iPt = 1;
  else if(pt < 25) iPt = 2;
  else if(pt < 30) iPt = 3;
  else  	   iPt = 4;

  int iEta = -1;
  if	 (TMath::Abs(eta) < 0.5) iEta = 0;
  else if(TMath::Abs(eta) < 1.0) iEta = 1;
  else if(TMath::Abs(eta) < 1.5) iEta = 2;
  else if(TMath::Abs(eta) < 2.0) iEta = 3;
  else  			 iEta = 4;

  if     (TMath::Abs(nsel) == 13 && period == 1 && (type== "medium" || type== "default"))  return fake_rate_m_25_medium[iPt][iEta]/(1.0-fake_rate_m_25_medium[iPt][iEta]);
  else if(TMath::Abs(nsel) == 11 && period == 1 &&  type== "medium")			   return fake_rate_e_25_medium[iPt][iEta]/(1.0-fake_rate_e_25_medium[iPt][iEta]);
  else if(TMath::Abs(nsel) == 13 && period == 1 &&  type== "tight")			   return fake_rate_m_25_tight [iPt][iEta]/(1.0-fake_rate_m_25_tight [iPt][iEta]);
  else if(TMath::Abs(nsel) == 11 && period == 1 && (type== "tight" || type== "default"))   return fake_rate_e_25_tight [iPt][iEta]/(1.0-fake_rate_e_25_tight [iPt][iEta]);

  assert(0);

  return 0.0;
}

double weightTruePileupFall15_74X(double ntrue){

  if(ntrue > 50) return 1.0;

  double w[50] = {
      1.0,
      56.8425,
      177.583,
      34.0851,
      16.3021,
      3.35763,
      1.95107,
      2.55352,
      3.43094,
      3.30801,
      3.00642,
      2.6491,
      2.08474,
      1.40075,
      0.790461,
      0.388935,
      0.180224,
      0.0992476,
      0.0706173,
      0.0586083,
      0.0526078,
      0.0504808,
      0.0516041,
      0.0536805,
      0.0577468,
      0.0618779,
      0.0658613,
      0.0714884,
      0.0765815,
      0.0789806,
      0.0737378,
      0.0596836,
      0.0366257,
      0.0227363,
      0.0106934,
      0.00547676,
      0.00289044,
      0.00138031,
      0.000810801,
      0.000408625,
      0.00024114,
      0.000115888,
      7.0553e-05,
      3.67137e-05,
      2.29998e-05,
      1.01618e-05,
      5.60541e-06,
      1.01912e-05,
      1.0745e-05,
      1.0
  };

 return w[(int)floor(ntrue)];

}
