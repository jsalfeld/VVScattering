char **strsplit(const char* str, const char* delim, size_t* numtokens);
double nPUScaleFactor(TH1D *fhDPU, float npu);
void InitializeJetIdCuts(Float_t fMVACut[4][4]);
bool passJetId(Float_t fMVACut[4][4], double mva, double pt, double eta);
double effScaleFactor(double pt, double eta, int nsel);
double fakeRateFactor(double pt, double eta, int nsel);
double selectIdIsoCut(TString type, int pdgId, double pt, double eta, double iso, int selBits);

double effSF_m_50[6] = {0.925917,1.008327,0.982077,1.023462,1.004476,0.997671};
double effSF_e_50[6] = {0.882519,0.969906,0.933259,0.999711,0.951164,0.960134};
double effSF_m_25[6] = {0.975344,1.287860,1.037558,1.015660,1.027284,1.026771};
double effSF_e_25[6] = {0.865479,1.056301,1.042451,1.113815,1.041915,0.954874};

double fake_rate_e_50[5][5] = {
0.196,0.206,0.174,0.108,0.030,
0.219,0.207,0.177,0.112,0.045,
0.272,0.254,0.218,0.142,0.096,
0.233,0.226,0.201,0.159,0.150,
0.226,0.243,0.232,0.201,0.198
};
double fake_rate_m_50[5][5] = {
0.316,0.231,0.200,0.183,0.120,
0.356,0.260,0.219,0.189,0.087,
0.383,0.308,0.264,0.274,0.155,
0.446,0.365,0.327,0.320,0.221,
0.456,0.366,0.322,0.333,0.200
};

double fake_rate_e_25[5][5] = {
0.123,0.118,0.101,0.016,0.062,
0.143,0.127,0.099,0.025,0.028,
0.159,0.160,0.111,0.061,0.029,
0.182,0.150,0.128,0.084,0.061,
0.147,0.142,0.127,0.106,0.087
};
double fake_rate_m_25[5][5] = {
0.140,0.110,0.122,0.097,0.030,
0.155,0.116,0.089,0.069,0.009,
0.178,0.133,0.125,0.111,0.020,
0.206,0.150,0.130,0.088,0.031,
0.206,0.179,0.101,0.078,0.286
};

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

double selectIdIsoCut(TString type, int pdgId, double pt, double eta, double iso, int selBits){
  bool isEB = TMath::Abs(eta) < 1.479;
  double isoCut = 0.;
  bool idCut = false;
  if     (TMath::Abs(pdgId) == 13) {
    isoCut = 0.12;
    if     (type == "medium") idCut = (selBits & BareLeptons::LepMediumIP) == BareLeptons::LepMediumIP;
    else if(type == "tight")  idCut = (selBits & BareLeptons::LepTightIP)  == BareLeptons::LepTightIP;
  }
  else if(TMath::Abs(pdgId) == 11) {
    if     (type == "medium") isoCut = (isEB ? 0.0766 : 0.0678);
    else if(type == "tight")  isoCut = (isEB ? 0.0354 : 0.0646);
    if     (type == "medium") idCut = (selBits & BareLeptons::LepMedium) == BareLeptons::LepMedium;
    else if(type == "tight")  idCut = (selBits & BareLeptons::LepTight)  == BareLeptons::LepTight;
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

double effScaleFactor(double pt, double eta, int nsel, int period){
  int iPt = -1;
  if	 (pt < 20) iPt = 0;
  else if(pt < 30) iPt = 1;
  else  	   iPt = 2;

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
  else assert(0);

  if     (TMath::Abs(nsel) == 13 && period == 0) return effSF_m_50[iPoint];
  else if(TMath::Abs(nsel) == 11 && period == 0) return effSF_e_50[iPoint];
  else if(TMath::Abs(nsel) == 13 && period == 1) return effSF_m_25[iPoint];
  else if(TMath::Abs(nsel) == 11 && period == 1) return effSF_e_25[iPoint];

  assert(0);

  return 0.0;
}

double fakeRateFactor(double pt, double eta, int nsel, int period){
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

  if     (TMath::Abs(nsel) == 13) return fake_rate_m_50[iPt][iEta]/(1.0-fake_rate_m_50[iPt][iEta]);
  else if(TMath::Abs(nsel) == 11) return fake_rate_e_50[iPt][iEta]/(1.0-fake_rate_e_50[iPt][iEta]);
  else if(TMath::Abs(nsel) == 11) return fake_rate_e_25[iPt][iEta]/(1.0-fake_rate_e_25[iPt][iEta]);
  else if(TMath::Abs(nsel) == 13) return fake_rate_m_25[iPt][iEta]/(1.0-fake_rate_m_25[iPt][iEta]);

  assert(0);

  return 0.0;
}
