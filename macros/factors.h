char **strsplit(const char* str, const char* delim, size_t* numtokens);
double nPUScaleFactor(TH1D *fhDPU, float npu);
void InitializeJetIdCuts(Float_t fMVACut[4][4]);
bool passJetId(Float_t fMVACut[4][4], double mva, double pt, double eta);
double effScaleFactor(double pt, double eta, int nsel);
double fakeRateFactor(double pt, double eta, int nsel);
double selectIdIsoCut(TString type, int pdgId, double pt, double eta, double iso, int selBits);

double effSF_m_50[6] = {0.906369,0.987096,0.961330,1.001347,0.983206,0.976484};
double effSF_e_50[6] = {0.863793,0.949225,0.913466,0.978506,0.930974,0.939785};
double effSF_m_25[6] = {0.910382,0.948537,0.969327,0.993834,0.990649,0.976285};
double effSF_e_25[6] = {0.903335,0.965685,0.933310,1.049067,0.963457,0.957635};

double fake_rate_e_50[5][5] = {
0.196,0.205,0.173,0.106,0.023,
0.219,0.206,0.176,0.109,0.040,
0.272,0.254,0.218,0.140,0.093,
0.233,0.225,0.201,0.159,0.149,
0.226,0.243,0.232,0.201,0.197
};
double fake_rate_m_50[5][5] = {
0.330,0.230,0.199,0.182,0.112,
0.370,0.260,0.218,0.187,0.079,
0.397,0.307,0.263,0.272,0.147,
0.460,0.364,0.326,0.318,0.214,
0.470,0.365,0.321,0.331,0.186
};

double fake_rate_e_25[5][5] = {
0.236,0.247,0.217,0.192,0.153,
0.256,0.255,0.223,0.200,0.159,
0.308,0.285,0.243,0.209,0.169,
0.262,0.231,0.211,0.194,0.173,
0.250,0.247,0.236,0.221,0.219
};
double fake_rate_m_25[5][5] = {
0.301,0.208,0.196,0.176,0.150,
0.305,0.232,0.214,0.198,0.164,
0.355,0.278,0.256,0.228,0.184,
0.381,0.320,0.300,0.273,0.232,
0.393,0.316,0.313,0.292,0.241
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
