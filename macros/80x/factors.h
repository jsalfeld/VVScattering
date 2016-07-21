double effSF_m_25_medium[10] = {0.836574,0.829962,0.877106,0.892587,0.935558,0.949318,0.946449,0.943031,0.973128,0.935147};
double effSF_e_25_medium[10] = {0.739427,0.690631,0.815957,0.822430,0.856509,0.859530,0.902980,0.923910,0.976901,0.942210};
double effSF_e_25_tight[10]  = {0.742330,0.661762,0.812364,0.820798,0.851298,0.846768,0.888956,0.915187,0.965834,0.927740};

double fake_rate_m_25_medium[5][5] = {
0.400,0.293,0.271,0.252,0.198,
0.406,0.307,0.285,0.262,0.221,
0.449,0.354,0.334,0.307,0.290,
0.490,0.409,0.397,0.371,0.330,
0.506,0.444,0.426,0.415,0.380
};
double fake_rate_e_25_medium[5][5] = {
0.205,0.184,0.178,0.167,0.198,
0.199,0.193,0.178,0.164,0.143,
0.217,0.212,0.186,0.195,0.168,
0.331,0.268,0.262,0.255,0.288,
0.359,0.340,0.348,0.351,0.381
};
double fake_rate_e_25_tight[5][5] = {
0.121,0.103,0.092,0.073,0.106,
0.113,0.105,0.094,0.077,0.048,
0.137,0.121,0.101,0.099,0.077,
0.232,0.177,0.163,0.139,0.131,
0.225,0.188,0.179,0.163,0.181
};
double fake_rate_e_25_medium_mva[5][5] = {
0.684,0.587,0.552,0.519,0.453,
0.517,0.465,0.430,0.408,0.369,
0.307,0.291,0.290,0.275,0.291,
0.329,0.321,0.354,0.372,0.378,
0.351,0.319,0.317,0.322,0.318
};
double fake_rate_e_25_tight_mva[5][5] = {
0.519,0.451,0.425,0.380,0.321,
0.386,0.332,0.306,0.296,0.245,
0.173,0.177,0.172,0.162,0.164,
0.204,0.216,0.235,0.240,0.223,
0.224,0.193,0.183,0.183,0.177
};

double jetEpsBtagB[5][5] = {
0.5587,0.6607,0.6934,0.7238,0.7229,
0.5645,0.6653,0.7002,0.7265,0.7322,
0.5047,0.6125,0.6523,0.6809,0.6797,
0.4743,0.5802,0.6267,0.6528,0.6405,
0.1597,0.2316,0.2773,0.3092,0.3507
};
double jetEpsBtagC[5][5] = { 
0.1653,0.2069,0.1983,0.2045,0.1981,
0.1692,0.2107,0.2046,0.2140,0.2124,
0.1426,0.1871,0.1832,0.1917,0.1911,
0.1313,0.1720,0.1781,0.1881,0.1903,
0.0372,0.0515,0.0579,0.0643,0.0784
};
double jetEpsBtagL[5][5] = {
0.0237,0.0327,0.0281,0.0266,0.0245,
0.0258,0.0356,0.0301,0.0288,0.0266,
0.0239,0.0336,0.0305,0.0283,0.0264,
0.0262,0.0366,0.0338,0.0308,0.0304,
0.0048,0.0110,0.0131,0.0131,0.0155
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
  double mynpu = TMath::Min(npu,(float)50.999);
  Int_t npuxbin = fhDPU->GetXaxis()->FindBin(mynpu);
  return fhDPU->GetBinContent(npuxbin);
}

double ratioFactor(TH1D *fhDVar, float var){
  Int_t nbin = fhDVar->GetXaxis()->FindBin(var);
  return fhDVar->GetBinContent(nbin);
}

double zpt_correction(double ptll, int type){
  
  if     (type == 0) return (0.95 - 0.1*TMath::Erf((ptll-14)/8.8));
  else if(type == 1) return (1.02852 - 0.0949640*TMath::Erf((ptll-19.0422)/10.4487) + 0.0758834*TMath::Erf((ptll-56.1146)/41.1653));

  return 1.0;
}

double selectIdIsoCut(TString type, int pdgId, double pt, double eta, double iso, int selBits, double mva){
  bool isEB = TMath::Abs(eta) < 1.479;
  double isoCut = 0.;
  bool idCut = false;
  if     (TMath::Abs(pdgId) == 13) {
    if (type=="loose" || type=="veto") isoCut = 0.25;
    else isoCut=0.15;
    if     (type == "default" || type == "medium" || type == "default_mva" || type == "medium_mva") idCut = (selBits & BareLeptons::LepTightIP) == BareLeptons::LepTightIP;
    else if(type == "loose") idCut= (selBits & BareLeptons::LepLoose) == BareLeptons::LepLoose;
    else if(type == "veto")  idCut= (selBits & BareLeptons::LepLoose) == BareLeptons::LepLoose;
    else printf("Problem with selectIsoCut!\n");

    return (idCut && iso/pt < isoCut);
  }
  else if(TMath::Abs(pdgId) == 11 && (type == "default" || type == "medium" || type == "veto" || type == "loose")) {
    if     (type == "medium")  isoCut = (isEB ? 0.0766 : 0.0678);
    else if(type == "default") isoCut = (isEB ? 0.0354 : 0.0646);
    else if(type == "loose")   isoCut = (isEB ? 0.0893 : 0.121 ); 
    else if(type == "veto")    isoCut = (isEB ? 0.126  : 0.144 ); 

    if     (type == "medium")  idCut = (selBits & BareLeptons::LepMedium) == BareLeptons::LepMedium;
    else if(type == "default") idCut = (selBits & BareLeptons::LepTight)  == BareLeptons::LepTight;
    else if(type == "loose"  ) idCut = (selBits & BareLeptons::LepLoose)  == BareLeptons::LepLoose;
    else if(type == "veto"   ) idCut = (selBits & BareLeptons::LepVeto)   == BareLeptons::LepVeto;

    return (idCut && iso/pt < isoCut);
  }
  else if(TMath::Abs(pdgId) == 11 && (type == "default_mva" || type == "medium_mva")) {
    idCut = (selBits & BareLeptons::LepFake) == BareLeptons::LepFake;
    double mvaCut = 3.0; 
    if     (type == "medium_mva"  && TMath::Abs(eta) < 0.800) mvaCut = 0.972153;
    else if(type == "medium_mva"  && TMath::Abs(eta) < 1.479) mvaCut = 0.922126;
    else if(type == "medium_mva")                             mvaCut = 0.610764;
    else if(type == "default_mva" && TMath::Abs(eta) < 0.800) mvaCut = 0.988153;
    else if(type == "default_mva" && TMath::Abs(eta) < 1.479) mvaCut = 0.967910;
    else if(type == "default_mva")                            mvaCut = 0.841729;
    else printf("Problem with selectMVACut!\n");

    return (idCut && mva > mvaCut);
  }
  else {
    printf("Problem with selectIsoCut!\n");
    assert(0);
  }
  return false;
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

double effhDScaleFactor(double pt, double eta, int nsel, TString type, TH2D *fhDMuMediumSF, TH2D *fhDElMediumSF, TH2D *fhDElTightSF, TH1D *fhDTrkSF){
  if(pt>=100) pt = +99.999;
  if     (eta>=+2.4) eta = +2.399;
  else if(eta<=-2.4) eta = -2.399;
  double trkSF = 1.0;
  if(TMath::Abs(nsel) == 13){
    Int_t binXT = fhDTrkSF->GetXaxis()->FindFixBin(eta);
    trkSF = fhDTrkSF->GetBinContent(binXT);
    if(trkSF <= 0) printf("trkSF <= 0! %f %d- %f %f\n",trkSF,binXT,pt,eta);
  }

  eta = abs(eta);
  Int_t binXA = 0;
  Int_t binYA = 0;
  Int_t binXB = 0;
  Int_t binYB = 0;

  if     (TMath::Abs(nsel) == 13 && (type== "medium" || type== "default"))  {binXA = fhDMuMediumSF->GetXaxis()->FindFixBin(eta);binYA = fhDMuMediumSF->GetYaxis()->FindFixBin(pt);}
  else if(TMath::Abs(nsel) == 11 && type== "medium")			    {binXA = fhDElMediumSF->GetXaxis()->FindFixBin(eta);binYA = fhDElMediumSF->GetYaxis()->FindFixBin(pt);}
  else if(TMath::Abs(nsel) == 11 && type== "default")			    {binXA = fhDElTightSF ->GetXaxis()->FindFixBin(eta);binYA = fhDElTightSF ->GetYaxis()->FindFixBin(pt);}
  else    printf("PROBLEM WITH BINS\n");

  double result = 0.0;
  if     (TMath::Abs(nsel) == 13 && (type== "medium" || type== "default")) result = fhDMuMediumSF->GetBinContent(binXA, binYA);
  else if(TMath::Abs(nsel) == 11 && type== "medium")	                   result = fhDElMediumSF->GetBinContent(binXA, binYA);
  else if(TMath::Abs(nsel) == 11 && type== "default")	                   result = fhDElTightSF ->GetBinContent(binXA, binYA);
  
  if(result <= 0) printf("Result <= 0! %f %d %d %d - %f %f\n",result,nsel,binXA,binYA,pt,eta);
  
  return result*trkSF;
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

  if     (TMath::Abs(nsel) == 13 && period == 1 && (type== "medium" || type== "default")) return effSF_m_25_medium[iPoint];
  else if(TMath::Abs(nsel) == 11 && period == 1 &&  type== "medium")                      return effSF_e_25_medium[iPoint];
  else if(TMath::Abs(nsel) == 11 && period == 1 &&  type== "default")                     return effSF_e_25_tight [iPoint];

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

  if     (TMath::Abs(nsel) == 13 && period == 1 && (type== "medium" || type== "default" || type== "medium_mva" || type== "default_mva")) return fake_rate_m_25_medium    [iPt][iEta]/(1.0-fake_rate_m_25_medium    [iPt][iEta]);
  else if(TMath::Abs(nsel) == 11 && period == 1 &&  type== "medium")			  return fake_rate_e_25_medium    [iPt][iEta]/(1.0-fake_rate_e_25_medium    [iPt][iEta]);
  else if(TMath::Abs(nsel) == 11 && period == 1 &&  type== "default")			  return fake_rate_e_25_tight     [iPt][iEta]/(1.0-fake_rate_e_25_tight     [iPt][iEta]);
  else if(TMath::Abs(nsel) == 11 && period == 1 &&  type== "medium_mva")		  return fake_rate_e_25_medium_mva[iPt][iEta]/(1.0-fake_rate_e_25_medium_mva[iPt][iEta]);
  else if(TMath::Abs(nsel) == 11 && period == 1 &&  type== "default_mva")		  return fake_rate_e_25_tight_mva [iPt][iEta]/(1.0-fake_rate_e_25_tight_mva [iPt][iEta]);
  else    printf("PROBLEM WITH FAKES\n");

  assert(0);

  return 0.0;
}

double weightTruePileupFall15_74X(double ntrue){

  if(ntrue > 50) return 1.0;

  double w[50] = {
126.337,
153.753,
108.762,
32.4929,
17.9949,
3.34053,
1.97929,
2.5719,
3.42167,
3.33126,
3.00845,
2.65362,
2.08864,
1.39939,
0.792047,
0.386897,
0.18109,
0.0990121,
0.0705014,
0.0587098,
0.0525807,
0.0503863,
0.051229,
0.0538919,
0.0574271,
0.0613761,
0.0657269,
0.0706353,
0.0757276,
0.078729,
0.0744712,
0.0591124,
0.038286,
0.02139,
0.0111463,
0.0057112,
0.00295235,
0.00155243,
0.000828797,
0.000446272,
0.000240566,
0.000129006,
6.85054e-05,
3.59166e-05,
1.85598e-05,
9.44456e-06,
1.16979e-05,
8.50725e-06,
1.49492e-05,
1
  };

 return w[(int)floor(ntrue)];

}


double weightTruePileupFall15_74X_wisconsin(double ntrue){

  if(ntrue > 50) return 1.0;

  double w[50] = {
58.1138048875,
93.2612967205,
98.6162081996,
28.7104534509,
15.4617381103,
 2.7153649567,
 1.4508190201,
 1.7113554671,
 2.3002683151,
 2.4788907956,
 2.5634937262,
 2.6060767979,
 2.3623139729,
 1.8314456721,
 1.1939779551,
 0.6482593491,
 0.3000319858,
 0.1272652275,
 0.0578769686,
 0.0312937347,
 0.0179469656,
 0.0093677189,
 0.0041699519,
 0.0015958001,
 0.0005656543,
 0.0002159007,
 0.0001028891,
 0.0000618602,
 0.0000433336,
 0.0000324305,
 0.0000233797,
 0.0000145054,
 0.0000073363,
 0.0000031279,
 0.0000011932,
 0.0000004231,
 0.0000001420,
 0.0000000454,
 0.0000000139,
 0.0000000041,
 0.0000000011,
 0.0000000003,
 0.0000000001,
 0.0000000000,
 0.0000000000,
 0.0000000000,
 0.0000000000,
 0.0000000000,
 0.0000000000,
 0.0000000000
};

 return w[(int)floor(ntrue)];

}

float kfactor_qqZZ_qcd_M(float GENmassZZ, int finalState = 2)
{

    // finalState=1 : 4e/4mu/4tau
    // finalState=2 : 2e2mu/2mutau/2e2tau

    float k=0.0;

    if (finalState==1) {
        k+=1.23613311013*(abs(GENmassZZ)>0.0&&abs(GENmassZZ)<=25.0);
        k+=1.17550314639*(abs(GENmassZZ)>25.0&&abs(GENmassZZ)<=50.0);
        k+=1.17044565911*(abs(GENmassZZ)>50.0&&abs(GENmassZZ)<=75.0);
        k+=1.03141209689*(abs(GENmassZZ)>75.0&&abs(GENmassZZ)<=100.0);
        k+=1.05285574912*(abs(GENmassZZ)>100.0&&abs(GENmassZZ)<=125.0);
        k+=1.11287217794*(abs(GENmassZZ)>125.0&&abs(GENmassZZ)<=150.0);
        k+=1.13361441158*(abs(GENmassZZ)>150.0&&abs(GENmassZZ)<=175.0);
        k+=1.10355603327*(abs(GENmassZZ)>175.0&&abs(GENmassZZ)<=200.0);
        k+=1.10053981637*(abs(GENmassZZ)>200.0&&abs(GENmassZZ)<=225.0);
        k+=1.10972676811*(abs(GENmassZZ)>225.0&&abs(GENmassZZ)<=250.0);
        k+=1.12069120525*(abs(GENmassZZ)>250.0&&abs(GENmassZZ)<=275.0);
        k+=1.11589101635*(abs(GENmassZZ)>275.0&&abs(GENmassZZ)<=300.0);
        k+=1.13906170314*(abs(GENmassZZ)>300.0&&abs(GENmassZZ)<=325.0);
        k+=1.14854594271*(abs(GENmassZZ)>325.0&&abs(GENmassZZ)<=350.0);
        k+=1.14616229031*(abs(GENmassZZ)>350.0&&abs(GENmassZZ)<=375.0);
        k+=1.14573157789*(abs(GENmassZZ)>375.0&&abs(GENmassZZ)<=400.0);
        k+=1.13829430515*(abs(GENmassZZ)>400.0&&abs(GENmassZZ)<=425.0);
        k+=1.15521193686*(abs(GENmassZZ)>425.0&&abs(GENmassZZ)<=450.0);
        k+=1.13679822698*(abs(GENmassZZ)>450.0&&abs(GENmassZZ)<=475.0);
        k+=1.13223956942*(abs(GENmassZZ)>475.0);
    }

    if (finalState==2) {
        k+=1.25094466582*(abs(GENmassZZ)>0.0&&abs(GENmassZZ)<=25.0);
        k+=1.22459455362*(abs(GENmassZZ)>25.0&&abs(GENmassZZ)<=50.0);
        k+=1.19287368979*(abs(GENmassZZ)>50.0&&abs(GENmassZZ)<=75.0);
        k+=1.04597506451*(abs(GENmassZZ)>75.0&&abs(GENmassZZ)<=100.0);
        k+=1.08323413771*(abs(GENmassZZ)>100.0&&abs(GENmassZZ)<=125.0);
        k+=1.09994968030*(abs(GENmassZZ)>125.0&&abs(GENmassZZ)<=150.0);
        k+=1.16698455800*(abs(GENmassZZ)>150.0&&abs(GENmassZZ)<=175.0);
        k+=1.10399053155*(abs(GENmassZZ)>175.0&&abs(GENmassZZ)<=200.0);
        k+=1.10592664340*(abs(GENmassZZ)>200.0&&abs(GENmassZZ)<=225.0);
        k+=1.10690381480*(abs(GENmassZZ)>225.0&&abs(GENmassZZ)<=250.0);
        k+=1.11194928918*(abs(GENmassZZ)>250.0&&abs(GENmassZZ)<=275.0);
        k+=1.13522586553*(abs(GENmassZZ)>275.0&&abs(GENmassZZ)<=300.0);
        k+=1.11895090244*(abs(GENmassZZ)>300.0&&abs(GENmassZZ)<=325.0);
        k+=1.13898508615*(abs(GENmassZZ)>325.0&&abs(GENmassZZ)<=350.0);
        k+=1.15463977506*(abs(GENmassZZ)>350.0&&abs(GENmassZZ)<=375.0);
        k+=1.17341664594*(abs(GENmassZZ)>375.0&&abs(GENmassZZ)<=400.0);
        k+=1.20093349763*(abs(GENmassZZ)>400.0&&abs(GENmassZZ)<=425.0);
        k+=1.18915554919*(abs(GENmassZZ)>425.0&&abs(GENmassZZ)<=450.0);
        k+=1.18546007375*(abs(GENmassZZ)>450.0&&abs(GENmassZZ)<=475.0);
        k+=1.12864505708*(abs(GENmassZZ)>475.0);
    }

    if (k==0.0) return 1.1;
    else return k; // if something goes wrong return inclusive k-factor

}
float kfactor_qqZZ_qcd_dPhi(float GENdPhiZZ, int finalState = 2)
{

    // finalState=1 : 4e/4mu/4tau
    // finalState=2 : 2e2mu/2mutau/2e2tau

    float k=0.0;    
    if (finalState==1) {        
        k+=1.515838921760*(abs(GENdPhiZZ)>0.0&&abs(GENdPhiZZ)<=0.1);
        k+=1.496256665410*(abs(GENdPhiZZ)>0.1&&abs(GENdPhiZZ)<=0.2);
        k+=1.495522061910*(abs(GENdPhiZZ)>0.2&&abs(GENdPhiZZ)<=0.3);
        k+=1.483273154250*(abs(GENdPhiZZ)>0.3&&abs(GENdPhiZZ)<=0.4);
        k+=1.465589701130*(abs(GENdPhiZZ)>0.4&&abs(GENdPhiZZ)<=0.5);
        k+=1.491500887510*(abs(GENdPhiZZ)>0.5&&abs(GENdPhiZZ)<=0.6);
        k+=1.441183580450*(abs(GENdPhiZZ)>0.6&&abs(GENdPhiZZ)<=0.7);
        k+=1.440830603990*(abs(GENdPhiZZ)>0.7&&abs(GENdPhiZZ)<=0.8);
        k+=1.414339019120*(abs(GENdPhiZZ)>0.8&&abs(GENdPhiZZ)<=0.9);
        k+=1.422534218560*(abs(GENdPhiZZ)>0.9&&abs(GENdPhiZZ)<=1.0);
        k+=1.401037066000*(abs(GENdPhiZZ)>1.0&&abs(GENdPhiZZ)<=1.1);
        k+=1.408539428810*(abs(GENdPhiZZ)>1.1&&abs(GENdPhiZZ)<=1.2);
        k+=1.381247744080*(abs(GENdPhiZZ)>1.2&&abs(GENdPhiZZ)<=1.3);
        k+=1.370553357430*(abs(GENdPhiZZ)>1.3&&abs(GENdPhiZZ)<=1.4);
        k+=1.347323316000*(abs(GENdPhiZZ)>1.4&&abs(GENdPhiZZ)<=1.5);
        k+=1.340113437450*(abs(GENdPhiZZ)>1.5&&abs(GENdPhiZZ)<=1.6);
        k+=1.312661036510*(abs(GENdPhiZZ)>1.6&&abs(GENdPhiZZ)<=1.7);
        k+=1.290055062010*(abs(GENdPhiZZ)>1.7&&abs(GENdPhiZZ)<=1.8);
        k+=1.255322614790*(abs(GENdPhiZZ)>1.8&&abs(GENdPhiZZ)<=1.9);
        k+=1.254455642450*(abs(GENdPhiZZ)>1.9&&abs(GENdPhiZZ)<=2.0);
        k+=1.224047664420*(abs(GENdPhiZZ)>2.0&&abs(GENdPhiZZ)<=2.1);
        k+=1.178816782670*(abs(GENdPhiZZ)>2.1&&abs(GENdPhiZZ)<=2.2);
        k+=1.162624827140*(abs(GENdPhiZZ)>2.2&&abs(GENdPhiZZ)<=2.3);
        k+=1.105401140940*(abs(GENdPhiZZ)>2.3&&abs(GENdPhiZZ)<=2.4);
        k+=1.074749265690*(abs(GENdPhiZZ)>2.4&&abs(GENdPhiZZ)<=2.5);
        k+=1.021864599380*(abs(GENdPhiZZ)>2.5&&abs(GENdPhiZZ)<=2.6);
        k+=0.946334793286*(abs(GENdPhiZZ)>2.6&&abs(GENdPhiZZ)<=2.7);
        k+=0.857458082628*(abs(GENdPhiZZ)>2.7&&abs(GENdPhiZZ)<=2.8);
        k+=0.716607670482*(abs(GENdPhiZZ)>2.8&&abs(GENdPhiZZ)<=2.9);
        k+=1.132841784840*(abs(GENdPhiZZ)>2.9&&abs(GENdPhiZZ)<=TMath::Pi());
    }

    if (finalState==2) {
       k+=1.513834489150*(abs(GENdPhiZZ)>0.0&&abs(GENdPhiZZ)<=0.1);
       k+=1.541738780180*(abs(GENdPhiZZ)>0.1&&abs(GENdPhiZZ)<=0.2);
       k+=1.497829632510*(abs(GENdPhiZZ)>0.2&&abs(GENdPhiZZ)<=0.3);
       k+=1.534956782920*(abs(GENdPhiZZ)>0.3&&abs(GENdPhiZZ)<=0.4);
       k+=1.478217033060*(abs(GENdPhiZZ)>0.4&&abs(GENdPhiZZ)<=0.5);
       k+=1.504330859290*(abs(GENdPhiZZ)>0.5&&abs(GENdPhiZZ)<=0.6);
       k+=1.520626246850*(abs(GENdPhiZZ)>0.6&&abs(GENdPhiZZ)<=0.7);
       k+=1.507013090030*(abs(GENdPhiZZ)>0.7&&abs(GENdPhiZZ)<=0.8);
       k+=1.494243156250*(abs(GENdPhiZZ)>0.8&&abs(GENdPhiZZ)<=0.9);
       k+=1.450536096150*(abs(GENdPhiZZ)>0.9&&abs(GENdPhiZZ)<=1.0);
       k+=1.460812521660*(abs(GENdPhiZZ)>1.0&&abs(GENdPhiZZ)<=1.1);
       k+=1.471603622200*(abs(GENdPhiZZ)>1.1&&abs(GENdPhiZZ)<=1.2);
       k+=1.467700038200*(abs(GENdPhiZZ)>1.2&&abs(GENdPhiZZ)<=1.3);
       k+=1.422408690640*(abs(GENdPhiZZ)>1.3&&abs(GENdPhiZZ)<=1.4);
       k+=1.397184022730*(abs(GENdPhiZZ)>1.4&&abs(GENdPhiZZ)<=1.5);
       k+=1.375593447520*(abs(GENdPhiZZ)>1.5&&abs(GENdPhiZZ)<=1.6);
       k+=1.391901318370*(abs(GENdPhiZZ)>1.6&&abs(GENdPhiZZ)<=1.7);
       k+=1.368564350560*(abs(GENdPhiZZ)>1.7&&abs(GENdPhiZZ)<=1.8);
       k+=1.317884804290*(abs(GENdPhiZZ)>1.8&&abs(GENdPhiZZ)<=1.9);
       k+=1.314019950800*(abs(GENdPhiZZ)>1.9&&abs(GENdPhiZZ)<=2.0);
       k+=1.274641749910*(abs(GENdPhiZZ)>2.0&&abs(GENdPhiZZ)<=2.1);
       k+=1.242346606820*(abs(GENdPhiZZ)>2.1&&abs(GENdPhiZZ)<=2.2);
       k+=1.244727403840*(abs(GENdPhiZZ)>2.2&&abs(GENdPhiZZ)<=2.3);
       k+=1.146259351670*(abs(GENdPhiZZ)>2.3&&abs(GENdPhiZZ)<=2.4);
       k+=1.107804993520*(abs(GENdPhiZZ)>2.4&&abs(GENdPhiZZ)<=2.5);
       k+=1.042053646740*(abs(GENdPhiZZ)>2.5&&abs(GENdPhiZZ)<=2.6);
       k+=0.973608545141*(abs(GENdPhiZZ)>2.6&&abs(GENdPhiZZ)<=2.7);
       k+=0.872169942668*(abs(GENdPhiZZ)>2.7&&abs(GENdPhiZZ)<=2.8);
       k+=0.734505279177*(abs(GENdPhiZZ)>2.8&&abs(GENdPhiZZ)<=2.9);
       k+=1.163152837230*(abs(GENdPhiZZ)>2.9&&abs(GENdPhiZZ)<=TMath::Pi());       
    }
    if (k==0.0) return 1.1; // if something goes wrong return inclusive k-factor
    else return k;

}

float kfactor_qqZZ_qcd_Pt(float GENpTZZ, int finalState = 2)
{

    // finalState=1 : 4e/4mu/4tau
    // finalState=2 : 2e2mu/2mutau/2e2tau

    float k=0.0;

    if (finalState==1) {
        k+=0.64155491983*(abs(GENpTZZ)>0.0&&abs(GENpTZZ)<=5.0);
        k+=1.09985240531*(abs(GENpTZZ)>5.0&&abs(GENpTZZ)<=10.0);
        k+=1.29390628654*(abs(GENpTZZ)>10.0&&abs(GENpTZZ)<=15.0);
        k+=1.37859998571*(abs(GENpTZZ)>15.0&&abs(GENpTZZ)<=20.0);
        k+=1.42430263312*(abs(GENpTZZ)>20.0&&abs(GENpTZZ)<=25.0);
        k+=1.45038493266*(abs(GENpTZZ)>25.0&&abs(GENpTZZ)<=30.0);
        k+=1.47015377651*(abs(GENpTZZ)>30.0&&abs(GENpTZZ)<=35.0);
        k+=1.48828685748*(abs(GENpTZZ)>35.0&&abs(GENpTZZ)<=40.0);
        k+=1.50573440448*(abs(GENpTZZ)>40.0&&abs(GENpTZZ)<=45.0);
        k+=1.50211655928*(abs(GENpTZZ)>45.0&&abs(GENpTZZ)<=50.0);
        k+=1.50918720827*(abs(GENpTZZ)>50.0&&abs(GENpTZZ)<=55.0);
        k+=1.52463089491*(abs(GENpTZZ)>55.0&&abs(GENpTZZ)<=60.0);
        k+=1.52400838378*(abs(GENpTZZ)>60.0&&abs(GENpTZZ)<=65.0);
        k+=1.52418067701*(abs(GENpTZZ)>65.0&&abs(GENpTZZ)<=70.0);
        k+=1.55424382578*(abs(GENpTZZ)>70.0&&abs(GENpTZZ)<=75.0);
        k+=1.52544284222*(abs(GENpTZZ)>75.0&&abs(GENpTZZ)<=80.0);
        k+=1.57896384602*(abs(GENpTZZ)>80.0&&abs(GENpTZZ)<=85.0);
        k+=1.53034682567*(abs(GENpTZZ)>85.0&&abs(GENpTZZ)<=90.0);
        k+=1.56147329708*(abs(GENpTZZ)>90.0&&abs(GENpTZZ)<=95.0);
        k+=1.54468169268*(abs(GENpTZZ)>95.0&&abs(GENpTZZ)<=100.0);
        k+=1.57222952415*(abs(GENpTZZ)>100.0);
    }

    if (finalState==2) {
        k+=0.743602533303*(abs(GENpTZZ)>0.0&&abs(GENpTZZ)<=5.0);
        k+=1.14789453219*(abs(GENpTZZ)>5.0&&abs(GENpTZZ)<=10.0);
        k+=1.33815867892*(abs(GENpTZZ)>10.0&&abs(GENpTZZ)<=15.0);
        k+=1.41420044104*(abs(GENpTZZ)>15.0&&abs(GENpTZZ)<=20.0);
        k+=1.45511318916*(abs(GENpTZZ)>20.0&&abs(GENpTZZ)<=25.0);
        k+=1.47569225244*(abs(GENpTZZ)>25.0&&abs(GENpTZZ)<=30.0);
        k+=1.49053003693*(abs(GENpTZZ)>30.0&&abs(GENpTZZ)<=35.0);
        k+=1.50622827695*(abs(GENpTZZ)>35.0&&abs(GENpTZZ)<=40.0);
        k+=1.50328889799*(abs(GENpTZZ)>40.0&&abs(GENpTZZ)<=45.0);
        k+=1.52186945281*(abs(GENpTZZ)>45.0&&abs(GENpTZZ)<=50.0);
        k+=1.52043468754*(abs(GENpTZZ)>50.0&&abs(GENpTZZ)<=55.0);
        k+=1.53977869986*(abs(GENpTZZ)>55.0&&abs(GENpTZZ)<=60.0);
        k+=1.53491994434*(abs(GENpTZZ)>60.0&&abs(GENpTZZ)<=65.0);
        k+=1.51772882172*(abs(GENpTZZ)>65.0&&abs(GENpTZZ)<=70.0);
        k+=1.54494489131*(abs(GENpTZZ)>70.0&&abs(GENpTZZ)<=75.0);
        k+=1.57762411697*(abs(GENpTZZ)>75.0&&abs(GENpTZZ)<=80.0);
        k+=1.55078339014*(abs(GENpTZZ)>80.0&&abs(GENpTZZ)<=85.0);
        k+=1.57078191891*(abs(GENpTZZ)>85.0&&abs(GENpTZZ)<=90.0);
        k+=1.56162666568*(abs(GENpTZZ)>90.0&&abs(GENpTZZ)<=95.0);
        k+=1.54183774627*(abs(GENpTZZ)>95.0&&abs(GENpTZZ)<=100.0);
        k+=1.58485762205*(abs(GENpTZZ)>100.0);
    }

    if (k==0.0) return 1.1;
    else return k; // if something goes wrong return inclusive k-factor

}

double ratioFactor_gjets_zll(double pt){
  double p[4] = {-1.04869,0.0293761,-5.86547e-05,3.60149e-08};
  if     (pt<= 60) pt = 60;
  else if(pt>=400) pt = 400;
  return (p[0]+p[1]*pt+p[2]*pt*pt+p[3]*pt*pt*pt);
}
