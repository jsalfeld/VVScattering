double effSF_m_25_medium_tight[10] = {0.893181,0.903896,0.918628,0.970606,0.961617,0.988883,0.981609,0.988133,0.997004,0.980483};
double effSF_m_25_medium[10] = {0.883965,0.888433,0.911088,0.962104,0.950514,0.981674,0.969118,0.978893,0.983652,0.968170};
double effSF_e_25_medium[10] = {0.882568,0.897607,0.907722,0.966725,0.947403,0.987567,0.966422,0.984508,0.980736,0.969033};
double effSF_m_25_tight[10]  = {0.874456,0.881048,0.900731,0.958797,0.912054,1.023287,0.929299,1.017120,0.951510,0.948481};
double effSF_e_25_tight[10]  = {0.870300,0.838181,0.890720,0.925357,0.901283,1.013620,0.924440,1.013939,0.943693,0.942360};

double fake_rate_m_25_medium_tight[5][5] = {
0.289,0.199,0.183,0.168,0.165,
0.311,0.218,0.199,0.189,0.191,
0.345,0.260,0.241,0.228,0.225,
0.373,0.299,0.284,0.266,0.261,
0.376,0.304,0.294,0.294,0.289
};
double fake_rate_m_25_medium[5][5] = {
0.352,0.259,0.242,0.224,0.223,
0.377,0.283,0.261,0.251,0.250,
0.412,0.329,0.309,0.294,0.290,
0.444,0.370,0.355,0.338,0.333,
0.443,0.376,0.370,0.367,0.365
};
double fake_rate_e_25_medium[5][5] = {
0.233,0.219,0.187,0.202,0.186,
0.236,0.221,0.199,0.165,0.183,
0.293,0.248,0.215,0.191,0.194,
0.281,0.242,0.236,0.256,0.278,
0.284,0.262,0.268,0.291,0.317
};
double fake_rate_m_25_tight[5][5] = {
0.360,0.266,0.250,0.233,0.241,
0.381,0.285,0.264,0.255,0.260,
0.416,0.330,0.310,0.296,0.299,
0.447,0.371,0.355,0.338,0.340,
0.448,0.377,0.370,0.368,0.375
};
double fake_rate_e_25_tight[5][5] = {
0.136,0.112,0.101,0.112,0.099,
0.131,0.118,0.113,0.087,0.103,
0.170,0.144,0.115,0.107,0.103,
0.205,0.177,0.154,0.151,0.155,
0.206,0.165,0.154,0.168,0.185
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
    isoCut = 0.15; if(type == "default_tight") isoCut = 0.12;
    if     (type == "default" || type == "default_tight" || type == "medium") idCut = (selBits & BareLeptons::LepMediumIP) == BareLeptons::LepMediumIP;
    else if(type == "tight")                                                  idCut = (selBits & BareLeptons::LepTightIP)  == BareLeptons::LepTightIP;
  }
  else if(TMath::Abs(pdgId) == 11) {
    if     (type == "medium")                                                isoCut = (isEB ? 0.0766 : 0.0678);
    else if(type == "default" || type == "default_tight" || type == "tight") isoCut = (isEB ? 0.0354 : 0.0646);

    if     (type == "medium")                                                idCut = (selBits & BareLeptons::LepMedium) == BareLeptons::LepMedium;
    else if(type == "default" || type == "default_tight" || type == "tight") idCut = (selBits & BareLeptons::LepTight)  == BareLeptons::LepTight;
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

double effhDScaleFactor(bool isAbs, double pt, double eta, int nsel, int period, TString type, TH2D *fhDMuMediumSF, TH2D *fhDElMediumSF, TH2D *fhDElTightSF){
  if(isAbs == true) eta = abs(eta);
  Int_t binX = 0;
  Int_t binY = 0;

  if     (TMath::Abs(nsel) == 13 && period == 1 && (type== "medium" || type== "default")) {binX = fhDMuMediumSF->GetXaxis()->FindFixBin(TMath::Min(eta,2.399));binY = fhDMuMediumSF->GetYaxis()->FindFixBin(TMath::Min(pt,99.999));}
  else if(TMath::Abs(nsel) == 11 && period == 1 &&  type== "medium")                      {binX = fhDElMediumSF->GetXaxis()->FindFixBin(TMath::Min(eta,2.399));binY = fhDElMediumSF->GetYaxis()->FindFixBin(TMath::Min(pt,99.999));}
  else if(TMath::Abs(nsel) == 11 && period == 1 && (type== "default"))                    {binX = fhDElTightSF ->GetXaxis()->FindFixBin(TMath::Min(eta,2.399));binY = fhDElTightSF ->GetYaxis()->FindFixBin(TMath::Min(pt,99.999));}
  else    printf("PROBLEM WITH BINS\n");

  if     (TMath::Abs(nsel) == 13 && period == 1 && (type== "medium" || type== "default")) return fhDMuMediumSF->GetBinContent(binX, binY);
  else if(TMath::Abs(nsel) == 11 && period == 1 &&  type== "medium")                      return fhDElMediumSF->GetBinContent(binX, binY);
  else if(TMath::Abs(nsel) == 11 && period == 1 && (type== "default"))                    return fhDElTightSF ->GetBinContent(binX, binY);
  return 0;
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

  if     (TMath::Abs(nsel) == 13 && period == 1 &&  type== "default_tight")                                         return effSF_m_25_medium_tight[iPoint];
  else if(TMath::Abs(nsel) == 13 && period == 1 &&  type== "tight")                                                 return effSF_m_25_tight [iPoint];
  else if(TMath::Abs(nsel) == 13 && period == 1 && (type== "medium" || type== "default"))                           return effSF_m_25_medium[iPoint];
  else if(TMath::Abs(nsel) == 11 && period == 1 &&  type== "medium")                                                return effSF_e_25_medium[iPoint];
  else if(TMath::Abs(nsel) == 11 && period == 1 && (type== "tight" || type== "default" || type== "default_tight" )) return effSF_e_25_tight [iPoint];

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

  if     (TMath::Abs(nsel) == 13 && period == 1 &&  type== "default_tight")					    return fake_rate_m_25_medium_tight[iPt][iEta]/(1.0-fake_rate_m_25_medium_tight[iPt][iEta]);
  else if(TMath::Abs(nsel) == 13 && period == 1 &&  type== "tight")						    return fake_rate_m_25_tight [iPt][iEta]/(1.0-fake_rate_m_25_tight [iPt][iEta]);
  else if(TMath::Abs(nsel) == 13 && period == 1 && (type== "medium" || type== "default"))			    return fake_rate_m_25_medium[iPt][iEta]/(1.0-fake_rate_m_25_medium[iPt][iEta]);
  else if(TMath::Abs(nsel) == 11 && period == 1 &&  type== "medium")						    return fake_rate_e_25_medium[iPt][iEta]/(1.0-fake_rate_e_25_medium[iPt][iEta]);
  else if(TMath::Abs(nsel) == 11 && period == 1 && (type== "tight" || type== "default" || type== "default_tight" )) return fake_rate_e_25_tight [iPt][iEta]/(1.0-fake_rate_e_25_tight [iPt][iEta]);

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
