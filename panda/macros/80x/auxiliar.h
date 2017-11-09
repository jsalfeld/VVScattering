#include "MitAnalysisRunII/panda/macros/80x/auxiliar_data.h"

double trigger_sf(double theTrgEff[2][2][nTrgBinPt1][nTrgBinPt2][nTrgBinEta1][nTrgBinEta2], double pt1, double eta1, int pdg1, double pt2, double eta2, int pdg2){
  double ptaux = pt1;
  if(pt1 < pt2) {printf("Changed order pt1(%f) pt2(%f)\n",pt1,pt2); pt1 = pt2; pt2 = ptaux;}

  int npdgId1 = TMath::Min(abs(pdg1)-11,1);
  int npdgId2 = TMath::Min(abs(pdg2)-11,1);

  int npt1  = hDTrgBinPt1 ->GetXaxis()->FindFixBin(pt1)-1;
  int npt2  = hDTrgBinPt2 ->GetXaxis()->FindFixBin(pt2)-1;
  int neta1 = hDTrgBinEta1->GetXaxis()->FindFixBin(eta1)-1;
  int neta2 = hDTrgBinEta2->GetXaxis()->FindFixBin(eta2)-1;

  return theTrgEff[npdgId1][npdgId2][npt1][npt2][neta1][neta2];
}
