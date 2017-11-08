double effMM[72] = {0.98565,0.99304,0.99197,0.97609,0.98977,1.00590,0.95855,0.99551,0.99481,0.99772,0.98271,0.97754,0.98284,0.99117,1.00027,0.99547,0.98448,0.99665,0.99928,1.00244,0.97493,1.00190,0.97300,1.00136,0.95925,0.99527,0.98787,0.96118,0.99874,0.99339,0.99473,0.99743,0.99660,0.98180,1.00116,0.99718,0.99542,0.99579,0.98471,0.98373,0.99352,0.99369,0.99485,0.98047,0.97130,0.99841,0.98924,0.97516,0.98459,0.98298,0.98789,0.99689,0.99771,0.99773,0.98406,0.98518,0.98858,0.99402,0.99640,0.99416,0.98512,0.98311,0.98516,0.99525,0.99706,0.99530,0.95997,0.95837,0.97495,0.99227,0.98727,0.98644};
double effEE[72] = {0.79324,0.60749,0.85393,0.41745,0.77167,0.86425,0.66309,0.44008,0.57977,1.00000,0.22690,0.94076,0.90107,0.92328,0.96178,0.96518,0.77869,0.94818,0.92218,1.02941,0.75272,0.77919,0.87122,0.89652,0.83144,0.84840,0.92460,0.85279,0.98179,0.97939,0.98123,0.99042,0.99012,0.88422,0.96993,0.97326,0.99765,0.98814,0.83832,0.91796,0.97447,0.96760,0.99216,0.96294,0.97574,0.89834,0.89290,0.98135,0.97314,0.98661,0.98810,0.99291,0.99272,0.98739,0.99516,0.98174,0.99577,0.99767,0.99502,0.99418,0.93587,0.92965,0.96184,0.98037,0.98749,0.99252,0.80292,0.92725,0.91587,0.97125,0.97560,0.94559};
double effEM[72] = {0.83261,0.88222,0.94980,0.89576,0.83451,0.85746,0.56934,0.68406,0.85664,0.56528,0.90707,0.82426,0.97277,0.96872,0.97794,0.99445,0.96856,1.00634,0.99374,1.00489,0.88966,0.94798,0.97441,0.97554,0.93136,0.88845,0.88053,1.03506,0.97962,0.99329,0.98951,0.99565,0.99749,0.96953,0.99444,0.99861,0.99787,0.99932,0.96410,0.98092,0.98562,0.98915,0.98577,0.91049,0.97990,0.94287,0.98525,0.97917,0.98403,0.98733,0.99436,0.99811,0.99727,0.99795,0.98319,0.99364,0.98339,0.98927,0.99181,0.99418,0.99291,0.98213,0.98643,0.98099,0.98963,0.99290,0.95509,0.96373,0.98176,0.97141,0.97938,0.97138};

double trigger_sf(double pt1, double eta1, int pdg1, double pt2, double eta2, int pdg2){
      double ptaux = pt1;
      if(pt1 < pt2) {pt1 = pt2; pt2 = ptaux; printf("Changed order\n");}

      int iPt[2] = {-1, -1};
      if     (pt1 < 25) iPt[0] = 0;
      else if(pt1 < 30) iPt[0] = 1;
      else if(pt1 < 50) iPt[0] = 2;
      else		iPt[0] = 3;
      if     (pt2 < 15) iPt[1] = 0;
      else if(pt2 < 20) iPt[1] = 1;
      else if(pt2 < 25) iPt[1] = 2;
      else if(pt2 < 30) iPt[1] = 3;
      else if(pt2 < 50) iPt[1] = 4;
      else		iPt[1] = 5;

      int iEta[2] = {-1, -1};
      if    (TMath::Abs(eta1) < 1.5) iEta[0] = 0;
      else			     iEta[0] = 1;
      if    (TMath::Abs(eta1) < 1.5) iEta[1] = 0;
      else			     iEta[1] = 1;

      int theBin = -1;
      if     (iPt[0] == 0 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 0) theBin = 0;
      else if(iPt[0] == 0 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 0) theBin = 1;
      else if(iPt[0] == 0 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 0) theBin = 2;
      else if(iPt[0] == 0 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 1) theBin = 3;
      else if(iPt[0] == 0 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 1) theBin = 4;
      else if(iPt[0] == 0 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 1) theBin = 5;
      else if(iPt[0] == 0 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 0) theBin = 6;
      else if(iPt[0] == 0 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 0) theBin = 7;
      else if(iPt[0] == 0 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 0) theBin = 8;
      else if(iPt[0] == 0 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 1) theBin = 9;
      else if(iPt[0] == 0 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 1) theBin =10;
      else if(iPt[0] == 0 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 1) theBin =11;

      else if(iPt[0] == 1 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 0) theBin =12;
      else if(iPt[0] == 1 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 0) theBin =13;
      else if(iPt[0] == 1 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 0) theBin =14;
      else if(iPt[0] == 1 && iPt[1] == 3 && iEta[0] == 0 && iEta[1] == 0) theBin =15;
      else if(iPt[0] == 1 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 1) theBin =16;
      else if(iPt[0] == 1 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 1) theBin =17;
      else if(iPt[0] == 1 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 1) theBin =18;
      else if(iPt[0] == 1 && iPt[1] == 3 && iEta[0] == 0 && iEta[1] == 1) theBin =19;
      else if(iPt[0] == 1 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 0) theBin =20;
      else if(iPt[0] == 1 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 0) theBin =21;
      else if(iPt[0] == 1 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 0) theBin =22;
      else if(iPt[0] == 1 && iPt[1] == 3 && iEta[0] == 1 && iEta[1] == 0) theBin =23;
      else if(iPt[0] == 1 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 1) theBin =24;
      else if(iPt[0] == 1 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 1) theBin =25;
      else if(iPt[0] == 1 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 1) theBin =26;
      else if(iPt[0] == 1 && iPt[1] == 3 && iEta[0] == 1 && iEta[1] == 1) theBin =27;

      else if(iPt[0] == 2 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 0) theBin =28;
      else if(iPt[0] == 2 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 0) theBin =29;
      else if(iPt[0] == 2 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 0) theBin =30;
      else if(iPt[0] == 2 && iPt[1] == 3 && iEta[0] == 0 && iEta[1] == 0) theBin =31;
      else if(iPt[0] == 2 && iPt[1] == 4 && iEta[0] == 0 && iEta[1] == 0) theBin =32;
      else if(iPt[0] == 2 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 1) theBin =33;
      else if(iPt[0] == 2 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 1) theBin =34;
      else if(iPt[0] == 2 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 1) theBin =35;
      else if(iPt[0] == 2 && iPt[1] == 3 && iEta[0] == 0 && iEta[1] == 1) theBin =36;
      else if(iPt[0] == 2 && iPt[1] == 4 && iEta[0] == 0 && iEta[1] == 1) theBin =37;
      else if(iPt[0] == 2 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 0) theBin =38;
      else if(iPt[0] == 2 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 0) theBin =39;
      else if(iPt[0] == 2 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 0) theBin =40;
      else if(iPt[0] == 2 && iPt[1] == 3 && iEta[0] == 1 && iEta[1] == 0) theBin =41;
      else if(iPt[0] == 2 && iPt[1] == 4 && iEta[0] == 1 && iEta[1] == 0) theBin =42;
      else if(iPt[0] == 2 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 1) theBin =43;
      else if(iPt[0] == 2 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 1) theBin =44;
      else if(iPt[0] == 2 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 1) theBin =45;
      else if(iPt[0] == 2 && iPt[1] == 3 && iEta[0] == 1 && iEta[1] == 1) theBin =46;
      else if(iPt[0] == 2 && iPt[1] == 4 && iEta[0] == 1 && iEta[1] == 1) theBin =47;

      else if(iPt[0] == 3 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 0) theBin =48;
      else if(iPt[0] == 3 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 0) theBin =49;
      else if(iPt[0] == 3 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 0) theBin =50;
      else if(iPt[0] == 3 && iPt[1] == 3 && iEta[0] == 0 && iEta[1] == 0) theBin =51;
      else if(iPt[0] == 3 && iPt[1] == 4 && iEta[0] == 0 && iEta[1] == 0) theBin =52;
      else if(iPt[0] == 3 && iPt[1] == 5 && iEta[0] == 0 && iEta[1] == 0) theBin =53;
      else if(iPt[0] == 3 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 1) theBin =54;
      else if(iPt[0] == 3 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 1) theBin =55;
      else if(iPt[0] == 3 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 1) theBin =56;
      else if(iPt[0] == 3 && iPt[1] == 3 && iEta[0] == 0 && iEta[1] == 1) theBin =57;
      else if(iPt[0] == 3 && iPt[1] == 4 && iEta[0] == 0 && iEta[1] == 1) theBin =58;
      else if(iPt[0] == 3 && iPt[1] == 5 && iEta[0] == 0 && iEta[1] == 1) theBin =59;
      else if(iPt[0] == 3 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 0) theBin =60;
      else if(iPt[0] == 3 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 0) theBin =61;
      else if(iPt[0] == 3 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 0) theBin =62;
      else if(iPt[0] == 3 && iPt[1] == 3 && iEta[0] == 1 && iEta[1] == 0) theBin =63;
      else if(iPt[0] == 3 && iPt[1] == 4 && iEta[0] == 1 && iEta[1] == 0) theBin =64;
      else if(iPt[0] == 3 && iPt[1] == 5 && iEta[0] == 1 && iEta[1] == 0) theBin =65;
      else if(iPt[0] == 3 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 1) theBin =66;
      else if(iPt[0] == 3 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 1) theBin =67;
      else if(iPt[0] == 3 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 1) theBin =68;
      else if(iPt[0] == 3 && iPt[1] == 3 && iEta[0] == 1 && iEta[1] == 1) theBin =69;
      else if(iPt[0] == 3 && iPt[1] == 4 && iEta[0] == 1 && iEta[1] == 1) theBin =70;
      else if(iPt[0] == 3 && iPt[1] == 5 && iEta[0] == 1 && iEta[1] == 1) theBin =71;
 
      else {printf("IMPOSSIBLE\n");}

      if     (abs(pdg1) == 13 && abs(pdg2) == 13) return effMM[theBin];
      else if(abs(pdg1) == 11 && abs(pdg2) == 11) return effEE[theBin];
      else                                        return effEM[theBin];
      
      return 1.0;

}
