import json
import ROOT as rt
from array import array

f=open('RunBCDEF_SF_ID.json','r')

results = json.load(f)

etaBins=array('d',[0,0.9,1.2,2.1,2.4])
ptBins=array('d',[20.,25.,30.,40.,50.,60.,120.])


histEffMuID=rt.TH2F("histEffMuID","histEffMuID",4,etaBins,6,ptBins)


eta=0
for etaKey, values in sorted(results["NUM_TightID_DEN_genTracks"]["abseta_pt"].iteritems()) :
    pt=1
    eta=eta+1
    for ptKey, result in sorted(values.iteritems()) :
	histEffMuID.SetBinContent(eta,pt,float(result["value"]))
	print result["value"]	
	pt=pt+1
        print "|eta| bin: %s  pT bin: %s\tdata/MC SF: %f +/- %f" % (etaKey, ptKey, result["value"], result["error"])
    

histEffMuID.SaveAs("histEffMuID.root")
