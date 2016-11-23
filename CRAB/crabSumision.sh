#/bin/sh

cp crabNero_static.py crabNero_dinamic.py;
awk '{split($1,a,"/");printf("    config.General.requestName = \"%s\"\n",a[2]);printf("    config.Data.inputDataset = \"%s\"\n",$1);printf("    submit(config)\n\n");}' miniaod_samples.txt >> crabNero_dinamic.py;
