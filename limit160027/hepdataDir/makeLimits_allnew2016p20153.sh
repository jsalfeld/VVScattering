#!/bin/sh
combine -M Asymptotic data_combFinal.txt -m 200  
combine -M Asymptotic data_combFinal.txt -m 300  
combine -M Asymptotic data_combFinal.txt -m 400  
combine -M Asymptotic data_combFinal.txt -m 500  
combine -M Asymptotic data_combFinal.txt -m 600  
combine -M Asymptotic data_combFinal.txt -m 700  
combine -M Asymptotic data_combFinal.txt -m 800  
combine -M Asymptotic data_combFinal.txt -m 900  
combine -M Asymptotic data_combFinal.txt -m 1000  
combine -M Asymptotic data_combFinal.txt -m 1500  
combine -M Asymptotic data_combFinal.txt -m 2000  
	

python makeLimit3.py