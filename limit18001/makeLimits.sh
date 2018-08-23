#!/bin/sh

#datacard={$1}

echo "M200"
combine -M AsymptoticLimits datacardH.text -m 200
echo "M300" 
combine -M AsymptoticLimits datacardH.text -m 300
echo "M400"  
combine -M AsymptoticLimits datacardH.text -m 400
echo "M500"
combine -M AsymptoticLimits datacardH.text -m 500
echo "M600"
combine -M AsymptoticLimits datacardH.text -m 600
echo "M700" 
combine -M AsymptoticLimits datacardH.text -m 700
echo "M800"
combine -M AsymptoticLimits datacardH.text -m 800
echo "M900"
combine -M AsymptoticLimits datacardH.text -m 900
echo "M1000"
combine -M AsymptoticLimits datacardH.text -m 1000
echo "M1500"
combine -M AsymptoticLimits datacardH.text -m 1500
echo "M2000"
combine -M AsymptoticLimits datacardH.text -m 2000

python makeLimits.py