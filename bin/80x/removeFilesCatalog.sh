#/!/bin/sh

rm -f list_bad_files;

cat > list_bad_files <<EOF

EOF

awk '{print"grep -vwE "$1" ~/catalog/t2mit/"$3"/"$2"/Files       > ooo;wc ooo ~/catalog/t2mit/"$3"/"$2"/Files;      mv ooo ~/catalog/t2mit/"$3"/"$2"/Files;      "}' list_bad_files  > zzz
awk '{print"grep -vwE "$1" ~/catalog/t2mit/"$3"/"$2"/RawFiles.00 > ooo;wc ooo ~/catalog/t2mit/"$3"/"$2"/RawFiles.00;mv ooo ~/catalog/t2mit/"$3"/"$2"/RawFiles.00;"}' list_bad_files >> zzz
chmod a+x zzz;./zzz;rm -f zzz;
rm -f list_bad_files;
