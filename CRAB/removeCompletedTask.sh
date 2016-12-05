grep -e Executing -e "Task status" log|awk '{if(NR%2==1)printf("%s ",$8);else printf("%s\n",$3)}'|grep COMPLETED|awk '{printf("rm -rf %s\n",$1)}' > ll;
wc ll;
chmod a+x ll;./ll;rm -f ll
