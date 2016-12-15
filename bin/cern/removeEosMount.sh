#!/bin/bash

NODE=0

while [ ${NODE} -lt 100 ] 
do
  NODE=$[$NODE+1]
  ssh lxplus00${NODE} 'hostname;pkill -u ceballos eosfsd;killall -9 eosfsd;ps -uceballos -f|grep eosfsd'
  #ssh lxplus00${NODE} 'hostname;ps -uceballos -f|grep root'
done

NODE=10
while [ ${NODE} -lt 100 ] 
do
  NODE=$[$NODE+1]
  ssh lxplus0${NODE} 'hostname;pkill -u ceballos eosfsd;killall -9 eosfsd;ps -uceballos -f|grep eosfsd'
  #ssh lxplus0${NODE} 'hostname;ps -uceballos -f|grep root'
done

while [ ${NODE} -lt 150 ] 
do
  NODE=$[$NODE+1]
  ssh lxplus${NODE} 'hostname;pkill -u ceballos eosfsd;killall -9 eosfsd;ps -uceballos -f|grep eosfsd'
  #ssh lxplus${NODE} 'hostname;ps -uceballos -f|grep root'
done
