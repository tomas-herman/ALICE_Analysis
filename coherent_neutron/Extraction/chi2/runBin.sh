#!/bin/bash

if [ -z $1 ] || [ $1 -lt 0 ] || [ $1 -gt 5 ]
then
    echo "parameter needed: type an integer in [0,5]"
    exit -1 
fi

case $1 in
    0 ) echo "copying DataSimone/inputData3540.h to inputData.h"
	cp DataSimone/inputData3540.h inputData.h ;;
    1 ) echo "copying DataSimone/inputData3035.h to inputData.h"
	cp DataSimone/inputData3035.h inputData.h ;;
    2 ) echo "copying DataSimone/inputData2530.h to inputData.h"
	cp DataSimone/inputData2530.h inputData.h ;;
    3 ) echo "copying DataMichal/inputData2080.h to inputData.h"
	cp DataMichal/inputData2080.h inputData.h ;;
    4 ) echo "copying DataMichal/inputData0020.h to inputData.h"
	cp DataMichal/inputData0020.h inputData.h ;;
    5 ) echo "copying DataMichal/inputData3580.h to inputData.h"
	cp DataMichal/inputData3580.h inputData.h ;;
    6 ) echo "copying DataMichal/inputData1535.h to inputData.h"
	cp DataMichal/inputData1535.h inputData.h ;;
    7 ) echo "copying DataMichal/inputData0015.h to inputData.h"
	cp DataMichal/inputData0015.h inputData.h ;;
    ? ) echo "parameter needed: type an integer in [0,5]"
	exit -1 ;;
esac

echo "Fitting the data"
root.exe -q -b fitgPb.C+

