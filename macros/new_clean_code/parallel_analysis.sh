#!/bin/bash
export DISPLAY=lxir123.gsi.de:0.0
# Our custom function
cust_func(){
	   INPUT="$1"
	   echo $INPUT
	   root -l "test_jt.C(\"$INPUT\")"
  }
echo "HELLO"   
  while IFS= read -r url
  do
	          cust_func "$url" &
		  echo  "calibrating root file  RUN 273,subrun $url"
	  done < list_combined.txt
	   
	  wait
	  echo "All files calibrated now"
