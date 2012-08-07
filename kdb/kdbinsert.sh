#!/usr/bin/env sh

curl -d "reactant=`cat $1`" -d "saddle=`cat $2`" -d "product=`cat $3`" -d "mode=`cat $4`" -d "nf=0.2" -d "dc=0.3" -d "mac=0.7" -d"b1=0" -d"b2=0" -d"p1=0" -d"p2=0" "http://localhost:8080/insert"

echo

