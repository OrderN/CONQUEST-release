#!/bin/bash

file=Conquest_out

if [ -f "$file" ]; then
   #grep 'DFT total energy' Conquest_out | awk '{print $6}'
   awk '/DFT total energy/ {print $6}' Conquest_out > tmp.energy
   awk '/force: Force Residual:/ {print $4}' Conquest_out > tmp.force
   awk '/enthalpy/ {print $8}' Conquest_out > tmp.enth
   awk '/Maximum stress:/ {print $3}' Conquest_out > tmp.stress
fi 
