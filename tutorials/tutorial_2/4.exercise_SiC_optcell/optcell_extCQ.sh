#!/bin/bash

file=Conquest_out

if [ -f "$file" ]; then
   #grep 'DFT total energy' Conquest_out | awk '{print $6}'
   awk '/enthalpy/ {print $8}' Conquest_out > tmp.energy
   awk '/Maximum stress:/ {print $3}' Conquest_out > tmp.stress
fi 
