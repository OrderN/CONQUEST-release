#!/bin/bash

file=Conquest_out

if [ -f "$file" ]; then
   #grep 'DFT total energy' Conquest_out | awk '{print $6}'
   awk '/DFT total energy/ {print $6}' Conquest_out > tmp.energy
   awk '/force: Force Residual:/ {print $4}' Conquest_out > tmp.force
fi 
