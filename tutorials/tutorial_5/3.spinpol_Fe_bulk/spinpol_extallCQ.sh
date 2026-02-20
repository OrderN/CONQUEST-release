dirs=(
"1.bcc_AF"
"2.bcc_FM"
"3.bcc_NM"
"4.fcc_FM"
"5.fcc_NM"
"6.hcp_FM"
"7.hcp_NM"
)

for d in "${dirs[@]}"; do
  (
    cd "$d" 
    if [ -f "Conquest_out" ]; then
       energy=$(awk '/DFT total energy/ {print $6}' Conquest_out)
       printf "$d energy = %.4f Ha \n" "$energy"
    fi 
  )
done
