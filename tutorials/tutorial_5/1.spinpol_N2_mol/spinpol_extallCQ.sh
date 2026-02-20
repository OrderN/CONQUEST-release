dirs=(
"1.singlet" 
"2.triplet"
"3.triplet_spinfix"
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
