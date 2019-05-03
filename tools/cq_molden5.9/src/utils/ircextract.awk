BEGIN{j=0
at_symbol[1]="H"
at_symbol[6]="C"
at_symbol[7]="N"
at_symbol[8]="O"
}
{
 if ($1=="Z-Matrix" && $2=="orientation:") {
   getline
   getline
   getline
   getline
   getline
   i=0
   while (NF!=1) {
      i++
      at[i]=$2
      if (NF==5) {
        x[i]=$3
        y[i]=$4
        z[i]=$5 
        }
      else {
        x[i]=$4
        y[i]=$5
        z[i]=$6 
        }
      getline
      }
   }
  if ($1=="SCF" && $2=="Done:") {
    energy=$5
    }
  if ($3=="Threshold") {
    getline
    if ($5=="YES") uno=1
      else uno=0
    getline 
    if ($5=="YES") dos=1
      else dos=0
    getline 
    if ($5=="YES") tres=1
      else tres=0
    getline 
    if ($5=="YES") cuatro=1
      else cuatro=0
    if (uno==1 && dos==1 && tres==1 && cuatro==1) {
      j++
      print i
      print "Point ",j, " Energy= ",energy
      for (k=1;k<=i;k++) {
        print at_symbol[at[k]], x[k], y[k], z[k]
        }
      }
    }
}
