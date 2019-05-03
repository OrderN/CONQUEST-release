#!/bin/csh
#
if (! $#argv) then
  echo "Usage: multi_vrml.csh vrmlfile1 vrmlfile2 ... vrmlfileN"
  echo " "
  echo "Output: multi_vrml_orbitals.wrl"
endif
#
set outfile = multi_vrml_orbitals.wrl
#
#
#
cat << EOF > $outfile
#VRML V2.0 utf8
 NavigationInfo { type "EXAMINE" }
 Viewpoint { position 0 0    3.0     }
EOF
#
#
#
set cnt = 0
foreach file ( $* )
@ cnt = $cnt + 1
cat $file | sed -n -e "/PROTO Surface/,/^  }/p" | sed -e "s/Surface/Surface$cnt/" | sed -e "/^     Group/,/^     DEF/d" >> $outfile
echo " }" >> $outfile
end
#
#
#
cat << EOF >> $outfile
 Group {
   children [
     DEF MOLDEN_TOUCH TouchSensor { },
     DEF MOLDENS Switch {
       whichChoice 0
       choice [
EOF
#
#
#
set cnt = 0
foreach file ( $* )
@ cnt = $cnt + 1
cat << EOF >> $outfile
          Group {
            children [
EOF
cat $file | sed -n -e '/^ *Surface/,$p' | sed -e "s/Surface/Surface$cnt/" >> $outfile
cat << EOF >> $outfile
            ]
          }
EOF
end
@ cnt = $cnt - 1
#
#
#
cat << EOF >> $outfile
       ]
     }
   ]
 }
  
 DEF MOLDEN_TIMER TimeSensor { cycleInterval 8 }
  
 PROTO SwitchInterpolator [
   eventIn SFFloat set_fraction
   eventOut SFInt32 value_changed
   field MFFloat key [ ]
   field MFFloat keyValue [ ]
 ]
 {
   DEF INTERP ScalarInterpolator {
     set_fraction IS set_fraction
     key IS key
     keyValue IS keyValue
   }
   DEF INTEGERIZE Script {
     eventIn SFFloat scalarValue
     eventOut SFInt32 value_changed IS value_changed
     url [
       "vrmlscript:
          function scalarValue(scalar) {
            value_changed = Math.floor(scalar);
          }",
          "Float2Int.class"
     ]
   }
   ROUTE INTERP.value_changed TO 
   INTEGERIZE.scalarValue
 }
 DEF MOLDEN_ANIMATOR SwitchInterpolator {
   key [0,1]
   keyValue [0, $cnt]
 }
 ROUTE MOLDEN_TOUCH.touchTime TO
 MOLDEN_TIMER.startTime
 ROUTE MOLDEN_TIMER.fraction_changed TO 
 MOLDEN_ANIMATOR.set_fraction
 ROUTE MOLDEN_ANIMATOR.value_changed TO 
 MOLDENS.whichChoice
EOF
#
