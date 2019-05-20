#!/bin/sh
#
file_envar='environ_cq_tool'
file_path='.'
#
# test and source for environment variables tool file
#
if test -r $file_path/$file_envar ; then
   echo $file_envar' file found ; source it'
   . $file_path/$file_envar

else
   echo $file_envar 'file not found: check your repository'
   exit 1

fi
#
cq_tool='cq_'$tool
#
echo ''
echo '######################################'
echo '#' 
echo '#        '$cq_tool' installer         '
echo '#  (as part of the Conquest package)  '
echo '#'
echo '######################################'
echo ''
#
# check if tool is already available in the current directory
#
echo '> need:' $tool' v'$tool_version ; echo
#
tool_tar=$tool$tool_version.tar
tool_targz=$tool_tar.gz
tools_list="$tool_tar $tool_targz"
#cq_base=$(basename "$cq_tool" .tar.gz)
#
for file in $tools_list ; do
   if test ! -r $file ; then
      echo "$file not existent or not readable"
      status=0
   else 
      echo "$file found"   
      status=1
   fi
done
#
# ...if not try to download
#

if [ $status -eq 1 ] ; then
   echo 'untar '$file'...\c'
   tar -xf $file
   check_failure $? ; echo 'done' 

else
   if test "$grab" = "" ; then
      echo 'cannot download ; need wget or curl'
      echo 'Aborting' 
      exit 1

   else
      echo $tool' v'$tool_version' will be download' 
      echo "attempt to grab $tool_targz"
      $grab $tool_targz $tool_url$tool_targz
   fi
fi
#
# untar and rename source directory 
#
src_dir='src/'
#
if test -d $src_dir ; then
    echo 
    echo "ERROR: $src_dir already exists ; delete it and re-run installer"
    echo "Aborting"
    exit 1
fi
#
echo 'untar '$tool_targz'...\c'
tar -xf $tool_targz
mv $tool$tool_version $src_dir
check_failure $? ; echo 'done'  
#
# test and copy patch file 
#
file_patch=$tool$tool_version'_cq.patch'
#
if test -r $file_path/$file_patch ; then
   echo $file_patch' file found ; copy it'
   cp $file_path/$file_patch $src_dir 

else
   echo $file_patch 'file not found: check your repository'
   exit 1

fi
#
# apply patch
#
if test "$cq_patch" = "" ; then
   echo 'cannot apply '$tool$tool_version'_cq patch'
   echo 'Aborting' 
   exit 1
else
   cd $src_dir 
   echo 'apply patch '$file_patch'...'
   patch < $file_patch 
fi 
#
cd ..
#
# cleaning...
#
cleaning_list="$tool_targz $tool_tar" 
for file in $cleaning_list ; do 
   if test -r $file ; then
      rm -f $file 
   fi
done

