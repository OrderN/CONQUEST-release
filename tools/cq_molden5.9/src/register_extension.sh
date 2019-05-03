#!/bin/bash
APP="gmolden"
MP="MimeType="
EXTS="pdb xyz mol2 ogl fdat arc ins molf irc zmt msi-car freq shelx gulp mopac-input"
COMMENT="$APP's data file"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Create directories if missing
mkdir -p ~/.local/share/mime/packages
mkdir -p ~/.local/share/applications
mkdir -p ~/.icons

if [ -f ~/.icons/gmolden.png ]; then 
   echo " "
   echo "File extensions already registered, run again (y/n) ?"
   echo " "
   read ANSWER
   if [ "$ANSWER" = "n" ]; then
      exit
   fi
fi

NEXTS=""
OKEXTS=""

for EXT in $EXTS; do
  touch /tmp/file.$EXT
  OF=`gvfs-info /tmp/file.$EXT | grep standard::content-type | awk '{ print $(NF) }'`
  if [ $OF != 'text/plain' ]
  then
     if ! gvfs-mime --query $OF | grep 'No recommended applications' 
     then
        NEXTS+=$EXT" "
     else
        OKEXTS+=$EXT" "
     fi
  else
     OKEXTS+=$EXT" "
  fi
done

echo " "
echo 'Overriding existing associations:'
echo  $NEXTS
echo " "
echo "Continue (y/n) ?"
echo " "
read ANSWER
if [ "$ANSWER" = "n" ]; then
   exit
fi

for EXT in $EXTS; do
# Create mime xml 
echo "<?xml version=\"1.0\" encoding=\"UTF-8\"?>
<mime-info xmlns=\"http://www.freedesktop.org/standards/shared-mime-info\">
    <mime-type type=\"application/x-$EXT\">
        <comment>$COMMENT</comment>
        <icon name=\"application-x-$EXT\"/>
        <glob pattern=\"*.$EXT\"/>
    </mime-type>
</mime-info>" > ~/.local/share/mime/packages/application-x-$EXT.xml
MP=$MP"application/x-$EXT;chemical/x-$EXT;"
cp gmolden.png ~/.icons/application-x-$EXT.png
cp gmolden.png ~/.icons/chemical-x-$EXT.png
done

# Create application desktop
echo "[Desktop Entry]
Name=$APP
Exec=$DIR/$APP %U
$MP
Icon=$HOME/.icons/$APP.png
Terminal=false
Type=Application
StartupNotify=true
Categories=
Comment=
"> ~/.local/share/applications/$APP.desktop

# update databases for both application and mime
update-desktop-database ~/.local/share/applications
update-mime-database    ~/.local/share/mime

sudo cp ~/.local/share/applications/$APP.desktop /usr/share/app-install/desktop/$APP:$APP.desktop
#sudo cp gmolden.png /usr/share/app-install/icons/
#sudo cp gmolden.png /usr/share/pixmaps
#sudo gtk-update-icon-cache /usr/share/app-install
#sudo update-app-install
# copy associated icons to pixmaps
