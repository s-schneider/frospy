#!/bin/bash
startup=0
lib='full'

wdir=`pwd`
path="`python -m site --user-site`"

while test $# -gt 0; do
  case "$1" in
    -h|--help)
      echo "frospy - attempt to capture input"
      echo " "
      echo "install.sh [options] application [arguments]"
      echo " "
      echo "options:"
      echo "-h, --help                show brief help"
      echo "-l, --lib=PATH            specify an subfolder path to install"
      echo "-startup                  install startup file"
      exit 0
      ;;
    -l)
      shift
      if test $# -gt 0; then
        lib=$1
        echo "Installing only: frospy/$lib"
        path=$path/frospy
      else
        echo "no lib specified"
        exit 1
      fi
      shift
      ;;
    --lib*)
      lib=`echo $1 | sed -e 's/^[^=]*=//g'`
      echo "Installing only: $lib"
      shift
      ;;
    -startup)
      shift
      startup=1
      shift
      ;;
    *)
      break
      ;;
  esac
done


echo '=============== Updating python =========================='
echo '                - Removing .pyc files of installation'


#Look for all .pyc files that may be left here and delete them
find . -type f -name "*.pyc" -delete
echo '                - Installing new files'

if [ "$lib" == "full" ]; then
    rm -rf $path
    echo "                  $wdir/frospy --> $path"
    rsync -r --delete frospy $path
    
else
    rm -rf $path/$lib
    echo "                  $wdir/$lib --> $path/frospy/$lib"
    rsync -r --delete frospy/$lib $path/$lib
fi

if (( startup == 1 )); then
    echo '=============== Updating startup file ===================='
    # rm ~/dev/.ipython/profile_default/startup/startup.py
    # cp startup.py ~/dev/.ipython/profile_default/startup/.
    rm ~/.ipython/profile_default/startup/startup.py
    cp startup.py ~/.ipython/profile_default/startup/.
fi

echo "=============== frospy succesful updated ==================="
echo ""
echo ""
