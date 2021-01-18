#!/bin/sh

# strict failur on errors
set -e
set -o pipefail

if [ $# -eq 0 ]
then
    gitBranch=$(git branch | sed -n '/^\* /s///p')
    echo "=============== Updating $gitBranch branch ==================="
    if [ $gitBranch == 'master' ]; then
        echo "                                __"
        echo "                               / _)"
        echo "                      _.----._/ /"
        echo "                     /         /"
        echo "                  __/ (  | (  |"
        echo "                 /__.-'|_|--|_|"
    elif [ $gitBranch == 'development' ]; then
        echo "                       ,''' "
        echo "                      /     \ "
        echo "                     :       :"
        echo "                     :       :"
        echo "                       .___,'"
    fi
	git fetch
	git pull
elif [ $# -eq 1 ]
then
    echo "=============== Updating $1 branch =========================="
    if [ $1 == 'master' ]; then
        echo "                                __"
        echo "                               / _)"
        echo "                      _.----._/ /"
        echo "                     /         /"
        echo "                  __/ (  | (  |"
        echo "                 /__.-'|_|--|_|"
    elif [ $1 == 'development' ]; then
        echo "                       ,''' "
        echo "                      /     \ "
        echo "                     :       :"
        echo "                     :       :"
        echo "                       .___,'"
    fi
	git checkout $1
	git fetch
	git pull
else
	printf "USAGE: ./update.sh Branch_Name\n"
	exit
fi
echo '=============== Updating libraries ======================='
./install.sh
