#!/bin/bash

#Analyse list of root files 
rootFilesList="./rootFileList.dat"
outHistF="./hist.root"

#Or analyse single root file 
#inRootFiles="../data_sim/lstTel_00000_sim.root"
inRootFiles="../data_sim/lstTel_00000_039718583_sim.root"
outHistSingleF="./histSingle.root"

make -f Makefilelstana clean; make -f Makefilelstana runlstana;

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d  : single root file"
    echo " [0] -l  : list of root files"
    echo " [0] -h  : print help"
}

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-d" ]; then
	./runlstana 1 $inRootFiles $outHistSingleF
    elif [ "$1" = "-l" ]; then
	./runlstana 0 $rootFilesList $outHistF
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi

#espeak "I have done"
