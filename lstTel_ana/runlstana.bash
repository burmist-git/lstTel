#!/bin/bash

#Analyse list of root files 
rootFilesList="./rootFileList.dat"
outHistF="./hist.root"

#Or analyse single root file 
inRootFiles="../data_sim/lstTel_00000_sim.root"
#inRootFiles="../data_sim/lstTel_00000_039718583_sim.root"
outHistSingleF="./histSingle.root"

make -f Makefilelstana clean; make -f Makefilelstana runlstana;

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d    : single root file"
    echo " [0] --d01 : single root file Loop01"
    echo " [0] -l    : list of root files"
    echo " [0] -s    : scan"
    echo " [0] -h    : print help"
}

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-d" ]; then
	./runlstana 1 $inRootFiles $outHistSingleF
    elif [ "$1" = "--d01" ]; then
	#./runlstana 2 ../lstTel-build/lstTel_179theta_0phi_7mx_0my.root histSingle_lstTel_179theta_0phi_7mx_0my.root
	#./runlstana 2 ../lstTel-build/lstTel_179theta_0phi_5mx_5my.root histSingle_lstTel_179theta_0phi_5mx_5my.root
	#./runlstana 2 ../lstTel-build/lstTel_179theta_0phi_0mx_7my.root histSingle_lstTel_179theta_0phi_0mx_7my.root
	#./runlstana 2 ../lstTel-build/lstTel_179theta_0phi_m5mx_5my.root histSingle_lstTel_179theta_0phi_m5mx_5my.root
	#./runlstana 2 ../lstTel-build/lstTel_179theta_0phi_m7mx_0my.root histSingle_lstTel_179theta_0phi_m7mx_0my.root
	#./runlstana 2 ../lstTel-build/lstTel_179theta_0phi_m5mx_m5my.root histSingle_lstTel_179theta_0phi_m5mx_m5my.root
	#./runlstana 2 ../lstTel-build/lstTel_179theta_0phi_0mx_m7my.root histSingle_lstTel_179theta_0phi_0mx_m7my.root
	#./runlstana 2 ../lstTel-build/lstTel_179theta_0phi_5mx_m5my.root histSingle_lstTel_179theta_0phi_5mx_m5my.root	
	#
	#./runlstana 2 ../lstTel-build/lstTel_180theta_0phi_7mx_0my.root histSingle_lstTel_180theta_0phi_7mx_0my.root
	#./runlstana 2 ../lstTel-build/lstTel_180theta_0phi_5mx_5my.root histSingle_lstTel_180theta_0phi_5mx_5my.root
	#./runlstana 2 ../lstTel-build/lstTel_180theta_0phi_0mx_7my.root histSingle_lstTel_180theta_0phi_0mx_7my.root
	#./runlstana 2 ../lstTel-build/lstTel_180theta_0phi_m5mx_5my.root histSingle_lstTel_180theta_0phi_m5mx_5my.root
	#./runlstana 2 ../lstTel-build/lstTel_180theta_0phi_m7mx_0my.root histSingle_lstTel_180theta_0phi_m7mx_0my.root
	#./runlstana 2 ../lstTel-build/lstTel_180theta_0phi_m5mx_m5my.root histSingle_lstTel_180theta_0phi_m5mx_m5my.root
	#./runlstana 2 ../lstTel-build/lstTel_180theta_0phi_0mx_m7my.root histSingle_lstTel_180theta_0phi_0mx_m7my.root
	#./runlstana 2 ../lstTel-build/lstTel_180theta_0phi_5mx_m5my.root histSingle_lstTel_180theta_0phi_5mx_m5my.root	
	#
	#./runlstana 2 ../lstTel-build/lstTel_180theta_0phi_0mx_0my.root histSingle_lstTel_180theta_0phi_0mx_0my.root
	#./runlstana 2 ../lstTel-build/lstTel_180theta_0phi_1mx_0my.root histSingle_lstTel_180theta_0phi_1mx_0my.root
	#./runlstana 2 ../lstTel-build/lstTel_180theta_0phi_2mx_0my.root histSingle_lstTel_180theta_0phi_2mx_0my.root
	#./runlstana 2 ../lstTel-build/lstTel_180theta_0phi_3mx_0my.root histSingle_lstTel_180theta_0phi_3mx_0my.root
	#./runlstana 2 ../lstTel-build/lstTel_180theta_0phi_4mx_0my.root histSingle_lstTel_180theta_0phi_4mx_0my.root
	#./runlstana 2 ../lstTel-build/lstTel_180theta_0phi_5mx_0my.root histSingle_lstTel_180theta_0phi_5mx_0my.root
	#./runlstana 2 ../lstTel-build/lstTel_180theta_0phi_6mx_0my.root histSingle_lstTel_180theta_0phi_6mx_0my.root
	#./runlstana 2 ../lstTel-build/lstTel_180theta_0phi_7mx_0my.root histSingle_lstTel_180theta_0phi_7mx_0my.root
	#
	#./runlstana 2 ../lstTel-build/lstTel_179theta_0phi_0mx_0my.root histSingle_lstTel_179theta_0phi_0mx_0my.root
	#./runlstana 2 ../lstTel-build/lstTel_179theta_0phi_1mx_0my.root histSingle_lstTel_179theta_0phi_1mx_0my.root
	#./runlstana 2 ../lstTel-build/lstTel_179theta_0phi_2mx_0my.root histSingle_lstTel_179theta_0phi_2mx_0my.root
	#./runlstana 2 ../lstTel-build/lstTel_179theta_0phi_3mx_0my.root histSingle_lstTel_179theta_0phi_3mx_0my.root
	#./runlstana 2 ../lstTel-build/lstTel_179theta_0phi_4mx_0my.root histSingle_lstTel_179theta_0phi_4mx_0my.root
	#./runlstana 2 ../lstTel-build/lstTel_179theta_0phi_5mx_0my.root histSingle_lstTel_179theta_0phi_5mx_0my.root
	#./runlstana 2 ../lstTel-build/lstTel_179theta_0phi_6mx_0my.root histSingle_lstTel_179theta_0phi_6mx_0my.root
	#./runlstana 2 ../lstTel-build/lstTel_179theta_0phi_7mx_0my.root histSingle_lstTel_179theta_0phi_7mx_0my.root
	#
	./runlstana 2 ../lstTel-build/lstTel_179.5theta_0phi_0mx_0my.root histSingle_lstTel_179.5theta_0phi_0mx_0my.root
	./runlstana 2 ../lstTel-build/lstTel_179.5theta_0phi_1mx_0my.root histSingle_lstTel_179.5theta_0phi_1mx_0my.root
	./runlstana 2 ../lstTel-build/lstTel_179.5theta_0phi_2mx_0my.root histSingle_lstTel_179.5theta_0phi_2mx_0my.root
	./runlstana 2 ../lstTel-build/lstTel_179.5theta_0phi_3mx_0my.root histSingle_lstTel_179.5theta_0phi_3mx_0my.root
	./runlstana 2 ../lstTel-build/lstTel_179.5theta_0phi_4mx_0my.root histSingle_lstTel_179.5theta_0phi_4mx_0my.root
	./runlstana 2 ../lstTel-build/lstTel_179.5theta_0phi_5mx_0my.root histSingle_lstTel_179.5theta_0phi_5mx_0my.root
	./runlstana 2 ../lstTel-build/lstTel_179.5theta_0phi_6mx_0my.root histSingle_lstTel_179.5theta_0phi_6mx_0my.root
	./runlstana 2 ../lstTel-build/lstTel_179.5theta_0phi_7mx_0my.root histSingle_lstTel_179.5theta_0phi_7mx_0my.root
    elif [ "$1" = "-l" ]; then
	./runlstana 0 $rootFilesList $outHistF
    elif [ "$1" = "-s" ]; then
	rm -rf info.log
	touch info.log
	for i in $(seq 1 100);
	do
	    root -l -b -q fit_muon_ring.C\($i\) | grep info: | tee -a info.log
	done
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi

#espeak "I have done"
