#!/bin/bash

function compressPdf_bash {
    #    -dPDFSETTINGS=/screen lower quality, smaller size. (72 dpi)
    #    -dPDFSETTINGS=/ebook for better quality, but slightly larger pdfs. (150 dpi)
    #    -dPDFSETTINGS=/prepress output similar to Acrobat Distiller "Prepress Optimized" setting (300 dpi)
    #    -dPDFSETTINGS=/printer selects output similar to the Acrobat Distiller "Print Optimized" setting (300 dpi)
    #    -dPDFSETTINGS=/default selects output intended to be useful across a wide variety of uses, possibly at the expense of a larger output file
    echo "$1 : $2"
    output=$2"_$1.pdf"
    gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/$1 -dNOPAUSE -dQUIET -dBATCH -sOutputFile=$output $2
}

function mergePdf {
    #gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=027_FbK_pre_ampli_01_all.pdf 027_FbK_pre_ampli_01.pdf FBKAmplifier-AmpliNEG-1.pdf ./FBK_ampl/001.pdf ./FBK_ampl/002.pdf ./FBK_ampl/image001.pdf ./FBK_ampl/image002.pdf 20230419164646169.pdf 
    #gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=027_FbK_pre_ampli_01_all_short.pdf 027_FbK_pre_ampli_01.pdf FBKAmplifier-AmpliNEG-1.pdf ./FBK_ampl/001.pdf 20230419164646169.pdf 
    gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=histSingle_179theta_0phi_xmx_xmy.rot.all.pdf histSingle_179theta_0phi_7mx_0my.root.pdf histSingle_179theta_0phi_0mx_7my.root.pdf histSingle_179theta_0phi_m5mx_5my.root.pdf histSingle_179theta_0phi_m7mx_m0my.root.pdf histSingle_179theta_0phi_m5mx_m5my.root.pdf histSingle_179theta_0phi_m0mx_m7my.root.pdf histSingle_179theta_0phi_5mx_5my.root.pdf
}

function printHelp { 
    echo " --> ERROR in input arguments <-- "
    echo " [0] -l  : Lower quality, smaller size (72 dpi)."
    echo " [0] -b  : Better quality, but slightly larger pdfs (150 dpi)."
    echo " [0] -p  : Output similar to Acrobat Distiller -Prepress Optimized- setting (300 dpi)."
    echo " [0] -o  : Output similar to Acrobat Distiller -Print Optimized- setting (300 dpi)."
    echo " [0] -d  : Output intended to be useful across a wide variety of uses, possibly at the expense of a larger output file."
    echo " [0] -a  : All."
    echo " [1]     : Input file."
    echo " [0] -m  : merge"    
}

if [ $# -eq 0 ] ; then   
    printHelp
elif [ $# -eq 1 ] ; then
    if [ "$1" = "-m" ] ; then
	mergePdf
    else [ "$1" = "-m" ]
	 printHelp
    fi
else
    if [ "$1" = "-l" ] ; then
	compressPdf_bash screen $2
    elif [ "$1" = "-b" ] ; then
	compressPdf_bash ebook $2
    elif [ "$1" = "-p" ] ; then
	compressPdf_bash prepress $2
    elif [ "$1" = "-o" ] ; then
	compressPdf_bash printer $2
    elif [ "$1" = "-d" ] ; then
	compressPdf_bash default $2
    elif [ "$1" = "-a" ] ; then
	compressPdf_bash screen $2
	compressPdf_bash ebook $2
	compressPdf_bash prepress $2
	compressPdf_bash printer $2
	compressPdf_bash default $2
    else	
	printHelp 
    fi
fi
