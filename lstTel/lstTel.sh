#!/bin/bash

n_jobs=14

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d               : default"
    echo " [0] -v               : vis"
    echo " [0] --screen         : run screen jobs"
    echo " [0] --screen_monitor : monitor screen jobs"
    echo " [0] --screen_job     : run single screen job"    
    echo "  [1] showerID        : shower ID"
    echo " [0] -h               : print help"
}

function run_screen_job_sh {
    rm -rf lstTel000*
    rm -rf screenlog.0
    for i in $(seq 0 $n_jobs); do
	jobID=`printf "%05d" $i`
	screenName='lst'$jobID
        echo "$jobID/$n_jobs $screenName"
	#$PWD/lstTel.sh --screen_job $jobID
	screen -S $screenName -L -d -m $PWD/lstTel.sh --screen_job $jobID
        sleep 1
    done
}


if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-d" ]; then
	#./lstTel run.mac 1232 lstTel_180theta_0phi_0mx_0my.root mu- 100 180.0 0.0    0.0    0.0
	#./lstTel run.mac 1233 lstTel_180theta_0phi_1mx_0my.root mu- 100 180.0 0.0  100.0    0.0
	#./lstTel run.mac 121232 lstTel_180theta_0phi_2mx_0my.root mu- 100 180.0 0.0  200.0    0.0
	#./lstTel run.mac 123 lstTel_180theta_0phi_3mx_0my.root mu- 100 180.0 0.0  300.0    0.0
	#./lstTel run.mac 123 lstTel_180theta_0phi_4mx_0my.root mu- 100 180.0 0.0  400.0    0.0
	#./lstTel run.mac 123 lstTel_180theta_0phi_5mx_0my.root mu- 100 180.0 0.0  500.0    0.0
	#./lstTel run.mac 12313 lstTel_180theta_0phi_6mx_0my.root mu- 100 180.0 0.0  600.0    0.0
	#./lstTel run.mac 123 lstTel_180theta_0phi_7mx_0my.root mu- 100 180.0 0.0  700.0    0.0
	#
	#./lstTel run.mac 1232 lstTel_179theta_0phi_0mx_0my.root mu- 100 179.0 0.0    0.0    0.0
	#./lstTel run.mac 1233 lstTel_179theta_0phi_1mx_0my.root mu- 100 179.0 0.0  100.0    0.0
	#./lstTel run.mac 121232 lstTel_179theta_0phi_2mx_0my.root mu- 100 179.0 0.0  200.0    0.0
	#./lstTel run.mac 123 lstTel_179theta_0phi_3mx_0my.root mu- 100 179.0 0.0  300.0    0.0
	#./lstTel run.mac 123 lstTel_179theta_0phi_4mx_0my.root mu- 100 179.0 0.0  400.0    0.0
	#./lstTel run.mac 123 lstTel_179theta_0phi_5mx_0my.root mu- 100 179.0 0.0  500.0    0.0
	#./lstTel run.mac 12313 lstTel_179theta_0phi_6mx_0my.root mu- 100 179.0 0.0  600.0    0.0
	#./lstTel run.mac 123 lstTel_179theta_0phi_7mx_0my.root mu- 100 179.0 0.0  700.0    0.0
	#
	#./lstTel run.mac 1232 lstTel_179.5theta_0phi_0mx_0my.root mu- 100 179.5 0.0    0.0    0.0
	#./lstTel run.mac 1233 lstTel_179.5theta_0phi_1mx_0my.root mu- 100 179.5 0.0  100.0    0.0
	#./lstTel run.mac 121232 lstTel_179.5theta_0phi_2mx_0my.root mu- 100 179.5 0.0  200.0    0.0
	#./lstTel run.mac 123 lstTel_179.5theta_0phi_3mx_0my.root mu- 100 179.5 0.0  300.0    0.0
	#./lstTel run.mac 123 lstTel_179.5theta_0phi_4mx_0my.root mu- 100 179.5 0.0  400.0    0.0
	#./lstTel run.mac 123345 lstTel_179.5theta_0phi_5mx_0my.root mu- 100 179.5 0.0  500.0    0.0
	./lstTel run.mac 123133 lstTel_179.5theta_0phi_6mx_0my.root mu- 100 179.5 0.0  600.0    0.0
	./lstTel run.mac 123 lstTel_179.5theta_0phi_7mx_0my.root mu- 100 179.5 0.0  700.0    0.0
    elif [ "$1" = "-v" ]; then
	./lstTel vis.mac 123 lstTel_vis_0mx_0my.root mu- 100 180.0 0.0 0.0 0.0
	#./lstTel vis.mac 123 lstTel_vis_1mx_0my.root mu- 100 180.0 0.0 100.0 0.0
	#./lstTel vis.mac 123 lstTel_vis_2mx_0my.root mu- 100 180.0 0.0 200.0 0.0
	#./lstTel vis.mac 123 lstTel_vis_3mx_0my.root mu- 100 180.0 0.0 300.0 0.0
	#./lstTel vis.mac 123 lstTel_vis_4mx_0my.root mu- 100 180.0 0.0 400.0 0.0
	#./lstTel vis.mac 123 lstTel_vis_5mx_0my.root mu- 100 180.0 0.0 500.0 0.0
	#./lstTel vis.mac 123 lstTel_vis_6mx_0my.root mu- 100 180.0 0.0 600.0 0.0
	#./lstTel vis.mac 123 lstTel_vis_7mx_0my.root mu- 100 180.0 0.0 700.0 0.0
	#./lstTel vis.mac 123 lstTel_vis_8mx_0my.root mu- 100 180.0 0.0 800.0 0.0
	#./lstTel vis.mac 123 lstTel_vis_9mx_0my.root mu- 100 180.0 0.0 900.0 0.0
	#./lstTel vis.mac 123 lstTel_vis_10mx_0my.root mu- 100 180.0 0.0 1000.0 0.0
	#./lstTel vis.mac 123 lstTel_vis_11mx_0my.root mu- 100 180.0 0.0 1100.0 0.0
	#./lstTel vis.mac 123 lstTel_vis_12mx_0my.root mu- 100 180.0 0.0 1200.0 0.0
	#./lstTel vis.mac 123 lstTel_vis_13mx_0my.root mu- 100 180.0 0.0 1300.0 0.0
	#./lstTel vis.mac 123 lstTel_vis_14mx_0my.root mu- 100 180.0 0.0 1400.0 0.0
	#./lstTel vis.mac 123 lstTel_vis_15mx_0my.root mu- 100 180.0 0.0 1500.0 0.0
    elif [ "$1" = "--screen" ]; then
	run_screen_job_sh
    elif [ "$1" = "--screen_monitor" ]; then
	grep -i event screenlog.0 | grep Collection | grep Optical | wc -l
	screen -ls
	du -hs lstTel_000*_sim.root
    elif [ "$1" = "--screen_job" ]; then
	if [ $# -eq 2 ]; then
	    rndseed=$(date +%N)
	    jobID=$2
	    outname="lstTel_"$jobID"_"$rndseed"_sim.root"
	    echo "rndseed = $rndseed"
	    echo "jobID   = $jobID"
	    echo "outname = $outname"
	    #
	    cp lstTel lstTel$jobID
	    ./lstTel$jobID run.mac $rndseed $outname mu- 12 180.0 0
        else
            printHelp
        fi
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi
