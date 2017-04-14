PATHTOCONV=/home/gian/Public/CERN/MG5_aMCcurr/ttbartoemu/Events
THISPATH=`pwd`
cd $PATHTOCONV
for numb in 1 2  # {11..17};
do
	rtu="run_"
	if [[ $numb -le 9 ]]
	then
		rtu=$rtu"0"
	fi
	rtu=${rtu}${numb};
	gunzip $rtu/unweighted_events.lhe.gz;
	echo "Converting events in "$rtu
	/home/gian/Public/CERN/MG5_aMCcurr/ExRootAnalysis/ExRootLHEFConverter $rtu/unweighted_events.lhe $rtu/unweighted_events.root;
	gzip  $rtu/unweighted_events.lhe
done
cd $THISPATH

