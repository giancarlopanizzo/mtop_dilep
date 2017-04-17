for lhaid in 13200 21000 25000 262000
do
	printf "Combining results for lhaid "${lhaid}" ..." 
	source ../../mtop_dilep/run_delphesmoments.sh /gpfs/atlas/gpanizzo 2 0 $lhaid &> moments.${lhaid}.log
	echo " done."
done
