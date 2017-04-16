lhaid=$1
for numi in {0..7}; 
do 
	bsub -q short -cwd `pwd` -e mom.${lhaid}.${numi}.log -oo mom.${lhaid}.${numi}.log source ../../mtop_dilep/run_delphesmoments.sh /gpfs/atlas/gpanizzo 1 $numi $lhaid; 
done
