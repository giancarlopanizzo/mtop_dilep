for lhaid in 13200 21000 25000 262000;
do
       	source ../../mtop_dilep/createmomentlists.sh $lhaid;
done
WRUN=$1
rm derivations.*; 
for files in `grep "$WRUN" mt.moments.*.list | awk -F : '{print $2}'`; 
do 
	bsub -q short -e derivations.err -o derivations.out -cwd `pwd` source ../../mtop_dilep/run_derivations.sh /gpfs/atlas/gpanizzo $files; 
done
