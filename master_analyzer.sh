for lhaid in 13200 21000 25000 262000; 
do
 	source ../../mtop_dilep/createmomentlists.sh $lhaid; 
 	source ../../mtop_dilep/par_run_delphmom.sh $lhaid; 
done
