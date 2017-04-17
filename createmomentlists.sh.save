WPTH=/gpfs/atlas/gpanizzo/MG5_aMCcurr/ttbarLOemu_
lhaid=$1
ls ${WPTH}${lhaid}_1*/Events/run_0?/testscripts_delphes_events.root | grep -v "run_01" > mt.moments.${lhaid}.list
NUMFILES=`ls ${WPTH}${lhaid}_1*/Events/run_0?/testscripts_delphes_events.root | grep -v -c "run_01"`
NUMLINES=$((NUMFILES / 8))
split -l ${NUMLINES} -d mt.moments.${lhaid}.list mt.moments.${lhaid}.list.
