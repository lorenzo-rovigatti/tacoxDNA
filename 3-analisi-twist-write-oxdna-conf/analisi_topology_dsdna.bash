dir_script="/scratch/asuma/Project_dsDNA_through_nanopore/script_analisi_oxdna"
traj="output/trajectoryMD.dat"


#analisi linking number# va messo nel comando dello script python il numero 0
obs="linking_number.dat"

nsnap=$(cat output/trajectoryMD.dat | awk '{if($1=="t")print}' | awk '{if($3 % 5000000==0)print}' | wc -l)


#vediamo quali conf sono state elaborate #if($1 % 5000000==0)
cat $traj | grep "t =" | awk '{print $3}' | awk '{print}'   > _tmptime

mkdir Observables
touch Observables/$obs

for snap in $(cat Observables/$obs | awk '{print($1)}' > _tmptime2 ; diff _tmptime _tmptime2 | awk '{if ($1=="<") print($2)}')   

do
	echo $snap
	#extract #snap conf
	
	cat output/trajectoryMD.dat | awk -v snap=$snap '{if($1=="t")i=$3;if(i==snap)print $1,$2,$3}'  | tail -n +4 > _tmp
	
	
        LK=$(python $dir_script/analisi_topology_dsdna.py _tmp 0)
	
	echo $snap   $LK >> Observables/$obs






done

rm -f _tmptime _tmptime2
