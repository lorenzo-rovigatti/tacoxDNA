#script equilibrazione

start=1
for i in $(seq $start 10)
do



    cat > pbs.script.$i << EOF
#!/bin/bash
#PBS -N "xxdirxx_xxxRxxx_${i}"
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=20
#PBS -q regular
#PBS -m n
#PBS -o cluster.dat
#PBS -j oe
# start from launch dir
cd \$PBS_O_WORKDIR

module load openmpi/1.8.3/gnu/4.9.2

if [ $i -gt 1 ]; then
  bash restart_script.x
fi

mpirun -np 20 /scratch/asuma/lammps_oxdna/long_runs/scripts/lmp_mpi -in ./oxdna.lammps  -log none -echo screen >> Log.dat &

id1=\$!
sleep 20m

flag=0

while [ "\$flag" -eq "0" ]
do
	ls -v oxdna.restart.* | head -n -2 |xargs rm
	sleep 20m
	if ps -p \$id1 >&-; then 
		flag=0
	else
		flag=1
	fi
	
done


EOF

if [ $i -eq $start ]; then
  job=`qsub pbs.script.$i`
else
    job_next=`qsub -W depend=afterany:$job pbs.script.$i`
    job=$job_next
fi

rm pbs.script.$i
done

