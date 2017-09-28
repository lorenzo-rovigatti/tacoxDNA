#scripts="/scratch/asuma/lammps_oxdna/long_runs/scripts"
scripts="$(pwd)"
cat $1 | awk 'NR>3{print}' > tmp
python $scripts/converter.py tmp 
rm tmp
