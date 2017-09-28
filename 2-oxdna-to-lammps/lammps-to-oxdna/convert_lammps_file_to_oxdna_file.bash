cat $1 | awk 'NR>=10{print $3,$4,$5,$6,$7,$8,$9}' > tmp
python converter_lammps_to_oxdna.py tmp 
