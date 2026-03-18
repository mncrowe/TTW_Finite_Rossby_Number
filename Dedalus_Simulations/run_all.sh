####################################################################
# Arguments:
# 1 - Ro        : Rossby Number (1)
# 2 - E         : Ekman Number (0.1)
# 3 - Pr        : Prandtl Number (1)
####################################################################

savename=snapshots
Ro=(0.01 0.025 0.05 0.075 0.1 0.25 0.5 0.75 1.0 2.5 5 7.5 10.0)
E=(0.001 0.01 0.1 1.0)
Pr=1

if [ ! -d "$data" ]; then
  mkdir data
fi

source ../activate_D3.sh

for i in ${!Ro[@]}; do
  for j in ${!E[@]};   do
    #echo "${savename}_${Ro[$i]}_${E[$j]}" ${Ro[$i]} ${E[$j]}
    mpiexec -np 32 python3 TTW.py "${savename}_${Ro[$i]}_${E[$j]}" ${Ro[$i]} ${E[$j]}
  done
done

mv $savename* data