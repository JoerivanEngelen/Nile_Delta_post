#!/bin/bash
#SBATCH -t 5-00:00:00
#SBATCH --output=stds/ND_%j.out
#SBATCH --error=stds/ND_%j.err
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 24

set -x
echo $TMPDIR

module load p7zip
module load gnu-parallel/2015.01.22

export PATH=$HOME/anaconda3/bin:$PATH
export PYTHONPATH=$HOME/anaconda3/lib/python3.6/site-packages:$PYTHONPATH

modname=insert_modname_here

runname=ND_${modname}_n${SLURM_JOB_NUM_NODES}m${SLURM_NTASKS}_${SLURM_JOB_ID}
outfol=$HOME/NileDelta/results/${runname}

mkdir -m 755 ${outfol}

cp -r $HOME/NileDelta/$modname/* $TMPDIR
cd $TMPDIR
begintime=$(date +%s%N)

echo 0 >> $TMPDIR/init_times.txt

max_nr=$(\ls *[0-9].RUN | sort -V | egrep -o [0-9]+ | tail -1)

7za a $TMPDIR/${runname}.7z $TMPDIR/* 
cp $TMPDIR/${runname}.7z ${outfol}/${runname}.7z

lastp=`expr ${SLURM_NTASKS} - 1`

srunexc="srun --exclusive -N1 -n1"
parallelov="parallel --delay .2 -j $SLURM_NTASKS"


for i in $( seq 1 $max_nr ); do
   starttime=$(date +%s%N)
   
   #Run simulation
   srun $HOME/seawat/bin/seawat_svn317 ${modname}_Runfile_${i}.RUN
   wait $!
   calctime=$(date +%s%N)
   echo "calctime: $(echo "scale=3;(${calctime} - ${starttime})/(1*10^09)" | bc) seconds" >> $HOME/NileDelta/stds/ND_${SLURM_JOB_ID}.log

   #Create initial conditions of next run as idfs
   mkdir -m 755 $TMPDIR/inits_$(($i+1))
   wait $!
   python $HOME/NileDelta/combine_subds_new.py $TMPDIR'/results_'${i}'/*' $TMPDIR/inits_$((i+1))
   wait $!
   initstime=$(date +%s%N)
   echo "initstime: $(echo "scale=3;(${initstime} - ${calctime})/(1*10^09)" | bc) seconds" >> $HOME/NileDelta/stds/ND_${SLURM_JOB_ID}.log   
   cp -r $TMPDIR/inits_$(($i+1)) ${outfol}   

   #Get initial time for postprocessing
   python $HOME/NileDelta/get_init_times.py "${TMPDIR}"'/results_*/head' $TMPDIR/init_times.txt ${i}

   #Creating some shortcuts as the next part is otherwise completely unreadable
   x=$(printf "%03d" $i)
   script=$HOME/NileDelta/output_to_netcdf_subdomain.py
   inpdata=$TMPDIR'/results_'${i}'/*/*_p'
   outdata=$TMPDIR/results_${x}_
   
   #Convert idfs to netcdf for each subdomain in parallel
   $parallelov 'set -o noglob ; p=$(printf "%03d" {}) ; $srunexc python '"${script} ${inpdata}"'${p}* '"${outdata}"'{}.nc '"$TMPDIR/init_times.txt ${i}" ::: $(seq 0 $lastp)
   wait $!
   converttime=$(date +%s%N)
   echo "converttime: $(echo "scale=3;(${converttime} - ${initstime})/(1*10^09)" | bc) seconds" >> $HOME/NileDelta/stds/ND_${SLURM_JOB_ID}.log
 
   #Combine these netcdfs to one big netcdf
   python $HOME/NileDelta/combine_netcdfs.py "$TMPDIR"'/results_'"${x}"'_*.nc' $TMPDIR/results_${x}.nc $TMPDIR/init_times.txt ${i}
   wait $!
   combinetime=$(date +%s%N)
   echo "combinetime: $(echo "scale=3;(${combinetime} - ${converttime})/(1*10^09)" | bc) seconds" >> $HOME/NileDelta/stds/ND_${SLURM_JOB_ID}.log
   
   #Remove files in parallel
   $parallelov 'p=$(printf "%03d" {}) ; $srunexc rm -rf '"$TMPDIR/results_${i}"'/*_p${p}*' ::: $(seq 0 $lastp)
   wait $!
   removetime=$(date +%s%N)
   echo "removetime: $(echo "scale=3;(${removetime} - ${combinetime})/(1*10^09)" | bc) seconds" >> $HOME/NileDelta/stds/ND_${SLURM_JOB_ID}.log

   cp $TMPDIR/results_${x}.nc ${outfol}
   wait $!
done

endtime=$(date +%s%N)

echo "total: $(echo "scale=3;(${endtime} - ${begintime})/(1*10^09)" | bc) seconds"

#max_init=$(\ls -d inits_+([0-9]) | sort -V | egrep -o [0-9]+ | tail -1) #This line works in interactive shell but not in bash script....
max_init=$(\ls -d inits_[0-9]* | sort -V | egrep -o [0-9]+ | tail -1)
rm -rf ${outfol}/inits_*

#cp $TMPDIR/*.list.* ${outfol}
cp -r $TMPDIR/inits_${max_init} ${outfol}

