queue=8nh

dir=$1 && shift

wildcard=\*
[[ -n $1 ]] && wildcard=$1 && shift

for f in ${dir}/${wildcard}.dat; do
    if [[ -n $1 ]]; then
	for i in $(seq 0 $(($1-1))); do
	    rm -f ${f}_${i}.log
	    bsub -q $queue -o ${f}_${i}.log run.sh -- ./reduce.py $PWD/$f $1 $i
	done
    else
	rm -f $f.log
	bsub -q $queue -o $f.log run.sh -- ./reduce.py $PWD/$f
    fi
done
