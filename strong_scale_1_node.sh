cmake --build .
module load contrib/1.0
module load upcxx/2023.3.0
export GASNET_BACKTRACE=1
for j in $(seq 1 1 4); do 
    for i in $(seq 0 1 6); do
        cmake --build .
        num=$((2**i))
        echo "ranks: $num" >> strong_scale_large_1_node_log
        salloc -N 1 -A mp309 -t 10:00 --qos=debug -C cpu srun -N 1 -n $num ./main_dist_pcsr ../tests/dist_pcsr_inserts_speed.txt >> strong_scale_large_1_node_log
    done
done