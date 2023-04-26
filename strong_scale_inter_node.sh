cmake --build .
module load contrib/1.0
module load upcxx/2023.3.0
export GASNET_BACKTRACE=1
export GASNET_OFI_RECEIVE_BUFF_SIZE=single
for i in $(seq 0 1 3); do
    for j in $(seq 1 1 4); do 
        num=$((2**i))
        echo "nodes: $num" >> strong_scale_large_inter_node_log
        # srun -N $num -n 64 ./main_dist_pcsr ../tests/dist_pcsr_inserts_speed.txt >> strong_scale_large_inter_node_log
        salloc -N $num -A mp309 -t 10:00 --qos=debug -C cpu srun -N $num -n 64 ./main_dist_pcsr ../tests/dist_pcsr_inserts_speed.txt >> strong_scale_large_inter_node_log
    done
done