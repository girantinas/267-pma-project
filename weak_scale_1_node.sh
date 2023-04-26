cmake --build .
module load contrib/1.0
module load upcxx/2023.3.0
export GASNET_BACKTRACE=1
for j in $(seq 1 1 4); do 
    echo "ranks: 1" >> weak_scale_1_node_log
    salloc -N 1 -A mp309 -t 10:00 --qos=debug -C cpu srun -N 1 -n 1 ./main_dist_pcsr ../tests/dist_pcsr_inserts_speed.txt >> weak_scale_1_node_log
    echo "ranks: 2" >> weak_scale_1_node_log
    salloc -N 1 -A mp309 -t 10:00 --qos=debug -C cpu srun -N 1 -n 2 ./main_dist_pcsr ../tests/dist_pcsr_inserts_speed.txt >> weak_scale_1_node_log
    echo "ranks: 4" >> weak_scale_1_node_log
    salloc -N 1 -A mp309 -t 10:00 --qos=debug -C cpu srun -N 1 -n 4 ./main_dist_pcsr ../tests/dist_pcsr_inserts_speed.txt >> weak_scale_1_node_log
    echo "ranks: 8" >> weak_scale_1_node_log
    salloc -N 1 -A mp309 -t 10:00 --qos=debug -C cpu srun -N 1 -n 8 ./main_dist_pcsr ../tests/dist_pcsr_inserts_speed.txt >> weak_scale_1_node_log
    echo "ranks: 16" >> weak_scale_1_node_log
    salloc -N 1 -A mp309 -t 10:00 --qos=debug -C cpu srun -N 1 -n 16 ./main_dist_pcsr ../tests/dist_pcsr_inserts_speed.txt >> weak_scale_1_node_log
    echo "ranks: 32" >> weak_scale_1_node_log
    salloc -N 1 -A mp309 -t 10:00 --qos=debug -C cpu srun -N 1 -n 32 ./main_dist_pcsr ../tests/dist_pcsr_inserts_speed.txt >> weak_scale_1_node_log
    echo "ranks: 64" >> weak_scale_1_node_log
    salloc -N 1 -A mp309 -t 10:00 --qos=debug -C cpu srun -N 1 -n 64 ./main_dist_pcsr ../tests/dist_pcsr_inserts_speed.txt >> weak_scale_1_node_log
done