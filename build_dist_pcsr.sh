# upcxx -std=c++17 -Wall -o a.out main_dist_pcsr.cpp && upcxx-run -n 2 ./a.out ./tests/pcsr_benchmark_inserts.txt ./tests/pcsr_benchmark_outserts.txt > log.txt && diff tests/pcsr_benchmark_outserts.txt tests/pcsr_benchmark_solution.txt > diff_log.txt
cmake --build .
module load contrib/1.0
module load upcxx/2023.3.0
export GASNET_BACKTRACE=1
export GASNET_OFI_RECEIVE_BUFF_SIZE=single
rm test*.dat
rm test*.bfs
rm test*.pr
salloc -N 4 -A mp309 -t 60:00 -q debug --qos=interactive -C cpu srun -N 4 -n 256 ./main_dist_pcsr ../tests/rmat-tests/tests/rmat-inserts >> log.txt #test
rm ../tests/dist_pcsr_outserts.txt
rm ../tests/dist_pcsr_bfs_outserts.txt
rm ../tests/dist_pcsr_pr_outserts.txt
rm ../tests/sorted_dist_pcsr_solution.txt
rm ../tests/sorted_dist_pcsr_bfs_solution.txt
cat test*.dat | sort > ../tests/dist_pcsr_outserts.txt
cat test*.bfs | sort > ../tests/dist_pcsr_bfs_outserts.txt
cat test*.pr | sort > ../tests/dist_pcsr_pr_outserts.txt
cat ../tests/dist_pcsr_solution.txt | sort > ../tests/sorted_dist_pcsr_solution.txt
cat ../tests/dist_pcsr_bfs_solution.txt | sort > ../tests/sorted_dist_pcsr_bfs_solution.txt
diff ../tests/sorted_dist_pcsr_solution.txt ../tests/dist_pcsr_outserts.txt > ../diff-log.txt
diff ../tests/sorted_dist_pcsr_bfs_solution.txt ../tests/dist_pcsr_bfs_outserts.txt > ../bfs-diff-log.txt