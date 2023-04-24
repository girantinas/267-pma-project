# UPCXX_CODEMODE=opt upcxx -std=c++17 -Wall -Wno-sign-compare -o a.out main_dist_pcsr.cpp && upcxx-run -n 2 ./a.out ./tests/pcsr_benchmark_inserts.txt ./tests/pcsr_benchmark_outserts.txt > log.txt && diff tests/pcsr_benchmark_outserts.txt tests/pcsr_benchmark_solution.txt > diff_log.txt
cmake --build .
export GASNET_BACKTRACE=1
export GASNET_OFI_RECEIVE_BUFF_SIZE=single
export GASNET_BACKTRACE_SIGNAL=SIGINT
rm test*.dat
rm final*.txt
rm redistribute*.txt
salloc -N 2 -A mp309 -t 10:00 -q debug --qos=interactive -C cpu srun -N 2 -n 8 ./main_dist_pcsr ../tests/dist_pcsr_inserts.txt test
rm ../tests/dist_pcsr_outserts.txt
rm ../tests/sorted_dist_pcsr_solution.txt
cat test*.dat | sort > ../tests/dist_pcsr_outserts.txt
cat ../tests/dist_pcsr_solution.txt | sort > ../tests/sorted_dist_pcsr_solution.txt
diff ../tests/sorted_dist_pcsr_solution.txt ../tests/dist_pcsr_outserts.txt > ../diff-log.txt
