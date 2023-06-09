# UPCXX_CODEMODE=opt upcxx -std=c++17 -Wall -Wno-sign-compare -o a.out main_dist_pcsr.cpp && upcxx-run -n 2 ./a.out ./tests/pcsr_benchmark_inserts.txt ./tests/pcsr_benchmark_outserts.txt > log.txt && diff tests/pcsr_benchmark_outserts.txt tests/pcsr_benchmark_solution.txt > diff_log.txt
module load contrib/1.0
module load upcxx/2023.3.0
cmake --build .
# export GASNET_BACKTRACE=1
export GASNET_OFI_RECEIVE_BUFF_SIZE=single
# export GASNET_BACKTRACE_SIGNAL=SIGUSR1
export UPCXX_SEGMENT_MB=512

# salloc -N 16 -A mp309 -t 5:00 -q regular -C cpu srun -N 16 -n 1024 ./main_dist_pcsr
srun -N 2 -n 128 ./main_dist_pcsr
