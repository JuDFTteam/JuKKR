# README
# compile with debug options and run valgrind with the following options

valgrind --trace-children=yes --log-file=valgrind_out.txt --max-stackframe=99999999 --main-stacksize=99999999 --valgrind-stacksize=10485760 --leak-check=full --leak-resolution=high --show-reachable=yes --track-origins=yes ./kkr.x

