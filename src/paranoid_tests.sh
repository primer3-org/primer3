echo; echo; echo 'STARTING 1ST SET OF TESTS (-O2)'; echo; echo;
make clean
make test
make clean
echo; echo; echo 'STARTING 2ND SET OF TESTS (-O3)'; echo; echo;
make O_OPTS=-O3
make test
make clean
echo; echo; echo 'STARTING 3RD SET OF TESTS (-O0) and valgrind'; echo; echo;
make O_OPTS=-O0
cd ../test; make VAL=--valgrind
