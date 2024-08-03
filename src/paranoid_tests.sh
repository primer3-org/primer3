export PATH=/usr/bin:$PATH

res1=$(date +%s.%N)

echo; echo; echo 'STARTING 1ST SET OF TESTS (-O2)'; echo; echo;
make clean
make test
make clean
echo; echo; echo 'STARTING 2ND SET OF TESTS (-O3)'; echo; echo;
make O_OPTS=-O3
make test
make clean
echo; echo; echo 'STARTING 3RD SET OF TESTS (-O0 and valgrind)'; echo; echo;
make O_OPTS=-O0
cd ../test
make TESTOPTS=--valgrind

res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)

LC_NUMERIC=C printf "Total runtime: %d:%02d:%02d:%02.4f\n" $dd $dh $dm $ds
