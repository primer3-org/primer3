export PATH=/usr/bin:$PATH

res1=$(date +%s)

# Prints OK for exit status 0, FAILED otherwise
result() {
    if [ "$1" -eq 0 ]; then
        echo OK
    else
        echo FAILED
    fi
}

echo; echo; echo 'STARTING 1ST SET OF TESTS (-O2)'; echo; echo;
make clean
make test
r1=$?
make clean
echo; echo; echo 'STARTING 2ND SET OF TESTS (-O3)'; echo; echo;
make O_OPTS=-O3 && make test
r2=$?
make clean
echo; echo; echo 'STARTING 3RD SET OF TESTS (-O0 and valgrind)'; echo; echo;
make O_OPTS=-O0
r3=$?
cd ../test
if [ $r3 -eq 0 ]; then
    make TESTOPTS=--valgrind
    r3=$?
fi

# Whole seconds, computed by the shell (no bc needed)
dt=$(( $(date +%s) - res1 ))

# Each set prints its own summary above; this only collects the results
echo; echo; echo 'PARANOID TEST RESULTS'; echo
echo "-O2:      $(result $r1)"
echo "-O3:      $(result $r2)"
echo "valgrind: $(result $r3)"
printf "Total runtime: %d:%02d:%02d:%02d\n" \
    $((dt / 86400)) $((dt % 86400 / 3600)) $((dt % 3600 / 60)) $((dt % 60))

if [ $r1 -ne 0 ] || [ $r2 -ne 0 ] || [ $r3 -ne 0 ]; then
    exit 1
fi
