Design PCR primers from DNA sequence. Widely used (190k Google hits for "primer3").
From mispriming libraries to sequence quality data to the generation of internal
oligos, primer3 does it. C&perl.

Building
--------

Install the necessary dependencies. For example, on Debian and Ubuntu,

```
sudo apt-get install -y build-essential g++ git-all
```

Then

```
git clone https://github.com/primer3-org/primer3.git primer3
cd primer3
make
make test
```

Run Primer3
-----------

```
./src/primer3_core example
```

Installing
----------

After successfully building and testing `primer3`, you can install the program
to destinary directory `$DESTDIR` (e.g. `/usr/local/`) by

```
DESTDIR=/usr/local make install
```

Read the complete Primer3 manual
--------------------------------

[Primer3 Manual](http://primer3.org/manual.html)

or see `./src/primer3_manual.htm`
