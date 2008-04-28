OS X Universal Binary Installation Instructions (bcf; 4/22/2008)


How to install this software
============================

1.  Double click on the .tar.gz file to extract the archive.

2.  The binary files are located in the 'bin' [for 'binary'] folder

3.  (Optional) To run the tests, cd to the new directory and then the test folder

4.  (Optional) Within this folder run:
	a. `perl p3test.pl`
	
5.  (Optional) You should not see 'FAILED' during the tests.

6.  (Optional) *NOTE*:  If your perl command is not called perl 
	(for example, if it is called perl5) you will have to modify 
	the internals of the test scripts).
	
7.	Copy the following files to a location of your choice:
	a.  bin/long_seq_tm_test
	b.  bin/ntdpal
	c.  bin/oligotm
	d.  bin/primer3_core

8.	(Optional) Make sure this location is within your $PATH (see below)


Where to put the binary files
=============================

A good place to put these is within ~/bin/ (this means in your home folder, within a folder named `bin` [for 'binary']).

You can also just drag the 'bin' folder to a location within your home directory.

You can certainly also copy the files within 'bin' to /usr/local/bin (if you are an administrator) or another similar location.  

You may need to adjust the permissions on the binaries if you get fancy.


Add the location to your $PATH
==============================

This is an optional step, but it will allow you to run primer3 in any directory on your machine as your user just by typing its name (primer3_core).

*** You should be very careful when altering your $PATH as things can go very wrong.  See below for an alternate method. ***  

If you added the binaries to /usr/local/bin, then you do not need to do this.

If you added the binaries to a local directory (let's say ~/bin/), do the following:

	1.  Edit your ~/.bash_profile.  You can edit this on the command line (Terminal) with:

		nano ~/.bash_profile
	
	2. Add the following line if it is not present (replacing '~/bin' if you used another directory):

		PATH=$PATH:~/bin/

	3. If a PATH line *is* present, ensure you add a colon to the end of what is there 
	and then the directory, so if you have something like:

		a) PATH=$PATH:/usr/local/genome/bin:/sw/bin
	
	make it look like:
	
		b) PATH=$PATH:/usr/local/genome/bin:/sw/bin:~/bin
	
	4. Quit and restart terminal for the changes to take effect.
		
If you don't add the location to your $PATH
===========================================

Assuming you don't want to modify your $PATH, you can still run the binaries.  Let's assume you put the files in '~/bin/.  You may run primer3_core by doing either of the following:

	1. ~/bin/primer3_core < yourInputFile
	2. /Users/<your username>/bin/primer3_core < yourInputFile
	
The first option is just a shortcut to the second.
