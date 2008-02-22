Windows Installation Instructions (bcf; 4/3/2006:0945)
---------------------

primer3 release 1.1.3
---------------------

Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008
Whitehead Institute for Biomedical Research, Steve Rozen
(http://jura.wi.mit.edu/rozen), and Helen Skaletsky
All rights reserved.

How to install this software
============================

1. Unzip the '.zip' file downloaded from SourceForge.net
2. You will create a primer3-1.1.3 folder in the location where the file was unzipped
3. You may copy the files from the 'bin' directory of the primer3-1.1.3 folder to a
location of your choice.  The exact same files are located within the 'src' folder so that
the tests may be run (windows does not allow relative paths in shortcuts).

Running the tests
=================
We are working on integrating the test suite to windows.  However, substantial differences
between windows and Unix/Linux require some differences in the test script.  

You must also install a perl distribution to run the windows tests.  

We *strongly* recommend you install ActiveState perl (http://www.activestate.com/products/activeperl/) 
as this was used to test our primer3 builds, and it is known to work.

***The perl test script for windows has a different name, at the moment, than that for the
unix/linux versions. See below.***

1. Click on 'Start > Run...'
2. Type 'cmd' into the space provided
3. Hit enter (or select 'OK')
4. Navigate to the location of the tests:
    
    A. if you put it in C:/Documents and Settings/YourName/primer3-1.1.3/test/,
    you would type 'cd c:/Documents and Settings/YourName/primer3-1.1.3/test/'
    
    B. you can also type 'cd ' (don't forget the space after cd) 
    and drag the primer3-1.1.3 folder onto the command-line window from 
    windows explorer, this will fill in the location for you

5. On the command line, run 'perl p3testz.pl -w' in this directory
6. You should see [OK] for all of the tests.

Running the software
====================

To run the program, you must do so from the MS-DOS command-line.  The intricacies of the 
DOS commandline are beyond the scope of this document.  Google for more information, if 
needed.  Here is a quick summary:

1. Click on 'Start > Run...'
2. Type 'cmd' into the space provided
3. Hit enter (or select 'OK')
4. Navigate to the location of the binary:
    
    A. if you put it in C:/Documents and Settings/YourName/Temp,
    you would type 'cd c:/Documents and Settings/YourName/Temp'
    
    B. you can also type 'cd ' (don't forget the space after cd) 
    and drag the primer3 folder onto the command-line window from 
    windows explorer, this will fill in the location for you
    
5. Run the example file by typing:

    primer3_core.exe < example
    
Other files may be run in a similar fashion.  If your input filename
is 'MyData.txt' you can run primer3 using this file (in the correct 
format; see README) with:

    primer3_core.exe < MyData.txt
    
If your file is not in the folder containing primer3_core.exe, 
you could run the program from the primer3_core folder using:
    
    primer3_core.exe < c:/someOtherFolder/someOtherFolder/MyData.txt
    
Finally, if you want to run the program without going to its folder, 
assuming primer3_core.exe is in c:/Temp, you could run:

    c:/Temp/primer3_core.exe < c:/someOtherFolder/someOtherFolder/MyData.txt
    