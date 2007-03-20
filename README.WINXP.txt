Windows Installation Instructions (bcf; 3/20/2006:1200)
---------------------

primer3 release 1.1.0
---------------------

Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006
Whitehead Institute for Biomedical Research, Steve Rozen
(http://jura.wi.mit.edu/rozen), and Helen Skaletsky
All rights reserved.

How to install this software
============================

1. Unzip the '.zip' file downloaded from SourceForge.net
2. You will create a primer3-1.1.0-beta folder in the location where the file was unzipped
3. You may move this folder to a location of your choice

Running the software
====================

To run the program, you must do so from the MS-DOS command-line.  The intricacies of the DOS commandline are beyond the scope of this document.  Google for more information, if needed.  Here is a quick summary:

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
    