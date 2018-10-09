# ScanningTools.py

A series of astronomical utilities developed to study scanning strategies of ground-based instruments. 

This software has been developed in the context of the LSPE/STRIP telescope that will be installed at Mount Teide, Tenerife. It provides functions to convert:

  - time period to mean sidereal or mean solar seconds
  - angles from sexagesimal to decimal and vice-versa
  - angles from degrees to hours and vice-versa
  - Greenwich Calendar Date (GCD) to Julian Date (JD)
  - Local Civil Time (LCT) to the GCD
  - LCT to JD

and a lot of others tools exploiting `astropy` to perform coordinates conversion. A module named `Quaternions` allows to perform all the operations with the Quaternions.
  
  
To install, clone the repository in your local folder with:

`git clone https://github.com/fincardona/ScanningTools.py.git`

then go inside the folder and run:

`python setup.py install`

Enjoy!
