# ScanningTools

[![Build Status](https://travis-ci.com/fincardona/ScanningTools.py.svg?branch=master)](https://travis-ci.com/fincardona/ScanningTools.py)

A series of astronomical utilities developed to study scanning strategies of ground-based instruments. 

This software has been developed in the context of the LSPE/STRIP telescope that will be installed at Mount Teide, Tenerife. It provides functions to convert:

  - time period to mean sidereal or mean solar seconds
  - angles from sexagesimal to decimal and vice-versa
  - angles from degrees to hours and vice-versa
  - Greenwich Calendar Date (GCD) to Julian Date (JD)
  - Local Civil Time (LCT) to the GCD
  - LCT to JD
  - Horizontal coordinates (Alt, Az) to ICRS coordinates (Dec, Ra) 

and a lot of others tools to perform coordinates conversion in a fast and vey accurte way. A module named `Quaternions` allows to perform all the operations with the quaternions.
  
To install, clone the repository in your local folder with:

`git clone https://github.com/fincardona/ScanningTools.py.git`

then go inside the folder and run:

`python setup.py install`

Enjoy!
