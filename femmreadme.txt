Transformer Finite Element Analysis using FEMM with Python
version 1.0, Jun-11-2022. Contact for any questions at pravirajpg@gmail.com

1. Install these before you start anything:

2. Install FEMM (ver=4.2), from https://www.femm.info/wiki/Download

3. Install Python 3.7 or later version

4. Go to Python3.x/Scripts/ folder and then open "cmd" and do
    pip install pyfemm
    pip install openpyxl
    pip install scipy
    pip install matplotlib
    pip install numpy

5. It is the assumption that you are familiar with Andersen FLD12 program. I am trying to mimic
    what FLD12 program does. So an extensive documentation on how to input the Excel file is
    difficult to prepare (basically I am lazy). I could not implement the short circuit force
    calculations yet, may need somemore time.

6. "femmapp.py" is the main program, "femmppg.py" hosts all the necessary function definitions, and
    "femmppg.xlsx" where all the input data sheets stored.

7. To run the program, open a "cmd.exe" (windows standard command prompt), cd to the path where all
    files are stored, then enter following command and then enter the excel sheet name (case insensitive):

    python femmapp.py

8. Connection codes (Row # 6)
    1 = Wye, 2 = Delta,
    3 = Auto connected (HV terminal), 4 = Auto connected (LV terminal),
    5 = ZigZag Neutral Connected winding, 6 = ZigZag Line Connected Winding
    7 = Polygon Delta Large Winding, 8 = Polygon Delta Small Winding
    9 = Extended Delta (inside Delta) winding, 10 = Extended Delta Terminal connected winding

9. Note 1: When defining connections (except for Wye or Delta connection), terminal numbers shall
    always be in increasing order. E.g, when defining Auto transformer terminals, define the HV terminal
    first (code = 3) and then LV terminal (code = 4). or say when defining a ZigZag connection, define
    the Neutral connected winding (code = 5) and then line connected winding (code = 6). For all these
    connections (except Wye, Delta, Auto), rated voltage is given as volts per turn times the number of turns.
    In all cases, MVA is given as total for the two windings for both terminals, the program will take care of
    the MVA based on their phase shifts.

10. Note 2: This code currently doesn't support sheet windings (especially the eddy losses will be incorrect
    if you tried to run it with sheet windings)

11. A file named "report.txt" will be saved in the same folder after running. This file is formatted to be
    viewed with proper indentations using Notepad++.

12. If Ampere Turns deviation is greater than 0.5%, it is most likely that the input is incorrect.

13. A plot will be created using matplotlib and saved in the same folder as "femmplot.png"

Good luck,
Pravi
