# This is main starting program, run this program.
# Transformer Leakage Flux Study using FEMM with Python
# Thanks to Dr. David Meeker, http://www.femm.info/wiki/DavidMeeker
# Please read FEMMREADME before using this code
# Code: Praviraj PG (pravirajpg@gmail.com), Rev 0, Date: Jun-11-2022

from femm import *
from femmppg import *


def femmapp():
    sheet_name = input('Enter sheet name: ').upper()  # Enter sheet name
    xl_fp = r'femmppg.xlsx'  # Excel file name
    start = time.time()  # check some running time
    d1, d2, d3, d4 = femmdata(sheet_name, xl_fp)  # Read excel data

    units = d1[0];
    xx = 1  # Units, 1=inches, 2=mm
    if units == 2: xx = 25.4  # Conversion, inch to mm

    MVA = d1[1]  # Base MVA
    f = d1[2]  # Frequency, Hz
    Di = d1[3] / xx  # Outer diameter of Core
    Z = d1[4] / xx  # Core Window Height
    numt = d1[5]  # Number of Terminals
    numl = d1[6]  # Number of Layers
    dt = d1[7] / xx  # Distance from outer winding to tank

    max_id = np.argmax(d3[1])  # Index of the layer with largest inner diameter
    # R = Radial Outer Boundary (basically the distance to tank from center of the core)
    R = dt + d3[1][max_id] / 2 / xx + d3[2][max_id] / xx
    Ri = Di / 2  # Radial Inner Boundary (Core outer radius)

    # Variables & Constants...
    unit = 'inches'
    Cond = 0  # Conductivity of Cu... MS/m (set 0 for AC, neglect skin effect)
    muo = pi * 4e-7  # Permeability of Air/Oil...
    mur_fe = 4e4  # Relative Permeability of Core..
    mur_ms = 100  # Relative Permeability of Mild Steel..
    cond_ms = 7  # Conductivity of Mild Steel, MS/m..
    c0i = 39.37 / (muo * mur_fe * Ri)  # Inner (Core) Boundary Conditions..
    c0o = 39.37 / (muo * mur_ms * R)  # Outer (Tank) Boundary Conditions..

    # Initialize...
    openfemm(1)  # Open the FEMM program, param=1 for hiding, enter 0 otherwise
    main_resize(600, 700)  # Resize FEMM window (WIDTH X HEIGHT)..
    newdocument(0)  # Create a new document...
    mi_probdef(f, unit, 'axi', 1e-8, 1, 30)  # Define Problem.. (f=0, Magnetostatic, f>0, AC)

    # Draw a rectangle to use as the Outer boundary and assign boundary condition..
    mi_drawrectangle(Ri, 0, R, Z)  # Boundary...
    mi_addboundprop('Core', 0, 0, 0, 0, 0, 0, c0i, 0, 2, 0, 0)  # Core Leg Boundary (inside)..
    mi_addboundprop('Tank', 0, 0, 0, 0, 0, 0, c0o, 0, 2, 0, 0)  # Tank/Air Boundary (outside)..
    mi_selectsegment(R, Z / 2);    mi_setsegmentprop('Tank', 1, 1, 0, 0)
    mi_clearselected();            mi_selectsegment(Ri, Z / 2)
    mi_selectsegment(R / 2, 0);    mi_selectsegment(R / 2, Z)
    mi_setsegmentprop('Core', 1, 1, 0, 0);    mi_clearselected()

    # Add medium, say Oil and Label it..
    mi_addmaterial('Oil', 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0)  # The medium
    mi_addblocklabel(R - 1, Z - 1)  # Label it as Oil/Air..
    mi_selectlabel(R - 1, Z - 1)
    mi_setblockprop('Oil', 0, 1, '<None>', 0, 0, 0)
    mi_clearselected();    mi_zoomnatural()

    # Add Winding Properties..
    mi_addmaterial('Coil', 1, 1, 0, 0, Cond, 0, 0, 1, 0, 0, 0)

    # Add each segments..
    d2 = np.append(np.array(d2), np.zeros((1, numt)), axis=0)  # Append Terminal data to add Angles..
    max_seg = np.max(d3[3])  # Maximum number of segments
    AT = np.zeros((max_seg, 4))  # [Term #, Lay #, RMS NI, Phase Angle]
    Js = []  # Currents in each segment..
    Ly = np.zeros((max_seg, 5))  # Segment details..
    sn = 1;     a_r = 1  # Auto Ratio for Auto-transformers
    for ly in range(numl):
        tr = d3[0][ly]  # Terminal Number..
        ir = d3[1][ly] / 2 / xx  # Inside Radius..
        rd = d3[2][ly] / xx  # Radial..
        ls = d3[3][ly]  # Last Segment..
        npg = d3[4][ly]  # No. of Parallel Groups..
        cd = d3[5][ly]  # Current Direction..

        con = d2[0][tr - 1]  # Connection Type..
        mva = d2[1][tr - 1]  # 3 - Phase MVA of Terminal..
        kV = d2[2][tr - 1]  # Voltage of Terminal(Line - Line in kV)..
        Ip = 0;
        ph = 0  # Phase Current & Angles

        if con == 1 or con == 3:  # 1 - Wye or 3 - Auto (High Voltage Terminal)
            Ip = mva * 1000 / sqrt(3) / kV  # RMS Coil Current..

        elif con == 2:  # 2 - Delta
            Ip = mva * 1000 / 3 / kV  # RMS Coil Current..

        elif con == 4:  # 4 - Auto (Low Voltage Terminal)
            Ip = mva * 1000 / sqrt(3) / kV  # RMS Coil Current..
            Ip = Ip - (d2[1][tr - 2] * 1000 / sqrt(3) / d2[2][tr - 2])
            a_r = 1 - kV/d2[2][tr - 2]       # Auto-ratio

        elif con == 5 or con == 7:  # Zig Neutral (5) or Polygon Delta Main (7)
            cvzi = kV  # Coil Volts..
            cvza = d2[2][tr]
            cvz = cvza + cvzi * exp(-1j * pi / 3)  # 120° or -60° b/w windings
            trm, tra = polar(cvz)  # Get magnitude and angle
            ph = -pi / 3 + abs(tra)  # ZIG Coil at angle w.r.t. 0° Winding
            Ip = mva * 1000 / sqrt(3) / trm  # RMS Coil Current..
            if con == 7:
                Ip = Ip / sqrt(3)

        elif con == 6 or con == 8:  # Zag Line (6) or Polygon Delta Terminal (8)
            cvza = kV  # Coil Volts..
            cvzi = d2[2][tr - 2]
            cvz = cvza + cvzi * exp(-1j * pi / 3)  # 120° or -60° b/w windings
            trm, tra = polar(cvz)  # Get magnitude and angle
            ph = abs(tra)  # ZAG Coil at angle
            Ip = mva * 1000 / sqrt(3) / trm  # RMS Coil Current..
            if con == 8:
                Ip = Ip / sqrt(3)

        elif con == 9:  # 9 - Extended Delta Main winding
            cvem = kV  # Coil Volts..
            cvet = d2[2][tr]
            cve = (cvem + cvet) + cvet * exp(-1j * pi / 3)  # 120° or -60° b/w windings
            trm, tra = polar(cve)  # Get magnitude and angle
            ph = tra  # Main Coil Angle
            Ip = mva * 1000 / 3 / trm  # RMS Coil Current..

        elif con == 10:  # 10 - Extended Delta Terminal
            cvet = kV  # Coil Volts..
            cvem = d2[2][tr - 2]
            cve = (cvem + cvet) + cvet * exp(-1j * pi / 3)  # 120° or -60° b/w windings
            trm, tra = polar(cve)  # Get magnitude and angle
            ph = pi / 6 + tra  # Extended Coil angle
            Ip = mva * 1000 / sqrt(3) / trm  # RMS Coil Current..

        d2[3][tr - 1] = ph * 180 / pi  # For reporting..
        Ip = Ip * exp(1j * ph)  # Convert to AC format..
        Ip = Ip / npg * sqrt(2)  # Peak Coil Current per Group..

        while sn <= ls:

            # Draw Segments & Label it..
            x1 = ir;     x2 = ir + rd
            y1 = d4[0][sn - 1] / xx;     y2 = d4[1][sn - 1] / xx
            mi_drawrectangle(x1, y1, x2, y2)  # Draw rectangle..
            mi_addblocklabel(x1 + rd / 2, (y1 + y2) / 2)  # Label..
            seg_name = 'S' + str(sn)

            # Assign of Active Ampere Turns..
            at = d4[3][sn - 1]
            if at == 0:  # If Active Turns are zero, make current zero / some how FEMM takes it
                mi_addcircprop(seg_name, 0, 1)  # as 1 Turn otherwise..
                J = 0; Js.append(J)
            else:

                mi_addcircprop(seg_name, Ip * at * cd, 1)  # Set Ampere Turns(AT)
                J = cd * abs(Ip); Js.append(J)  # Peak current in each segment..

            mi_selectlabel(x1 + rd / 2, (y1 + y2) / 2)
            mi_setblockprop('Coil', 0, 1, seg_name, 0, 1, 1)  # Set Turns to 1, total AT is spec'd in Line 159
            mi_clearselected()
            AT[sn - 1] = [tr, ly + 1, (J*at/sqrt(2)), d2[3][tr - 1]]  # Segment data..
            Ly[sn - 1] = [x1, y1, x2, y2, J]  # Layer data for plotting..
            sn = sn + 1

    # Name the geometry, analyze, & load solution..
    mi_saveas('femmxfr.fem')
    mi_analyze()
    # mi_loadsolution(); mo_hidepoints()                 # Use this to see output in FEMM App
    # mo_showvectorplot(5, 1)                            # Use this to see output in FEMM App
    # mo_showcontourplot(0,0.003,-0.01,'real')

    # Call FLUXBRZ to calculate flux densities (B) from vector potential (A)
    E, Rr, Zz, Fr, Fz, Fa = femmflux()

    # Leakage Reactance Calculation from total energy
    X = (1 / 2) * 2 * E * (2 * pi * f) / (MVA * 1e6 / 3) * 100

    # Estimate eddy loss
    Pe, I2R = femmeddy(d1, d3, d4, Js, Fr, Fz)

    # Create report and save in same folder under name "report.txt"
    # Open the report.txt with Notepad++ for view in proper formatting
    femmreport(sheet_name, d1, d2, d3, d4, AT, X, E, I2R, Pe, a_r)

    # Plotting function
    plt = femmplot(sheet_name, Fr, Fz, Fa, Ri, R, Z, Ly, xx)

    # Print total time taken
    end = time.time()
    print(f'Time taken: {(end - start):.1f} sec')

    # Show the plot
    plt.show()

    # Hold FEMM until you press Enter
    # prompt('Press Enter')
    closefemm()


# Start
if __name__ == "__main__":
    femmapp()
