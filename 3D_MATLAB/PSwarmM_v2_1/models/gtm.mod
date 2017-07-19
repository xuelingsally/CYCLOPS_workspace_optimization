#  NLP written by GAMS Convert at 06/20/02 12:21:47
#  
#  Equation counts
#     Total       E       G       L       N       X
#        25       1      14      10       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        64      64       0       0       0       0       0       0
#  FX     4       4       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#       162     142      20       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 >= 0;
var x2 >= 0, <= 0.067;
var x3 >= 0, <= 0.067;
var x4 >= 0, <= 0.067;
var x5 >= 0, <= 0.067;
var x6 >= 0, <= 0.033;
var x7 >= 0;
var x8 >= 0;
var x9 >= 0, <= 0.3;
var x10 >= 0, <= 0.15;
var x11 >= 0, <= 0.1;
var x12 >= 0;
var x13 >= 0;
var x14 >= 0;
var x15 >= 0;
var x16 >= 0;
var x17 >= 0;
var x18 >= 0;
var x19 >= 0, <= 0.34;
var x20 >= 0, <= 0.35;
var x21 >= 0;
var x22 >= 0, <= 1.39;
var x23 >= 0, <= 1.06;
var x24 >= 0, <= 2;
var x25 >= 0, <= 2.62;
var x26 >= 0, <= 3.73;
var x27 >= 0, <= 0.62;
var x28 >= 0, <= 2.3;
var x29 >= 0, <= 1.03;
var x30 >= 0, <= 0.12;
var x31 >= 0, <= 1.45;
var x32 >= 0, <= 1.46;
var x33 >= 0, <= 0.48;
var x34 >= 0, <= 0.14;
var x35 >= 0;
var x36 >= 0, <= 0.1;
var x37 >= 0;
var x38 >= 0, <= 0.48;
var x39 >= 0, <= 0.8;
var x40 >= 0, <= 2.475;
var x41 >= 0, <= 3.7125;
var x42 >= 0, <= 0.297;
var x43 >= 0, <= 0.7128;
var x44 >= 0, <= 9.6525;
var x45 >= 0, <= 2.5245;
var x46 >= 0, <= 1.7028;
var x47 >= 0, <= 1.4256;
var x48 >= 0, <= 0.5148;
var x49 >= 0, <= 99;
var x50 := 2.2, >= 2.2, <= 2.2;
var x51 := 0.2, >= 0.2, <= 0.2;
var x52 := 1.47, >= 1.47, <= 1.47;
var x53 := 1.38, >= 1.38, <= 1.38;
var x54 := 0.2, >= 0.2;
var x55 := 0.2, >= 0.2;
var x56 := 0.2, >= 0.2;
var x57 := 0.2, >= 0.2;
var x58 := 0.2, >= 0.2;
var x59 := 0.2, >= 0.2;
var x60 := 0.2, >= 0.2;
var x61 := 0.2, >= 0.2;
var x62 := 0.2, >= 0.2;
var x63 := 0.2, >= 0.2;

minimize obj:  - (-4.84/x50 - 0.14/x51 - 6.4827/x52 - 6.6654/x53 - 
              8.89583741831423*x54^(-0.666666666666667) - 20.7788808225955*x55^
              (-0.515151515151515) - 12.8222379289592*x56^(-0.538461538461538)
               - 112.274462577384*x57^(-0.123595505617978) - 78.984522912416*
              x58^(-0.538461538461538) - 325.606233858943*x59^(-
              0.19047619047619) - 19.9925533406708*x60^(-0.492537313432836) - 
              20.2959676146409*x61^(-0.851851851851852) - 34.6492709112034*x62^
              (-1.32558139534884) - 2.07326743881507*x63^(-0.754385964912281)
               - (0.0372*x44 - 6.47537234042553*log(1 - 0.102564102564103*x44)
               - 0.489999999999999*log(1 - 1.38888888888889*x43) - 1.68*log(1
               - 0.392156862745098*x45) - 1.2271875*log(1 - 0.581395348837209*
              x46) - 0.2187*x46 - 0.979999999999999*log(1 - 0.694444444444444*
              x47) - 0.35*log(1 - 1.92307692307692*x48))) + 0.25*x1 + 2.29*x2
               + 2.22*x3 + 2.03*x4 + 1.96*x5 + 2.13*x6 + 0.4*x7 + 0.9*x8
               + 1.15*x9 + 1.1*x10 + 1.1*x11 + 0.8*x12 + 0.8*x13 + 0.65*x14
               + 0.7*x15 + 0.65*x16 + 1.5*x18 + 0.72*x19 + 0.46*x20 + 2.12*x21
               + 1.08*x22 + 1.01*x23 + 0.82*x24 + 0.75*x25 + 0.04*x26
               + 0.86*x27 + 0.14*x28 + 0.64*x29 + 0.77*x30 + 0.05*x31
               + 0.94*x32 + 0.53*x33 + 0.31*x34 + 0.58*x35 + 0.7*x36 + 1.91*x37
               + 0.43*x38 + 6*x39 + 2*x49;

subject to

e1:    x1 + x2 + x3 + x4 + x5 + x6 - x40 <= 0;

e2:    x7 + x8 + x9 + x10 + x11 + x12 + x13 + x14 + x15 + x16 - x41 <= 0;

e3:    x17 + x18 - x42 <= 0;

e4:    x19 + x20 - x43 <= 0;

e5:    x21 + x22 + x23 + x24 + x25 + x26 - x44 <= 0;

e6:    x27 + x28 + x29 - x45 <= 0;

e7:    x30 + x31 + x32 - x46 <= 0;

e8:    x33 + x34 + x35 + x36 + x37 - x47 <= 0;

e9:    x38 - x48 <= 0;

e10:    x39 - x49 <= 0;

e11:    x1 - x50 >= 0;

e12:    x17 - x51 >= 0;

e13:    x7 - x52 >= 0;

e14:    x8 - x53 >= 0;

e15:    x9 + x18 + x21 - x54 >= 0;

e16:    x2 + x10 + x22 - x55 >= 0;

e17:    x3 + x11 + x19 + x23 - x56 >= 0;

e18:    x4 + x24 - x57 >= 0;

e19:    x5 + x12 + x20 + x25 + x27 + x30 + x33 + x39 - x58 >= 0;

e20:    x26 + x28 + x31 - x59 >= 0;

e21:    x13 + x29 + x34 - x60 >= 0;

e22:    x14 + x35 - x61 >= 0;

e23:    x6 + x15 + x32 + x36 + x38 - x62 >= 0;

e24:    x16 + x37 - x63 >= 0;
