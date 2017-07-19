param a {1..10,1..5};
param b {1..10};
param c {1..5,1..5};
param d {1..5};
param e {1..5};

var x{1..5} >=0;

minimize obj:
   sum{i in 1..5} sum {j in 1..5} 
   c[i,j]*x[i]*x[j]
   +
   sum{j in 1..5} (e[j]*x[j]+d[j]*x[j]^3);

subject to constr1 {i in 1..10}:
   sum {j in 1..5}
   a[i,j]*x[j]-b[i]>=0;

data;
 param a :=
   1  1  -16
   1  2    2
   1  3    0
   1  4    1
   1  5    0
   2  1    0
   2  2   -2 
   2  3    0
   2  4    4
   2  5    2
   3  1   -3.5
   3  2    0
   3  3    2
   3  4    0
   3  5    0
   4  1    0
   4  2   -2
   4  3    0
   4  4   -4
   4  5   -1
   5  1    0
   5  2   -9
   5  3   -2
   5  4    1
   5  5   -2.8
   6  1    2
   6  2    0
   6  3   -4
   6  4    0
   6  5    0
   7  1   -1
   7  2   -1
   7  3   -1
   7  4   -1
   7  5   -1
   8  1   -1
   8  2   -2
   8  3   -3 
   8  4   -2
   8  5   -1
   9  1    1
   9  2    2
   9  3    3
   9  4    4
   9  5    5
  10  1    1
  10  2    1
  10  3    1
  10  4    1
  10  5    1
   ;

param c :=
   1  1   30
   1  2  -20
   1  3  -10
   1  4   32
   1  5  -10
   2  1  -20
   2  2   39
   2  3   -6
   2  4  -31
   2  5   32
   3  1  -10
   3  2   -6
   3  3   10
   3  4   -6
   3  5  -10
   4  1   32
   4  2  -31
   4  3   -6
   4  4   39
   4  5  -20
   5  1  -10
   5  2   32
   5  3  -10
   5  4  -20 
   5  5   30
   ; 

param e :=
   1  -15
   2  -27
   3  -36
   4  -18
   5  -12
   ;

param d := 
   1    4
   2    8
   3   10
   4    6
   5    2 
   ;

param b :=
   1  -40
   2   -2
   3    -.25
   4   -4
   5   -4
   6   -1
   7  -40
   8  -60
   9    5
  10    1
   ;

let x[1] := 0;
let x[2] := 0;
let x[3] := 0;
let x[4] := 0;
let x[5] := 1;

#let x[1] := 0.3;
#let x[2] := 0.33346761;
#let x[3] := 0.4;
#let x[4] := 0.42831010;
#let x[5] := 0.22396487;



#display obj+32.34867897;

