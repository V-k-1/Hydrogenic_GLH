(* ::Package:: *)

BeginPackage["dipolemoments`"];
g0::usage = "g0[n] calculates density normalized dipole transition moments for n,n-1->0,n.";
gup::usage = "gup[i][n,\[Kappa]] calculates density normalized dipole transition moments for n,n-i->\[Kappa],n-i+1.";
gdown::usage = "gdown[i][n,\[Kappa]] calculates density normalized dipole transition moments for n,n-i->\[Kappa],n-i-1.";
Rho::usage = "Rho[\[Kappa]] gives the density of states at \[Kappa]. g**[*][*,\[Kappa]]/Rho[\[Kappa]] gives dipole transitions.";

Begin["`Private`"];
Z=1;
a0 = 1;
Eh=1;
Ry =Eh/2;
(* Input validation *)
validateInput[n_, \[Kappa]_] := (
    If[!IntegerQ[n] || n < 1, Message[g0::intnm, n]; Return[$Failed]];
    If[!NumericQ[\[Kappa]] || \[Kappa] < 0, Message[gup::numk, \[Kappa]]; Return[$Failed]]
);
(* Remember : eqn. 24 already has Sqrt(1/2\[Kappa]) thats why we only need 1/Sqrt[Ry] to invert the density normalization gd0*)
(* Base case *)
g0[n_] := (
    validateInput[n, 0];
    n^2/Z^2 *Sqrt[1/(2*(2n-1)!)]*4*(4n)^n*Exp[-2n]*a0/Sqrt[Ry]
);

(* Transition cases *)
gup[1][n_, \[Kappa]_] := (
    validateInput[n, \[Kappa]];
    Module[{product},
        product = Product[(1 + s^2 \[Kappa]^2), {s, 1, n}];
        If[\[Kappa]==0,1,Sqrt[product/(1 - Exp[-2\[Pi]/\[Kappa]])] * Exp[2n - 2/\[Kappa] ArcTan[n \[Kappa]]] * (1 + n^2 \[Kappa]^2)^(-n-2)] * g0[n]
    ]
);
gup[2][n_, \[Kappa]_] := (
    validateInput[n, \[Kappa]];
    (1/2)*Sqrt[(2n-1)*(1+n^2 \[Kappa]^2)]*gup[1][n, \[Kappa]]
);

gdown[1][n_, \[Kappa]_] := (
    validateInput[n, \[Kappa]];
    1/(2n)Sqrt[(1+n^2\[Kappa]^2)/(1 + (n-1)^2\[Kappa]^2)] * gup[1][n,\[Kappa]]
);
gdown[2][n_, \[Kappa]_] := (
    validateInput[n, \[Kappa]];
    ((4+(n-1)(1+n^2\[Kappa]^2))/(2n)) * Sqrt[(2n-1)/(1 + (n-2)^2 \[Kappa]^2)] * gdown[1][n, \[Kappa]]
    );
    
(* Recursive case *)
gup[i_][n_, \[Kappa]_] := (
    validateInput[n, \[Kappa]];
    If[i < 3, Return[gup[i][n, \[Kappa]]]];
    Module[{l, term1, term2, coeff},
        l = n - i + 2;
        term1 = (4n^2 - 4l^2 + l(2l-1)(1+n^2 \[Kappa]^2))*gup[i-1][n, \[Kappa]];
        term2 = 2n*Sqrt[(n^2-l^2)*(1+(l+1)^2 \[Kappa]^2)]*gup[i-2][n, \[Kappa]];
        coeff = 2n*Sqrt[(n^2-(l-1)^2)*(1+l^2 \[Kappa]^2)];
        (term1 - term2)/coeff
    ]
);
  

gdown[i_][n_,  \[Kappa]_] /;i>2 := Module[{l,term1, term2, coeff},
\[NonBreakingSpace]\[NonBreakingSpace]l=n-i+1;
    term1 = (4n^2 - 4l^2 + l(2l+1)(1 + n^2 \[Kappa]^2)) * gdown[i-1][n, \[Kappa]];
    term2 = 2n * Sqrt[(n^2 - (l+1)^2)*(1 + (l)^2 \[Kappa]^2)] * gdown[i-2][n, \[Kappa]];
    coeff = 2n*Sqrt[(n^2 - l^2)*(1 + (l-1)^2 \[Kappa]^2)]; 
    (* Full solution *)
    (term1 - term2)/coeff
];
Rho[\[Kappa]_] := 1/(2*Ry*\[Kappa])
End[];
EndPackage[];
