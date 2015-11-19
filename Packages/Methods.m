(* ::Package:: *)

BeginPackage["Methods`"]
Needs["Developer`"]

\[Epsilon]1=10.;
\[Epsilon]2=5.;
H=12.;(*Height of first layer in angstroms*)
a= 3.904;(*unit cell size in angstroms*)
t1=.25;(*hopping parameters in eV*)
t2=.025;

MM=8;
NN=4;
\[Beta]=.95;(*percent of the old distribution to keep when calculating the new distribution*)
numPartitions=20;
numStates=20;

Ry=13.60569253;(*Rydberg unit of energy, in eV. 8\[Pi]Ry=e^2/(Subscript[a, 0]Subscript[\[Epsilon], 0]) in eV*)

vals={};
vecs={};
degeneracies={};
baseMatrix={};
dk=N[Pi]/numPartitions;
initShift=0.;
mm=Floor[(MM+1)/2];
EFermi= 0.;
errorTrackingList={};
SetSharedVariable[vals,vecs,degeneracies];

\[Alpha]=.529/a;(*Subscript[a, 0]/a, where Subscript[a, 0] is the Bohr radius*)
densityToNumStatesScaling=(NN(MM-1)dk)^-1;
lambdaToNumStatesScaling = ((\[Epsilon]2/\[Epsilon]1)(1/2 ArcTan[(MM a)/(2H)]-1/2 ArcTan[-((MM-2) a)/(2H)])NN (MM-1) 1/numPartitions)^-1;
\[Lambda] = lambdaToNumStatesScaling numStates;
totalCharge= -(\[Epsilon]2/\[Epsilon]1)(\[Lambda]/N[\[Pi]])(1/2 ArcTan[(MM a)/(2H)]-1/2 ArcTan[-(((MM-2) a)/(2H))]);

chargeDist[x_,y_]:=Sin[(\[Pi] y)/MM]Exp[-(2\[Pi] x)/MM];



TotalEnergySortedList::usage="TotalEnergySortedList[] returns {{vals}.{vecs},{degeneracies}} for every site over every available k value, sorted by order of increasing energy";
PoissonSolver::usage="PoissonSolver[chargeList_List] returns the solution to Poisson's equation over a NNxMM grid when given the charge distribution at each point as a list.";
PoissonPlot::usage="PoissonPlot[chargeList_List] plots the solution to Poisson's equation.";
eigensystem::usage="eigensystem[kz_] returns the eigensystem of the matrix operator for in the tight binding model in the form {{vals},{vecs},{degeneracies}}. This gives the eigenvalues of each atom for a particular value of crystal momentum, followed by the corresponding eigenvectors, followed by the corresponding degeneracies.";
NecessaryStates::usage="NecessaryStates[numStatesTotal_] find the lowest energies with associated eigenvalues and eigenvectors that are necessary to fill the desired number of states, accounting for fractional filling. This is equivalent to a zero temperature Fermi calculation.";
totalDensity::usage="totalDensity[n_] returns the density of the \!\(\*SuperscriptBox[\(n\), \(th\)]\) lowest eigenvalue as a list.";
orbitalDensity::usage="orbitalDensity[n_,\[Alpha]_] returns the orbitally projected density of the \!\(\*SuperscriptBox[\(n\), \(th\)]\) lowest eigenvalue as a list. \[Alpha]=1 for yz, 2 for zx, 3 for xy.";
totalOrbitallyProjectedDensity::usage="totalOrbitallyProjectedDensity[n_,\[Alpha]_] returns the density of the \!\(\*SuperscriptBox[\(n\), \(th\)]\) lowest eigenvalue due to a single orbital as a list. \[Alpha]=1 for yz, 2 for zx, 3 for xy.";
totalDensityPlot::usage="totalDensityPlot[n_] plots the density of the \!\(\*SuperscriptBox[\(n\), \(th\)]\) lowest eigenvalue over the full space.";
normalizedSumList::usage="normalizedSumList[numStatesTotal_] returns the charge distribution due to the lowest energies that fill the first numStatesTotal states, over the full space.";
normalizedProjectedSumList::usage="normalizedSumList[numStatesTotal_,\[Alpha]_] returns the orbitally projected charge distribution due to the lowest energies that fill the first numStatesTotal states, over the full space. \[Alpha]=1 for yz, 2 for zx, 3 for xy.";
normalizedSumPlot::usage="normalizedSumPlot[numStatesTotal_] plots the charge distribution due to the lowest energies that fill the first numStatesTotal states, over the full space.";
BandStructureList::usage="BandStructureList[listvals_] plots the band structure resulting from a charge distribution of listvals.";
BandStructurePlot::usage="BandStructurePlot[minE_,maxE_,list_List] plots the band structure given a list of {k,\!\(\*SubscriptBox[\(E\), \(n\)]\)(k)} values, with an energy range from minE to maxE.";
SelfConsistentDistribution::usage="SelfConsistentDistribution[firstGuess_,errorTolerance_] returns a self consistent charge distribution starting with firstGuess as the initial charge distribution and exits when successive iterations have a average sum squared differnece of less than errorTolerance.";
PlotFullSpaceDistribution::usage="PlotFullSpaceDistribution[listvals_List] plots a graph of listvals over a NNx(MM-1) grid, missing the periodic edge.";
Initialize::usage="Initialize[] initializes variables that are necessary for operation.";
SetParameters::usage="SetParameters[new\[Epsilon]1_,new\[Epsilon]2_,newH_, newa_,  newt1_,newt2_, newMM_, newNN_, new\[Beta]_, newnumPartitions_, newnumStates_ ] set the new parameters based on the new inputsinput";


Begin["Private`"]

Initialize[]:=Module[{},
listvals=Flatten[Table[
 N[chargeDist[i ,j ]],
{i,1,NN},
{j,1,MM-1}
],1];
baseMatrix=Developer`ToPackedArray[BaseMatrixTakesChargeList[listvals]];
{vals,vecs,degeneracies}=eigensystem[0];
DistributeDefinitions["Methods`"];
];

SetParameters[new\[Epsilon]1_,new\[Epsilon]2_,newH_, newa_,  newt1_,newt2_, newMM_, newNN_, new\[Beta]_, newnumPartitions_, newnumStates_ ]:=Module[{},
\[Epsilon]1=new\[Epsilon]1;
\[Epsilon]2=new\[Epsilon]2;
H=newH;
a=newa;
t1=newt1;
t2=newt2;
MM=newMM;
NN=newNN;
\[Beta]=new\[Beta];
numPartitions=newnumPartitions;
numStates=newnumStates;
densityToNumStatesScaling=(NN(MM-1)dk)^-1;
lambdaToNumStatesScaling = ((\[Epsilon]2/\[Epsilon]1)(1/2 ArcTan[(MM a)/(2H)]-1/2 ArcTan[-((MM-2) a)/(2H)])NN (MM-1) 1/numPartitions)^-1;
\[Alpha]=.529/a;(*Subscript[a, 0]/a, where Subscript[a, 0] is the Bohr radius*)
\[Lambda] = lambdaToNumStatesScaling numStates;
totalCharge= -(\[Epsilon]2/\[Epsilon]1)(\[Lambda]/N[\[Pi]])(1/2 ArcTan[(MM a)/(2H)]-1/2 ArcTan[-(((MM-2) a)/(2H))]);
mm=Floor[(MM+1)/2];
]

(*Poisson solver that implements the first part of the algorithm described below, solving over full space, including periodic edge - it represents the Laplacian as a matrix operator then inverts it. Any solution to Poisson's equation is now found by multiplying the inverse matrix by the (forcing function \[Rho]? idk if you can call it that)*)
(*This is a symmetric finite difference method to solve the Poisson Equation*) 
(*en.wikipedia.org/wiki/Discrete_Poisson_equation*)
PoissonSolver[chargeList_List ]:=Module[
{m,m1,m2,m3,m4,m5,m6,m7,m8,m9,newm,\[Gamma],potential,exteriorPoints,chargeVector},

(*construct D, the matrix operator that corresponds to the Laplacian. This is a 9 point stencil. It's periodic with respect to two edges: Born-von Karman boundaries. One edge is fixed to 0 everywhere: Dirichlet boundaries. The last edge has the normal component of the flux at the boundary: Neumann boundaries.*)
periodicblock[c_,b_]:=c IdentityMatrix[MM-1]+DiagonalMatrix[Table[b,{i,0,MM-3}],1]+DiagonalMatrix[Table[b,{i,0,MM-3}],-1]+DiagonalMatrix[{b},-MM+2]+DiagonalMatrix[{b},MM-2];
m1 = periodicblock[-2./3,-1/6];
m2 = periodicblock[10/3,-2/3];
m3 =KroneckerProduct[DiagonalMatrix[Join[Table[1,{i,1,NN}],{0}],0],m2];
m4 =KroneckerProduct[DiagonalMatrix[Join[Table[1,{i,1,NN-1}],{0}],1],m1];
m5 =KroneckerProduct[DiagonalMatrix[Join[Table[1,{i,1,NN-1}],{0}],-1],m1];
m6 =KroneckerProduct[DiagonalMatrix[{1},NN],m1];
m7=1/2 IdentityMatrix[MM-1];
(*The matrix has been augmented so that the bottom blocks correspond to a row of "ghost" points above the top row. THey are used to establish a centered second order approximation for the first derivative on V normal to the boundary. This *)
m8 =KroneckerProduct[DiagonalMatrix[{0,1},-NN+1],m7];
m9=KroneckerProduct[DiagonalMatrix[Join[Table[0,{i,1,NN}],{1}],0],-m7];
m = m3+m4+m5+m6+m8+m9;
newm=m3+m4+m5;

(*coefficient in front of n on the right hand side of the poisson equation*)
\[Gamma]=-((8N[\[Pi]]Ry \[Alpha])/\[Epsilon]2);
e1[j_]:=(8N[\[Pi]]Ry \[Alpha] \[Lambda] H a)/(2 N[\[Pi]]\[Epsilon]1 (H^2+(j-Floor[(MM+1)/2])^2 a^2));(*Normal component of electric field at top edge of the grid*)

(*Defines the list with the charge distribution Subscript[\[Rho], ij] generated from the user, augmented with the boundary conditions *)
chargeVector  = \[Gamma] chargeList;
exteriorPoints = 
Table[
  N[e1[j ]],
{j,1,MM-1}
]//Flatten;
chargeVector = Join[chargeVector, exteriorPoints];
(*correction to get 4th order error*)
chargeVector=chargeVector-1/2 newm.chargeVector;
(*calculates the potential and formats it as a matrix, while dropping the row of ghost points*)
potential=Drop[Partition[LinearSolve[m,chargeVector],MM-1],-1];

(*inserts periodic edge*)
potential=Insert[potential//Transpose,potential[[All,1]],-1]//Transpose//Flatten;

(*Reindexes the potential and outputs a list of voltage formatted as {x,y,V(x,y)}*)
Flatten[Table[
{i,j,potential[[j + (MM)(i-1)]]},
{i,1,NN},
{j,1,MM}
],1]
]
PoissonPlot[chargeList_List]:=Module[{plotList=PoissonSolver[chargeList ]},
ListPlot3D[plotList,PlotRange->All,PlotLabel->"Calculated Electric Potential",PlotRange->All]
]


(*gets rid of periodic edge, formats list for hamiltonian*)
TripleVoltage[chargeList_List ]:=Module[{list=Developer`ToPackedArray[PoissonSolver[chargeList][[All,3]]]},
list=Drop[list,{MM,-1,MM}];
Partition[Join[list,list,list],list//Length]//Transpose//Flatten
]

BaseMatrixTakesChargeList[chargeList_List ]:=Module[{list=TripleVoltage[chargeList],voltageMatrix,tx,ty,m1,m2,m3,m4},
voltageMatrix=DiagonalMatrix[list];
tx = DiagonalMatrix[{t2,t1,t1},0];
ty = DiagonalMatrix[{t1,t2,t1},0];

(*These 4 matrices are interatomic bonding in each of the 4 directions in the x-y plane*)
m1 = KroneckerProduct[DiagonalMatrix[Table[1,{NN(MM-1)-1}],1],ty];
m2 = KroneckerProduct[DiagonalMatrix[Table[1,{NN(MM-1)-1}],-1],ty];
m3 = KroneckerProduct[DiagonalMatrix[Table[1,{NN(MM-1)-(MM-1)}],(MM-1)],tx];
m4 = KroneckerProduct[DiagonalMatrix[Table[1,{NN(MM-1)-(MM-1)}],-(MM-1)],tx];
m1+m2+m3+m4+voltageMatrix
]

(*Returns eigensystem for matrix for a single k state*)
eigensystem[kz_]:=Module[
{Eyz,Ezx,Exy, hamiltonian,m1,vals,vecs,degeneracyList,matrix},
(*Eyz = 2t1 -2 t1 Cos[kz]+2t2+2t1;
Ezx = 2t1 -2 t1 Cos[kz]+2t1+2t2;
Exy = 2t2 -2 t2 Cos[kz]+2t1+2t1;*)
Eyz = 2t1 -2 t1 Cos[kz];
Ezx = 2t1 -2 t1 Cos[kz];
Exy = 2t2 -2 t2 Cos[kz];
hamiltonian = DiagonalMatrix[{Eyz,Ezx,Exy},0];
KroneckerProduct[IdentityMatrix[NN*(MM-1)],hamiltonian];
m1 = KroneckerProduct[IdentityMatrix[NN*(MM-1)],hamiltonian];
matrix=baseMatrix+m1;
{vals,vecs}=Eigensystem[matrix];
degeneracyList=Developer`ToPackedArray[spinKDegeneracyFactor[kz] Table[1,{i,1,vals//Length}]];
{vals,vecs,degeneracyList}
];

(*returns sorted list of energies for ALL k states*)
TotalEnergySortedList[]:=Module[{vecs,vals,degeneracies,sys,orderList,temp},

sys=ParallelTable[eigensystem[i dk],{i,0,numPartitions}];
vals=Developer`ToPackedArray[Flatten[sys[[All,1]]]];
vecs=Flatten[sys[[All,2]],1];
degeneracies=Developer`ToPackedArray[Flatten[sys[[All,3]]]];

orderList=Ordering[vals];
vals=vals[[orderList]];
vecs=vecs[[orderList]];
degeneracies=degeneracies[[orderList]];

{vals,vecs,degeneracies}
]

(*returns portion of the energy list to use and computers the fractional filling*)
NecessaryStates[numStatesTotal_]:=Module[{list=TotalEnergySortedList[],sum=0,count=0,vals,vecs,degeneracies},
degeneracies=list[[3]];
(*avoids computing all partial sums, for increased efficiency with long lists*)
Do[
sum+=degeneracies[[i]];
If[sum>numStatesTotal,
count=i;
Break[];
];
,{i,1,degeneracies//Length}
];

vals=list[[1]][[1;;count]];
vecs=list[[2]][[1;;count]];
degeneracies=degeneracies[[1;;count-1]];

degeneracies=Join[degeneracies,{numStatesTotal-(sum-list[[3,count]])}];

If[sum-numStatesTotal== list[[3,count]],vals=Drop[vals,-1];vecs=Drop[vecs,-1];degeneracies=Drop[degeneracies,-1]];
EFermi=vals[[-1]];
{vals,vecs,degeneracies}

];

(*\[Alpha]=1 for yz, 2 for zx, 3 for xy*)
orbitalDensity[n_,\[Alpha]_]:=Module[{list=vecs[[n]][[\[Alpha];;;;3]]},
Flatten[Table[
{i ,j ,Norm[list[[j + (MM-1)(i-1)]]]^2},
{i,1,NN},
{j,1,MM-1}
],1]
];

totalOrbitallyProjectedDensity[n_,\[Alpha]_]:=orbitalDensity[n,\[Alpha]][[All,3]];
totalDensity[n_]:=totalOrbitallyProjectedDensity[n,1]+totalOrbitallyProjectedDensity[n,2]+totalOrbitallyProjectedDensity[n,3];
(*\[Alpha]=1 for yz, 2 for zx, 3 for xy*)
totalDensityPlot[n_]:=
PlotFullSpaceDistribution[totalDensity[n]];


spinKDegeneracyFactor[kz_]:=2 If[kz==0 ||kz==N[Pi],1,2];

(*returns CHARGE distribution, not density of states*)
normalizedSumList[numStatesTotal_]:=Module[{totalList},
{vals,vecs,degeneracies}=NecessaryStates[numStatesTotal];
densityToNumStatesScaling ParallelSum[degeneracies[[n]]totalDensity[n],{n,1,Length[vals]}]
];

(*returns orbitally projected charge distribution*)
normalizedProjectedSumList[numStatesTotal_,\[Alpha]_]:=Module[{totalList},
{vals,vecs,degeneracies}=NecessaryStates[numStatesTotal];
densityToNumStatesScaling ParallelSum[degeneracies[[n]]totalOrbitallyProjectedDensity[n,\[Alpha]],{n,1,Length[vals]}]
];

normalizedSumPlot[numStatesTotal_]:=PlotFullSpaceDistribution[normalizedSumList[numStatesTotal]];

BandStructureList[chargeList_List]:=Module[{table,length},
baseMatrix=BaseMatrixTakesChargeList[chargeList];
table=ParallelTable[eigensystem[i dk][[1]],{i,0,numPartitions}];
table=table-table[[1,1]];
length=table[[1]]//Length;
Table[
{i dk, table[[i+1]][[j]]},
{j,1,length},
{i,0,numPartitions}
]
]

BandStructurePlot[minE_,maxE_,list_List]:=Module[{FermiLine,f1,f2},
FermiLine = Table[{i dk,EFermi},{i,0,numPartitions}];
f1=ListPlot[list,Joined-> True,Mesh-> All,PlotRange->{minE,maxE}];
f2 = ListPlot[FermiLine,Joined-> True,PlotStyle-> {Thick,Black}];
Show[f1,f2]
]

SelfConsistentDistribution[firstGuess_,errorTolerance_]:=Module[{oldDist,error,errorList},
listvals=firstGuess;
oldDist=listvals;
errorList={};
error=0;
errorTrackingList={};
While[True,
baseMatrix=BaseMatrixTakesChargeList[listvals];
oldDist=listvals; 
listvals = normalizedSumList[numStates];
errorList = oldDist-listvals;
error=errorList.errorList/(NN (MM-1));
listvals=(1-\[Beta])listvals +\[Beta] oldDist; 
listvals=(firstGuess//Abs//Total )/(listvals//Abs//Total) listvals;
AppendTo[errorTrackingList,error];
If[error<errorTolerance,Break[]];
];
{listvals,errorTrackingList}
]

PlotFullSpaceDistribution[listvals_List]:=Module[{interpList,testInterp},
interpList=Flatten[Table[
{{i ,j} ,listvals[[j + (MM-1)(i-1)]]},
{i,1,NN},
{j,1,MM-1}
],1];
testInterp = Interpolation[interpList,InterpolationOrder->1];
Print[Plot3D[testInterp[x,y],{x,1,NN},{y,1,MM-1},PlotRange-> All]]
]

End[]

EndPackage[]
