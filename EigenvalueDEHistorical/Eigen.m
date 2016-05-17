(* ::Package:: *)

BeginPackage["Eigen`"]


(*Public functions*)
ClearAll[EigenNDSolve,PlotSpectrum,funcFromCheby]

EigenNDSolve::usage="EigenNDsolve[\!\(\*
StyleBox[\"eqns\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"{\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[
StyleBox[\"u\",\nFontSlant->\"Italic\"], \"1\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"u\", \"2\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[Ellipsis]\",\nFontSlant->\"Italic\"]\)}, {\!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"min\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"max\"],\nFontSlant->\"Italic\"]\)}, \!\(\*
StyleBox[\"\[Omega]\",\nFontSlant->\"Italic\"]\), Options]
Numerically solves for the eigenfunctions and eigenvalues of \!\(\*
StyleBox[\"eqns\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)using a spectral 
expansion in Chebyshev polynomials.
\!\(\*
StyleBox[\"{\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[
StyleBox[\"u\",\nFontSlant->\"Italic\"], \"1\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"u\", \"2\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[Ellipsis]\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"}\",\nFontSlant->\"Italic\"]\) is a list of the dependent variables;
\!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\) is the independent variable, \!\(\*
StyleBox[SubscriptBox[\"x\", \"min\"],\nFontSlant->\"Italic\"]\) and \!\(\*
StyleBox[SubscriptBox[\"x\", \"max\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)is the range;
\!\(\*
StyleBox[\"\[Omega]\",\nFontSlant->\"Italic\"]\) is the eigenvalue parameter (\!\(\*
StyleBox[\"eqns\",\nFontSlant->\"Italic\"]\) must be linear in \[Omega]).

Options:
Npoly\[Rule]80 - Number of polynomials to use in Chebyshev expansion (actually Npoly+1 are used 
	due to \!\(\*SubscriptBox[\(T\), \(0\)]\)).
returnInterp\[Rule]True - If True, returns lists of eigenfunctions as InterpolatingFunction objects.
	If False, returns a list of the Chebyshev coefficients for each eigenstate.
precision\[Rule]MachinePrecision - Precision of the calculation. If not MachinePrecision,
	calculation will be much slower; however this is necessary for accurate calculation
	with certain types of differential equations.
useParallel\[Rule]False - Use parallelization in computing the operator matrices. Will only 
	substantially improve speed with particulatly complex differential equations.
noRowReduce\[Rule]False - If True, over-rides solving for the boundary conditions and forming a
	reduced matrix (see Dongarra et al. 1997). The default option, False, will remove singular
	matrices from the generalized eigenvalue problem and remove spurious eigenvalues in
	certain cases.
	After some experimentation, it seems that this step generally does not substantially 
	improve results aside from removing spurious eigenvalues.";



PlotSpectrum::usage="PlotSpectrum[\!\(\*
StyleBox[\"esystem\",\nFontSlant->\"Italic\"]\), {{\!\(\*
StyleBox[SubscriptBox[\"x\", \"min\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"max\"],\nFontSlant->\"Italic\"]\)},{\!\(\*
StyleBox[SubscriptBox[\"y\", \"min\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"y\", \"max\"],\nFontSlant->\"Italic\"]\)}}, Options]
Plots the eigenvalue spectrum \!\(\*
StyleBox[\"esystem\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)(as output from EigenNDSolve) over the region 
(\!\(\*
StyleBox[SubscriptBox[\"x\", \"min\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[Rule]\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"max\"],\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[SubscriptBox[\"y\", \"min\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[Rule]\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"y\", \"max\"],\nFontSlant->\"Italic\"]\)), where \!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)and \!\(\*
StyleBox[\"y\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)refer to the real and imaginary axes respectively.
Clicking on eigenvalues will display the corresponding eigenfunction in the inset and 
with mouseover the value of the eigenvalue is displayed.

Options:
insetSize\[Rule]0.5 - Size of the inset eigenfunction plot in units of the full plot size.
insetPosition->{0,0} - Position of the lower LH corner of the inset plot in units scaled
	to the full plot size.";




funcFromCheby::usage="funcFromCheby[\!\(\*
StyleBox[\"cheblist\",\nFontSlant->\"Italic\"]\), Options]
Outputs an InterpolatingFunction object given by the list of Chebyshev coefficients
on the range {-1,1}.
funcFromCheby[\!\(\*
StyleBox[\"cheblist\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"{\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[
StyleBox[\"x\",\nFontSlant->\"Italic\"], \"min\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"max\"],\nFontSlant->\"Italic\"]\)}, Options]
Outputs an InterpolatingFunction object given by the list of Chebyshev coefficients
on the range {\!\(\*
StyleBox[SubscriptBox[
StyleBox[\"x\",\nFontSlant->\"Italic\"], \"min\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"x\", \"max\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"}\",\nFontSlant->\"Italic\"]\).

Options:
Npoints\[Rule]101 - Number of points from which to construct the InterpolatingFunction.
returnList\[Rule]False - Returns a list of real space points rather than an InterpolatingFunction.";



(* ::Section:: *)
(*Main Package*)


Begin["`Private`"]




Options[EigenNDSolve]={Npoly->80,returnInterp->True,noRowReduce->False,
	precision->MachinePrecision,useParallel->False};
EigenNDSolve::eqns="Supplied equations not in the correct form!";
EigenNDSolve::inhomog="Upgrade needed to deal with inhomogenous equations!";
EigenNDSolve::numeq="Number of variables does not match number of equations!";
EigenNDSolve::largeEl="Some derivative matrix elements are very large, `1`, and may be ill represented by machine precision numbers. Could increase Option, precision.";
EigenNDSolve::bounds="Unexpected boundary conditions supplied: `1`";
EigenNDSolve::noVBC="Not performing row reduction on boundary conditions due to unexpected boundary
	condition properties. Spurious eigenvalues will probably be produced";
EigenNDSolve::findzeros="Unable to find correct combination of non-zero elements to allow 
row-reduction. Proceeding calculation with full matrix";

EigenNDSolve[eqns_List,varsin_,range_List,evalvar_,OptionsPattern[]]:=
	Module[{vars,z(*z is internal independent variable to use*),eqnBC,eqn,BC,dord,varQ,
eqnmat(*Terms multiplying each derivative of each variable*),BCnum,rowRedQ,inhomog,
T,dt,brules(*Abstract Chebyshev polynomials*),createOpMat,Npol=OptionValue[Npoly],Nprec=OptionValue[precision],sizelist,
OPnBC,BCrows,lv,chebvs,OPfull,
$dcon,myMap,
bcolind,browind,zeropos,choice,rowRinds,
BCfinal,BCSolfunc,OPfinal,partsfun,
eigsys},

(*Add list for single variable*)
If[Head[varsin]=!=List,vars={varsin},vars=varsin];
(*Put on a {-1,1} region*)
varQ=(MemberQ[vars,#]&);(*Test is a variable is a var*)
eqnBC=transCo[eqns,varQ,range,{z,-1,1}];
(*Split into equations and boundary conditions*)
eqn=Select[eqnBC,MemberQ[#,f_?varQ[z]|Derivative[_][f_?varQ][z],Infinity]&];
BC=Complement[eqnBC,eqn];

If[Length[eqn]!= Length[vars],Message[EigenNDSolve::numeq];Abort[]];
(*Maximum derivative order list for each variable*)
dord=Max[Cases[eqns,Derivative[n_][#][_]:>n,Infinity]]&/@vars;
(*Convert eqns into single expression that equals zero*)
{eqn,BC}=If[And@@(Head[#]===Equal&/@#),Expand/@(Subtract@@@#),Message[EigenNDSolve::eqns];Abort[];]&/@{eqn,BC};
(*Check boundary conditions and calculate number of extra polynomials for each variable*)
{BCnum,rowRedQ}=checkBCform[BC,vars,dord];
If[OptionValue[noRowReduce],rowRedQ=False;];

(*Put equations into "matrix" form*)
eqnmat=Function[{expr},Function[{var},Coefficient[expr,Derivative[#][var][z]]&/@Range[0,Max@@dord]]/@vars]/@eqn;
inhomog=Join[eqn,BC]/.{f_?varQ[_]->0,Derivative[_][f_?varQ][_]->0};
If[Or@@(0=!=#&/@Chop[inhomog]),Message[EigenNDSolve::inhomog];Print[inhomog];Abort[]];

(*Define abstract Chebyshev properties to use in constructing matrix*)
brules=defineAbCheb[T,dt];


$dcon=If[OptionValue[useParallel],DistributedContexts -> {"Eigen`Private`"},{1}];
myMap=If[OptionValue[useParallel],ParallelMap,Map];
(*Function to make matrices*)
createOpMat[coef_,size_List]:=Block[{do=Length[coef],consFromT,nozmats,zfuncs,chebexpans,phivec,diffeq},
	(*Size is {N,N+BC} so actual matix size will be {N+1,N+BC+1}*)
	consFromT={Complement[Range[do],#],#}&@Select[Range[do],MemberQ[coef[[#]],z,{0,Infinity}]&];
	(*Calculate mat for those that don't contain z, this is seperate for speed*)
	nozmats=Total[coef[[#]]constructDmat[#-1,size+1]&/@consFromT[[1]]];
	If[nozmats===0,nozmats=SparseArray[{1,1}->0,size+1]];
	If[Max@Abs[Flatten[nozmats]]>10^17&&OptionValue[precision]===MachinePrecision,
			Message[EigenNDSolve::largeEl,Max@Abs[Flatten[nozmats]]>10^17]];

	(*Calculate for those that do contain z*)
	(*First, get Chebyshev expansions*)
	zfuncs=Function[Evaluate@{z},Evaluate[#]]&/@coef[[consFromT[[2]]]];
	phivec=(T[#]&/@Range[0,size[[2]]]);
	chebexpans=Map[Dot[#,phivec]&,Chop[CCapprox[#,size[[2]],evalvar,Nprec],10^-14]&/@zfuncs];
	(*Plug this into Chebyshev expansion*)
	
	diffeq=Total[myMap[ExpandAll,
		chebexpans (Nest[myMap[Expand[dt[#]]&,#,
			$dcon]&,phivec,#]&/@(consFromT[[2]]-1)),$dcon]];
	nozmats+Join[{diffeq/.T[_]->0},
		myMap[Coefficient[diffeq,T[#]]&, Range[size[[1]]],$dcon]]
];

(*Create full matrix without Boundary conditions*)
sizelist=ConstantArray[{Npol,#}&/@(Npol+BCnum),Length[BCnum]];

(*This is just an implementation of MapThread[f,l,2] using only Map, as Parallelizing was somewhat broken*)
OPnBC=N[Chop[ArrayFlatten[Map[createOpMat@@##&,Transpose[{eqnmat,sizelist},{3,1,2}],{2}
		]],10^-14],Nprec];


(*Formulate boundary conditions and add into matrix*)
lv=Npol+BCnum;
chebvs=Map[Flatten@MapThread[#1(T[#]&/@Range[0,#2])&,{vars,lv}]/.#->1/.Thread[vars->0]&,vars];
BCrows=(BC/.Thread[vars->chebvs])/.{f_List[pm1_]:>(f/.brules[[-pm1]]),Derivative[n_][f_List][pm1_]:>
	(Nest[Map[Expand[dt[#]]&,#]&,f,n]/.brules[[-pm1]])};
OPfull=Join[OPnBC,N[BCrows,Nprec]];

If[rowRedQ,
(*First need to find {row,column} boundary points to be able to do a reduction*)
(*Indices of "tau" polynomials*)
bcolind=Sort@Flatten@MapThread[#1-Range[#2]&,{Rest@Accumulate[Prepend[lv+1,1]],BCnum}];
browind=Length[vars](Npol+1)+Range[Total[BCnum]];
(*Positions of zeros in each column*)
zeropos=Flatten@Position[OPfull[[browind,#]],a_/;a!=0]&/@bcolind;
(*Choose a combination that gives all 4 numbers without zero*)
choice=Flatten@Select[Tuples[zeropos],Sort[#]==Range[Length[zeropos]]&,1];


If[Length@choice==0,
(*In case it doesn't work, keep the program working*)
Message[EigenNDSolve::findzeros];OPfinal=OPfull,

rowRinds={browind[[choice]],bcolind}\[Transpose];
OPfull=Fold[RowR,OPfull,rowRinds];

BCfinal=OPfull[[browind]];
(*Function to solve for BCs given solution vector*)
BCSolfunc[solvec_]:=Block[{aBC,nbc=Range@Length[bcolind],fullsolvec},
fullsolvec=Fold[Insert[#1,#2[[2]],#2[[1]]]&,solvec,{bcolind,aBC/@nbc}\[Transpose]];
fullsolvec/.First@Solve[Chop[BCfinal.{fullsolvec}\[Transpose]]==0,aBC/@nbc]
];

partsfun=Complement[Range[Total[lv+1]],#]&;
OPfinal=OPfull[[partsfun[rowRinds[[All,1]]],partsfun[rowRinds[[All,2]]]]]
],

OPfinal=OPfull
];



(*Normal operation*)
eigsys=Eigensystem[{-OPfinal/.evalvar->0,Coefficient[OPfinal,evalvar]}];
If[rowRedQ,
eigsys[[2]]=Map[BCSolfunc,eigsys[[2]]]];
eigsys[[2]]=Internal`PartitionRagged[#,lv+1]&/@eigsys[[2]];
(*Order according to size of imaginary part*)
eigsys=Join[{eigsys[[1]]},eigsys[[2]]\[Transpose]][[All,#]]&@
	Ordering[eigsys[[1]],All,
		If[Head[#1]===Complex&&Head[#2]===Complex,Im[#1]>Im[#2],Abs[#2]>Abs[#1]]&];

If[OptionValue[returnInterp],
eigsToInterp[eigsys,Rest@range],
eigsys]

]



(* ::Subsection:: *)
(*Plotting*)


(*Plotting e-values*)
plotnum=1;
Options[PlotSpectrum]={insetPosition->{0,0},insetSize->0.5};
PlotSpectrum::estateform="E-states must be given in the form of an interpolating function";
PlotSpectrum[evallist_List,range_,OptionsPattern[]]:=DynamicModule[{imre,evals,efuncs,funrange},
imre={Re[#],Im[#]}&;

efuncs=Rest@evallist;
evals=First@evallist;
If[Head[efuncs[[1,1]]]=!=InterpolatingFunction,Message[PlotSpectrum::estateform],
funrange=efuncs[[1,1,1]]];

Show[
ListPlot[MapIndexed[Button[Tooltip[#1,ToString[#2[[1]]]<>": "<>ToString[#1[[1]]+I #1[[2]]]],
	plotnum=#2[[1]]]&,imre/@evals],PlotStyle->Directive[PointSize[Medium],Black],
	PlotRange->range],
Epilog->Inset[  
	Dynamic[intFuncsPlot[efuncs[[All,plotnum]]]],
		ImageScaled[OptionValue[insetPosition]],ImageScaled[{0,0}],
			OptionValue[insetSize](range[[1,2]]-range[[1,1]])]
]


]


(* ::Subsection:: *)
(*Auxilary functions*)


ClearAll@constructDmat;
constructDmat[ord_,size_List]:=Block[{lsiz,arrrules,tempmat},
If[ord==0,
SparseArray[Band[{1,1}]->ConstantArray[1,Min[size]],size],
lsiz=Max[size];
(*Since Dmat will be very common, construct using formula for speed*)
arrrules=Flatten[Table[{{1,2j}->2j-1,{i+1,i+2j}->2(i+2j-1)},{i,1,lsiz},{j,1,lsiz}]]/.({_,n_/;n>lsiz}->_):>Sequence[];
tempmat=SparseArray[arrrules,{lsiz,lsiz}];
MatrixPower[tempmat,ord][[1;;size[[1]],1;;size[[2]]]]]
]



ClearAll@CCapprox;
CCapprox[func_,n_,eigv_,prec_:MachinePrecision]:=Block[{cnodes,funclist,cc},
cnodes=N[Cos[Pi Range[0,n]/n],prec];
funclist={#/.eigv->0,Coefficient[#,eigv]}&[func/@cnodes];
cc=FourierDCT[#,1]Sqrt[2/n]&/@funclist;
cc[[All,{1,-1}]]/=2;
cc=SetPrecision[cc,prec];
cc[[1]]+eigv cc[[2]]]


ClearAll@transCo;
transCo[eqn_,varQ_,xrange_,zrange_]:=
	Module[{x=xrange[[1]],xr=xrange[[2;;]],z=zrange[[1]],zr=zrange[[2;;]],trans},
		trans={Subtract@@xr/Subtract@@zr,(xr[[2]]zr[[1]]-xr[[1]]zr[[2]])/(zr[[1]]-zr[[2]])+Subtract@@xr/Subtract@@zr #&,
				(zr[[2]]xr[[1]]-zr[[1]]xr[[2]])/(xr[[1]]-xr[[2]])+Subtract@@zr/Subtract@@xr #&};
		eqn/.Derivative[n_][f_?varQ][v_]:>1/ (trans[[1]]^n) Derivative[n][f][trans[[3]][v]]/.
			{f_?varQ[v_]:>f[trans[[3]][v]]}/.x->trans[[2]][z]//ExpandAll
]



ClearAll[calc$chebyshevLists];
(*Produces interpolating function from a list of Chebyshev coefficients*)
calc$chebyshevLists[npol_,npoints_]:=$chebyshevLists=Table[N[ChebyshevT[i,Range[-1,1,2/(npoints-1)]]],{i,0,npol}];
calc$chebyshevLists[120,101];

Options[funcFromCheby]={Npoints->101,returnList->False};
funcFromCheby[chebco_,rang_:{-1,1},OptionsPattern[]]:=Module[{},
	If[Length[chebco]>Length@$chebyshevLists||OptionValue[Npoints]!=Length@First@$chebyshevLists,
	calc$chebyshevLists[Length[chebco],OptionValue[Npoints]]];
	If[OptionValue[returnList],
		chebco.$chebyshevLists[[1;;Length[chebco]]],
		ListInterpolation[chebco.$chebyshevLists[[1;;Length[chebco]]],{rang}]
	]
]




ClearAll@intFuncsPlot;
(*Plots a list of interpolating functions in a nicish way, showing real and imaginary parts*)
intFuncsPlot[intfuncs_List]:=Module[{plotlist,range,collist,x},
plotlist=Flatten[{Re@#,Im@#}&/@Through[intfuncs[x]]];
range=intfuncs[[1,1,1]];

collist=Flatten[{#,Directive[#,Dashed]}&/@(ColorData[1,#]&/@Range[Length[intfuncs]])];

Plot[plotlist,Evaluate[{x,Sequence@@range}],PlotRange->All,PlotStyle->collist]
]




ClearAll[RowR]
RowR[mat_,{n_,m_}]:=Block[{matt,i},
matt=Table[-(mat[[i,m]]/mat[[n,m]])mat[[n]]+mat[[i]],Evaluate[{i,Complement[Range[Length[mat]],{n}]}]];
Chop[Insert[matt,mat[[n]],n],10^-14]]


ClearAll[checkBCform]
checkBCform[BCs_,vars_,dord_]:=
Block[{totalBCs=Length[BCs],BClist,BCmatch,fornone,RRQ,nBCs,lens,vl,vcl,overlaps,numvars,numBCfin},
BClist=Function[{var},Map[Cases[#,var[_]|Derivative[_][var][_],{0,Infinity}]&,BCs]]/@vars;
If[Or@@Thread[Abs[Flatten[BClist][[All,1]]]!=1],
	Message[EigenNDSolve::bounds,Thread[BCs==0]];Abort[]];
(*Don't want two of the same variable in each BC list, then double count*)
fornone=If[Length[#]!=0,{First@#},#]&;
BClist=Map[fornone,BClist,{2}];
(*Check to see if BCs can be "row reduced" from the matrix*)
nBCs=Length/@(Flatten/@BClist);
BCmatch=Thread[dord<=nBCs];

(*The True/False returned refers to whether the row reduction step should be performed*)
RRQ=If[Length[BClist[[1]]]==Total[dord] && And@@BCmatch,True,
	Message[EigenNDSolve::noVBC];False];

(*Choose the number of extra polynomials needed for each variable*)
lens=Map[Length,BClist,{2}];
vl=Range[Length[vars]];
(*Finds the number of overlapping boundary conditions for each variable*)
overlaps=Function[{vnum},Total[Flatten[If[lens[[vnum,#]]!=0,
	Boole[#!=0]&/@lens[[Complement[vl,{vnum}],#]],0]&/@Range[totalBCs]]]]/@vl;

If[Length[BClist[[1]]]!=Total[dord],
(*If number of boundary conditions doesn't match dord, just reduce some choice of those that overlap*)
	numvars[nbc_,maxred_]:=Block[{i=nbc,j=maxred},While[j>0,--i;--j];i];
	numBCfin=MapThread[numvars,{nBCs,overlaps}],
(*Else work out possible reduction from dord*)
(*number of BCs,possible reduction,dord*)
	numvars[nbc_,maxred_,do_]:=Block[{i=nbc,j=maxred},While[i>do&&j>0,--i;--j];i];
	numBCfin=MapThread[numvars,{nBCs,overlaps,dord}]
];
{numBCfin,RRQ}

]



ClearAll[defineAbCheb]
defineAbCheb[T_,dt_]:=Module[{},
(*Call this new each time to stop variable overlap*)
ClearAll[dt,T];
(*derivative properties (needed if we have a(x)\!\(
\*SubscriptBox[\(\[PartialD]\), \(x\)]\[Xi]\) since only \!\(
\*SubscriptBox[\(\[PartialD]\), \({x, n}\)]\[Xi]\) is stored*)
dt[a_+b_]:=dt[a]+dt[b];
dt[a_ b_]:=a dt[b]+b dt[a];
dt[a_?NumericQ]=0;
(*Derivatives of polynomials*)
T[0]:=1;
dt[T[n_]]:=If[EvenQ[n-1],Expand[n (2Sum[T[j],{j,0,n-1,2}]-1)],Expand[n(2Sum[T[j],{j,1,n-1,2}])]];
(*Products of polynomials*)
T/:T[i_Integer]T[j_Integer]:=1/2(T[i+j]+T[Abs[i-j]]);
T/:T[j_Integer]^2:=1/2(T[2j]+T[0]);
(*Rules for boundaries*)
{T[j_]:>(-1)^j,T[j_]:>1}
]



ClearAll[calc$chebyshevLists]
(*Produces interpolating function from a list of Chebyshev coefficients*)
calc$chebyshevLists[npol_,npoints_]:=$chebyshevLists=Table[N[ChebyshevT[i,Range[-1,1,2/(npoints-1)]]],{i,0,npol}];
calc$chebyshevLists[120,101];
Options[funcFromCheby]={Npoints->101,returnList->False};
funcFromCheby[chebco_,rang_:{-1,1},OptionsPattern[]]:=Module[{},
	If[Length[chebco]>Length@$chebyshevLists||OptionValue[Npoints]!=Length@First@$chebyshevLists,
	calc$chebyshevLists[Length[chebco],OptionValue[Npoints]]];
	If[OptionValue[returnList],
		chebco.$chebyshevLists[[1;;Length[chebco]]],
		ListInterpolation[chebco.$chebyshevLists[[1;;Length[chebco]]],{rang}]
	]
]

eigsToInterp[eigsys_,range_:{-1,1}]:=Block[{esys=eigsys},
esys[[2;;]]=Map[funcFromCheby[#,range]&,esys[[2;;]],{2}];
esys]






ClearAll[remEnds,addEnds]
remEnds=Most@Rest@#&;
addEnds=Join[{0},#,{0}]&;



End[]
EndPackage[]
