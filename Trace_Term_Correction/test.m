(* ::Package:: *)

Print[Integrate[Log[1-x],{x,0,1}]]


Print[Integrate[Log[1-x]^2,{x,0,1}]]


Export["dat1.wdx",NIntegrate[Log[1-x]^2,{x,0,1}]]


Export["dat2.wdx",Integrate[Log[1-x]^2,x]]
