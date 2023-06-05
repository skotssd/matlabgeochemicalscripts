function [HFO,FeIII,HFOP,AAboundFe,MASSERR]=FeClNaPMgCaSC_AA_tableau(pH,pe,T,flag1,flag2,flag3,flag4,flag5,database)

% input tableau.  change this part % ----------------------------------------------

% constants from PHREEQC database version 3.7.7-1
% struvite Ksp from Aseel's mfile
% strengite from minteq.v4.dat
% HFO surface complexes from minteq.v4.dat (just repeats strong and weak
% twice)
% AA from NIST.  redox potential from Harris 
%(convert to e based log by  Eo/59.2 *2
% Aseel looked up the AA and Ox values. including the Ksp of MG and Fe(II)
% oxalate pptes

Tableau=[...
    %H       %e      %FeIII      %CO3      %PO4  %SO4     %Mg2+ %Ca2+ %Cl- %Na+     %DAA  %HFOs  %HFOw  %logK    %phase   %name
    % component species identity matrix 
    1        0         0         0         0      0       0      0     0    0        0     0     0        0           0     {'H'}
    0        1         0         0         0      0       0      0     0    0        0     0     0        0           0     {'e'}
    0        0         1         0         0      0       0      0     0    0        0     0     0        0           0     {'FeIII'}
    0        0         0         1         0      0       0      0     0    0        0     0     0        0           0     {'CO3'}
    0        0         0         0         1      0       0      0     0    0        0     0     0        0           0     {'PO4'}
    0        0         0         0         0      1       0      0     0    0        0     0     0        0           0     {'NH3'}
    0        0         0         0         0      0       1      0     0    0        0     0     0        0           0     {'Mg'}
    0        0         0         0         0      0       0      1     0    0        0     0     0        0           0     {'Ca'}
    0        0         0         0         0      0       0      0     1    0        0     0     0        0           0     {'Cl'}
    0        0         0         0         0      0       0      0     0    1        0     0     0        0           0     {'Na'}
    0        0         0         0         0      0       0      0     0    0        1     0     0        0           0     {'DAA'}
    0        0         0         0         0      0       0      0     0    0        0     1     0        0           0     {'HFOs'}
    0        0         0         0         0      0       0      0     0    0        0     0     1        0           0     {'HFOw'}
    % Kw
    -1       0         0         0         0      0       0      0     0    0        0     0     0        -14         0     {'OH'}
    % hydrolsis products of Fe
    %-3       0         1         0         0      0       0      0     0    0        0     0     0        -4.891      1     {'HFO'}
    -3       0         1         0         0      0       0      0     0    0        0     0     0        -3.191      1     {'HFO'}
    0        1         1         0         0      0       0      0     0    0        0     0     0        13.02       0     {'FeII'}
    -1       0         1         0         0      0       0      0     0    0        0     0     0        -2.19       0     {'FeIIIOH'}
    -2       0         1         0         0      0       0      0     0    0        0     0     0        -5.67       0     {'FeIIIOH2'}
    -3       0         1         0         0      0       0      0     0    0        0     0     0        -12.56      0     {'FeIIIOH3'}
    -4       0         1         0         0      0       0      0     0    0        0     0     0        -21.6       0     {'FeIIIOH4'}
    -2       0         2         0         0      0       0      0     0    0        0     0     0        -2.95       0     {'FeIII2OH2'}
    -4       0         3         0         0      0       0      0     0    0        0     0     0        -6.3        0     {'FeIII3OH4'}
    -1       1         1         0         0      0       0      0     0    0        0     0     0        3.52        0     {'FeIIOH'}
    %-2       1         1         0         0      0       0      0     0    0        0     0     0        -7.55       0     {'FeIIOH2'}
    %-3       1         1         0         0      0       0      0     0    0        0     0     0        -17.98      0     {'FeIIOH3'}
   % %phosphoric acid and strengite and vivianite
    1        0         0         0         1      0       0      0     0    0        0     0     0        12.346      0     {'HPO4'}
    2        0         0         0         1      0       0      0     0    0        0     0     0        19.553      0     {'H2PO4'}
    3        0         0         0         1      0       0      0     0    0        0     0     0        21.721      0     {'H3PO4'}
    0        0         1         0         1      0       0      0     0    0        0     0     0        26.4        1     {'Strengite'}
    0        3         3         0         2      0       0      0     0    0        0     0     0        75.06       1     {'Vivianite'}
   %  %carbonic acid and FeCO3(s)
    1        0         0         1         0      0       0      0     0    0        0     0     0        10.329      0     {'HCO3'}
    2        0         0         1         0      0       0      0     0    0        0     0     0        16.681      0     {'H2CO3'}
    0        1         1         1         0      0       0      0     0    0        0     0     0        23.91       1     {'Siderite'}
    0        1         1         1         0      0       0      0     0    0        0     0     0        17.4        1     {'FeIICO3'}
    1        1         1         1         0      0       0      0     0    0        0     0     0        25.349      1     {'FeIIHCO3'}
    0        0         0         1         0      0       0      0     0    1        0     0     0        1.270       0     {'NaCO3'}
    0        0         0         1         0      0       1      0     0    0        0     0     0        2.98        0     {'MgCO3'}
    0        0         0         1         0      0       0      1     0    0        0     0     0        3.224       0     {'CaCO3'}
    1        0         0         1         0      0       0      0     0    1        0     0     0        10.079      0     {'NaHCO3'}
    1        0         0         1         0      0       1      0     0    0        0     0     0        11.399      0     {'MgHCO3'}
    1        0         0         1         0      0       0      1     0    0        0     0     0        11.435      0     {'CaHCO3'}
   %  % Fe-P soluble complexes
    1        0         1         0         1      0       0      0     0   0         0     0     0        17.776      0     {'FeIIIHPO4'}
    2        0         1         0         1      0       0      0     0   0         0     0     0        24.983      0     {'FeIIIH2PO4'}
    1        1         1         0         1      0       0      0     0   0         0     0     0        28.966      0     {'FeIIHPO4'}
    2        1         1         0         1      0       0      0     0   0         0     0     0        35.273      0     {'FeIIH2PO4'}
   %  % protonation SO4
    1        0         0         0         0      1       0      0     0   0        0      0     0        1.99        0     {'HSO4'}
    0        1         1         0         0      1       0      0     0   0        0      0     0        15.27       0     {'FeIISO4'}
    1        1         1         0         0      1       0      0     0   0        0      0     0        16.09       0     {'FeIIHSO4'}
    0        0         1         0         0      1       0      0     0   0        0      0     0        4.04        0     {'FeIIISO4'}
    0        0         1         0         0      2       0      0     0   0        0      0     0        5.38        0     {'FeIIISO42'}
   %  % Mg soluble species
    -1       0         0         0         0      0       1      0     0   0        0      0     0        -11.44      0     {'MgOH'}
    0        0         0         1         0      0       1      0     0   0        0      0     0        2.98        0     {'MgCO3'}
    1        0         0         1         0      0       1      0     0   0        0      0     0        11.399      0     {'MgHCO3'}
    0        0         0         0         1      0       1      0     0   0        0      0     0        6.589       0     {'MgPO4'}
    1        0         0         0         1      0       1      0     0   0        0      0     0        15.216      0     {'MgHPO4'}
    2        0         0         0         1      0       1      0     0   0        0      0     0        21.066      0     {'MgH2PO4'}
   %  % ascorbic acid reactions
    2        2         0         0         0      0       0      0     0   0        1      0     0        14.03       0     {'AAH2'}
    1        2         0         0         0      0       0      0     0   0        1      0     0        9.72        0     {'AAH'}
    0        2         0         0         0      0       0      0     0   0        1      0     0        -2.1        0     {'AA'}
    0        2         1         0         0      0       0      0     0   0        1      0     0        7.4         0     {'FeIIIAA'}
    -1       2         1         0         0      0       0      0     0   0        1      0     0        5           0     {'FeIIIOHAA'}
    0        4         1         0         0      0       0      0     0   0        2      0     0        14.6        0     {'FeIIIAA2'}
    -2       2         1         0         0      0       0      0     0   0        1      0     0        1.53        0     {'FeIIIOH2AA'}
    0        3         1         0         0      0       0      0     0   0        1      0     0        16.16       0     {'FeIIAA'}
    1        3         1         0         0      0       0      0     0   0        1      0     0        22.95       0     {'FeIIHAA'}
    % surface complexes PHREEQC
    1        0         0         0         0      0       0      0     0   0        0      1     0        8.93        0     {'HFOsH'}
    2        0         0         0         0      0       0      0     0   0        0      1     0        16.22       0     {'HFOsH2'}
    1        0         0         0         0      0       0      0     0   0        0      0     1        8.93        0     {'HFOwH'}
    2        0         0         0         0      0       0      0     0   0        0      0     1        16.22       0     {'HFOwH2'}
    2        0         0         0         1      0       0      0     0   0        0      0     1        26.65       0     {'HFOwPO4'}
    3        0         0         0         1      0       0      0     0   0        0      0     1        34.32       0     {'HFOwHPO4'}
    4        0         0         0         1      0       0      0     0   0        0      0     1        40.22       0     {'HFOwH2PO4'}
    2        0         0         0         1      0       0      0     0   0        0      1     0        26.65       0     {'HFOsPO4'}
    3        0         0         0         1      0       0      0     0   0        0      1     0        34.32       0     {'HFOsHPO4'}
    4        0         0         0         1      0       0      0     0   0        0      1     0        40.22       0     {'HFOsH2PO4'}
    % assorted weak complexes
    %H       %e      %FeIII      %CO3      %PO4  %SO4     %Mg2+ %Ca2+ %Cl- %Na+     %DAA  %HFOs  %HFOw  %logK    %phase   %name
    0        0         0         0         0      1       0      1     0    0        0     0     0        2.3         0     {'CaSO4'}
    1        0         0         0         1      0       0      1     0    0        0     0     0        15.085      0     {'CaHPO4'}
    2        0         0         0         1      0       0      1     0    0        0     0     0        20.961      0     {'CaH2PO4'}
    0        0         0         0         1      0       0      1     0    0        0     0     0        6.459       0     {'CaPO4'}
    -1       0         0         0         0      0       0      1     0    0        0     0     0        -12.78      0     {'CaOH'}
    1        0         0         0         0      1       0      1     0    0        0     0     0        3.07        0     {'CaHSO4'}
    0        1         1         0         0      0       0      0     1    0        0     0     0        13.16       0     {'FeIICl'}
    0        0         1         0         0      0       0      0     1    0        0     0     0        1.48        0     {'FeIIICl'}
    0        0         1         0         0      0       0      0     2    0        0     0     0        2.13        0     {'FeIIICl2'}
    0        0         1         0         0      0       0      0     3    0        0     0     0        1.13        0     {'FeIIICl3'}
  ];

% end of tableau.  ------------------ % ----------------------------------------------

[KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES]=processtableau(Tableau,pH,pe);

[SPECIESCONCS,SPECIATIONNAMES,MASSERR,X]=returnspeciationRE(KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES,T,flag1,flag2,flag3,flag4,flag5,database);

for k=1:size(SPECIESCONCS,1)
      txt=[SPECIATIONNAMES(k,:),'=SPECIESCONCS(k);'];
      eval(txt)
end

MASSERR=max(MASSERR);

HFOP=HFOwPO4+HFOwH2PO4+HFOsPO4+HFOsHPO4+HFOsH2PO4;

%AAboundFe=FeIIAAH+FeIIAAH2+FeIIIAAH2;
%AAboundFe=FeIIIAAH2;
%AAboundFe=FeIIIAA+FeIIAA;
AAboundFe=FeIIIAA+FeIIIOHAA+FeIIIAA2+FeIIAA+FeIIHAA;
end