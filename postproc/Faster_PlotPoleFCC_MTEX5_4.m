%% Import Script for ODF Data
%
% This script was automatically created by the import wizard. You should
% run the whoole script or parts of it in order to import your data. There
% is no problem in making any changes to this script.

%% Specify Crystal and Specimen Symmetries
clear ;
close all ;

addpath('~/Documents/research/software/MATLABlibrary/mtex-5.4.0/')
startup_mtex

CS = crystalSymmetry('cubic');
%% Specify sample symmetry
SS = specimenSymmetry('triclinic');
% plotting convention
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','outOfPlane');

%% Specify File Names
% path to files
%pname = '/Users/kteferra/Documents/research/projects/AMICME/codes/CA/results/onyx/8/';
pname = '/Users/kteferra/Documents/research/projects/LLNL/';
%pname = '/Users/kteferra/Documents/research/projects/AMICME/codes/CA/tmp/';
% which files to be imported
%fname = [pname 'CA3DsimNRLdata1_1547_vf27_DREAM3d_ODF.csv'];
%fname = [pname 'CA3DsimNRLdata1_1547_ODF.csv'];
%fname = [pname 'CA3Dsim6_719_ODF.csv'];
%fname = [pname 'CA3Dnlv7_730_ODF.csv'];
fname = [pname 'CA3Dnlv7_2040_ODF.csv'];
%fname = [pname 'CA3Dnlv8_612_ODF.csv'];
%fname = [pname 'Orient_Samples.csv'];
%% Import the Data
% specify kernel
psi = deLaValleePoussinKernel('halfwidth',5*degree); % this is obsolete
%but works with MTEX version 5.4

% load the ODF into the variable odf
odf = ODF.load(fname,CS,SS,'density','kernel',psi,'resolution',5*degree,...
  'interface','generic',...
  'ColumnNames', { 'phi1' 'Phi' 'phi2' 'Weight'}, 'Columns', [1 2 3 4], 'Bunge');
h = [Miller(0,0,1,CS,SS),Miller(0,1,1,CS,SS),Miller(1,1,1,CS,SS)];
%plotPDF(odf,h,'antipodal','silent',CS);
pf1 = calcPoleFigure(odf,Miller({0,0,1},h.CS),'resolution',5*degree,'complete') ;
pf2 = calcPoleFigure(odf,Miller({0,1,1},h.CS),'resolution',5*degree,'complete') ;
pf3 = calcPoleFigure(odf,Miller({1,1,1},h.CS),'resolution',5*degree,'complete') ;
figure ;
plot(pf1','smooth','colorrange','equal') ;
colorbar ;
figure ;
plot(pf2,'smooth','colorrange','equal') ;
colorbar  ;
figure ;
plot(pf3,'smooth','colorrange','equal') ;
colorbar ;


fodf = FourierODF(odf,256);
% f_hat = calcFourier(fodf,2) ;
plotFourier(fodf);

