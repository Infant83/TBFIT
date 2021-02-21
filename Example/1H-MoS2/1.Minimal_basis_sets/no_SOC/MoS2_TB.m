#!~/usr/bin/local/octave
clc
clear all

%program started...
setenv ("GNUTERM", "x11");
%   global verbose;verbose = 1;
	ifit=0
    pkg load tbfit
	fname_TBA_out="band_structure_TBA_tbfit.dat"

%set unitcell information
	a0 = 3.1716343032172727 ;
	a1=[a0, 0, 0];  a2=[a0/2 , a0/2*sqrt(3), 0 ];  a3=[0, 0, 15. ];
	[b1,b2,b3]=get_reci(a1,a2,a3);

%set k-path
	K(:,1)=[1/2, 1/2,  0]; %M
	K(:,2)=[1/3, 2/3,  0]; %K
	K(:,3)=[0  ,   0,  0]; %G
	K(:,4)=[2/3, 1/3,  0]; %K'
	K(:,5)=[1/2, 1/2,  0]; %M
	nline=size(K)(2)-1;

    p_name=importdata('param_tbfit.dat','=').rowheaders'; p=importdata ('param_tbfit.dat','=').data';

	[E,V,C,NE,NK,kdist]=get_ham_en_line(K,b1,b2,b3,p,nline,ndiv=50,'ham'); 

%write to the file
	print_energy(NE,fname_TBA_out,E',V,kdist(:)',IORB=1) ; % IORB=1 print orbital composition
    system ("gnuplot gnuBAND-TB_.gp") % run gnuplot
