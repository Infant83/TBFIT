%% Tight binding program suited for 1T'-MoTe2 system: function for Hamiltonian matrix
%% Written by Dr. Hyun-Jung Kim (Korea Institute for Advanced Study)
%% Infant@kias.re.kr

 function [H]=ham(k_,p)

%predifine orbital index (should not be modified!!)
s=6;   px=11; py=12; pz=13;   dz2=21; dx2=22; dxy=23; dxz=24; dyz=25;

%define parameter set and scaling (decaying) factor sets (valid for when scaling function type 1 to 3 is called)
% the parameter set sequence should be sigma(s), pi(p), and delta(d)
%          SK-parameters                decaying factor
%           sig,pi,del   sig,pi,del  
 edz2 =  p( 1) ; 
 edx2 =  p( 2) ;
 dds  =  p( 3) ;
 ddp  =  p( 4) ;
 ddd  =  p( 5) ;

 a0 = 3.1716343032172727;
 R1 = [ a0 0 0];
 R2 = -R1 ;
 R3 = [1.5858171516086363    2.7467158781003151    0.0000000000000001];
 R4 = -R3 ;
 R5 = R3 - R1 ;
 R6 = -R5 ;

 t11 = 0.75*ddd + 0.25*dds ;
 t12 = 0.375*ddd - 0.375*dds ;
 t13 = sqrt(3)/4 * ( ddd - dds ) ;
 t22 = ddp ;
 t22_= 3/16*ddd + 1/4*ddp + 9/16*dds ;
 t23 =-sqrt(3)/16*ddd + sqrt(3)/4*ddp -sqrt(3)*3/16*dds ; 
 t33 = 0.25*ddd + 0.75*dds ;
 t33_= 0.0625*ddd + 0.75*ddp + 0.1875*dds ;


    H11 = edz2+ t11*( fk(k_,R1) + fk(k_,R2) + fk(k_,R3) + fk(k_,R4) + fk(k_,R5) + fk(k_,R6)) ;
    H12 =       t12*( fk(k_,R3) + fk(k_,R4) - fk(k_,R5) - fk(k_,R6) ) ;
    H13 =       t13*( fk(k_,R1) + fk(k_,R2) - fk(k_,R3)/2 - fk(k_,R4)/2 - fk(k_,R5)/2 - fk(k_,R6)/2 ) ;

    H22 = edx2+(t22*( fk(k_,R1) + fk(k_,R2))+ t22_*(fk(k_,R3) + fk(k_,R4) +fk(k_,R5) + fk(k_,R6)) ) ;
    H23 =       t23*( fk(k_,R3) + fk(k_,R4) - fk(k_,R5) - fk(k_,R6) ) ;

    H33 = edx2+ t33*( fk(k_,R1) + fk(k_,R2))+ t33_*(fk(k_,R3) + fk(k_,R4) +fk(k_,R5) + fk(k_,R6)) ;

%      1    2    3  
%      A    B    C  
%     dz2  dxy  dx2 
 H = [H11  H12  H13 ,
      H12' H22  H23 ,
      H13' H23' H33 ];

%kx = k_(1)*a0/2; ky = k_(2)*a0/2*sqrt(3);
%t0 = -0.184 ; t1 = 0.401; t2 = 0.507 ; t11 = 0.218 ; t12 = 0.338 ; t22 = 0.057 ;
%e1 = 1.046 ; e2 = 2.104;
%H0 = e1 + 2* t0 *( cos(2*kx) + 2 * cos(kx) * cos(ky) ) ;
%H1 =    -2*sqrt(3)*t2*sin(kx)*sin(ky) + 2*1i*t1*(sin(2*kx) + sin(kx)*cos(ky) ) ;
%H2 =     2*t2*(cos(2*kx) - cos(kx)*cos(ky)) + 2*sqrt(3)*1i*t1*cos(kx)*sin(ky) ;
%H11= e2+ 2*t11*cos(2*kx) + (t11+ 3*t22) * cos(kx)*cos(ky) ;
%H22= e2+ 2*t22*cos(2*kx) + (3*t11+t22) * cos(kx)*cos(ky) ;
%H12=     sqrt(3)*(t22 - t11) * sin(kx)*sin(ky) + 4*1i*t12*sin(kx)*(cos(kx) - cos(ky)) ;

%H = [H0  H1   H2;
%     H1' H11  H12;
%     H2' H12' H22];


endfunction
