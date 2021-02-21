set term post portrait  enhanced color "Helvetica,14"
set output 'BAND.eps'
#set title "Energy band BiSi110 model" font "Curier Bold,18,"
#generated from BAND-VASP.sh 

pi2= 4*atan(1) * 2
    KINIT=       0.000000 ; KNAME_INIT="{/Symbol \G}"
       K2=       1.480653 ; KNAME_2   ="M"
       K3=       2.335508 ; KNAME_3   ="K"
       K4=       4.045218 ; KNAME_4   ="{/Symbol \G}"
     KEND=       4.254657 ; KNAME_END ="A"
 set xtics (KNAME_INIT KINIT, KNAME_2 K2, KNAME_3 K3, KNAME_4 K4, KNAME_END KEND) nomirror
 set xrange[KINIT:KEND]
 set yrange[-6.:6.]
set ytics -6,3,6 nomirror
set ylabel "Energy (eV)" font "Helvetica,24" offset 1,0
 set arrow from 0,0 to KEND,0 nohead lt 3 lw .1 lc rgb 'black'
cup1='#1e90ff'    ##dogger blue
cup2='#5f9ea0'    ##Cadet Blue
cup3='#b22222'    ##Firebrick
cup4='#e6e6fa'    ##Lavender
cup5='#2e8b57'    ##Sea Green
cup6='#f4a460'    ##Sandy Brown
cdn1='#ff4500'    ##Orange Red
cdn2='#bdb76b'    ##dark khaki
cdn3='#ffc0cb'    ##Pink
cdn4='#8a2be2'    ##Blue Violet
cdn5='#ffd700'    ##Gold
cdn6='#cdc9c9'    ##Snow 3
zecolor='#00ced1'
set key right bot samplen 1  spacing 1.5  font "Helvetica, 15"
set palette defined (-10 cup1,0 'white', 10 cdn1 )
set cbrange [-10:10]
 set multiplot
 set origin 0,0
 set size nosquare 0.6,0.4 
 s=1

 result_dft='DOS_atom_projected.dat'
 result_tba='band_structure_TBA.dat'
 result_eff='band_structure_TBA_EFF.dat'
 result_eff='band_structure_TBA_EFF.dat'
 target_dft='band_structure_DFT.dat'

  plot result_tba          u 1:2        w l         lw  1            lc rgb 'black'   ti "TBA",\
       result_tba          u 1:2:($6)*s     w p  ps vari  lw  .1 pt 7      lc rgb 'blue'   ti "p_z",\
       result_eff          u 1:2            w p  ps  .5   lw  1  pt 6      lc rgb 'red'   ti "effective",\
