set term post portrait  enhanced color "Helvetica,15"
set output 'BAND.eps'
#set title "Energy band BiSi110 model" font "Curier Bold,15,"
#generated from BAND-VASP.sh 

pi2= 4*atan(1) * 2
    KINIT=       0.000000 ; KNAME_INIT="-Y"
       K2=       0.393219 ; KNAME_2   ="{/Symbol \G}"
     KEND=       0.788434 ; KNAME_END ="Y"
 set xtics (KNAME_INIT KINIT,  KNAME_2 K2, KNAME_END KEND) nomirror
 set size nosquare 1,0.55 ;
 set xrange[KINIT:KEND]
 set yrange[-0.5:0.5]
set ytics -0.5,.5,0.5 nomirror
set ylabel "Energy (eV)" font "Helvetica,24" offset 1,0
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
set key right top samplen 1  spacing 1.   font "Helvetica, 15"
set palette defined (-10 cup1,0 'gray', 10 cdn1 )
set cbrange [-10:10]

 set multiplot
 set origin 0,0
 set size nosquare 0.6,0.4
 s=.6

 result_dft='DOS_atom_projected.dat'
 result_tba='band_structure_TBA.dat'
 result_brc='BERRYCURV_TBA.total.dat'
 target_dft='band_structure_DFT.dat'
 result_wcc='Z2.WCC.0.0-B3.B1_B2-PLANE.dat'
 result_gap='Z2.GAP.0.0-B3.B1_B2-PLANE.dat'
 result_crn='Z2.GAP.0.0-B3.B1_B2-PLANE.dat'

# plot result_tba          u 1:2:($3  +$4  +$5  +$6  )*.5  w p lw .5 ps vari pt 6 lc rgb 'red'    ti "L_A",\
#      result_tba          u 1:2:($51 +$52 +$53 +$54 )*.5  w p lw .5 ps vari pt 6 lc rgb 'blue'   ti "L_B",\
#      result_tba          u 1:2:($47 +$48 +$49 +$50 )*.5  w p lw .5 ps vari pt 4 lc rgb 'orange' ti "R_A",\
#      result_tba          u 1:2:($31 +$32 +$33 +$34 )*.5  w p lw .5 ps vari pt 4 lc rgb 'green'  ti "R_B",\
#      result_tba          u 1:2        w l       lw 0.1          lc rgb 'black' noti

  plot result_tba          u 1:2        w l       lw 0.1          lc rgb 'black' noti

