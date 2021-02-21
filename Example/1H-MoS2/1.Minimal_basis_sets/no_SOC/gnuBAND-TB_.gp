set term post portrait  enhanced color "Helvetica,18"
set output 'BAND_tbfit.eps'
#set title "Energy band BiSi110 model" font "Curier Bold,18,"
#generated from BAND-VASP.sh 

pi2= 4*atan(1) * 2
    KINIT=       0.000000 ; KNAME_INIT="M"
       K2=       0.660352 ; KNAME_2   ="K"
       K3=       1.981056 ; KNAME_3   ="{/Symbol \G}"
       K4=       3.301760 ; KNAME_4   ="K'"
     KEND=       3.962112 ; KNAME_END ="M"
#set xtics (KNAME_INIT KINIT,  KNAME_2 K2, KNAME_3 K3, KNAME_4 K4, KNAME_END KEND) nomirror
 set xtics (KNAME_3 K3, KNAME_4 K4, KNAME_END KEND) nomirror
 set size nosquare 1,0.55 ;
 set xrange[K3:KEND]
 set yrange[-6.:5.]
set ytics -6,1,6 nomirror
set ylabel "Energy (eV)" font "Helvetica,24" offset 1,0
 set arrow from K3,0 to KEND,0 nohead lt 3 lw .1 lc rgb 'black'
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
set key right top samplen 1  spacing 1.5  font "Helvetica, 15"

 set multiplot
 set origin 0,0
 set size nosquare 0.8,0.55
 s=.4

 result_tba='band_structure_TBA_tbfit.dat'

 plot  result_tba          u 1:2 w p  pt 7     ps vari  lc rgb 'blue'    ti "TBA"
