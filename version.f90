module version
  use time
  use mpi_setup

contains 

  subroutine version_stamp(t_start)
    implicit none
    character*132  ver_tag
    real*8         t_start

    ver_tag='# TBFIT 2019. Mar. 21. (build ' // __DATE__// ' ' //__TIME__// ') '

    write(6,*) ver_tag
    write(6,*)" "
    write(6,*)"#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(6,*)"#! Copyright 2018. Hyun-Jung Kim All rights reserved.!"
    write(6,*)"#!           (KIAS, Infant@kias.re.kr)               !"
    write(6,*)"#!           (angpangmokjang@hanmail.net)            !"
    write(6,*)"#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(6,*)" "
    call timestamp('Program start on', t_start)

    return
  endsubroutine

endmodule
