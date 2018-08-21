module phase_factor
contains
  function F_IJ(k,R_IJ)
    use parameters, only : zi
    implicit none
    real*8, intent(in) :: k(3)
    real*8, intent(in) :: R_IJ(3)
    complex*16            F_IJ 

    F_IJ = exp( zi * dot_product( k(1:3), R_IJ(1:3) ) )

    return
  endfunction

  function dxF_IJ(k,R_IJ)
    use parameters, only : zi
    implicit none
    real*8, intent(in) :: k(3)
    real*8, intent(in) :: R_IJ(3)
    complex*16            dxF_IJ

    dxF_IJ = exp( zi * dot_product( k(1:3), R_IJ(1:3) ) ) * R_IJ(1) * zi

    return
  endfunction

  function dyF_IJ(k,R_IJ)
    use parameters, only : zi
    implicit none
    real*8, intent(in) :: k(3)
    real*8, intent(in) :: R_IJ(3)
    complex*16            dyF_IJ

    dyF_IJ = exp( zi * dot_product( k(1:3), R_IJ(1:3) ) ) * R_IJ(2) * zi

    return
  endfunction

  function dzF_IJ(k,R_IJ)
    use parameters, only : zi
    implicit none
    real*8, intent(in) :: k(3)
    real*8, intent(in) :: R_IJ(3)
    complex*16            dzF_IJ

    dzF_IJ = exp( zi * dot_product( k(1:3), R_IJ(1:3) ) ) * R_IJ(3) * zi

    return
  endfunction

endmodule
