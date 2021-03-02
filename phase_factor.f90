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

  function dRxF_IJ(k,R_IJ)
    use parameters, only : zi
    implicit none
    real*8, intent(in) :: k(3)
    real*8, intent(in) :: R_IJ(3)
    complex*16            dRxF_IJ

    dRxF_IJ = exp( zi * dot_product( k(1:3), R_IJ(1:3) ) ) * k(1) * zi

    return
  endfunction

  function dRyF_IJ(k,R_IJ)
    use parameters, only : zi
    implicit none
    real*8, intent(in) :: k(3)
    real*8, intent(in) :: R_IJ(3)
    complex*16            dRyF_IJ

    dRyF_IJ = exp( zi * dot_product( k(1:3), R_IJ(1:3) ) ) * k(2) * zi 

    return
  endfunction

  function dRzF_IJ(k,R_IJ)
    use parameters, only : zi
    implicit none
    real*8, intent(in) :: k(3)
    real*8, intent(in) :: R_IJ(3)
    complex*16            dRzF_IJ

    dRzF_IJ = exp( zi * dot_product( k(1:3), R_IJ(1:3) ) ) * k(3) * zi

    return
  endfunction

endmodule
