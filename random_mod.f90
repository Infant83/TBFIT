module random_mod
    use mykind

contains

function random() result(r)
    implicit none
    real(kind=dp) :: r
    call random_number(r)
endfunction

subroutine random_init(iseed)
    implicit none
    integer(kind=sp), intent(in) :: iseed
    integer(kind=sp),allocatable :: seed(:)
    integer(kind=sp)             :: n
    call random_seed(size=n)
    allocate(seed(n))
    seed = iseed
    call random_seed(put=seed)
    deallocate(seed)
endsubroutine

endmodule
