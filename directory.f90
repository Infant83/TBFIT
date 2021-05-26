module directory
    implicit none

contains

    subroutine makedir(path)
        character(len=*), intent(in) :: path
        
        call system( 'mkdir -p '//adjustl(trim(path)) )
        
        return
    endsubroutine

endmodule
