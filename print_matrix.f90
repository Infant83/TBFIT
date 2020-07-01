#include "alias.inc"
module print_matrix
contains

subroutine print_matrix_i(H,msize_row,msize_col, title, iflag, fmt_)
 use parameters, only : pid_matrix, eta, zi
 implicit none
 integer*4 i,j,iflag,msize_row,msize_col
 integer*4  H(msize_row,msize_col)
 character (len = *) title
 character*80 fname
 character(len = *), optional :: fmt_ ! format, ex) f12.5
 character*80 fmt

 if(present(fmt_)) then
   write(fmt,'(3A)')'(4x,*(',trim(fmt_),'))'
 else
   write(fmt,'(A)')"(4x,*(I10))"
 endif

 ! iflag=0 : print to monitor, 1: print to file with name 'title'
 if (iflag .eq. 0) then
  write(6,'(A,A)')trim(title),'= ['
  do i=1,msize_row
    write(6,fmt,ADVANCE='NO')H(i,1:msize_col)
    if(i .eq. msize_row) then
      write(6,'(A)',ADVANCE='YES')'];'
    elseif(i .lt. msize_row) then
      write(6,'(A)',ADVANCE='YES')';'
    endif
  enddo

 elseif(iflag .eq. 1) then
   write(fname,'(A,A)')trim(title),'.matrix.dat'
   open(pid_matrix,file=fname,status='unknown')
   write(pid_matrix,'(A,A,A)')"# ",trim(title),'='
   do i=1,msize_row
     write(pid_matrix,fmt,ADVANCE='YES')H(i,1:msize_col)
   enddo
   close(pid_matrix)
 endif
return

end subroutine

subroutine print_matrix_c(H_,msize_row,msize_col, title, iflag, fmt_)
 use parameters, only : pid_matrix, eta, zi
 implicit none
 integer*4 i,j,iflag,msize_row,msize_col
 integer*4 ii
 complex*16 H(msize_row,msize_col)
 complex*16 H_(msize_row,msize_col)
 character (len = *) title
 character*80 fname
 character(len = *), optional :: fmt_ ! format, ex) f12.5
 character*80 fmt, fmt_f, fmtp, fmtm
 real*8       a
 if(present(fmt_)) then
   write(fmt,'(5A)')'(*(1x,',trim(fmt_),",'+',",trim(fmt_),",'i'))"
   write(fmtp,'(5A)')'(*(1x,',trim(fmt_),",'+',",trim(fmt_),",'i'))"
   write(fmtm,'(5A)')'(*(1x,',trim(fmt_),",'-',",trim(fmt_),",'i'))"
   write(fmt_f,'(5A)')'(*(1x,',trim(fmt_),",'  ',",trim(fmt_),"    ))"
 else
   write(fmt,'(A)')"(*(2x,F7.3,'+',F7.3,'i'))"
   write(fmt_f,'(A)')"(*(2x,F7.3,' ',F7.3,' '))"
 endif

 
 H = H_

 ! iflag=0 : print to monitor, 1: print to file with name 'title'
 if (iflag .eq. 0) then
  write(6,'(A,A)')trim(title),'= ['
  do i=1,msize_row
   do j = 1, msize_col
     if(abs( real(H(i,j))) .lt. 10d0*eta) then 
       a = aimag(H(i,j))
       H(i,j) = 0d0 + a * zi  ! be careful !!
     endif
     if(abs(aimag(H(i,j))) .lt. 10d0*eta) then 
        a = real(H(i,j))
        H(i,j) =  a+ 0d0 * zi ! be careful !!
     endif
   enddo
  !write(6,fmt,ADVANCE='NO')H(i,1:msize_col)
   do ii = 1, msize_col
     if( aimag(H(i,ii)) .lt. 0 ) then
       write(6,fmtm,ADVANCE='NO') real(H(i,ii)), -aimag(H(i,ii))
     else
       write(6,fmtp,ADVANCE='NO')H(i,ii)
     endif
   enddo
   if(i .eq. msize_row) then
     write(6,'(A)',ADVANCE='YES')'];'
   elseif(i .lt. msize_row) then
     write(6,'(A)',ADVANCE='YES')';'
   endif
  enddo
  write(6,*)''
 elseif(iflag .eq. 1) then
  write(fname,'(A,A,A)')'Matrix.', trim(title),'.dat'
  open(pid_matrix,file=fname,status='unknown')
  write(pid_matrix,'(A,A,A)')"  ",trim(title),'= ['
  do i=1,msize_row
   do j = 1, msize_col
     if(abs( real(H(i,j))) .lt. 10d0*eta) then 
       a = aimag(H(i,j))
       H(i,j) = 0d0 + a * zi  ! be careful !!
     endif
     if(abs(aimag(H(i,j))) .lt. 10d0*eta) then
        a = real(H(i,j))
        H(i,j) =  a+ 0d0 * zi ! be careful !!
     endif
   enddo
  !write(pid_matrix,fmt,ADVANCE='NO' )H(i,1:msize_col)
   do ii = 1, msize_col
     if( aimag(H(i,ii)) .lt. 0 ) then
       write(pid_matrix,fmtm,ADVANCE='NO') real(H(i,ii)), -aimag(H(i,ii))
     else
       write(pid_matrix,fmtp,ADVANCE='NO')H(i,ii)
     endif
   enddo

   if(i .eq. msize_row) then
     write(pid_matrix,'(A)',ADVANCE='YES')'];'
   elseif(i .lt. msize_row) then
     write(pid_matrix,'(A)',ADVANCE='YES')';'
   endif

  enddo
  write(pid_matrix,*)''
  close(pid_matrix)
 endif
999 format( *(2x,F7.3,'+',F7.3,'i') )
998 format( *(2x,F15.8,' + ',F15.8,' i') )
return
end subroutine

subroutine print_matrix_r(H,msize_row,msize_col, title, iflag,fmt_)
 use parameters, only : pid_matrix, eta
 implicit none
 integer*4 i,j,iflag,msize_row,msize_col
 real*8 H(msize_row,msize_col)
 character (len = *) title
 character*80 fname
 character(len = *), optional :: fmt_
 character*80 fmt

 if(present(fmt_)) then
   write(fmt,'(3A)')'(4x,*(',trim(fmt_),'))'
 else
   write(fmt,'(A)')"(4x,*(F10.5))"
 endif

 ! iflag=0 : print to monitor, 1: print to file with name 'title'
 if (iflag .eq. 0) then
  write(6,'(A,A)')trim(title),'= ['
  do i=1,msize_row
   do j = 1, msize_col
     if(abs(H(i,j)) .lt.100d0*eta) H(i,j) = 0d0 ! be careful !!
   enddo
   write(6,fmt,ADVANCE='NO')H(i,1:msize_col)
   if(i .eq. msize_row) then
     write(6,'(A)',ADVANCE='YES')'];'
   elseif(i .lt. msize_row) then
     write(6,'(A)',ADVANCE='YES')';'
   endif

  enddo
 elseif(iflag .eq. 1) then
  write(fname,'(A,A)')trim(title),'.matrix.dat'
  open(pid_matrix,file=fname,status='unknown')
  write(pid_matrix,'(A,A,A)')"# ",trim(title),'='
  do i=1,msize_row
   do j = 1, msize_col
     if(abs(H(i,j)) .lt.100d0*eta) H(i,j) = 0d0 ! be careful !!
   enddo
   write(pid_matrix,fmt,ADVANCE='YES')H(i,1:msize_col)
  enddo
  close(pid_matrix)
 endif
return
end subroutine
endmodule
