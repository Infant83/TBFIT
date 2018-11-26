module sorting

contains
subroutine get_sort_index_1D(isort_index, target_variable, nvariable, sort_mode)
   implicit none
   integer*4    nvariable
   integer*4    i, j, itemp
   integer*4    isort_index(nvariable)
   real*8       target_variable(nvariable)
   real*8       target_variable_(nvariable)
   real*8       temp
   character*10 sort_mode

   isort_index(1:nvariable) = (/1:nvariable/)
   target_variable_ = target_variable

   select case (sort_mode(1:3))

     case('inc','asc')
       do i = nvariable - 1,1, -1
         do j = 1, i
           if(target_variable_(j+1) .lt. target_variable_(j)) then
             temp = target_variable_(j)
             target_variable_(j) = target_variable_(j+1)
             target_variable_(j+1) = temp
      
             itemp = isort_index(j)
             isort_index(j) = isort_index(j+1)
             isort_index(j+1) = itemp
           endif
         enddo
       enddo
      
      case('dec','des')
       do i = nvariable - 1,1, -1
         do j = 1, i
           if(target_variable_(j+1) .gt. target_variable_(j)) then
             temp = target_variable_(j)
             target_variable_(j) = target_variable_(j+1)
             target_variable_(j+1) = temp
      
             itemp = isort_index(j)
             isort_index(j) = isort_index(j+1)
             isort_index(j+1) = itemp
           endif
         enddo
       enddo

   endselect

   return
endsubroutine

subroutine get_sort_index_1D_int(isort_index, target_variable, nvariable, sort_mode)
   implicit none
   integer*4    nvariable
   integer*4    i, j, itemp
   integer*4    isort_index(nvariable)
   integer*4       target_variable(nvariable)
   integer*4       target_variable_(nvariable)
   integer*4       temp
   character(*), intent(in) :: sort_mode

   isort_index(1:nvariable) = (/1:nvariable/)
   target_variable_ = target_variable

   select case (sort_mode(1:3))

     case('inc','asc')
       do i = nvariable - 1,1, -1
         do j = 1, i
           if(target_variable_(j+1) .lt. target_variable_(j)) then
             temp = target_variable_(j)
             target_variable_(j) = target_variable_(j+1)
             target_variable_(j+1) = temp

             itemp = isort_index(j)
             isort_index(j) = isort_index(j+1)
             isort_index(j+1) = itemp
           endif
         enddo
       enddo

      case('dec','des')
       do i = nvariable - 1,1, -1
         do j = 1, i
           if(target_variable_(j+1) .gt. target_variable_(j)) then
             temp = target_variable_(j)
             target_variable_(j) = target_variable_(j+1)
             target_variable_(j+1) = temp

             itemp = isort_index(j)
             isort_index(j) = isort_index(j+1)
             isort_index(j+1) = itemp
           endif
         enddo
       enddo

   endselect

   return
endsubroutine

subroutine get_sort_variable_1D(target_variable, nvariable, sort_mode)
   implicit none
   integer*4    i
   integer*4    nvariable
   integer*4    isort_index(nvariable)
   real*8       target_variable(nvariable)
   real*8       target_variable_(nvariable)
   character(*), intent(in) :: sort_mode

   target_variable_ = target_variable

   call get_sort_index_1D(isort_index, target_variable_, nvariable, trim(sort_mode))

   do i = 1, nvariable
     target_variable(i) = target_variable_(isort_index(i))
   enddo

   return
endsubroutine

subroutine get_sort_variable_1D_int(target_variable, nvariable, sort_mode)
   implicit none
   integer*4    i
   integer*4    nvariable
   integer*4    isort_index(nvariable)
   integer*4    target_variable(nvariable)
   integer*4    target_variable_(nvariable)
   character(*), intent(in) :: sort_mode

   target_variable_ = target_variable

   call get_sort_index_1D_int(isort_index, target_variable_, nvariable, trim(sort_mode))

   do i = 1, nvariable
     target_variable(i) = target_variable_(isort_index(i))
   enddo

   return
endsubroutine

subroutine get_sort_index(isort_index, target_variable, target_dimension, nvariable)
   implicit none
   integer*4             i, j
   integer*4             idimension
   integer*4             nvariable, target_dimension
   integer*4             isort_index(nvariable)
   real*8, intent(in) :: target_variable(target_dimension, nvariable)
   real*8                target_variable_(target_dimension, nvariable)
   real*8                temp(target_dimension)
   integer*4             itemp
   real*8 t1,t2

   isort_index(1:nvariable) = (/1:nvariable/)
   target_variable_ = target_variable
   call cpu_time(t1)

   do idimension = target_dimension, 1, -1
     if( target_dimension - idimension .ge. 1 ) then
       do i = nvariable - 1,1, -1
         do j = 1, i
          if( all(target_variable_(idimension+1:target_dimension, j+1) .eq. &
                  target_variable_(idimension+1:target_dimension, j)) ) then
             if(target_variable_(idimension, j+1) .gt. target_variable_(idimension, j)) then
!              call get_reversed_value_real(target_variable_(:,j), target_variable_(:,j+1), target_dimension)
!              call get_reversed_value_integer(isort_index(j), isort_index(j+1) ,1)
               temp(:) = target_variable_(:,j)
               target_variable_(:,j) = target_variable_(:,j+1)
               target_variable_(:,j+1) = temp(:)
            
               itemp = isort_index(j)
               isort_index(j) = isort_index(j+1)
               isort_index(j+1) = itemp
             endif
           endif
         enddo
       enddo
     elseif( target_dimension - idimension .eq. 0 ) then
       do i = nvariable - 1,1, -1
         do j = 1, i
           if(target_variable_(idimension, j+1) .gt. target_variable_(idimension, j)) then
!            call get_reversed_value_real(target_variable_(:,j), target_variable_(:,j+1), target_dimension)
!            call get_reversed_value_integer(isort_index(j), isort_index(j+1) ,1)
             temp(:) = target_variable_(:,j)
             target_variable_(:,j) = target_variable_(:,j+1)
             target_variable_(:,j+1) = temp(:)

             itemp = isort_index(j)
             isort_index(j) = isort_index(j+1)
             isort_index(j+1) = itemp
           endif
         enddo
       enddo
     endif
   enddo

  call cpu_time(t2)
  write(6,'(A,F12.4)')"  TIME SORT: ", t2 - t1

   return
endsubroutine

subroutine get_reversed_value_real(A, B, N)
   implicit none
   integer*4    n
   real*8, intent(inout) ::       A(n)
   real*8, intent(inout) ::       B(n)
   real*8       TEMP(n)

   TEMP(:) = A(:)
   A(:)    = B(:)
   B(:)    = TEMP(:)

   return
endsubroutine

subroutine get_reversed_value_integer(A, B, N)
   implicit none
   integer*4    n
   integer*4, intent(inout) ::    A(n)
   integer*4, intent(inout) ::    B(n)
   integer*4    TEMP(n)
   
   TEMP(:) = A(:)
   A(:)    = B(:)
   B(:)    = TEMP(:)

   return
endsubroutine

endmodule
