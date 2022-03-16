! Module pyfit defined in file pyfit.f90

subroutine f90wrap_incar_py__get__flag_python_module(this, f90wrap_flag_python_module)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_python_module
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_python_module = this_ptr%p%flag_python_module
end subroutine f90wrap_incar_py__get__flag_python_module

subroutine f90wrap_incar_py__set__flag_python_module(this, f90wrap_flag_python_module)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_python_module
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_python_module = f90wrap_flag_python_module
end subroutine f90wrap_incar_py__set__flag_python_module

subroutine f90wrap_incar_py_initialise(this)
    use pyfit, only: incar_py
    implicit none
    
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    type(incar_py_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_incar_py_initialise

subroutine f90wrap_incar_py_finalise(this)
    use pyfit, only: incar_py
    implicit none
    
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    type(incar_py_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_incar_py_finalise

subroutine f90wrap_init(comm)
    use pyfit, only: init
    implicit none
    
    integer, intent(in) :: comm
    call init(comm=comm)
end subroutine f90wrap_init

! End of module pyfit defined in file pyfit.f90

