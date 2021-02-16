! Module pyfit defined in file tbfitpy_mod.f90

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

subroutine f90wrap_incar_py__get__flag_get_band(this, f90wrap_flag_get_band)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_get_band
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_get_band = this_ptr%p%flag_get_band
end subroutine f90wrap_incar_py__get__flag_get_band

subroutine f90wrap_incar_py__set__flag_get_band(this, f90wrap_flag_get_band)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_get_band
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_get_band = f90wrap_flag_get_band
end subroutine f90wrap_incar_py__set__flag_get_band

subroutine f90wrap_incar_py__get__flag_fit_degeneracy(this, f90wrap_flag_fit_degeneracy)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_fit_degeneracy
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_fit_degeneracy = this_ptr%p%flag_fit_degeneracy
end subroutine f90wrap_incar_py__get__flag_fit_degeneracy

subroutine f90wrap_incar_py__set__flag_fit_degeneracy(this, f90wrap_flag_fit_degeneracy)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_fit_degeneracy
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_fit_degeneracy = f90wrap_flag_fit_degeneracy
end subroutine f90wrap_incar_py__set__flag_fit_degeneracy

subroutine f90wrap_incar_py__get__flag_report_geom(this, f90wrap_flag_report_geom)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_report_geom
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_report_geom = this_ptr%p%flag_report_geom
end subroutine f90wrap_incar_py__get__flag_report_geom

subroutine f90wrap_incar_py__set__flag_report_geom(this, f90wrap_flag_report_geom)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_report_geom
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_report_geom = f90wrap_flag_report_geom
end subroutine f90wrap_incar_py__set__flag_report_geom

subroutine f90wrap_incar_py__get__flag_phase(this, f90wrap_flag_phase)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_phase
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_phase = this_ptr%p%flag_phase
end subroutine f90wrap_incar_py__get__flag_phase

subroutine f90wrap_incar_py__set__flag_phase(this, f90wrap_flag_phase)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_phase
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_phase = f90wrap_flag_phase
end subroutine f90wrap_incar_py__set__flag_phase

subroutine f90wrap_incar_py__get__flag_tbfit(this, f90wrap_flag_tbfit)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_tbfit
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_tbfit = this_ptr%p%flag_tbfit
end subroutine f90wrap_incar_py__get__flag_tbfit

subroutine f90wrap_incar_py__set__flag_tbfit(this, f90wrap_flag_tbfit)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_tbfit
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_tbfit = f90wrap_flag_tbfit
end subroutine f90wrap_incar_py__set__flag_tbfit

subroutine f90wrap_incar_py__get__flag_tbfit_finish(this, f90wrap_flag_tbfit_finish)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_tbfit_finish
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_tbfit_finish = this_ptr%p%flag_tbfit_finish
end subroutine f90wrap_incar_py__get__flag_tbfit_finish

subroutine f90wrap_incar_py__set__flag_tbfit_finish(this, f90wrap_flag_tbfit_finish)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_tbfit_finish
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_tbfit_finish = f90wrap_flag_tbfit_finish
end subroutine f90wrap_incar_py__set__flag_tbfit_finish

subroutine f90wrap_incar_py__get__flag_print_only_target(this, f90wrap_flag_print_only_target)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_print_only_target
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_print_only_target = this_ptr%p%flag_print_only_target
end subroutine f90wrap_incar_py__get__flag_print_only_target

subroutine f90wrap_incar_py__set__flag_print_only_target(this, f90wrap_flag_print_only_target)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_print_only_target
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_print_only_target = f90wrap_flag_print_only_target
end subroutine f90wrap_incar_py__set__flag_print_only_target

subroutine f90wrap_incar_py__get__flag_print_energy_diff(this, f90wrap_flag_print_energy_diff)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_print_energy_diff
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_print_energy_diff = this_ptr%p%flag_print_energy_diff
end subroutine f90wrap_incar_py__get__flag_print_energy_diff

subroutine f90wrap_incar_py__set__flag_print_energy_diff(this, f90wrap_flag_print_energy_diff)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_print_energy_diff
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_print_energy_diff = f90wrap_flag_print_energy_diff
end subroutine f90wrap_incar_py__set__flag_print_energy_diff

subroutine f90wrap_incar_py__get__flag_print_orbital(this, f90wrap_flag_print_orbital)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_print_orbital
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_print_orbital = this_ptr%p%flag_print_orbital
end subroutine f90wrap_incar_py__get__flag_print_orbital

subroutine f90wrap_incar_py__set__flag_print_orbital(this, f90wrap_flag_print_orbital)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_print_orbital
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_print_orbital = f90wrap_flag_print_orbital
end subroutine f90wrap_incar_py__set__flag_print_orbital

subroutine f90wrap_incar_py__get__flag_get_orbital(this, f90wrap_flag_get_orbital)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_get_orbital
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_get_orbital = this_ptr%p%flag_get_orbital
end subroutine f90wrap_incar_py__get__flag_get_orbital

subroutine f90wrap_incar_py__set__flag_get_orbital(this, f90wrap_flag_get_orbital)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_get_orbital
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_get_orbital = f90wrap_flag_get_orbital
end subroutine f90wrap_incar_py__set__flag_get_orbital

subroutine f90wrap_incar_py__get__flag_print_mag(this, f90wrap_flag_print_mag)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_print_mag
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_print_mag = this_ptr%p%flag_print_mag
end subroutine f90wrap_incar_py__get__flag_print_mag

subroutine f90wrap_incar_py__set__flag_print_mag(this, f90wrap_flag_print_mag)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_print_mag
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_print_mag = f90wrap_flag_print_mag
end subroutine f90wrap_incar_py__set__flag_print_mag

subroutine f90wrap_incar_py__get__flag_print_single(this, f90wrap_flag_print_single)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_print_single
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_print_single = this_ptr%p%flag_print_single
end subroutine f90wrap_incar_py__get__flag_print_single

subroutine f90wrap_incar_py__set__flag_print_single(this, f90wrap_flag_print_single)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_print_single
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_print_single = f90wrap_flag_print_single
end subroutine f90wrap_incar_py__set__flag_print_single

subroutine f90wrap_incar_py__get__flag_local_charge(this, f90wrap_flag_local_charge)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_local_charge
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_local_charge = this_ptr%p%flag_local_charge
end subroutine f90wrap_incar_py__get__flag_local_charge

subroutine f90wrap_incar_py__set__flag_local_charge(this, f90wrap_flag_local_charge)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_local_charge
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_local_charge = f90wrap_flag_local_charge
end subroutine f90wrap_incar_py__set__flag_local_charge

subroutine f90wrap_incar_py__get__flag_collinear(this, f90wrap_flag_collinear)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_collinear
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_collinear = this_ptr%p%flag_collinear
end subroutine f90wrap_incar_py__get__flag_collinear

subroutine f90wrap_incar_py__set__flag_collinear(this, f90wrap_flag_collinear)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_collinear
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_collinear = f90wrap_flag_collinear
end subroutine f90wrap_incar_py__set__flag_collinear

subroutine f90wrap_incar_py__get__flag_noncollinear(this, f90wrap_flag_noncollinear)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_noncollinear
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_noncollinear = this_ptr%p%flag_noncollinear
end subroutine f90wrap_incar_py__get__flag_noncollinear

subroutine f90wrap_incar_py__set__flag_noncollinear(this, f90wrap_flag_noncollinear)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_noncollinear
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_noncollinear = f90wrap_flag_noncollinear
end subroutine f90wrap_incar_py__set__flag_noncollinear

subroutine f90wrap_incar_py__get__flag_soc(this, f90wrap_flag_soc)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_soc
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_soc = this_ptr%p%flag_soc
end subroutine f90wrap_incar_py__get__flag_soc

subroutine f90wrap_incar_py__set__flag_soc(this, f90wrap_flag_soc)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_soc
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_soc = f90wrap_flag_soc
end subroutine f90wrap_incar_py__set__flag_soc

subroutine f90wrap_incar_py__get__flag_plus_U(this, f90wrap_flag_plus_U)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_plus_U
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_plus_U = this_ptr%p%flag_plus_U
end subroutine f90wrap_incar_py__get__flag_plus_U

subroutine f90wrap_incar_py__set__flag_plus_U(this, f90wrap_flag_plus_U)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_plus_U
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_plus_U = f90wrap_flag_plus_U
end subroutine f90wrap_incar_py__set__flag_plus_U

subroutine f90wrap_incar_py__get__flag_scissor(this, f90wrap_flag_scissor)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_scissor
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_scissor = this_ptr%p%flag_scissor
end subroutine f90wrap_incar_py__get__flag_scissor

subroutine f90wrap_incar_py__set__flag_scissor(this, f90wrap_flag_scissor)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_scissor
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_scissor = f90wrap_flag_scissor
end subroutine f90wrap_incar_py__set__flag_scissor

subroutine f90wrap_incar_py__get__flag_print_proj(this, f90wrap_flag_print_proj)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_print_proj
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_print_proj = this_ptr%p%flag_print_proj
end subroutine f90wrap_incar_py__get__flag_print_proj

subroutine f90wrap_incar_py__set__flag_print_proj(this, f90wrap_flag_print_proj)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_print_proj
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_print_proj = f90wrap_flag_print_proj
end subroutine f90wrap_incar_py__set__flag_print_proj

subroutine f90wrap_incar_py__get__flag_print_proj_sum(this, f90wrap_flag_print_proj_sum)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_print_proj_sum
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_print_proj_sum = this_ptr%p%flag_print_proj_sum
end subroutine f90wrap_incar_py__get__flag_print_proj_sum

subroutine f90wrap_incar_py__set__flag_print_proj_sum(this, f90wrap_flag_print_proj_sum)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_print_proj_sum
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_print_proj_sum = f90wrap_flag_print_proj_sum
end subroutine f90wrap_incar_py__set__flag_print_proj_sum

subroutine f90wrap_incar_py__get__flag_plot_fit(this, f90wrap_flag_plot_fit)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_plot_fit
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_plot_fit = this_ptr%p%flag_plot_fit
end subroutine f90wrap_incar_py__get__flag_plot_fit

subroutine f90wrap_incar_py__set__flag_plot_fit(this, f90wrap_flag_plot_fit)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_plot_fit
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_plot_fit = f90wrap_flag_plot_fit
end subroutine f90wrap_incar_py__set__flag_plot_fit

subroutine f90wrap_incar_py__get__flag_plot(this, f90wrap_flag_plot)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_plot
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_plot = this_ptr%p%flag_plot
end subroutine f90wrap_incar_py__get__flag_plot

subroutine f90wrap_incar_py__set__flag_plot(this, f90wrap_flag_plot)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_plot
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_plot = f90wrap_flag_plot
end subroutine f90wrap_incar_py__set__flag_plot

subroutine f90wrap_incar_py__get__flag_get_band_order(this, f90wrap_flag_get_band_order)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_get_band_order
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_get_band_order = this_ptr%p%flag_get_band_order
end subroutine f90wrap_incar_py__get__flag_get_band_order

subroutine f90wrap_incar_py__set__flag_get_band_order(this, f90wrap_flag_get_band_order)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_get_band_order
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_get_band_order = f90wrap_flag_get_band_order
end subroutine f90wrap_incar_py__set__flag_get_band_order

subroutine f90wrap_incar_py__get__flag_get_band_order_print_only(this, f90wrap_flag_get_band_order_print_only)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_get_band_order_print_only
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_get_band_order_print_only = this_ptr%p%flag_get_band_order_print_only
end subroutine f90wrap_incar_py__get__flag_get_band_order_print_only

subroutine f90wrap_incar_py__set__flag_get_band_order_print_only(this, f90wrap_flag_get_band_order_print_only)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_get_band_order_print_only
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_get_band_order_print_only = f90wrap_flag_get_band_order_print_only
end subroutine f90wrap_incar_py__set__flag_get_band_order_print_only

subroutine f90wrap_incar_py__get__flag_use_weight(this, f90wrap_flag_use_weight)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_use_weight
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_use_weight = this_ptr%p%flag_use_weight
end subroutine f90wrap_incar_py__get__flag_use_weight

subroutine f90wrap_incar_py__set__flag_use_weight(this, f90wrap_flag_use_weight)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_use_weight
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_use_weight = f90wrap_flag_use_weight
end subroutine f90wrap_incar_py__set__flag_use_weight

subroutine f90wrap_incar_py__get__flag_fit_orbital(this, f90wrap_flag_fit_orbital)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_fit_orbital
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_fit_orbital = this_ptr%p%flag_fit_orbital
end subroutine f90wrap_incar_py__get__flag_fit_orbital

subroutine f90wrap_incar_py__set__flag_fit_orbital(this, f90wrap_flag_fit_orbital)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_fit_orbital
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_fit_orbital = f90wrap_flag_fit_orbital
end subroutine f90wrap_incar_py__set__flag_fit_orbital

subroutine f90wrap_incar_py__get__ptol(this, f90wrap_ptol)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_ptol
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ptol = this_ptr%p%ptol
end subroutine f90wrap_incar_py__get__ptol

subroutine f90wrap_incar_py__set__ptol(this, f90wrap_ptol)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_ptol
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ptol = f90wrap_ptol
end subroutine f90wrap_incar_py__set__ptol

subroutine f90wrap_incar_py__get__ftol(this, f90wrap_ftol)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_ftol
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ftol = this_ptr%p%ftol
end subroutine f90wrap_incar_py__get__ftol

subroutine f90wrap_incar_py__set__ftol(this, f90wrap_ftol)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_ftol
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ftol = f90wrap_ftol
end subroutine f90wrap_incar_py__set__ftol

subroutine f90wrap_incar_py__get__fdiff(this, f90wrap_fdiff)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_fdiff
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_fdiff = this_ptr%p%fdiff
end subroutine f90wrap_incar_py__get__fdiff

subroutine f90wrap_incar_py__set__fdiff(this, f90wrap_fdiff)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_fdiff
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%fdiff = f90wrap_fdiff
end subroutine f90wrap_incar_py__set__fdiff

subroutine f90wrap_incar_py__get__r_scissor(this, f90wrap_r_scissor)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_r_scissor
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_r_scissor = this_ptr%p%r_scissor
end subroutine f90wrap_incar_py__get__r_scissor

subroutine f90wrap_incar_py__set__r_scissor(this, f90wrap_r_scissor)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_r_scissor
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%r_scissor = f90wrap_r_scissor
end subroutine f90wrap_incar_py__set__r_scissor

subroutine f90wrap_incar_py__get__nsystem(this, f90wrap_nsystem)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nsystem
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nsystem = this_ptr%p%nsystem
end subroutine f90wrap_incar_py__get__nsystem

subroutine f90wrap_incar_py__set__nsystem(this, f90wrap_nsystem)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nsystem
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nsystem = f90wrap_nsystem
end subroutine f90wrap_incar_py__set__nsystem

subroutine f90wrap_incar_py__get__lmmax(this, f90wrap_lmmax)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_lmmax
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_lmmax = this_ptr%p%lmmax
end subroutine f90wrap_incar_py__get__lmmax

subroutine f90wrap_incar_py__set__lmmax(this, f90wrap_lmmax)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_lmmax
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%lmmax = f90wrap_lmmax
end subroutine f90wrap_incar_py__set__lmmax

subroutine f90wrap_incar_py__get__miter(this, f90wrap_miter)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_miter
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_miter = this_ptr%p%miter
end subroutine f90wrap_incar_py__get__miter

subroutine f90wrap_incar_py__set__miter(this, f90wrap_miter)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_miter
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%miter = f90wrap_miter
end subroutine f90wrap_incar_py__set__miter

subroutine f90wrap_incar_py__get__mxfit(this, f90wrap_mxfit)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_mxfit
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_mxfit = this_ptr%p%mxfit
end subroutine f90wrap_incar_py__get__mxfit

subroutine f90wrap_incar_py__set__mxfit(this, f90wrap_mxfit)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_mxfit
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%mxfit = f90wrap_mxfit
end subroutine f90wrap_incar_py__set__mxfit

subroutine f90wrap_incar_py__array__nn_max(this, nd, dtype, dshape, dloc)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in) :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%nn_max)
    dloc = loc(this_ptr%p%nn_max)
end subroutine f90wrap_incar_py__array__nn_max

subroutine f90wrap_incar_py__get__ispin(this, f90wrap_ispin)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_ispin
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ispin = this_ptr%p%ispin
end subroutine f90wrap_incar_py__get__ispin

subroutine f90wrap_incar_py__set__ispin(this, f90wrap_ispin)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_ispin
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ispin = f90wrap_ispin
end subroutine f90wrap_incar_py__set__ispin

subroutine f90wrap_incar_py__get__ispinor(this, f90wrap_ispinor)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_ispinor
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ispinor = this_ptr%p%ispinor
end subroutine f90wrap_incar_py__get__ispinor

subroutine f90wrap_incar_py__set__ispinor(this, f90wrap_ispinor)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_ispinor
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ispinor = f90wrap_ispinor
end subroutine f90wrap_incar_py__set__ispinor

subroutine f90wrap_incar_py__get__nspin(this, f90wrap_nspin)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nspin
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nspin = this_ptr%p%nspin
end subroutine f90wrap_incar_py__get__nspin

subroutine f90wrap_incar_py__set__nspin(this, f90wrap_nspin)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nspin
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nspin = f90wrap_nspin
end subroutine f90wrap_incar_py__set__nspin

subroutine f90wrap_incar_py__get__nproj_sum(this, f90wrap_nproj_sum)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nproj_sum
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nproj_sum = this_ptr%p%nproj_sum
end subroutine f90wrap_incar_py__get__nproj_sum

subroutine f90wrap_incar_py__set__nproj_sum(this, f90wrap_nproj_sum)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nproj_sum
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nproj_sum = f90wrap_nproj_sum
end subroutine f90wrap_incar_py__set__nproj_sum

subroutine f90wrap_incar_py__get__iverbose(this, f90wrap_iverbose)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_iverbose
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_iverbose = this_ptr%p%iverbose
end subroutine f90wrap_incar_py__get__iverbose

subroutine f90wrap_incar_py__set__iverbose(this, f90wrap_iverbose)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_iverbose
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%iverbose = f90wrap_iverbose
end subroutine f90wrap_incar_py__set__iverbose

subroutine f90wrap_incar_py__get__i_scissor(this, f90wrap_i_scissor)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_i_scissor
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_i_scissor = this_ptr%p%i_scissor
end subroutine f90wrap_incar_py__get__i_scissor

subroutine f90wrap_incar_py__set__i_scissor(this, f90wrap_i_scissor)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_i_scissor
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%i_scissor = f90wrap_i_scissor
end subroutine f90wrap_incar_py__set__i_scissor

subroutine f90wrap_incar_py__get__ls_type(this, f90wrap_ls_type)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    character(8), intent(out) :: f90wrap_ls_type
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ls_type = this_ptr%p%ls_type
end subroutine f90wrap_incar_py__get__ls_type

subroutine f90wrap_incar_py__set__ls_type(this, f90wrap_ls_type)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    character(8), intent(in) :: f90wrap_ls_type
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ls_type = f90wrap_ls_type
end subroutine f90wrap_incar_py__set__ls_type

subroutine f90wrap_incar_py__get__fnamelog(this, f90wrap_fnamelog)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    character(40), intent(out) :: f90wrap_fnamelog
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_fnamelog = this_ptr%p%fnamelog
end subroutine f90wrap_incar_py__get__fnamelog

subroutine f90wrap_incar_py__set__fnamelog(this, f90wrap_fnamelog)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    character(40), intent(in) :: f90wrap_fnamelog
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%fnamelog = f90wrap_fnamelog
end subroutine f90wrap_incar_py__set__fnamelog

subroutine f90wrap_incar_py__array__proj_atom(this, nd, dtype, dshape, dloc)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in) :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%proj_atom)) then
        dshape(1:2) = shape(this_ptr%p%proj_atom)
        dloc = loc(this_ptr%p%proj_atom)
    else
        dloc = 0
    end if
end subroutine f90wrap_incar_py__array__proj_atom

subroutine f90wrap_incar_py__array__proj_natom(this, nd, dtype, dshape, dloc)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in) :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%proj_natom)) then
        dshape(1:1) = shape(this_ptr%p%proj_natom)
        dloc = loc(this_ptr%p%proj_natom)
    else
        dloc = 0
    end if
end subroutine f90wrap_incar_py__array__proj_natom

subroutine f90wrap_incar_py__array__ifilenm(this, nd, dtype, dshape, dloc)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in) :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%ifilenm)) then
        dshape(1:2) = (/len(this_ptr%p%ifilenm(1)), shape(this_ptr%p%ifilenm)/)
        dloc = loc(this_ptr%p%ifilenm)
    else
        dloc = 0
    end if
end subroutine f90wrap_incar_py__array__ifilenm

subroutine f90wrap_incar_py__array__title(this, nd, dtype, dshape, dloc)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in) :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%title)) then
        dshape(1:2) = (/len(this_ptr%p%title(1)), shape(this_ptr%p%title)/)
        dloc = loc(this_ptr%p%title)
    else
        dloc = 0
    end if
end subroutine f90wrap_incar_py__array__title

subroutine f90wrap_incar_py__get__flag_get_total_energy(this, f90wrap_flag_get_total_energy)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_get_total_energy
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_get_total_energy = this_ptr%p%flag_get_total_energy
end subroutine f90wrap_incar_py__get__flag_get_total_energy

subroutine f90wrap_incar_py__set__flag_get_total_energy(this, f90wrap_flag_get_total_energy)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_get_total_energy
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_get_total_energy = f90wrap_flag_get_total_energy
end subroutine f90wrap_incar_py__set__flag_get_total_energy

subroutine f90wrap_incar_py__get__electronic_temperature(this, f90wrap_electronic_temperature)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_electronic_temperature
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_electronic_temperature = this_ptr%p%electronic_temperature
end subroutine f90wrap_incar_py__get__electronic_temperature

subroutine f90wrap_incar_py__set__electronic_temperature(this, f90wrap_electronic_temperature)
    use pyfit, only: incar_py
    implicit none
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(incar_py_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_electronic_temperature
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%electronic_temperature = f90wrap_electronic_temperature
end subroutine f90wrap_incar_py__set__electronic_temperature

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

subroutine f90wrap_params_py__get__flag_set_param_const(this, f90wrap_flag_set_param_const)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_set_param_const
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_set_param_const = this_ptr%p%flag_set_param_const
end subroutine f90wrap_params_py__get__flag_set_param_const

subroutine f90wrap_params_py__set__flag_set_param_const(this, f90wrap_flag_set_param_const)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_set_param_const
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_set_param_const = f90wrap_flag_set_param_const
end subroutine f90wrap_params_py__set__flag_set_param_const

subroutine f90wrap_params_py__get__flag_pfile_index(this, f90wrap_flag_pfile_index)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_pfile_index
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_pfile_index = this_ptr%p%flag_pfile_index
end subroutine f90wrap_params_py__get__flag_pfile_index

subroutine f90wrap_params_py__set__flag_pfile_index(this, f90wrap_flag_pfile_index)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_pfile_index
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_pfile_index = f90wrap_flag_pfile_index
end subroutine f90wrap_params_py__set__flag_pfile_index

subroutine f90wrap_params_py__get__flag_use_overlap(this, f90wrap_flag_use_overlap)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_use_overlap
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_use_overlap = this_ptr%p%flag_use_overlap
end subroutine f90wrap_params_py__get__flag_use_overlap

subroutine f90wrap_params_py__set__flag_use_overlap(this, f90wrap_flag_use_overlap)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_use_overlap
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_use_overlap = f90wrap_flag_use_overlap
end subroutine f90wrap_params_py__set__flag_use_overlap

subroutine f90wrap_params_py__get__flag_slater_koster(this, f90wrap_flag_slater_koster)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_slater_koster
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_slater_koster = this_ptr%p%flag_slater_koster
end subroutine f90wrap_params_py__get__flag_slater_koster

subroutine f90wrap_params_py__set__flag_slater_koster(this, f90wrap_flag_slater_koster)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_slater_koster
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_slater_koster = f90wrap_flag_slater_koster
end subroutine f90wrap_params_py__set__flag_slater_koster

subroutine f90wrap_params_py__get__flag_nrl_slater_koster(this, f90wrap_flag_nrl_slater_koster)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_nrl_slater_koster
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_nrl_slater_koster = this_ptr%p%flag_nrl_slater_koster
end subroutine f90wrap_params_py__get__flag_nrl_slater_koster

subroutine f90wrap_params_py__set__flag_nrl_slater_koster(this, f90wrap_flag_nrl_slater_koster)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_nrl_slater_koster
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_nrl_slater_koster = f90wrap_flag_nrl_slater_koster
end subroutine f90wrap_params_py__set__flag_nrl_slater_koster

subroutine f90wrap_params_py__get__l_broaden(this, f90wrap_l_broaden)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_l_broaden
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_l_broaden = this_ptr%p%l_broaden
end subroutine f90wrap_params_py__get__l_broaden

subroutine f90wrap_params_py__set__l_broaden(this, f90wrap_l_broaden)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_l_broaden
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%l_broaden = f90wrap_l_broaden
end subroutine f90wrap_params_py__set__l_broaden

subroutine f90wrap_params_py__get__slater_koster_type(this, f90wrap_slater_koster_type)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_slater_koster_type
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_slater_koster_type = this_ptr%p%slater_koster_type
end subroutine f90wrap_params_py__get__slater_koster_type

subroutine f90wrap_params_py__set__slater_koster_type(this, f90wrap_slater_koster_type)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_slater_koster_type
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%slater_koster_type = f90wrap_slater_koster_type
end subroutine f90wrap_params_py__set__slater_koster_type

subroutine f90wrap_params_py__get__nparam(this, f90wrap_nparam)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nparam
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nparam = this_ptr%p%nparam
end subroutine f90wrap_params_py__get__nparam

subroutine f90wrap_params_py__set__nparam(this, f90wrap_nparam)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nparam
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nparam = f90wrap_nparam
end subroutine f90wrap_params_py__set__nparam

subroutine f90wrap_params_py__get__nparam_const(this, f90wrap_nparam_const)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nparam_const
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nparam_const = this_ptr%p%nparam_const
end subroutine f90wrap_params_py__get__nparam_const

subroutine f90wrap_params_py__set__nparam_const(this, f90wrap_nparam_const)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nparam_const
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nparam_const = f90wrap_nparam_const
end subroutine f90wrap_params_py__set__nparam_const

subroutine f90wrap_params_py__get__nparam_nrl(this, f90wrap_nparam_nrl)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nparam_nrl
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nparam_nrl = this_ptr%p%nparam_nrl
end subroutine f90wrap_params_py__get__nparam_nrl

subroutine f90wrap_params_py__set__nparam_nrl(this, f90wrap_nparam_nrl)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nparam_nrl
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nparam_nrl = f90wrap_nparam_nrl
end subroutine f90wrap_params_py__set__nparam_nrl

subroutine f90wrap_params_py__get__nparam_nrl_free(this, f90wrap_nparam_nrl_free)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nparam_nrl_free
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nparam_nrl_free = this_ptr%p%nparam_nrl_free
end subroutine f90wrap_params_py__get__nparam_nrl_free

subroutine f90wrap_params_py__set__nparam_nrl_free(this, f90wrap_nparam_nrl_free)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nparam_nrl_free
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nparam_nrl_free = f90wrap_nparam_nrl_free
end subroutine f90wrap_params_py__set__nparam_nrl_free

subroutine f90wrap_params_py__get__param_nsub_max(this, f90wrap_param_nsub_max)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_param_nsub_max
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_param_nsub_max = this_ptr%p%param_nsub_max
end subroutine f90wrap_params_py__get__param_nsub_max

subroutine f90wrap_params_py__set__param_nsub_max(this, f90wrap_param_nsub_max)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_param_nsub_max
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%param_nsub_max = f90wrap_param_nsub_max
end subroutine f90wrap_params_py__set__param_nsub_max

subroutine f90wrap_params_py__get__nparam_free(this, f90wrap_nparam_free)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nparam_free
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nparam_free = this_ptr%p%nparam_free
end subroutine f90wrap_params_py__get__nparam_free

subroutine f90wrap_params_py__set__nparam_free(this, f90wrap_nparam_free)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nparam_free
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nparam_free = f90wrap_nparam_free
end subroutine f90wrap_params_py__set__nparam_free

subroutine f90wrap_params_py__get__pfilenm(this, f90wrap_pfilenm)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    character(132), intent(out) :: f90wrap_pfilenm
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_pfilenm = this_ptr%p%pfilenm
end subroutine f90wrap_params_py__get__pfilenm

subroutine f90wrap_params_py__set__pfilenm(this, f90wrap_pfilenm)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    character(132), intent(in) :: f90wrap_pfilenm
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%pfilenm = f90wrap_pfilenm
end subroutine f90wrap_params_py__set__pfilenm

subroutine f90wrap_params_py__get__pfileoutnm(this, f90wrap_pfileoutnm)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    character(132), intent(out) :: f90wrap_pfileoutnm
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_pfileoutnm = this_ptr%p%pfileoutnm
end subroutine f90wrap_params_py__get__pfileoutnm

subroutine f90wrap_params_py__set__pfileoutnm(this, f90wrap_pfileoutnm)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in)   :: this(2)
    type(params_py_ptr_type) :: this_ptr
    character(132), intent(in) :: f90wrap_pfileoutnm
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%pfileoutnm = f90wrap_pfileoutnm
end subroutine f90wrap_params_py__set__pfileoutnm

subroutine f90wrap_params_py__array__param(this, nd, dtype, dshape, dloc)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in) :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%param)) then
        dshape(1:1) = shape(this_ptr%p%param)
        dloc = loc(this_ptr%p%param)
    else
        dloc = 0
    end if
end subroutine f90wrap_params_py__array__param

subroutine f90wrap_params_py__array__param_nrl(this, nd, dtype, dshape, dloc)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in) :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%param_nrl)) then
        dshape(1:2) = shape(this_ptr%p%param_nrl)
        dloc = loc(this_ptr%p%param_nrl)
    else
        dloc = 0
    end if
end subroutine f90wrap_params_py__array__param_nrl

subroutine f90wrap_params_py__array__param_const(this, nd, dtype, dshape, dloc)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in) :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%param_const)) then
        dshape(1:2) = shape(this_ptr%p%param_const)
        dloc = loc(this_ptr%p%param_const)
    else
        dloc = 0
    end if
end subroutine f90wrap_params_py__array__param_const

subroutine f90wrap_params_py__array__param_const_nrl(this, nd, dtype, dshape, dloc)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in) :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%param_const_nrl)) then
        dshape(1:3) = shape(this_ptr%p%param_const_nrl)
        dloc = loc(this_ptr%p%param_const_nrl)
    else
        dloc = 0
    end if
end subroutine f90wrap_params_py__array__param_const_nrl

subroutine f90wrap_params_py__array__param_nsub(this, nd, dtype, dshape, dloc)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in) :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%param_nsub)) then
        dshape(1:1) = shape(this_ptr%p%param_nsub)
        dloc = loc(this_ptr%p%param_nsub)
    else
        dloc = 0
    end if
end subroutine f90wrap_params_py__array__param_nsub

subroutine f90wrap_params_py__array__iparam_free(this, nd, dtype, dshape, dloc)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in) :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%iparam_free)) then
        dshape(1:1) = shape(this_ptr%p%iparam_free)
        dloc = loc(this_ptr%p%iparam_free)
    else
        dloc = 0
    end if
end subroutine f90wrap_params_py__array__iparam_free

subroutine f90wrap_params_py__array__iparam_free_nrl(this, nd, dtype, dshape, dloc)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in) :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%iparam_free_nrl)) then
        dshape(1:1) = shape(this_ptr%p%iparam_free_nrl)
        dloc = loc(this_ptr%p%iparam_free_nrl)
    else
        dloc = 0
    end if
end subroutine f90wrap_params_py__array__iparam_free_nrl

subroutine f90wrap_params_py__array__param_name(this, nd, dtype, dshape, dloc)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in) :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%param_name)) then
        dshape(1:2) = (/len(this_ptr%p%param_name(1)), shape(this_ptr%p%param_name)/)
        dloc = loc(this_ptr%p%param_name)
    else
        dloc = 0
    end if
end subroutine f90wrap_params_py__array__param_name

subroutine f90wrap_params_py__array__c_const(this, nd, dtype, dshape, dloc)
    use pyfit, only: params_py
    implicit none
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    integer, intent(in) :: this(2)
    type(params_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%c_const)) then
        dshape(1:3) = (/len(this_ptr%p%c_const(1,1)), shape(this_ptr%p%c_const)/)
        dloc = loc(this_ptr%p%c_const)
    else
        dloc = 0
    end if
end subroutine f90wrap_params_py__array__c_const

subroutine f90wrap_params_py_initialise(this)
    use pyfit, only: params_py
    implicit none
    
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    type(params_py_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_params_py_initialise

subroutine f90wrap_params_py_finalise(this)
    use pyfit, only: params_py
    implicit none
    
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    type(params_py_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_params_py_finalise

subroutine f90wrap_kpoints_py__get__flag_klinemode(this, f90wrap_flag_klinemode)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_klinemode
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_klinemode = this_ptr%p%flag_klinemode
end subroutine f90wrap_kpoints_py__get__flag_klinemode

subroutine f90wrap_kpoints_py__set__flag_klinemode(this, f90wrap_flag_klinemode)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_klinemode
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_klinemode = f90wrap_flag_klinemode
end subroutine f90wrap_kpoints_py__set__flag_klinemode

subroutine f90wrap_kpoints_py__get__flag_kgridmode(this, f90wrap_flag_kgridmode)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_kgridmode
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_kgridmode = this_ptr%p%flag_kgridmode
end subroutine f90wrap_kpoints_py__get__flag_kgridmode

subroutine f90wrap_kpoints_py__set__flag_kgridmode(this, f90wrap_flag_kgridmode)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_kgridmode
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_kgridmode = f90wrap_flag_kgridmode
end subroutine f90wrap_kpoints_py__set__flag_kgridmode

subroutine f90wrap_kpoints_py__get__flag_gamma(this, f90wrap_flag_gamma)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_gamma
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_gamma = this_ptr%p%flag_gamma
end subroutine f90wrap_kpoints_py__get__flag_gamma

subroutine f90wrap_kpoints_py__set__flag_gamma(this, f90wrap_flag_gamma)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_gamma
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_gamma = f90wrap_flag_gamma
end subroutine f90wrap_kpoints_py__set__flag_gamma

subroutine f90wrap_kpoints_py__get__flag_reciprocal(this, f90wrap_flag_reciprocal)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_reciprocal
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_reciprocal = this_ptr%p%flag_reciprocal
end subroutine f90wrap_kpoints_py__get__flag_reciprocal

subroutine f90wrap_kpoints_py__set__flag_reciprocal(this, f90wrap_flag_reciprocal)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_reciprocal
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_reciprocal = f90wrap_flag_reciprocal
end subroutine f90wrap_kpoints_py__set__flag_reciprocal

subroutine f90wrap_kpoints_py__get__flag_cartesianK(this, f90wrap_flag_cartesianK)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_cartesianK
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_cartesianK = this_ptr%p%flag_cartesianK
end subroutine f90wrap_kpoints_py__get__flag_cartesianK

subroutine f90wrap_kpoints_py__set__flag_cartesianK(this, f90wrap_flag_cartesianK)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_cartesianK
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_cartesianK = f90wrap_flag_cartesianK
end subroutine f90wrap_kpoints_py__set__flag_cartesianK

subroutine f90wrap_kpoints_py__array__k_shift(this, nd, dtype, dshape, dloc)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in) :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%k_shift)
    dloc = loc(this_ptr%p%k_shift)
end subroutine f90wrap_kpoints_py__array__k_shift

subroutine f90wrap_kpoints_py__get__mysystem(this, f90wrap_mysystem)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_mysystem
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_mysystem = this_ptr%p%mysystem
end subroutine f90wrap_kpoints_py__get__mysystem

subroutine f90wrap_kpoints_py__set__mysystem(this, f90wrap_mysystem)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_mysystem
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%mysystem = f90wrap_mysystem
end subroutine f90wrap_kpoints_py__set__mysystem

subroutine f90wrap_kpoints_py__get__nkpoint(this, f90wrap_nkpoint)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nkpoint
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nkpoint = this_ptr%p%nkpoint
end subroutine f90wrap_kpoints_py__get__nkpoint

subroutine f90wrap_kpoints_py__set__nkpoint(this, f90wrap_nkpoint)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nkpoint
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nkpoint = f90wrap_nkpoint
end subroutine f90wrap_kpoints_py__set__nkpoint

subroutine f90wrap_kpoints_py__get__nline(this, f90wrap_nline)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nline
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nline = this_ptr%p%nline
end subroutine f90wrap_kpoints_py__get__nline

subroutine f90wrap_kpoints_py__set__nline(this, f90wrap_nline)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nline
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nline = f90wrap_nline
end subroutine f90wrap_kpoints_py__set__nline

subroutine f90wrap_kpoints_py__get__n_ndiv(this, f90wrap_n_ndiv)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_n_ndiv
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_n_ndiv = this_ptr%p%n_ndiv
end subroutine f90wrap_kpoints_py__get__n_ndiv

subroutine f90wrap_kpoints_py__set__n_ndiv(this, f90wrap_n_ndiv)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_n_ndiv
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%n_ndiv = f90wrap_n_ndiv
end subroutine f90wrap_kpoints_py__set__n_ndiv

subroutine f90wrap_kpoints_py__get__idiv_mode(this, f90wrap_idiv_mode)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_idiv_mode
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_idiv_mode = this_ptr%p%idiv_mode
end subroutine f90wrap_kpoints_py__get__idiv_mode

subroutine f90wrap_kpoints_py__set__idiv_mode(this, f90wrap_idiv_mode)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_idiv_mode
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%idiv_mode = f90wrap_idiv_mode
end subroutine f90wrap_kpoints_py__set__idiv_mode

subroutine f90wrap_kpoints_py__get__kreduce(this, f90wrap_kreduce)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_kreduce
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_kreduce = this_ptr%p%kreduce
end subroutine f90wrap_kpoints_py__get__kreduce

subroutine f90wrap_kpoints_py__set__kreduce(this, f90wrap_kreduce)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_kreduce
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%kreduce = f90wrap_kreduce
end subroutine f90wrap_kpoints_py__set__kreduce

subroutine f90wrap_kpoints_py__get__kfilenm(this, f90wrap_kfilenm)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    character(132), intent(out) :: f90wrap_kfilenm
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_kfilenm = this_ptr%p%kfilenm
end subroutine f90wrap_kpoints_py__get__kfilenm

subroutine f90wrap_kpoints_py__set__kfilenm(this, f90wrap_kfilenm)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    character(132), intent(in) :: f90wrap_kfilenm
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%kfilenm = f90wrap_kfilenm
end subroutine f90wrap_kpoints_py__set__kfilenm

subroutine f90wrap_kpoints_py__get__ribbon_kfilenm(this, f90wrap_ribbon_kfilenm)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    character(132), intent(out) :: f90wrap_ribbon_kfilenm
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ribbon_kfilenm = this_ptr%p%ribbon_kfilenm
end subroutine f90wrap_kpoints_py__get__ribbon_kfilenm

subroutine f90wrap_kpoints_py__set__ribbon_kfilenm(this, f90wrap_ribbon_kfilenm)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    character(132), intent(in) :: f90wrap_ribbon_kfilenm
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ribbon_kfilenm = f90wrap_ribbon_kfilenm
end subroutine f90wrap_kpoints_py__set__ribbon_kfilenm

subroutine f90wrap_kpoints_py__get__kline_type(this, f90wrap_kline_type)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    character(1024), intent(out) :: f90wrap_kline_type
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_kline_type = this_ptr%p%kline_type
end subroutine f90wrap_kpoints_py__get__kline_type

subroutine f90wrap_kpoints_py__set__kline_type(this, f90wrap_kline_type)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    character(1024), intent(in) :: f90wrap_kline_type
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%kline_type = f90wrap_kline_type
end subroutine f90wrap_kpoints_py__set__kline_type

subroutine f90wrap_kpoints_py__get__kunit(this, f90wrap_kunit)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    character(1024), intent(out) :: f90wrap_kunit
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_kunit = this_ptr%p%kunit
end subroutine f90wrap_kpoints_py__get__kunit

subroutine f90wrap_kpoints_py__set__kunit(this, f90wrap_kunit)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in)   :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    character(1024), intent(in) :: f90wrap_kunit
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%kunit = f90wrap_kunit
end subroutine f90wrap_kpoints_py__set__kunit

subroutine f90wrap_kpoints_py__array__kpoint(this, nd, dtype, dshape, dloc)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in) :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%kpoint)) then
        dshape(1:2) = shape(this_ptr%p%kpoint)
        dloc = loc(this_ptr%p%kpoint)
    else
        dloc = 0
    end if
end subroutine f90wrap_kpoints_py__array__kpoint

subroutine f90wrap_kpoints_py__array__kline(this, nd, dtype, dshape, dloc)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in) :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%kline)) then
        dshape(1:2) = shape(this_ptr%p%kline)
        dloc = loc(this_ptr%p%kline)
    else
        dloc = 0
    end if
end subroutine f90wrap_kpoints_py__array__kline

subroutine f90wrap_kpoints_py__array__kpoint_reci(this, nd, dtype, dshape, dloc)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in) :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%kpoint_reci)) then
        dshape(1:2) = shape(this_ptr%p%kpoint_reci)
        dloc = loc(this_ptr%p%kpoint_reci)
    else
        dloc = 0
    end if
end subroutine f90wrap_kpoints_py__array__kpoint_reci

subroutine f90wrap_kpoints_py__array__ndiv(this, nd, dtype, dshape, dloc)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in) :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%ndiv)) then
        dshape(1:1) = shape(this_ptr%p%ndiv)
        dloc = loc(this_ptr%p%ndiv)
    else
        dloc = 0
    end if
end subroutine f90wrap_kpoints_py__array__ndiv

subroutine f90wrap_kpoints_py__array__k_name(this, nd, dtype, dshape, dloc)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in) :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%k_name)) then
        dshape(1:2) = (/len(this_ptr%p%k_name(1)), shape(this_ptr%p%k_name)/)
        dloc = loc(this_ptr%p%k_name)
    else
        dloc = 0
    end if
end subroutine f90wrap_kpoints_py__array__k_name

subroutine f90wrap_kpoints_py__array__k_name2(this, nd, dtype, dshape, dloc)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in) :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%k_name2)) then
        dshape(1:2) = (/len(this_ptr%p%k_name2(1)), shape(this_ptr%p%k_name2)/)
        dloc = loc(this_ptr%p%k_name2)
    else
        dloc = 0
    end if
end subroutine f90wrap_kpoints_py__array__k_name2

subroutine f90wrap_kpoints_py__array__k_name_index(this, nd, dtype, dshape, dloc)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in) :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%k_name_index)) then
        dshape(1:1) = shape(this_ptr%p%k_name_index)
        dloc = loc(this_ptr%p%k_name_index)
    else
        dloc = 0
    end if
end subroutine f90wrap_kpoints_py__array__k_name_index

subroutine f90wrap_kpoints_py__array__kdist(this, nd, dtype, dshape, dloc)
    use pyfit, only: kpoints_py
    implicit none
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    integer, intent(in) :: this(2)
    type(kpoints_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%kdist)) then
        dshape(1:1) = shape(this_ptr%p%kdist)
        dloc = loc(this_ptr%p%kdist)
    else
        dloc = 0
    end if
end subroutine f90wrap_kpoints_py__array__kdist

subroutine f90wrap_kpoints_py_initialise(this)
    use pyfit, only: kpoints_py
    implicit none
    
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    type(kpoints_py_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_kpoints_py_initialise

subroutine f90wrap_kpoints_py_finalise(this)
    use pyfit, only: kpoints_py
    implicit none
    
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    type(kpoints_py_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_kpoints_py_finalise

subroutine f90wrap_weight_py__get__flag_weight_default(this, f90wrap_flag_weight_default)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_weight_default
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_weight_default = this_ptr%p%flag_weight_default
end subroutine f90wrap_weight_py__get__flag_weight_default

subroutine f90wrap_weight_py__set__flag_weight_default(this, f90wrap_flag_weight_default)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_weight_default
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_weight_default = f90wrap_flag_weight_default
end subroutine f90wrap_weight_py__set__flag_weight_default

subroutine f90wrap_weight_py__get__flag_weight_orb(this, f90wrap_flag_weight_orb)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_weight_orb
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_weight_orb = this_ptr%p%flag_weight_orb
end subroutine f90wrap_weight_py__get__flag_weight_orb

subroutine f90wrap_weight_py__set__flag_weight_orb(this, f90wrap_flag_weight_orb)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_weight_orb
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_weight_orb = f90wrap_flag_weight_orb
end subroutine f90wrap_weight_py__set__flag_weight_orb

subroutine f90wrap_weight_py__get__efile_ef(this, f90wrap_efile_ef)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_efile_ef
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_efile_ef = this_ptr%p%efile_ef
end subroutine f90wrap_weight_py__get__efile_ef

subroutine f90wrap_weight_py__set__efile_ef(this, f90wrap_efile_ef)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_efile_ef
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%efile_ef = f90wrap_efile_ef
end subroutine f90wrap_weight_py__set__efile_ef

subroutine f90wrap_weight_py__get__mysystem(this, f90wrap_mysystem)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_mysystem
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_mysystem = this_ptr%p%mysystem
end subroutine f90wrap_weight_py__get__mysystem

subroutine f90wrap_weight_py__set__mysystem(this, f90wrap_mysystem)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_mysystem
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%mysystem = f90wrap_mysystem
end subroutine f90wrap_weight_py__set__mysystem

subroutine f90wrap_weight_py__get__itarget_e_start(this, f90wrap_itarget_e_start)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_itarget_e_start
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_itarget_e_start = this_ptr%p%itarget_e_start
end subroutine f90wrap_weight_py__get__itarget_e_start

subroutine f90wrap_weight_py__set__itarget_e_start(this, f90wrap_itarget_e_start)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_itarget_e_start
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%itarget_e_start = f90wrap_itarget_e_start
end subroutine f90wrap_weight_py__set__itarget_e_start

subroutine f90wrap_weight_py__get__read_energy_column_index(this, f90wrap_read_energy_column_index)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_read_energy_column_index
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_read_energy_column_index = this_ptr%p%read_energy_column_index
end subroutine f90wrap_weight_py__get__read_energy_column_index

subroutine f90wrap_weight_py__set__read_energy_column_index(this, f90wrap_read_energy_column_index)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_read_energy_column_index
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%read_energy_column_index = f90wrap_read_energy_column_index
end subroutine f90wrap_weight_py__set__read_energy_column_index

subroutine f90wrap_weight_py__get__read_energy_column_index_dn(this, f90wrap_read_energy_column_index_dn)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_read_energy_column_index_dn
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_read_energy_column_index_dn = this_ptr%p%read_energy_column_index_dn
end subroutine f90wrap_weight_py__get__read_energy_column_index_dn

subroutine f90wrap_weight_py__set__read_energy_column_index_dn(this, f90wrap_read_energy_column_index_dn)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_read_energy_column_index_dn
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%read_energy_column_index_dn = f90wrap_read_energy_column_index_dn
end subroutine f90wrap_weight_py__set__read_energy_column_index_dn

subroutine f90wrap_weight_py__get__ie_cutoff(this, f90wrap_ie_cutoff)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_ie_cutoff
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ie_cutoff = this_ptr%p%ie_cutoff
end subroutine f90wrap_weight_py__get__ie_cutoff

subroutine f90wrap_weight_py__set__ie_cutoff(this, f90wrap_ie_cutoff)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_ie_cutoff
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ie_cutoff = f90wrap_ie_cutoff
end subroutine f90wrap_weight_py__set__ie_cutoff

subroutine f90wrap_weight_py__get__efilenmu(this, f90wrap_efilenmu)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    character(132), intent(out) :: f90wrap_efilenmu
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_efilenmu = this_ptr%p%efilenmu
end subroutine f90wrap_weight_py__get__efilenmu

subroutine f90wrap_weight_py__set__efilenmu(this, f90wrap_efilenmu)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    character(132), intent(in) :: f90wrap_efilenmu
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%efilenmu = f90wrap_efilenmu
end subroutine f90wrap_weight_py__set__efilenmu

subroutine f90wrap_weight_py__get__efilenmd(this, f90wrap_efilenmd)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    character(132), intent(out) :: f90wrap_efilenmd
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_efilenmd = this_ptr%p%efilenmd
end subroutine f90wrap_weight_py__get__efilenmd

subroutine f90wrap_weight_py__set__efilenmd(this, f90wrap_efilenmd)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    character(132), intent(in) :: f90wrap_efilenmd
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%efilenmd = f90wrap_efilenmd
end subroutine f90wrap_weight_py__set__efilenmd

subroutine f90wrap_weight_py__get__efile_type(this, f90wrap_efile_type)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    character(1024), intent(out) :: f90wrap_efile_type
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_efile_type = this_ptr%p%efile_type
end subroutine f90wrap_weight_py__get__efile_type

subroutine f90wrap_weight_py__set__efile_type(this, f90wrap_efile_type)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in)   :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    character(1024), intent(in) :: f90wrap_efile_type
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%efile_type = f90wrap_efile_type
end subroutine f90wrap_weight_py__set__efile_type

subroutine f90wrap_weight_py__array__WT(this, nd, dtype, dshape, dloc)
    use pyfit, only: weight_py
    implicit none
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    integer, intent(in) :: this(2)
    type(weight_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%WT)) then
        dshape(1:2) = shape(this_ptr%p%WT)
        dloc = loc(this_ptr%p%WT)
    else
        dloc = 0
    end if
end subroutine f90wrap_weight_py__array__WT

subroutine f90wrap_weight_py_initialise(this)
    use pyfit, only: weight_py
    implicit none
    
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    type(weight_py_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_weight_py_initialise

subroutine f90wrap_weight_py_finalise(this)
    use pyfit, only: weight_py
    implicit none
    
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    type(weight_py_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_weight_py_finalise

subroutine f90wrap_poscar_py__get__flag_selective(this, f90wrap_flag_selective)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_selective
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_selective = this_ptr%p%flag_selective
end subroutine f90wrap_poscar_py__get__flag_selective

subroutine f90wrap_poscar_py__set__flag_selective(this, f90wrap_flag_selective)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_selective
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_selective = f90wrap_flag_selective
end subroutine f90wrap_poscar_py__set__flag_selective

subroutine f90wrap_poscar_py__get__flag_direct(this, f90wrap_flag_direct)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_direct
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_direct = this_ptr%p%flag_direct
end subroutine f90wrap_poscar_py__get__flag_direct

subroutine f90wrap_poscar_py__set__flag_direct(this, f90wrap_flag_direct)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_direct
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_direct = f90wrap_flag_direct
end subroutine f90wrap_poscar_py__set__flag_direct

subroutine f90wrap_poscar_py__get__flag_cartesian(this, f90wrap_flag_cartesian)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_cartesian
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_cartesian = this_ptr%p%flag_cartesian
end subroutine f90wrap_poscar_py__get__flag_cartesian

subroutine f90wrap_poscar_py__set__flag_cartesian(this, f90wrap_flag_cartesian)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_cartesian
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_cartesian = f90wrap_flag_cartesian
end subroutine f90wrap_poscar_py__set__flag_cartesian

subroutine f90wrap_poscar_py__get__mysystem(this, f90wrap_mysystem)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_mysystem
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_mysystem = this_ptr%p%mysystem
end subroutine f90wrap_poscar_py__get__mysystem

subroutine f90wrap_poscar_py__set__mysystem(this, f90wrap_mysystem)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_mysystem
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%mysystem = f90wrap_mysystem
end subroutine f90wrap_poscar_py__set__mysystem

subroutine f90wrap_poscar_py__get__neig(this, f90wrap_neig)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_neig
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_neig = this_ptr%p%neig
end subroutine f90wrap_poscar_py__get__neig

subroutine f90wrap_poscar_py__set__neig(this, f90wrap_neig)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_neig
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%neig = f90wrap_neig
end subroutine f90wrap_poscar_py__set__neig

subroutine f90wrap_poscar_py__get__neig_total(this, f90wrap_neig_total)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_neig_total
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_neig_total = this_ptr%p%neig_total
end subroutine f90wrap_poscar_py__get__neig_total

subroutine f90wrap_poscar_py__set__neig_total(this, f90wrap_neig_total)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_neig_total
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%neig_total = f90wrap_neig_total
end subroutine f90wrap_poscar_py__set__neig_total

subroutine f90wrap_poscar_py__get__neig_target(this, f90wrap_neig_target)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_neig_target
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_neig_target = this_ptr%p%neig_target
end subroutine f90wrap_poscar_py__get__neig_target

subroutine f90wrap_poscar_py__set__neig_target(this, f90wrap_neig_target)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_neig_target
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%neig_target = f90wrap_neig_target
end subroutine f90wrap_poscar_py__set__neig_target

subroutine f90wrap_poscar_py__get__nbasis(this, f90wrap_nbasis)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nbasis
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nbasis = this_ptr%p%nbasis
end subroutine f90wrap_poscar_py__get__nbasis

subroutine f90wrap_poscar_py__set__nbasis(this, f90wrap_nbasis)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nbasis
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nbasis = f90wrap_nbasis
end subroutine f90wrap_poscar_py__set__nbasis

subroutine f90wrap_poscar_py__get__neig_eff(this, f90wrap_neig_eff)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_neig_eff
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_neig_eff = this_ptr%p%neig_eff
end subroutine f90wrap_poscar_py__get__neig_eff

subroutine f90wrap_poscar_py__set__neig_eff(this, f90wrap_neig_eff)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_neig_eff
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%neig_eff = f90wrap_neig_eff
end subroutine f90wrap_poscar_py__set__neig_eff

subroutine f90wrap_poscar_py__get__init_erange(this, f90wrap_init_erange)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_init_erange
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_init_erange = this_ptr%p%init_erange
end subroutine f90wrap_poscar_py__get__init_erange

subroutine f90wrap_poscar_py__set__init_erange(this, f90wrap_init_erange)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_init_erange
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%init_erange = f90wrap_init_erange
end subroutine f90wrap_poscar_py__set__init_erange

subroutine f90wrap_poscar_py__get__fina_erange(this, f90wrap_fina_erange)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_fina_erange
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_fina_erange = this_ptr%p%fina_erange
end subroutine f90wrap_poscar_py__get__fina_erange

subroutine f90wrap_poscar_py__set__fina_erange(this, f90wrap_fina_erange)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_fina_erange
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%fina_erange = f90wrap_fina_erange
end subroutine f90wrap_poscar_py__set__fina_erange

subroutine f90wrap_poscar_py__get__nband(this, f90wrap_nband)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nband
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nband = this_ptr%p%nband
end subroutine f90wrap_poscar_py__get__nband

subroutine f90wrap_poscar_py__set__nband(this, f90wrap_nband)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nband
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nband = f90wrap_nband
end subroutine f90wrap_poscar_py__set__nband

subroutine f90wrap_poscar_py__get__n_spec(this, f90wrap_n_spec)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_n_spec
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_n_spec = this_ptr%p%n_spec
end subroutine f90wrap_poscar_py__get__n_spec

subroutine f90wrap_poscar_py__set__n_spec(this, f90wrap_n_spec)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_n_spec
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%n_spec = f90wrap_n_spec
end subroutine f90wrap_poscar_py__set__n_spec

subroutine f90wrap_poscar_py__get__n_atom(this, f90wrap_n_atom)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_n_atom
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_n_atom = this_ptr%p%n_atom
end subroutine f90wrap_poscar_py__get__n_atom

subroutine f90wrap_poscar_py__set__n_atom(this, f90wrap_n_atom)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_n_atom
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%n_atom = f90wrap_n_atom
end subroutine f90wrap_poscar_py__set__n_atom

subroutine f90wrap_poscar_py__get__max_orb(this, f90wrap_max_orb)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_max_orb
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_max_orb = this_ptr%p%max_orb
end subroutine f90wrap_poscar_py__get__max_orb

subroutine f90wrap_poscar_py__set__max_orb(this, f90wrap_max_orb)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_max_orb
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%max_orb = f90wrap_max_orb
end subroutine f90wrap_poscar_py__set__max_orb

subroutine f90wrap_poscar_py__get__title(this, f90wrap_title)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    character(132), intent(out) :: f90wrap_title
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_title = this_ptr%p%title
end subroutine f90wrap_poscar_py__get__title

subroutine f90wrap_poscar_py__set__title(this, f90wrap_title)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    character(132), intent(in) :: f90wrap_title
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%title = f90wrap_title
end subroutine f90wrap_poscar_py__set__title

subroutine f90wrap_poscar_py__get__gfilenm(this, f90wrap_gfilenm)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    character(132), intent(out) :: f90wrap_gfilenm
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_gfilenm = this_ptr%p%gfilenm
end subroutine f90wrap_poscar_py__get__gfilenm

subroutine f90wrap_poscar_py__set__gfilenm(this, f90wrap_gfilenm)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    character(132), intent(in) :: f90wrap_gfilenm
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%gfilenm = f90wrap_gfilenm
end subroutine f90wrap_poscar_py__set__gfilenm

subroutine f90wrap_poscar_py__get__system_name(this, f90wrap_system_name)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    character(40), intent(out) :: f90wrap_system_name
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_system_name = this_ptr%p%system_name
end subroutine f90wrap_poscar_py__get__system_name

subroutine f90wrap_poscar_py__set__system_name(this, f90wrap_system_name)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    character(40), intent(in) :: f90wrap_system_name
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%system_name = f90wrap_system_name
end subroutine f90wrap_poscar_py__set__system_name

subroutine f90wrap_poscar_py__get__a_scale(this, f90wrap_a_scale)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_a_scale
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_a_scale = this_ptr%p%a_scale
end subroutine f90wrap_poscar_py__get__a_scale

subroutine f90wrap_poscar_py__set__a_scale(this, f90wrap_a_scale)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in)   :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_a_scale
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%a_scale = f90wrap_a_scale
end subroutine f90wrap_poscar_py__set__a_scale

subroutine f90wrap_poscar_py__array__a_latt(this, nd, dtype, dshape, dloc)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in) :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%a_latt)
    dloc = loc(this_ptr%p%a_latt)
end subroutine f90wrap_poscar_py__array__a_latt

subroutine f90wrap_poscar_py__array__b_latt(this, nd, dtype, dshape, dloc)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in) :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%b_latt)
    dloc = loc(this_ptr%p%b_latt)
end subroutine f90wrap_poscar_py__array__b_latt

subroutine f90wrap_poscar_py__array__nelect(this, nd, dtype, dshape, dloc)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in) :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%nelect)) then
        dshape(1:1) = shape(this_ptr%p%nelect)
        dloc = loc(this_ptr%p%nelect)
    else
        dloc = 0
    end if
end subroutine f90wrap_poscar_py__array__nelect

subroutine f90wrap_poscar_py__array__n_orbital(this, nd, dtype, dshape, dloc)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in) :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%n_orbital)) then
        dshape(1:1) = shape(this_ptr%p%n_orbital)
        dloc = loc(this_ptr%p%n_orbital)
    else
        dloc = 0
    end if
end subroutine f90wrap_poscar_py__array__n_orbital

subroutine f90wrap_poscar_py__array__orb_index(this, nd, dtype, dshape, dloc)
    use pyfit, only: poscar_py
    implicit none
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    integer, intent(in) :: this(2)
    type(poscar_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%orb_index)) then
        dshape(1:1) = shape(this_ptr%p%orb_index)
        dloc = loc(this_ptr%p%orb_index)
    else
        dloc = 0
    end if
end subroutine f90wrap_poscar_py__array__orb_index

subroutine f90wrap_poscar_py_initialise(this)
    use pyfit, only: poscar_py
    implicit none
    
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    type(poscar_py_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_poscar_py_initialise

subroutine f90wrap_poscar_py_finalise(this)
    use pyfit, only: poscar_py
    implicit none
    
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    type(poscar_py_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_poscar_py_finalise

subroutine f90wrap_hopping_py__get__flag_efield(this, f90wrap_flag_efield)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in)   :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_efield
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_efield = this_ptr%p%flag_efield
end subroutine f90wrap_hopping_py__get__flag_efield

subroutine f90wrap_hopping_py__set__flag_efield(this, f90wrap_flag_efield)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in)   :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_efield
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_efield = f90wrap_flag_efield
end subroutine f90wrap_hopping_py__set__flag_efield

subroutine f90wrap_hopping_py__get__flag_efield_frac(this, f90wrap_flag_efield_frac)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in)   :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_efield_frac
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_efield_frac = this_ptr%p%flag_efield_frac
end subroutine f90wrap_hopping_py__get__flag_efield_frac

subroutine f90wrap_hopping_py__set__flag_efield_frac(this, f90wrap_flag_efield_frac)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in)   :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_efield_frac
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_efield_frac = f90wrap_flag_efield_frac
end subroutine f90wrap_hopping_py__set__flag_efield_frac

subroutine f90wrap_hopping_py__get__flag_efield_cart(this, f90wrap_flag_efield_cart)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in)   :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_flag_efield_cart
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_flag_efield_cart = this_ptr%p%flag_efield_cart
end subroutine f90wrap_hopping_py__get__flag_efield_cart

subroutine f90wrap_hopping_py__set__flag_efield_cart(this, f90wrap_flag_efield_cart)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in)   :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_flag_efield_cart
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%flag_efield_cart = f90wrap_flag_efield_cart
end subroutine f90wrap_hopping_py__set__flag_efield_cart

subroutine f90wrap_hopping_py__get__onsite_tolerance(this, f90wrap_onsite_tolerance)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in)   :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_onsite_tolerance
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_onsite_tolerance = this_ptr%p%onsite_tolerance
end subroutine f90wrap_hopping_py__get__onsite_tolerance

subroutine f90wrap_hopping_py__set__onsite_tolerance(this, f90wrap_onsite_tolerance)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in)   :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_onsite_tolerance
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%onsite_tolerance = f90wrap_onsite_tolerance
end subroutine f90wrap_hopping_py__set__onsite_tolerance

subroutine f90wrap_hopping_py__array__efield(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%efield)
    dloc = loc(this_ptr%p%efield)
end subroutine f90wrap_hopping_py__array__efield

subroutine f90wrap_hopping_py__array__efield_origin(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%efield_origin)
    dloc = loc(this_ptr%p%efield_origin)
end subroutine f90wrap_hopping_py__array__efield_origin

subroutine f90wrap_hopping_py__array__efield_origin_cart(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%efield_origin_cart)
    dloc = loc(this_ptr%p%efield_origin_cart)
end subroutine f90wrap_hopping_py__array__efield_origin_cart

subroutine f90wrap_hopping_py__get__mysystem(this, f90wrap_mysystem)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in)   :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_mysystem
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_mysystem = this_ptr%p%mysystem
end subroutine f90wrap_hopping_py__get__mysystem

subroutine f90wrap_hopping_py__set__mysystem(this, f90wrap_mysystem)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in)   :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_mysystem
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%mysystem = f90wrap_mysystem
end subroutine f90wrap_hopping_py__set__mysystem

subroutine f90wrap_hopping_py__get__n_neighbor(this, f90wrap_n_neighbor)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in)   :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_n_neighbor
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_n_neighbor = this_ptr%p%n_neighbor
end subroutine f90wrap_hopping_py__get__n_neighbor

subroutine f90wrap_hopping_py__set__n_neighbor(this, f90wrap_n_neighbor)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in)   :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_n_neighbor
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%n_neighbor = f90wrap_n_neighbor
end subroutine f90wrap_hopping_py__set__n_neighbor

subroutine f90wrap_hopping_py__get__max_nn_pair(this, f90wrap_max_nn_pair)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in)   :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_max_nn_pair
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_max_nn_pair = this_ptr%p%max_nn_pair
end subroutine f90wrap_hopping_py__get__max_nn_pair

subroutine f90wrap_hopping_py__set__max_nn_pair(this, f90wrap_max_nn_pair)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in)   :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_max_nn_pair
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%max_nn_pair = f90wrap_max_nn_pair
end subroutine f90wrap_hopping_py__set__max_nn_pair

subroutine f90wrap_hopping_py__get__nnfilenm(this, f90wrap_nnfilenm)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in)   :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    character(132), intent(out) :: f90wrap_nnfilenm
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nnfilenm = this_ptr%p%nnfilenm
end subroutine f90wrap_hopping_py__get__nnfilenm

subroutine f90wrap_hopping_py__set__nnfilenm(this, f90wrap_nnfilenm)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in)   :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    character(132), intent(in) :: f90wrap_nnfilenm
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nnfilenm = f90wrap_nnfilenm
end subroutine f90wrap_hopping_py__set__nnfilenm

subroutine f90wrap_hopping_py__get__nnfilenmo(this, f90wrap_nnfilenmo)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in)   :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    character(132), intent(out) :: f90wrap_nnfilenmo
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nnfilenmo = this_ptr%p%nnfilenmo
end subroutine f90wrap_hopping_py__get__nnfilenmo

subroutine f90wrap_hopping_py__set__nnfilenmo(this, f90wrap_nnfilenmo)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in)   :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    character(132), intent(in) :: f90wrap_nnfilenmo
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nnfilenmo = f90wrap_nnfilenmo
end subroutine f90wrap_hopping_py__set__nnfilenmo

subroutine f90wrap_hopping_py__array__flag_site_cindex(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%flag_site_cindex)) then
        dshape(1:1) = shape(this_ptr%p%flag_site_cindex)
        dloc = loc(this_ptr%p%flag_site_cindex)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__flag_site_cindex

subroutine f90wrap_hopping_py__array__i_coord(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%i_coord)) then
        dshape(1:2) = shape(this_ptr%p%i_coord)
        dloc = loc(this_ptr%p%i_coord)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__i_coord

subroutine f90wrap_hopping_py__array__j_coord(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%j_coord)) then
        dshape(1:2) = shape(this_ptr%p%j_coord)
        dloc = loc(this_ptr%p%j_coord)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__j_coord

subroutine f90wrap_hopping_py__array__Rij(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%Rij)) then
        dshape(1:2) = shape(this_ptr%p%Rij)
        dloc = loc(this_ptr%p%Rij)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__Rij

subroutine f90wrap_hopping_py__array__R(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%R)) then
        dshape(1:2) = shape(this_ptr%p%R)
        dloc = loc(this_ptr%p%R)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__R

subroutine f90wrap_hopping_py__array__Dij(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%Dij)) then
        dshape(1:1) = shape(this_ptr%p%Dij)
        dloc = loc(this_ptr%p%Dij)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__Dij

subroutine f90wrap_hopping_py__array__Dij0(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%Dij0)) then
        dshape(1:1) = shape(this_ptr%p%Dij0)
        dloc = loc(this_ptr%p%Dij0)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__Dij0

subroutine f90wrap_hopping_py__array__i_sign(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%i_sign)) then
        dshape(1:1) = shape(this_ptr%p%i_sign)
        dloc = loc(this_ptr%p%i_sign)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__i_sign

subroutine f90wrap_hopping_py__array__j_sign(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%j_sign)) then
        dshape(1:1) = shape(this_ptr%p%j_sign)
        dloc = loc(this_ptr%p%j_sign)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__j_sign

subroutine f90wrap_hopping_py__array__R_nn(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%R_nn)) then
        dshape(1:2) = shape(this_ptr%p%R_nn)
        dloc = loc(this_ptr%p%R_nn)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__R_nn

subroutine f90wrap_hopping_py__array__R0_nn(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%R0_nn)) then
        dshape(1:2) = shape(this_ptr%p%R0_nn)
        dloc = loc(this_ptr%p%R0_nn)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__R0_nn

subroutine f90wrap_hopping_py__array__local_charge(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%local_charge)) then
        dshape(1:1) = shape(this_ptr%p%local_charge)
        dloc = loc(this_ptr%p%local_charge)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__local_charge

subroutine f90wrap_hopping_py__array__local_moment(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%local_moment)) then
        dshape(1:2) = shape(this_ptr%p%local_moment)
        dloc = loc(this_ptr%p%local_moment)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__local_moment

subroutine f90wrap_hopping_py__array__local_moment_rot(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%local_moment_rot)) then
        dshape(1:2) = shape(this_ptr%p%local_moment_rot)
        dloc = loc(this_ptr%p%local_moment_rot)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__local_moment_rot

subroutine f90wrap_hopping_py__array__i_atom(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%i_atom)) then
        dshape(1:1) = shape(this_ptr%p%i_atom)
        dloc = loc(this_ptr%p%i_atom)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__i_atom

subroutine f90wrap_hopping_py__array__j_atom(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%j_atom)) then
        dshape(1:1) = shape(this_ptr%p%j_atom)
        dloc = loc(this_ptr%p%j_atom)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__j_atom

subroutine f90wrap_hopping_py__array__i_matrix(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%i_matrix)) then
        dshape(1:1) = shape(this_ptr%p%i_matrix)
        dloc = loc(this_ptr%p%i_matrix)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__i_matrix

subroutine f90wrap_hopping_py__array__j_matrix(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%j_matrix)) then
        dshape(1:1) = shape(this_ptr%p%j_matrix)
        dloc = loc(this_ptr%p%j_matrix)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__j_matrix

subroutine f90wrap_hopping_py__array__n_class(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%n_class)) then
        dshape(1:1) = shape(this_ptr%p%n_class)
        dloc = loc(this_ptr%p%n_class)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__n_class

subroutine f90wrap_hopping_py__array__sk_index_set(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%sk_index_set)) then
        dshape(1:2) = shape(this_ptr%p%sk_index_set)
        dloc = loc(this_ptr%p%sk_index_set)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__sk_index_set

subroutine f90wrap_hopping_py__array__cc_index_set(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%cc_index_set)) then
        dshape(1:2) = shape(this_ptr%p%cc_index_set)
        dloc = loc(this_ptr%p%cc_index_set)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__cc_index_set

subroutine f90wrap_hopping_py__array__n_nn(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%n_nn)) then
        dshape(1:1) = shape(this_ptr%p%n_nn)
        dloc = loc(this_ptr%p%n_nn)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__n_nn

subroutine f90wrap_hopping_py__array__j_nn(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%j_nn)) then
        dshape(1:2) = shape(this_ptr%p%j_nn)
        dloc = loc(this_ptr%p%j_nn)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__j_nn

subroutine f90wrap_hopping_py__array__l_onsite_param_index(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%l_onsite_param_index)) then
        dshape(1:1) = shape(this_ptr%p%l_onsite_param_index)
        dloc = loc(this_ptr%p%l_onsite_param_index)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__l_onsite_param_index

subroutine f90wrap_hopping_py__array__stoner_I_param_index(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%stoner_I_param_index)) then
        dshape(1:1) = shape(this_ptr%p%stoner_I_param_index)
        dloc = loc(this_ptr%p%stoner_I_param_index)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__stoner_I_param_index

subroutine f90wrap_hopping_py__array__local_U_param_index(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%local_U_param_index)) then
        dshape(1:1) = shape(this_ptr%p%local_U_param_index)
        dloc = loc(this_ptr%p%local_U_param_index)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__local_U_param_index

subroutine f90wrap_hopping_py__array__plus_U_param_index(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%plus_U_param_index)) then
        dshape(1:1) = shape(this_ptr%p%plus_U_param_index)
        dloc = loc(this_ptr%p%plus_U_param_index)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__plus_U_param_index

subroutine f90wrap_hopping_py__array__soc_param_index(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%soc_param_index)) then
        dshape(1:1) = shape(this_ptr%p%soc_param_index)
        dloc = loc(this_ptr%p%soc_param_index)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__soc_param_index

subroutine f90wrap_hopping_py__array__ci_orb(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%ci_orb)) then
        dshape(1:2) = (/len(this_ptr%p%ci_orb(1)), shape(this_ptr%p%ci_orb)/)
        dloc = loc(this_ptr%p%ci_orb)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__ci_orb

subroutine f90wrap_hopping_py__array__cj_orb(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%cj_orb)) then
        dshape(1:2) = (/len(this_ptr%p%cj_orb(1)), shape(this_ptr%p%cj_orb)/)
        dloc = loc(this_ptr%p%cj_orb)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__cj_orb

subroutine f90wrap_hopping_py__array__p_class(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%p_class)) then
        dshape(1:2) = (/len(this_ptr%p%p_class(1)), shape(this_ptr%p%p_class)/)
        dloc = loc(this_ptr%p%p_class)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__p_class

subroutine f90wrap_hopping_py__array__site_cindex(this, nd, dtype, dshape, dloc)
    use pyfit, only: hopping_py
    implicit none
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    integer, intent(in) :: this(2)
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%site_cindex)) then
        dshape(1:2) = (/len(this_ptr%p%site_cindex(1)), shape(this_ptr%p%site_cindex)/)
        dloc = loc(this_ptr%p%site_cindex)
    else
        dloc = 0
    end if
end subroutine f90wrap_hopping_py__array__site_cindex

subroutine f90wrap_hopping_py_initialise(this)
    use pyfit, only: hopping_py
    implicit none
    
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_hopping_py_initialise

subroutine f90wrap_hopping_py_finalise(this)
    use pyfit, only: hopping_py
    implicit none
    
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    type(hopping_py_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_hopping_py_finalise

subroutine f90wrap_energy_py__get__mysystem(this, f90wrap_mysystem)
    use pyfit, only: energy_py
    implicit none
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    integer, intent(in)   :: this(2)
    type(energy_py_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_mysystem
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_mysystem = this_ptr%p%mysystem
end subroutine f90wrap_energy_py__get__mysystem

subroutine f90wrap_energy_py__set__mysystem(this, f90wrap_mysystem)
    use pyfit, only: energy_py
    implicit none
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    integer, intent(in)   :: this(2)
    type(energy_py_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_mysystem
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%mysystem = f90wrap_mysystem
end subroutine f90wrap_energy_py__set__mysystem

subroutine f90wrap_energy_py__array__E(this, nd, dtype, dshape, dloc)
    use pyfit, only: energy_py
    implicit none
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    integer, intent(in) :: this(2)
    type(energy_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%E)) then
        dshape(1:2) = shape(this_ptr%p%E)
        dloc = loc(this_ptr%p%E)
    else
        dloc = 0
    end if
end subroutine f90wrap_energy_py__array__E

subroutine f90wrap_energy_py__array__dE(this, nd, dtype, dshape, dloc)
    use pyfit, only: energy_py
    implicit none
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    integer, intent(in) :: this(2)
    type(energy_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%dE)) then
        dshape(1:2) = shape(this_ptr%p%dE)
        dloc = loc(this_ptr%p%dE)
    else
        dloc = 0
    end if
end subroutine f90wrap_energy_py__array__dE

subroutine f90wrap_energy_py__array__ORB(this, nd, dtype, dshape, dloc)
    use pyfit, only: energy_py
    implicit none
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    integer, intent(in) :: this(2)
    type(energy_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%ORB)) then
        dshape(1:3) = shape(this_ptr%p%ORB)
        dloc = loc(this_ptr%p%ORB)
    else
        dloc = 0
    end if
end subroutine f90wrap_energy_py__array__ORB

subroutine f90wrap_energy_py__array__V(this, nd, dtype, dshape, dloc)
    use pyfit, only: energy_py
    implicit none
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    integer, intent(in) :: this(2)
    type(energy_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 15
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%V)) then
        dshape(1:3) = shape(this_ptr%p%V)
        dloc = loc(this_ptr%p%V)
    else
        dloc = 0
    end if
end subroutine f90wrap_energy_py__array__V

subroutine f90wrap_energy_py__array__SV(this, nd, dtype, dshape, dloc)
    use pyfit, only: energy_py
    implicit none
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    integer, intent(in) :: this(2)
    type(energy_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 15
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%SV)) then
        dshape(1:3) = shape(this_ptr%p%SV)
        dloc = loc(this_ptr%p%SV)
    else
        dloc = 0
    end if
end subroutine f90wrap_energy_py__array__SV

subroutine f90wrap_energy_py__array__E_BAND(this, nd, dtype, dshape, dloc)
    use pyfit, only: energy_py
    implicit none
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    integer, intent(in) :: this(2)
    type(energy_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%E_BAND)) then
        dshape(1:1) = shape(this_ptr%p%E_BAND)
        dloc = loc(this_ptr%p%E_BAND)
    else
        dloc = 0
    end if
end subroutine f90wrap_energy_py__array__E_BAND

subroutine f90wrap_energy_py__array__E_TOT(this, nd, dtype, dshape, dloc)
    use pyfit, only: energy_py
    implicit none
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    integer, intent(in) :: this(2)
    type(energy_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%E_TOT)) then
        dshape(1:1) = shape(this_ptr%p%E_TOT)
        dloc = loc(this_ptr%p%E_TOT)
    else
        dloc = 0
    end if
end subroutine f90wrap_energy_py__array__E_TOT

subroutine f90wrap_energy_py__array__F_OCC(this, nd, dtype, dshape, dloc)
    use pyfit, only: energy_py
    implicit none
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    integer, intent(in) :: this(2)
    type(energy_py_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%F_OCC)) then
        dshape(1:2) = shape(this_ptr%p%F_OCC)
        dloc = loc(this_ptr%p%F_OCC)
    else
        dloc = 0
    end if
end subroutine f90wrap_energy_py__array__F_OCC

subroutine f90wrap_energy_py__get__E_F(this, f90wrap_E_F)
    use pyfit, only: energy_py
    implicit none
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    integer, intent(in)   :: this(2)
    type(energy_py_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_E_F
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_E_F = this_ptr%p%E_F
end subroutine f90wrap_energy_py__get__E_F

subroutine f90wrap_energy_py__set__E_F(this, f90wrap_E_F)
    use pyfit, only: energy_py
    implicit none
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    integer, intent(in)   :: this(2)
    type(energy_py_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_E_F
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%E_F = f90wrap_E_F
end subroutine f90wrap_energy_py__set__E_F

subroutine f90wrap_energy_py_initialise(this)
    use pyfit, only: energy_py
    implicit none
    
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    type(energy_py_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_energy_py_initialise

subroutine f90wrap_energy_py_finalise(this)
    use pyfit, only: energy_py
    implicit none
    
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    type(energy_py_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_energy_py_finalise

subroutine f90wrap_init(comm, pinpt_py, ppram_py, pkpts_py, pwght_py, pgeom_py, nn_table_py, edft_py, etba_py)
    use pyfit, only: params_py, hopping_py, weight_py, kpoints_py, incar_py, poscar_py, init, energy_py
    implicit none
    
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    integer, intent(in) :: comm
    type(incar_py_ptr_type) :: pinpt_py_ptr
    integer, intent(in), dimension(2) :: pinpt_py
    type(params_py_ptr_type) :: ppram_py_ptr
    integer, intent(in), dimension(2) :: ppram_py
    type(kpoints_py_ptr_type) :: pkpts_py_ptr
    integer, intent(in), dimension(2) :: pkpts_py
    type(weight_py_ptr_type) :: pwght_py_ptr
    integer, intent(in), dimension(2) :: pwght_py
    type(poscar_py_ptr_type) :: pgeom_py_ptr
    integer, intent(in), dimension(2) :: pgeom_py
    type(hopping_py_ptr_type) :: nn_table_py_ptr
    integer, intent(in), dimension(2) :: nn_table_py
    type(energy_py_ptr_type) :: edft_py_ptr
    integer, intent(in), dimension(2) :: edft_py
    type(energy_py_ptr_type) :: etba_py_ptr
    integer, intent(in), dimension(2) :: etba_py
    pinpt_py_ptr = transfer(pinpt_py, pinpt_py_ptr)
    ppram_py_ptr = transfer(ppram_py, ppram_py_ptr)
    pkpts_py_ptr = transfer(pkpts_py, pkpts_py_ptr)
    pwght_py_ptr = transfer(pwght_py, pwght_py_ptr)
    pgeom_py_ptr = transfer(pgeom_py, pgeom_py_ptr)
    nn_table_py_ptr = transfer(nn_table_py, nn_table_py_ptr)
    edft_py_ptr = transfer(edft_py, edft_py_ptr)
    etba_py_ptr = transfer(etba_py, etba_py_ptr)
    call init(comm=comm, PINPT_PY=pinpt_py_ptr%p, PPRAM_PY=ppram_py_ptr%p, PKPTS_PY=pkpts_py_ptr%p, PWGHT_PY=pwght_py_ptr%p, &
        PGEOM_PY=pgeom_py_ptr%p, NN_TABLE_PY=nn_table_py_ptr%p, EDFT_PY=edft_py_ptr%p, ETBA_PY=etba_py_ptr%p)
end subroutine f90wrap_init

subroutine f90wrap_eig(comm, pinpt_py, ppram_py, pkpts_py, pgeom_py, nn_table_py, etba_py)
    use pyfit, only: params_py, hopping_py, eig, kpoints_py, incar_py, poscar_py, energy_py
    implicit none
    
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    integer, intent(inout) :: comm
    type(incar_py_ptr_type) :: pinpt_py_ptr
    integer, intent(in), dimension(2) :: pinpt_py
    type(params_py_ptr_type) :: ppram_py_ptr
    integer, intent(in), dimension(2) :: ppram_py
    type(kpoints_py_ptr_type) :: pkpts_py_ptr
    integer, intent(in), dimension(2) :: pkpts_py
    type(poscar_py_ptr_type) :: pgeom_py_ptr
    integer, intent(in), dimension(2) :: pgeom_py
    type(hopping_py_ptr_type) :: nn_table_py_ptr
    integer, intent(in), dimension(2) :: nn_table_py
    type(energy_py_ptr_type) :: etba_py_ptr
    integer, intent(in), dimension(2) :: etba_py
    pinpt_py_ptr = transfer(pinpt_py, pinpt_py_ptr)
    ppram_py_ptr = transfer(ppram_py, ppram_py_ptr)
    pkpts_py_ptr = transfer(pkpts_py, pkpts_py_ptr)
    pgeom_py_ptr = transfer(pgeom_py, pgeom_py_ptr)
    nn_table_py_ptr = transfer(nn_table_py, nn_table_py_ptr)
    etba_py_ptr = transfer(etba_py, etba_py_ptr)
    call eig(comm=comm, PINPT_PY=pinpt_py_ptr%p, PPRAM_PY=ppram_py_ptr%p, PKPTS_PY=pkpts_py_ptr%p, PGEOM_PY=pgeom_py_ptr%p, &
        NN_TABLE_PY=nn_table_py_ptr%p, ETBA_PY=etba_py_ptr%p)
end subroutine f90wrap_eig

subroutine f90wrap_toten(comm, pinpt_py, pkpts_py, pgeom_py, etba_py)
    use pyfit, only: kpoints_py, incar_py, toten, poscar_py, energy_py
    implicit none
    
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    integer, intent(in) :: comm
    type(incar_py_ptr_type) :: pinpt_py_ptr
    integer, intent(in), dimension(2) :: pinpt_py
    type(kpoints_py_ptr_type) :: pkpts_py_ptr
    integer, intent(in), dimension(2) :: pkpts_py
    type(poscar_py_ptr_type) :: pgeom_py_ptr
    integer, intent(in), dimension(2) :: pgeom_py
    type(energy_py_ptr_type) :: etba_py_ptr
    integer, intent(in), dimension(2) :: etba_py
    pinpt_py_ptr = transfer(pinpt_py, pinpt_py_ptr)
    pkpts_py_ptr = transfer(pkpts_py, pkpts_py_ptr)
    pgeom_py_ptr = transfer(pgeom_py, pgeom_py_ptr)
    etba_py_ptr = transfer(etba_py, etba_py_ptr)
    call toten(comm=comm, PINPT_PY=pinpt_py_ptr%p, PKPTS_PY=pkpts_py_ptr%p, PGEOM_PY=pgeom_py_ptr%p, ETBA_PY=etba_py_ptr%p)
end subroutine f90wrap_toten

subroutine f90wrap_fit(comm, pinpt_py, ppram_py, pkpts_py, pwght_py, pgeom_py, nn_table_py, edft_py, etba_py)
    use pyfit, only: params_py, hopping_py, weight_py, fit, kpoints_py, incar_py, poscar_py, energy_py
    implicit none
    
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    integer, intent(in) :: comm
    type(incar_py_ptr_type) :: pinpt_py_ptr
    integer, intent(in), dimension(2) :: pinpt_py
    type(params_py_ptr_type) :: ppram_py_ptr
    integer, intent(in), dimension(2) :: ppram_py
    type(kpoints_py_ptr_type) :: pkpts_py_ptr
    integer, intent(in), dimension(2) :: pkpts_py
    type(weight_py_ptr_type) :: pwght_py_ptr
    integer, intent(in), dimension(2) :: pwght_py
    type(poscar_py_ptr_type) :: pgeom_py_ptr
    integer, intent(in), dimension(2) :: pgeom_py
    type(hopping_py_ptr_type) :: nn_table_py_ptr
    integer, intent(in), dimension(2) :: nn_table_py
    type(energy_py_ptr_type) :: edft_py_ptr
    integer, intent(in), dimension(2) :: edft_py
    type(energy_py_ptr_type) :: etba_py_ptr
    integer, intent(in), dimension(2) :: etba_py
    pinpt_py_ptr = transfer(pinpt_py, pinpt_py_ptr)
    ppram_py_ptr = transfer(ppram_py, ppram_py_ptr)
    pkpts_py_ptr = transfer(pkpts_py, pkpts_py_ptr)
    pwght_py_ptr = transfer(pwght_py, pwght_py_ptr)
    pgeom_py_ptr = transfer(pgeom_py, pgeom_py_ptr)
    nn_table_py_ptr = transfer(nn_table_py, nn_table_py_ptr)
    edft_py_ptr = transfer(edft_py, edft_py_ptr)
    etba_py_ptr = transfer(etba_py, etba_py_ptr)
    call fit(comm=comm, PINPT_PY=pinpt_py_ptr%p, PPRAM_PY=ppram_py_ptr%p, PKPTS_PY=pkpts_py_ptr%p, PWGHT_PY=pwght_py_ptr%p, &
        PGEOM_PY=pgeom_py_ptr%p, NN_TABLE_PY=nn_table_py_ptr%p, EDFT_PY=edft_py_ptr%p, ETBA_PY=etba_py_ptr%p)
end subroutine f90wrap_fit

subroutine f90wrap_init_incar_py(ifilenm, ret_pinpt_py, nsystem, n0)
    use pyfit, only: init_incar_py, incar_py
    implicit none
    
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    character(1024), intent(in), dimension(n0) :: ifilenm
    type(incar_py_ptr_type) :: ret_pinpt_py_ptr
    integer, intent(out), dimension(2) :: ret_pinpt_py
    integer(4), intent(in) :: nsystem
    integer :: n0
    !f2py intent(hide), depend(ifilenm) :: n0 = shape(ifilenm,0)
    allocate(ret_pinpt_py_ptr%p)
    ret_pinpt_py_ptr%p = init_incar_py(ifilenm=ifilenm, nsystem=nsystem)
    ret_pinpt_py = transfer(ret_pinpt_py_ptr, ret_pinpt_py)
end subroutine f90wrap_init_incar_py

subroutine f90wrap_init_params_py(ret_ppram_py)
    use pyfit, only: params_py, init_params_py
    implicit none
    
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    type(params_py_ptr_type) :: ret_ppram_py_ptr
    integer, intent(out), dimension(2) :: ret_ppram_py
    allocate(ret_ppram_py_ptr%p)
    ret_ppram_py_ptr%p = init_params_py()
    ret_ppram_py = transfer(ret_ppram_py_ptr, ret_ppram_py)
end subroutine f90wrap_init_params_py

subroutine f90wrap_init_kpoints_py(ret_pkpts_py)
    use pyfit, only: init_kpoints_py, kpoints_py
    implicit none
    
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    type(kpoints_py_ptr_type) :: ret_pkpts_py_ptr
    integer, intent(out), dimension(2) :: ret_pkpts_py
    allocate(ret_pkpts_py_ptr%p)
    ret_pkpts_py_ptr%p = init_kpoints_py()
    ret_pkpts_py = transfer(ret_pkpts_py_ptr, ret_pkpts_py)
end subroutine f90wrap_init_kpoints_py

subroutine f90wrap_init_weight_py(ret_pwght_py)
    use pyfit, only: init_weight_py, weight_py
    implicit none
    
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    type(weight_py_ptr_type) :: ret_pwght_py_ptr
    integer, intent(out), dimension(2) :: ret_pwght_py
    allocate(ret_pwght_py_ptr%p)
    ret_pwght_py_ptr%p = init_weight_py()
    ret_pwght_py = transfer(ret_pwght_py_ptr, ret_pwght_py)
end subroutine f90wrap_init_weight_py

subroutine f90wrap_init_poscar_py(ret_pgeom_py)
    use pyfit, only: poscar_py, init_poscar_py
    implicit none
    
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    type(poscar_py_ptr_type) :: ret_pgeom_py_ptr
    integer, intent(out), dimension(2) :: ret_pgeom_py
    allocate(ret_pgeom_py_ptr%p)
    ret_pgeom_py_ptr%p = init_poscar_py()
    ret_pgeom_py = transfer(ret_pgeom_py_ptr, ret_pgeom_py)
end subroutine f90wrap_init_poscar_py

subroutine f90wrap_init_hopping_py(ret_nn_table_py)
    use pyfit, only: init_hopping_py, hopping_py
    implicit none
    
    type hopping_py_ptr_type
        type(hopping_py), pointer :: p => NULL()
    end type hopping_py_ptr_type
    type(hopping_py_ptr_type) :: ret_nn_table_py_ptr
    integer, intent(out), dimension(2) :: ret_nn_table_py
    allocate(ret_nn_table_py_ptr%p)
    ret_nn_table_py_ptr%p = init_hopping_py()
    ret_nn_table_py = transfer(ret_nn_table_py_ptr, ret_nn_table_py)
end subroutine f90wrap_init_hopping_py

subroutine f90wrap_init_energy_py(ret_e_py)
    use pyfit, only: energy_py, init_energy_py
    implicit none
    
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    type(energy_py_ptr_type) :: ret_e_py_ptr
    integer, intent(out), dimension(2) :: ret_e_py
    allocate(ret_e_py_ptr%p)
    ret_e_py_ptr%p = init_energy_py()
    ret_e_py = transfer(ret_e_py_ptr, ret_e_py)
end subroutine f90wrap_init_energy_py

subroutine f90wrap_print_param_py(pinpt_py, ppram_py, pfileoutnm)
    use pyfit, only: print_param_py, params_py, incar_py
    implicit none
    
    type params_py_ptr_type
        type(params_py), pointer :: p => NULL()
    end type params_py_ptr_type
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    type(incar_py_ptr_type) :: pinpt_py_ptr
    integer, intent(in), dimension(2) :: pinpt_py
    type(params_py_ptr_type) :: ppram_py_ptr
    integer, intent(in), dimension(2) :: ppram_py
    character(132), intent(inout) :: pfileoutnm
    pinpt_py_ptr = transfer(pinpt_py, pinpt_py_ptr)
    ppram_py_ptr = transfer(ppram_py, ppram_py_ptr)
    call print_param_py(PINPT_PY=pinpt_py_ptr%p, PPRAM_PY=ppram_py_ptr%p, pfileoutnm=pfileoutnm)
end subroutine f90wrap_print_param_py

subroutine f90wrap_print_weight(pwght_py, wfileoutnm)
    use pyfit, only: weight_py, print_weight
    implicit none
    
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    type(weight_py_ptr_type) :: pwght_py_ptr
    integer, intent(in), dimension(2) :: pwght_py
    character(132), intent(inout) :: wfileoutnm
    pwght_py_ptr = transfer(pwght_py, pwght_py_ptr)
    call print_weight(PWGHT_PY=pwght_py_ptr%p, wfileoutnm=wfileoutnm)
end subroutine f90wrap_print_weight

subroutine f90wrap_print_target(pinpt_py, pkpts_py, edft_py, pwght_py, pgeom_py, tfileoutnm)
    use pyfit, only: weight_py, print_target, kpoints_py, incar_py, poscar_py, energy_py
    implicit none
    
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    type(incar_py_ptr_type) :: pinpt_py_ptr
    integer, intent(in), dimension(2) :: pinpt_py
    type(kpoints_py_ptr_type) :: pkpts_py_ptr
    integer, intent(in), dimension(2) :: pkpts_py
    type(energy_py_ptr_type) :: edft_py_ptr
    integer, intent(in), dimension(2) :: edft_py
    type(weight_py_ptr_type) :: pwght_py_ptr
    integer, intent(in), dimension(2) :: pwght_py
    type(poscar_py_ptr_type) :: pgeom_py_ptr
    integer, intent(in), dimension(2) :: pgeom_py
    character(132), intent(inout) :: tfileoutnm
    pinpt_py_ptr = transfer(pinpt_py, pinpt_py_ptr)
    pkpts_py_ptr = transfer(pkpts_py, pkpts_py_ptr)
    edft_py_ptr = transfer(edft_py, edft_py_ptr)
    pwght_py_ptr = transfer(pwght_py, pwght_py_ptr)
    pgeom_py_ptr = transfer(pgeom_py, pgeom_py_ptr)
    call print_target(PINPT_PY=pinpt_py_ptr%p, PKPTS_PY=pkpts_py_ptr%p, EDFT_PY=edft_py_ptr%p, PWGHT_PY=pwght_py_ptr%p, &
        PGEOM_PY=pgeom_py_ptr%p, tfileoutnm=tfileoutnm)
end subroutine f90wrap_print_target

subroutine f90wrap_print_etba(pinpt_py, pkpts_py, etba_py, edft_py, pwght_py, pgeom_py, suffix, flag_use_overlap)
    use pyfit, only: weight_py, print_etba, kpoints_py, incar_py, poscar_py, energy_py
    implicit none
    
    type weight_py_ptr_type
        type(weight_py), pointer :: p => NULL()
    end type weight_py_ptr_type
    type kpoints_py_ptr_type
        type(kpoints_py), pointer :: p => NULL()
    end type kpoints_py_ptr_type
    type incar_py_ptr_type
        type(incar_py), pointer :: p => NULL()
    end type incar_py_ptr_type
    type poscar_py_ptr_type
        type(poscar_py), pointer :: p => NULL()
    end type poscar_py_ptr_type
    type energy_py_ptr_type
        type(energy_py), pointer :: p => NULL()
    end type energy_py_ptr_type
    type(incar_py_ptr_type) :: pinpt_py_ptr
    integer, intent(in), dimension(2) :: pinpt_py
    type(kpoints_py_ptr_type) :: pkpts_py_ptr
    integer, intent(in), dimension(2) :: pkpts_py
    type(energy_py_ptr_type) :: etba_py_ptr
    integer, intent(in), dimension(2) :: etba_py
    type(energy_py_ptr_type) :: edft_py_ptr
    integer, intent(in), dimension(2) :: edft_py
    type(weight_py_ptr_type) :: pwght_py_ptr
    integer, intent(in), dimension(2) :: pwght_py
    type(poscar_py_ptr_type) :: pgeom_py_ptr
    integer, intent(in), dimension(2) :: pgeom_py
    character(132), intent(inout) :: suffix
    logical :: flag_use_overlap
    pinpt_py_ptr = transfer(pinpt_py, pinpt_py_ptr)
    pkpts_py_ptr = transfer(pkpts_py, pkpts_py_ptr)
    etba_py_ptr = transfer(etba_py, etba_py_ptr)
    edft_py_ptr = transfer(edft_py, edft_py_ptr)
    pwght_py_ptr = transfer(pwght_py, pwght_py_ptr)
    pgeom_py_ptr = transfer(pgeom_py, pgeom_py_ptr)
    call print_etba(PINPT_PY=pinpt_py_ptr%p, PKPTS_PY=pkpts_py_ptr%p, ETBA_PY=etba_py_ptr%p, EDFT_PY=edft_py_ptr%p, &
        PWGHT_PY=pwght_py_ptr%p, PGEOM_PY=pgeom_py_ptr%p, suffix=suffix, flag_use_overlap=flag_use_overlap)
end subroutine f90wrap_print_etba

subroutine f90wrap_set_verbose(iverbose_mod)
    use pyfit, only: set_verbose
    implicit none
    
    integer(4) :: iverbose_mod
    call set_verbose(iverbose_mod=iverbose_mod)
end subroutine f90wrap_set_verbose

! End of module pyfit defined in file tbfitpy_mod.f90

