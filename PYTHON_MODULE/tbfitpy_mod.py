from __future__ import print_function, absolute_import, division
import _tbfitpy_mod
import f90wrap.runtime
import logging

class Pyfit(f90wrap.runtime.FortranModule):
    """
    Module pyfit
    
    
    Defined at tbfitpy_mod.f90 lines 2-1991
    
    """
    @f90wrap.runtime.register_class("tbfitpy_mod.incar_py")
    class incar_py(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=incar_py)
        
        
        Defined at tbfitpy_mod.f90 lines 10-66
        
        """
        def __init__(self, handle=None):
            """
            self = Incar_Py()
            
            
            Defined at tbfitpy_mod.f90 lines 10-66
            
            
            Returns
            -------
            this : Incar_Py
            	Object to be constructed
            
            
            Automatically generated constructor for incar_py
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _tbfitpy_mod.f90wrap_incar_py_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Incar_Py
            
            
            Defined at tbfitpy_mod.f90 lines 10-66
            
            Parameters
            ----------
            this : Incar_Py
            	Object to be destructed
            
            
            Automatically generated destructor for incar_py
            """
            if self._alloc:
                _tbfitpy_mod.f90wrap_incar_py_finalise(this=self._handle)
        
        @property
        def flag_python_module(self):
            """
            Element flag_python_module ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 12
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_python_module(self._handle)
        
        @flag_python_module.setter
        def flag_python_module(self, flag_python_module):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_python_module(self._handle, \
                flag_python_module)
        
        @property
        def flag_get_band(self):
            """
            Element flag_get_band ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 13
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_get_band(self._handle)
        
        @flag_get_band.setter
        def flag_get_band(self, flag_get_band):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_get_band(self._handle, flag_get_band)
        
        @property
        def flag_fit_degeneracy(self):
            """
            Element flag_fit_degeneracy ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 14
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_fit_degeneracy(self._handle)
        
        @flag_fit_degeneracy.setter
        def flag_fit_degeneracy(self, flag_fit_degeneracy):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_fit_degeneracy(self._handle, \
                flag_fit_degeneracy)
        
        @property
        def flag_report_geom(self):
            """
            Element flag_report_geom ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 15
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_report_geom(self._handle)
        
        @flag_report_geom.setter
        def flag_report_geom(self, flag_report_geom):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_report_geom(self._handle, \
                flag_report_geom)
        
        @property
        def flag_phase(self):
            """
            Element flag_phase ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 16
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_phase(self._handle)
        
        @flag_phase.setter
        def flag_phase(self, flag_phase):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_phase(self._handle, flag_phase)
        
        @property
        def flag_tbfit(self):
            """
            Element flag_tbfit ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 17
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_tbfit(self._handle)
        
        @flag_tbfit.setter
        def flag_tbfit(self, flag_tbfit):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_tbfit(self._handle, flag_tbfit)
        
        @property
        def flag_tbfit_finish(self):
            """
            Element flag_tbfit_finish ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 18
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_tbfit_finish(self._handle)
        
        @flag_tbfit_finish.setter
        def flag_tbfit_finish(self, flag_tbfit_finish):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_tbfit_finish(self._handle, \
                flag_tbfit_finish)
        
        @property
        def flag_print_only_target(self):
            """
            Element flag_print_only_target ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 19
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_print_only_target(self._handle)
        
        @flag_print_only_target.setter
        def flag_print_only_target(self, flag_print_only_target):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_print_only_target(self._handle, \
                flag_print_only_target)
        
        @property
        def flag_print_energy_diff(self):
            """
            Element flag_print_energy_diff ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 20
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_print_energy_diff(self._handle)
        
        @flag_print_energy_diff.setter
        def flag_print_energy_diff(self, flag_print_energy_diff):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_print_energy_diff(self._handle, \
                flag_print_energy_diff)
        
        @property
        def flag_print_orbital(self):
            """
            Element flag_print_orbital ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 21
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_print_orbital(self._handle)
        
        @flag_print_orbital.setter
        def flag_print_orbital(self, flag_print_orbital):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_print_orbital(self._handle, \
                flag_print_orbital)
        
        @property
        def flag_get_orbital(self):
            """
            Element flag_get_orbital ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 22
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_get_orbital(self._handle)
        
        @flag_get_orbital.setter
        def flag_get_orbital(self, flag_get_orbital):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_get_orbital(self._handle, \
                flag_get_orbital)
        
        @property
        def flag_print_mag(self):
            """
            Element flag_print_mag ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 23
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_print_mag(self._handle)
        
        @flag_print_mag.setter
        def flag_print_mag(self, flag_print_mag):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_print_mag(self._handle, flag_print_mag)
        
        @property
        def flag_print_single(self):
            """
            Element flag_print_single ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 24
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_print_single(self._handle)
        
        @flag_print_single.setter
        def flag_print_single(self, flag_print_single):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_print_single(self._handle, \
                flag_print_single)
        
        @property
        def flag_local_charge(self):
            """
            Element flag_local_charge ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 25
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_local_charge(self._handle)
        
        @flag_local_charge.setter
        def flag_local_charge(self, flag_local_charge):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_local_charge(self._handle, \
                flag_local_charge)
        
        @property
        def flag_collinear(self):
            """
            Element flag_collinear ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 26
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_collinear(self._handle)
        
        @flag_collinear.setter
        def flag_collinear(self, flag_collinear):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_collinear(self._handle, flag_collinear)
        
        @property
        def flag_noncollinear(self):
            """
            Element flag_noncollinear ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 27
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_noncollinear(self._handle)
        
        @flag_noncollinear.setter
        def flag_noncollinear(self, flag_noncollinear):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_noncollinear(self._handle, \
                flag_noncollinear)
        
        @property
        def flag_soc(self):
            """
            Element flag_soc ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 28
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_soc(self._handle)
        
        @flag_soc.setter
        def flag_soc(self, flag_soc):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_soc(self._handle, flag_soc)
        
        @property
        def flag_plus_u(self):
            """
            Element flag_plus_u ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 29
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_plus_u(self._handle)
        
        @flag_plus_u.setter
        def flag_plus_u(self, flag_plus_u):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_plus_u(self._handle, flag_plus_u)
        
        @property
        def flag_scissor(self):
            """
            Element flag_scissor ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 30
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_scissor(self._handle)
        
        @flag_scissor.setter
        def flag_scissor(self, flag_scissor):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_scissor(self._handle, flag_scissor)
        
        @property
        def flag_print_proj(self):
            """
            Element flag_print_proj ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 31
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_print_proj(self._handle)
        
        @flag_print_proj.setter
        def flag_print_proj(self, flag_print_proj):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_print_proj(self._handle, \
                flag_print_proj)
        
        @property
        def flag_print_proj_sum(self):
            """
            Element flag_print_proj_sum ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 32
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_print_proj_sum(self._handle)
        
        @flag_print_proj_sum.setter
        def flag_print_proj_sum(self, flag_print_proj_sum):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_print_proj_sum(self._handle, \
                flag_print_proj_sum)
        
        @property
        def flag_print_proj_atom(self):
            """
            Element flag_print_proj_atom ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 33
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_print_proj_atom(self._handle)
        
        @flag_print_proj_atom.setter
        def flag_print_proj_atom(self, flag_print_proj_atom):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_print_proj_atom(self._handle, \
                flag_print_proj_atom)
        
        @property
        def flag_plot_fit(self):
            """
            Element flag_plot_fit ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 34
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_plot_fit(self._handle)
        
        @flag_plot_fit.setter
        def flag_plot_fit(self, flag_plot_fit):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_plot_fit(self._handle, flag_plot_fit)
        
        @property
        def flag_plot(self):
            """
            Element flag_plot ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 35
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_plot(self._handle)
        
        @flag_plot.setter
        def flag_plot(self, flag_plot):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_plot(self._handle, flag_plot)
        
        @property
        def flag_get_band_order(self):
            """
            Element flag_get_band_order ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 36
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_get_band_order(self._handle)
        
        @flag_get_band_order.setter
        def flag_get_band_order(self, flag_get_band_order):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_get_band_order(self._handle, \
                flag_get_band_order)
        
        @property
        def flag_get_band_order_print_only(self):
            """
            Element flag_get_band_order_print_only ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 37
            
            """
            return \
                _tbfitpy_mod.f90wrap_incar_py__get__flag_get_band_order_print_only(self._handle)
        
        @flag_get_band_order_print_only.setter
        def flag_get_band_order_print_only(self, flag_get_band_order_print_only):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_get_band_order_print_only(self._handle, \
                flag_get_band_order_print_only)
        
        @property
        def flag_use_weight(self):
            """
            Element flag_use_weight ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 38
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_use_weight(self._handle)
        
        @flag_use_weight.setter
        def flag_use_weight(self, flag_use_weight):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_use_weight(self._handle, \
                flag_use_weight)
        
        @property
        def flag_fit_orbital(self):
            """
            Element flag_fit_orbital ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 39
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_fit_orbital(self._handle)
        
        @flag_fit_orbital.setter
        def flag_fit_orbital(self, flag_fit_orbital):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_fit_orbital(self._handle, \
                flag_fit_orbital)
        
        @property
        def flag_fit_orbital_parse(self):
            """
            Element flag_fit_orbital_parse ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 40
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_fit_orbital_parse(self._handle)
        
        @flag_fit_orbital_parse.setter
        def flag_fit_orbital_parse(self, flag_fit_orbital_parse):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_fit_orbital_parse(self._handle, \
                flag_fit_orbital_parse)
        
        @property
        def flag_pso_with_lmdif(self):
            """
            Element flag_pso_with_lmdif ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 41
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_pso_with_lmdif(self._handle)
        
        @flag_pso_with_lmdif.setter
        def flag_pso_with_lmdif(self, flag_pso_with_lmdif):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_pso_with_lmdif(self._handle, \
                flag_pso_with_lmdif)
        
        @property
        def flag_get_unfold(self):
            """
            Element flag_get_unfold ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 42
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_get_unfold(self._handle)
        
        @flag_get_unfold.setter
        def flag_get_unfold(self, flag_get_unfold):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_get_unfold(self._handle, \
                flag_get_unfold)
        
        @property
        def ptol(self):
            """
            Element ptol ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 43
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__ptol(self._handle)
        
        @ptol.setter
        def ptol(self, ptol):
            _tbfitpy_mod.f90wrap_incar_py__set__ptol(self._handle, ptol)
        
        @property
        def ftol(self):
            """
            Element ftol ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 44
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__ftol(self._handle)
        
        @ftol.setter
        def ftol(self, ftol):
            _tbfitpy_mod.f90wrap_incar_py__set__ftol(self._handle, ftol)
        
        @property
        def fdiff(self):
            """
            Element fdiff ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 45
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__fdiff(self._handle)
        
        @fdiff.setter
        def fdiff(self, fdiff):
            _tbfitpy_mod.f90wrap_incar_py__set__fdiff(self._handle, fdiff)
        
        @property
        def r_scissor(self):
            """
            Element r_scissor ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 46
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__r_scissor(self._handle)
        
        @r_scissor.setter
        def r_scissor(self, r_scissor):
            _tbfitpy_mod.f90wrap_incar_py__set__r_scissor(self._handle, r_scissor)
        
        @property
        def nsystem(self):
            """
            Element nsystem ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 47
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__nsystem(self._handle)
        
        @nsystem.setter
        def nsystem(self, nsystem):
            _tbfitpy_mod.f90wrap_incar_py__set__nsystem(self._handle, nsystem)
        
        @property
        def lmmax(self):
            """
            Element lmmax ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 48
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__lmmax(self._handle)
        
        @lmmax.setter
        def lmmax(self, lmmax):
            _tbfitpy_mod.f90wrap_incar_py__set__lmmax(self._handle, lmmax)
        
        @property
        def miter(self):
            """
            Element miter ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 49
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__miter(self._handle)
        
        @miter.setter
        def miter(self, miter):
            _tbfitpy_mod.f90wrap_incar_py__set__miter(self._handle, miter)
        
        @property
        def mxfit(self):
            """
            Element mxfit ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 50
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__mxfit(self._handle)
        
        @mxfit.setter
        def mxfit(self, mxfit):
            _tbfitpy_mod.f90wrap_incar_py__set__mxfit(self._handle, mxfit)
        
        @property
        def nn_max(self):
            """
            Element nn_max ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 51
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_incar_py__array__nn_max(self._handle)
            if array_handle in self._arrays:
                nn_max = self._arrays[array_handle]
            else:
                nn_max = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_incar_py__array__nn_max)
                self._arrays[array_handle] = nn_max
            return nn_max
        
        @nn_max.setter
        def nn_max(self, nn_max):
            self.nn_max[...] = nn_max
        
        @property
        def ispin(self):
            """
            Element ispin ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 52
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__ispin(self._handle)
        
        @ispin.setter
        def ispin(self, ispin):
            _tbfitpy_mod.f90wrap_incar_py__set__ispin(self._handle, ispin)
        
        @property
        def ispinor(self):
            """
            Element ispinor ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 53
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__ispinor(self._handle)
        
        @ispinor.setter
        def ispinor(self, ispinor):
            _tbfitpy_mod.f90wrap_incar_py__set__ispinor(self._handle, ispinor)
        
        @property
        def nspin(self):
            """
            Element nspin ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 54
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__nspin(self._handle)
        
        @nspin.setter
        def nspin(self, nspin):
            _tbfitpy_mod.f90wrap_incar_py__set__nspin(self._handle, nspin)
        
        @property
        def nproj_sum(self):
            """
            Element nproj_sum ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 55
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__nproj_sum(self._handle)
        
        @nproj_sum.setter
        def nproj_sum(self, nproj_sum):
            _tbfitpy_mod.f90wrap_incar_py__set__nproj_sum(self._handle, nproj_sum)
        
        @property
        def iverbose(self):
            """
            Element iverbose ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 56
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__iverbose(self._handle)
        
        @iverbose.setter
        def iverbose(self, iverbose):
            _tbfitpy_mod.f90wrap_incar_py__set__iverbose(self._handle, iverbose)
        
        @property
        def i_scissor(self):
            """
            Element i_scissor ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 57
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__i_scissor(self._handle)
        
        @i_scissor.setter
        def i_scissor(self, i_scissor):
            _tbfitpy_mod.f90wrap_incar_py__set__i_scissor(self._handle, i_scissor)
        
        @property
        def ls_type(self):
            """
            Element ls_type ftype=character(len=10) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 58
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__ls_type(self._handle)
        
        @ls_type.setter
        def ls_type(self, ls_type):
            _tbfitpy_mod.f90wrap_incar_py__set__ls_type(self._handle, ls_type)
        
        @property
        def fnamelog(self):
            """
            Element fnamelog ftype=character(len=40) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 59
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__fnamelog(self._handle)
        
        @fnamelog.setter
        def fnamelog(self, fnamelog):
            _tbfitpy_mod.f90wrap_incar_py__set__fnamelog(self._handle, fnamelog)
        
        @property
        def proj_atom(self):
            """
            Element proj_atom ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 60
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_incar_py__array__proj_atom(self._handle)
            if array_handle in self._arrays:
                proj_atom = self._arrays[array_handle]
            else:
                proj_atom = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_incar_py__array__proj_atom)
                self._arrays[array_handle] = proj_atom
            return proj_atom
        
        @proj_atom.setter
        def proj_atom(self, proj_atom):
            self.proj_atom[...] = proj_atom
        
        @property
        def proj_natom(self):
            """
            Element proj_natom ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 61
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_incar_py__array__proj_natom(self._handle)
            if array_handle in self._arrays:
                proj_natom = self._arrays[array_handle]
            else:
                proj_natom = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_incar_py__array__proj_natom)
                self._arrays[array_handle] = proj_natom
            return proj_natom
        
        @proj_natom.setter
        def proj_natom(self, proj_natom):
            self.proj_natom[...] = proj_natom
        
        @property
        def ifilenm(self):
            """
            Element ifilenm ftype=character(len=132) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 62
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_incar_py__array__ifilenm(self._handle)
            if array_handle in self._arrays:
                ifilenm = self._arrays[array_handle]
            else:
                ifilenm = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_incar_py__array__ifilenm)
                self._arrays[array_handle] = ifilenm
            return ifilenm
        
        @ifilenm.setter
        def ifilenm(self, ifilenm):
            self.ifilenm[...] = ifilenm
        
        @property
        def title(self):
            """
            Element title ftype=character(len=132) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 63
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_incar_py__array__title(self._handle)
            if array_handle in self._arrays:
                title = self._arrays[array_handle]
            else:
                title = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_incar_py__array__title)
                self._arrays[array_handle] = title
            return title
        
        @title.setter
        def title(self, title):
            self.title[...] = title
        
        @property
        def flag_get_total_energy(self):
            """
            Element flag_get_total_energy ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 64
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__flag_get_total_energy(self._handle)
        
        @flag_get_total_energy.setter
        def flag_get_total_energy(self, flag_get_total_energy):
            _tbfitpy_mod.f90wrap_incar_py__set__flag_get_total_energy(self._handle, \
                flag_get_total_energy)
        
        @property
        def electronic_temperature(self):
            """
            Element electronic_temperature ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 65
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__electronic_temperature(self._handle)
        
        @electronic_temperature.setter
        def electronic_temperature(self, electronic_temperature):
            _tbfitpy_mod.f90wrap_incar_py__set__electronic_temperature(self._handle, \
                electronic_temperature)
        
        @property
        def orbital_fit_smearing(self):
            """
            Element orbital_fit_smearing ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 66
            
            """
            return _tbfitpy_mod.f90wrap_incar_py__get__orbital_fit_smearing(self._handle)
        
        @orbital_fit_smearing.setter
        def orbital_fit_smearing(self, orbital_fit_smearing):
            _tbfitpy_mod.f90wrap_incar_py__set__orbital_fit_smearing(self._handle, \
                orbital_fit_smearing)
        
        def __str__(self):
            ret = ['<incar_py>{\n']
            ret.append('    flag_python_module : ')
            ret.append(repr(self.flag_python_module))
            ret.append(',\n    flag_get_band : ')
            ret.append(repr(self.flag_get_band))
            ret.append(',\n    flag_fit_degeneracy : ')
            ret.append(repr(self.flag_fit_degeneracy))
            ret.append(',\n    flag_report_geom : ')
            ret.append(repr(self.flag_report_geom))
            ret.append(',\n    flag_phase : ')
            ret.append(repr(self.flag_phase))
            ret.append(',\n    flag_tbfit : ')
            ret.append(repr(self.flag_tbfit))
            ret.append(',\n    flag_tbfit_finish : ')
            ret.append(repr(self.flag_tbfit_finish))
            ret.append(',\n    flag_print_only_target : ')
            ret.append(repr(self.flag_print_only_target))
            ret.append(',\n    flag_print_energy_diff : ')
            ret.append(repr(self.flag_print_energy_diff))
            ret.append(',\n    flag_print_orbital : ')
            ret.append(repr(self.flag_print_orbital))
            ret.append(',\n    flag_get_orbital : ')
            ret.append(repr(self.flag_get_orbital))
            ret.append(',\n    flag_print_mag : ')
            ret.append(repr(self.flag_print_mag))
            ret.append(',\n    flag_print_single : ')
            ret.append(repr(self.flag_print_single))
            ret.append(',\n    flag_local_charge : ')
            ret.append(repr(self.flag_local_charge))
            ret.append(',\n    flag_collinear : ')
            ret.append(repr(self.flag_collinear))
            ret.append(',\n    flag_noncollinear : ')
            ret.append(repr(self.flag_noncollinear))
            ret.append(',\n    flag_soc : ')
            ret.append(repr(self.flag_soc))
            ret.append(',\n    flag_plus_u : ')
            ret.append(repr(self.flag_plus_u))
            ret.append(',\n    flag_scissor : ')
            ret.append(repr(self.flag_scissor))
            ret.append(',\n    flag_print_proj : ')
            ret.append(repr(self.flag_print_proj))
            ret.append(',\n    flag_print_proj_sum : ')
            ret.append(repr(self.flag_print_proj_sum))
            ret.append(',\n    flag_print_proj_atom : ')
            ret.append(repr(self.flag_print_proj_atom))
            ret.append(',\n    flag_plot_fit : ')
            ret.append(repr(self.flag_plot_fit))
            ret.append(',\n    flag_plot : ')
            ret.append(repr(self.flag_plot))
            ret.append(',\n    flag_get_band_order : ')
            ret.append(repr(self.flag_get_band_order))
            ret.append(',\n    flag_get_band_order_print_only : ')
            ret.append(repr(self.flag_get_band_order_print_only))
            ret.append(',\n    flag_use_weight : ')
            ret.append(repr(self.flag_use_weight))
            ret.append(',\n    flag_fit_orbital : ')
            ret.append(repr(self.flag_fit_orbital))
            ret.append(',\n    flag_fit_orbital_parse : ')
            ret.append(repr(self.flag_fit_orbital_parse))
            ret.append(',\n    flag_pso_with_lmdif : ')
            ret.append(repr(self.flag_pso_with_lmdif))
            ret.append(',\n    flag_get_unfold : ')
            ret.append(repr(self.flag_get_unfold))
            ret.append(',\n    ptol : ')
            ret.append(repr(self.ptol))
            ret.append(',\n    ftol : ')
            ret.append(repr(self.ftol))
            ret.append(',\n    fdiff : ')
            ret.append(repr(self.fdiff))
            ret.append(',\n    r_scissor : ')
            ret.append(repr(self.r_scissor))
            ret.append(',\n    nsystem : ')
            ret.append(repr(self.nsystem))
            ret.append(',\n    lmmax : ')
            ret.append(repr(self.lmmax))
            ret.append(',\n    miter : ')
            ret.append(repr(self.miter))
            ret.append(',\n    mxfit : ')
            ret.append(repr(self.mxfit))
            ret.append(',\n    nn_max : ')
            ret.append(repr(self.nn_max))
            ret.append(',\n    ispin : ')
            ret.append(repr(self.ispin))
            ret.append(',\n    ispinor : ')
            ret.append(repr(self.ispinor))
            ret.append(',\n    nspin : ')
            ret.append(repr(self.nspin))
            ret.append(',\n    nproj_sum : ')
            ret.append(repr(self.nproj_sum))
            ret.append(',\n    iverbose : ')
            ret.append(repr(self.iverbose))
            ret.append(',\n    i_scissor : ')
            ret.append(repr(self.i_scissor))
            ret.append(',\n    ls_type : ')
            ret.append(repr(self.ls_type))
            ret.append(',\n    fnamelog : ')
            ret.append(repr(self.fnamelog))
            ret.append(',\n    proj_atom : ')
            ret.append(repr(self.proj_atom))
            ret.append(',\n    proj_natom : ')
            ret.append(repr(self.proj_natom))
            ret.append(',\n    ifilenm : ')
            ret.append(repr(self.ifilenm))
            ret.append(',\n    title : ')
            ret.append(repr(self.title))
            ret.append(',\n    flag_get_total_energy : ')
            ret.append(repr(self.flag_get_total_energy))
            ret.append(',\n    electronic_temperature : ')
            ret.append(repr(self.electronic_temperature))
            ret.append(',\n    orbital_fit_smearing : ')
            ret.append(repr(self.orbital_fit_smearing))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("tbfitpy_mod.params_py")
    class params_py(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=params_py)
        
        
        Defined at tbfitpy_mod.f90 lines 68-106
        
        """
        def __init__(self, handle=None):
            """
            self = Params_Py()
            
            
            Defined at tbfitpy_mod.f90 lines 68-106
            
            
            Returns
            -------
            this : Params_Py
            	Object to be constructed
            
            
            Automatically generated constructor for params_py
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _tbfitpy_mod.f90wrap_params_py_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Params_Py
            
            
            Defined at tbfitpy_mod.f90 lines 68-106
            
            Parameters
            ----------
            this : Params_Py
            	Object to be destructed
            
            
            Automatically generated destructor for params_py
            """
            if self._alloc:
                _tbfitpy_mod.f90wrap_params_py_finalise(this=self._handle)
        
        @property
        def flag_set_param_const(self):
            """
            Element flag_set_param_const ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 69
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__flag_set_param_const(self._handle)
        
        @flag_set_param_const.setter
        def flag_set_param_const(self, flag_set_param_const):
            _tbfitpy_mod.f90wrap_params_py__set__flag_set_param_const(self._handle, \
                flag_set_param_const)
        
        @property
        def flag_pfile_index(self):
            """
            Element flag_pfile_index ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 70
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__flag_pfile_index(self._handle)
        
        @flag_pfile_index.setter
        def flag_pfile_index(self, flag_pfile_index):
            _tbfitpy_mod.f90wrap_params_py__set__flag_pfile_index(self._handle, \
                flag_pfile_index)
        
        @property
        def flag_use_overlap(self):
            """
            Element flag_use_overlap ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 71
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__flag_use_overlap(self._handle)
        
        @flag_use_overlap.setter
        def flag_use_overlap(self, flag_use_overlap):
            _tbfitpy_mod.f90wrap_params_py__set__flag_use_overlap(self._handle, \
                flag_use_overlap)
        
        @property
        def flag_slater_koster(self):
            """
            Element flag_slater_koster ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 72
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__flag_slater_koster(self._handle)
        
        @flag_slater_koster.setter
        def flag_slater_koster(self, flag_slater_koster):
            _tbfitpy_mod.f90wrap_params_py__set__flag_slater_koster(self._handle, \
                flag_slater_koster)
        
        @property
        def flag_nrl_slater_koster(self):
            """
            Element flag_nrl_slater_koster ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 73
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__flag_nrl_slater_koster(self._handle)
        
        @flag_nrl_slater_koster.setter
        def flag_nrl_slater_koster(self, flag_nrl_slater_koster):
            _tbfitpy_mod.f90wrap_params_py__set__flag_nrl_slater_koster(self._handle, \
                flag_nrl_slater_koster)
        
        @property
        def flag_fit_plain(self):
            """
            Element flag_fit_plain ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 74
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__flag_fit_plain(self._handle)
        
        @flag_fit_plain.setter
        def flag_fit_plain(self, flag_fit_plain):
            _tbfitpy_mod.f90wrap_params_py__set__flag_fit_plain(self._handle, \
                flag_fit_plain)
        
        @property
        def l_broaden(self):
            """
            Element l_broaden ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 75
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__l_broaden(self._handle)
        
        @l_broaden.setter
        def l_broaden(self, l_broaden):
            _tbfitpy_mod.f90wrap_params_py__set__l_broaden(self._handle, l_broaden)
        
        @property
        def pso_c1(self):
            """
            Element pso_c1 ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 76
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__pso_c1(self._handle)
        
        @pso_c1.setter
        def pso_c1(self, pso_c1):
            _tbfitpy_mod.f90wrap_params_py__set__pso_c1(self._handle, pso_c1)
        
        @property
        def pso_c2(self):
            """
            Element pso_c2 ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 77
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__pso_c2(self._handle)
        
        @pso_c2.setter
        def pso_c2(self, pso_c2):
            _tbfitpy_mod.f90wrap_params_py__set__pso_c2(self._handle, pso_c2)
        
        @property
        def pso_w(self):
            """
            Element pso_w ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 78
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__pso_w(self._handle)
        
        @pso_w.setter
        def pso_w(self, pso_w):
            _tbfitpy_mod.f90wrap_params_py__set__pso_w(self._handle, pso_w)
        
        @property
        def pso_max_noise_amplitude(self):
            """
            Element pso_max_noise_amplitude ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 79
            
            """
            return \
                _tbfitpy_mod.f90wrap_params_py__get__pso_max_noise_amplitude(self._handle)
        
        @pso_max_noise_amplitude.setter
        def pso_max_noise_amplitude(self, pso_max_noise_amplitude):
            _tbfitpy_mod.f90wrap_params_py__set__pso_max_noise_amplitude(self._handle, \
                pso_max_noise_amplitude)
        
        @property
        def pso_miter(self):
            """
            Element pso_miter ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 80
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__pso_miter(self._handle)
        
        @pso_miter.setter
        def pso_miter(self, pso_miter):
            _tbfitpy_mod.f90wrap_params_py__set__pso_miter(self._handle, pso_miter)
        
        @property
        def slater_koster_type(self):
            """
            Element slater_koster_type ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 81
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__slater_koster_type(self._handle)
        
        @slater_koster_type.setter
        def slater_koster_type(self, slater_koster_type):
            _tbfitpy_mod.f90wrap_params_py__set__slater_koster_type(self._handle, \
                slater_koster_type)
        
        @property
        def nparam(self):
            """
            Element nparam ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 82
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__nparam(self._handle)
        
        @nparam.setter
        def nparam(self, nparam):
            _tbfitpy_mod.f90wrap_params_py__set__nparam(self._handle, nparam)
        
        @property
        def nparam_const(self):
            """
            Element nparam_const ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 83
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__nparam_const(self._handle)
        
        @nparam_const.setter
        def nparam_const(self, nparam_const):
            _tbfitpy_mod.f90wrap_params_py__set__nparam_const(self._handle, nparam_const)
        
        @property
        def nparam_nrl(self):
            """
            Element nparam_nrl ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 84
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__nparam_nrl(self._handle)
        
        @nparam_nrl.setter
        def nparam_nrl(self, nparam_nrl):
            _tbfitpy_mod.f90wrap_params_py__set__nparam_nrl(self._handle, nparam_nrl)
        
        @property
        def nparam_nrl_free(self):
            """
            Element nparam_nrl_free ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 85
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__nparam_nrl_free(self._handle)
        
        @nparam_nrl_free.setter
        def nparam_nrl_free(self, nparam_nrl_free):
            _tbfitpy_mod.f90wrap_params_py__set__nparam_nrl_free(self._handle, \
                nparam_nrl_free)
        
        @property
        def param_nsub_max(self):
            """
            Element param_nsub_max ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 86
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__param_nsub_max(self._handle)
        
        @param_nsub_max.setter
        def param_nsub_max(self, param_nsub_max):
            _tbfitpy_mod.f90wrap_params_py__set__param_nsub_max(self._handle, \
                param_nsub_max)
        
        @property
        def nparam_free(self):
            """
            Element nparam_free ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 87
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__nparam_free(self._handle)
        
        @nparam_free.setter
        def nparam_free(self, nparam_free):
            _tbfitpy_mod.f90wrap_params_py__set__nparam_free(self._handle, nparam_free)
        
        @property
        def pso_nparticles(self):
            """
            Element pso_nparticles ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 88
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__pso_nparticles(self._handle)
        
        @pso_nparticles.setter
        def pso_nparticles(self, pso_nparticles):
            _tbfitpy_mod.f90wrap_params_py__set__pso_nparticles(self._handle, \
                pso_nparticles)
        
        @property
        def pfilenm(self):
            """
            Element pfilenm ftype=character(len=132) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 89
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__pfilenm(self._handle)
        
        @pfilenm.setter
        def pfilenm(self, pfilenm):
            _tbfitpy_mod.f90wrap_params_py__set__pfilenm(self._handle, pfilenm)
        
        @property
        def pfileoutnm(self):
            """
            Element pfileoutnm ftype=character(len=132) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 90
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__pfileoutnm(self._handle)
        
        @pfileoutnm.setter
        def pfileoutnm(self, pfileoutnm):
            _tbfitpy_mod.f90wrap_params_py__set__pfileoutnm(self._handle, pfileoutnm)
        
        @property
        def param(self):
            """
            Element param ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 91
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_params_py__array__param(self._handle)
            if array_handle in self._arrays:
                param = self._arrays[array_handle]
            else:
                param = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_params_py__array__param)
                self._arrays[array_handle] = param
            return param
        
        @param.setter
        def param(self, param):
            self.param[...] = param
        
        @property
        def param_best(self):
            """
            Element param_best ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 92
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_params_py__array__param_best(self._handle)
            if array_handle in self._arrays:
                param_best = self._arrays[array_handle]
            else:
                param_best = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_params_py__array__param_best)
                self._arrays[array_handle] = param_best
            return param_best
        
        @param_best.setter
        def param_best(self, param_best):
            self.param_best[...] = param_best
        
        @property
        def param_nrl(self):
            """
            Element param_nrl ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 93
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_params_py__array__param_nrl(self._handle)
            if array_handle in self._arrays:
                param_nrl = self._arrays[array_handle]
            else:
                param_nrl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_params_py__array__param_nrl)
                self._arrays[array_handle] = param_nrl
            return param_nrl
        
        @param_nrl.setter
        def param_nrl(self, param_nrl):
            self.param_nrl[...] = param_nrl
        
        @property
        def param_const(self):
            """
            Element param_const ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 94
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_params_py__array__param_const(self._handle)
            if array_handle in self._arrays:
                param_const = self._arrays[array_handle]
            else:
                param_const = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_params_py__array__param_const)
                self._arrays[array_handle] = param_const
            return param_const
        
        @param_const.setter
        def param_const(self, param_const):
            self.param_const[...] = param_const
        
        @property
        def param_const_nrl(self):
            """
            Element param_const_nrl ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 96
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_params_py__array__param_const_nrl(self._handle)
            if array_handle in self._arrays:
                param_const_nrl = self._arrays[array_handle]
            else:
                param_const_nrl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_params_py__array__param_const_nrl)
                self._arrays[array_handle] = param_const_nrl
            return param_const_nrl
        
        @param_const_nrl.setter
        def param_const_nrl(self, param_const_nrl):
            self.param_const_nrl[...] = param_const_nrl
        
        @property
        def param_nsub(self):
            """
            Element param_nsub ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 97
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_params_py__array__param_nsub(self._handle)
            if array_handle in self._arrays:
                param_nsub = self._arrays[array_handle]
            else:
                param_nsub = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_params_py__array__param_nsub)
                self._arrays[array_handle] = param_nsub
            return param_nsub
        
        @param_nsub.setter
        def param_nsub(self, param_nsub):
            self.param_nsub[...] = param_nsub
        
        @property
        def iparam_free(self):
            """
            Element iparam_free ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 98
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_params_py__array__iparam_free(self._handle)
            if array_handle in self._arrays:
                iparam_free = self._arrays[array_handle]
            else:
                iparam_free = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_params_py__array__iparam_free)
                self._arrays[array_handle] = iparam_free
            return iparam_free
        
        @iparam_free.setter
        def iparam_free(self, iparam_free):
            self.iparam_free[...] = iparam_free
        
        @property
        def iparam_free_nrl(self):
            """
            Element iparam_free_nrl ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 99
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_params_py__array__iparam_free_nrl(self._handle)
            if array_handle in self._arrays:
                iparam_free_nrl = self._arrays[array_handle]
            else:
                iparam_free_nrl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_params_py__array__iparam_free_nrl)
                self._arrays[array_handle] = iparam_free_nrl
            return iparam_free_nrl
        
        @iparam_free_nrl.setter
        def iparam_free_nrl(self, iparam_free_nrl):
            self.iparam_free_nrl[...] = iparam_free_nrl
        
        @property
        def param_name(self):
            """
            Element param_name ftype=character(len=40) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 100
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_params_py__array__param_name(self._handle)
            if array_handle in self._arrays:
                param_name = self._arrays[array_handle]
            else:
                param_name = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_params_py__array__param_name)
                self._arrays[array_handle] = param_name
            return param_name
        
        @param_name.setter
        def param_name(self, param_name):
            self.param_name[...] = param_name
        
        @property
        def c_const(self):
            """
            Element c_const ftype=character(len=40) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 101
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_params_py__array__c_const(self._handle)
            if array_handle in self._arrays:
                c_const = self._arrays[array_handle]
            else:
                c_const = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_params_py__array__c_const)
                self._arrays[array_handle] = c_const
            return c_const
        
        @c_const.setter
        def c_const(self, c_const):
            self.c_const[...] = c_const
        
        @property
        def pso_cost_history(self):
            """
            Element pso_cost_history ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 102
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_params_py__array__pso_cost_history(self._handle)
            if array_handle in self._arrays:
                pso_cost_history = self._arrays[array_handle]
            else:
                pso_cost_history = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_params_py__array__pso_cost_history)
                self._arrays[array_handle] = pso_cost_history
            return pso_cost_history
        
        @pso_cost_history.setter
        def pso_cost_history(self, pso_cost_history):
            self.pso_cost_history[...] = pso_cost_history
        
        @property
        def pso_cost_history_i(self):
            """
            Element pso_cost_history_i ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 103
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_params_py__array__pso_cost_history_i(self._handle)
            if array_handle in self._arrays:
                pso_cost_history_i = self._arrays[array_handle]
            else:
                pso_cost_history_i = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_params_py__array__pso_cost_history_i)
                self._arrays[array_handle] = pso_cost_history_i
            return pso_cost_history_i
        
        @pso_cost_history_i.setter
        def pso_cost_history_i(self, pso_cost_history_i):
            self.pso_cost_history_i[...] = pso_cost_history_i
        
        @property
        def pso_pbest_history(self):
            """
            Element pso_pbest_history ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 104
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_params_py__array__pso_pbest_history(self._handle)
            if array_handle in self._arrays:
                pso_pbest_history = self._arrays[array_handle]
            else:
                pso_pbest_history = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_params_py__array__pso_pbest_history)
                self._arrays[array_handle] = pso_pbest_history
            return pso_pbest_history
        
        @pso_pbest_history.setter
        def pso_pbest_history(self, pso_pbest_history):
            self.pso_pbest_history[...] = pso_pbest_history
        
        @property
        def cost_history(self):
            """
            Element cost_history ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 105
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_params_py__array__cost_history(self._handle)
            if array_handle in self._arrays:
                cost_history = self._arrays[array_handle]
            else:
                cost_history = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_params_py__array__cost_history)
                self._arrays[array_handle] = cost_history
            return cost_history
        
        @cost_history.setter
        def cost_history(self, cost_history):
            self.cost_history[...] = cost_history
        
        @property
        def niter(self):
            """
            Element niter ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 106
            
            """
            return _tbfitpy_mod.f90wrap_params_py__get__niter(self._handle)
        
        @niter.setter
        def niter(self, niter):
            _tbfitpy_mod.f90wrap_params_py__set__niter(self._handle, niter)
        
        def __str__(self):
            ret = ['<params_py>{\n']
            ret.append('    flag_set_param_const : ')
            ret.append(repr(self.flag_set_param_const))
            ret.append(',\n    flag_pfile_index : ')
            ret.append(repr(self.flag_pfile_index))
            ret.append(',\n    flag_use_overlap : ')
            ret.append(repr(self.flag_use_overlap))
            ret.append(',\n    flag_slater_koster : ')
            ret.append(repr(self.flag_slater_koster))
            ret.append(',\n    flag_nrl_slater_koster : ')
            ret.append(repr(self.flag_nrl_slater_koster))
            ret.append(',\n    flag_fit_plain : ')
            ret.append(repr(self.flag_fit_plain))
            ret.append(',\n    l_broaden : ')
            ret.append(repr(self.l_broaden))
            ret.append(',\n    pso_c1 : ')
            ret.append(repr(self.pso_c1))
            ret.append(',\n    pso_c2 : ')
            ret.append(repr(self.pso_c2))
            ret.append(',\n    pso_w : ')
            ret.append(repr(self.pso_w))
            ret.append(',\n    pso_max_noise_amplitude : ')
            ret.append(repr(self.pso_max_noise_amplitude))
            ret.append(',\n    pso_miter : ')
            ret.append(repr(self.pso_miter))
            ret.append(',\n    slater_koster_type : ')
            ret.append(repr(self.slater_koster_type))
            ret.append(',\n    nparam : ')
            ret.append(repr(self.nparam))
            ret.append(',\n    nparam_const : ')
            ret.append(repr(self.nparam_const))
            ret.append(',\n    nparam_nrl : ')
            ret.append(repr(self.nparam_nrl))
            ret.append(',\n    nparam_nrl_free : ')
            ret.append(repr(self.nparam_nrl_free))
            ret.append(',\n    param_nsub_max : ')
            ret.append(repr(self.param_nsub_max))
            ret.append(',\n    nparam_free : ')
            ret.append(repr(self.nparam_free))
            ret.append(',\n    pso_nparticles : ')
            ret.append(repr(self.pso_nparticles))
            ret.append(',\n    pfilenm : ')
            ret.append(repr(self.pfilenm))
            ret.append(',\n    pfileoutnm : ')
            ret.append(repr(self.pfileoutnm))
            ret.append(',\n    param : ')
            ret.append(repr(self.param))
            ret.append(',\n    param_best : ')
            ret.append(repr(self.param_best))
            ret.append(',\n    param_nrl : ')
            ret.append(repr(self.param_nrl))
            ret.append(',\n    param_const : ')
            ret.append(repr(self.param_const))
            ret.append(',\n    param_const_nrl : ')
            ret.append(repr(self.param_const_nrl))
            ret.append(',\n    param_nsub : ')
            ret.append(repr(self.param_nsub))
            ret.append(',\n    iparam_free : ')
            ret.append(repr(self.iparam_free))
            ret.append(',\n    iparam_free_nrl : ')
            ret.append(repr(self.iparam_free_nrl))
            ret.append(',\n    param_name : ')
            ret.append(repr(self.param_name))
            ret.append(',\n    c_const : ')
            ret.append(repr(self.c_const))
            ret.append(',\n    pso_cost_history : ')
            ret.append(repr(self.pso_cost_history))
            ret.append(',\n    pso_cost_history_i : ')
            ret.append(repr(self.pso_cost_history_i))
            ret.append(',\n    pso_pbest_history : ')
            ret.append(repr(self.pso_pbest_history))
            ret.append(',\n    cost_history : ')
            ret.append(repr(self.cost_history))
            ret.append(',\n    niter : ')
            ret.append(repr(self.niter))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("tbfitpy_mod.kpoints_py")
    class kpoints_py(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=kpoints_py)
        
        
        Defined at tbfitpy_mod.f90 lines 108-131
        
        """
        def __init__(self, handle=None):
            """
            self = Kpoints_Py()
            
            
            Defined at tbfitpy_mod.f90 lines 108-131
            
            
            Returns
            -------
            this : Kpoints_Py
            	Object to be constructed
            
            
            Automatically generated constructor for kpoints_py
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _tbfitpy_mod.f90wrap_kpoints_py_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Kpoints_Py
            
            
            Defined at tbfitpy_mod.f90 lines 108-131
            
            Parameters
            ----------
            this : Kpoints_Py
            	Object to be destructed
            
            
            Automatically generated destructor for kpoints_py
            """
            if self._alloc:
                _tbfitpy_mod.f90wrap_kpoints_py_finalise(this=self._handle)
        
        @property
        def flag_klinemode(self):
            """
            Element flag_klinemode ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 109
            
            """
            return _tbfitpy_mod.f90wrap_kpoints_py__get__flag_klinemode(self._handle)
        
        @flag_klinemode.setter
        def flag_klinemode(self, flag_klinemode):
            _tbfitpy_mod.f90wrap_kpoints_py__set__flag_klinemode(self._handle, \
                flag_klinemode)
        
        @property
        def flag_kgridmode(self):
            """
            Element flag_kgridmode ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 110
            
            """
            return _tbfitpy_mod.f90wrap_kpoints_py__get__flag_kgridmode(self._handle)
        
        @flag_kgridmode.setter
        def flag_kgridmode(self, flag_kgridmode):
            _tbfitpy_mod.f90wrap_kpoints_py__set__flag_kgridmode(self._handle, \
                flag_kgridmode)
        
        @property
        def flag_gamma(self):
            """
            Element flag_gamma ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 111
            
            """
            return _tbfitpy_mod.f90wrap_kpoints_py__get__flag_gamma(self._handle)
        
        @flag_gamma.setter
        def flag_gamma(self, flag_gamma):
            _tbfitpy_mod.f90wrap_kpoints_py__set__flag_gamma(self._handle, flag_gamma)
        
        @property
        def flag_reciprocal(self):
            """
            Element flag_reciprocal ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 112
            
            """
            return _tbfitpy_mod.f90wrap_kpoints_py__get__flag_reciprocal(self._handle)
        
        @flag_reciprocal.setter
        def flag_reciprocal(self, flag_reciprocal):
            _tbfitpy_mod.f90wrap_kpoints_py__set__flag_reciprocal(self._handle, \
                flag_reciprocal)
        
        @property
        def flag_cartesiank(self):
            """
            Element flag_cartesiank ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 113
            
            """
            return _tbfitpy_mod.f90wrap_kpoints_py__get__flag_cartesiank(self._handle)
        
        @flag_cartesiank.setter
        def flag_cartesiank(self, flag_cartesiank):
            _tbfitpy_mod.f90wrap_kpoints_py__set__flag_cartesiank(self._handle, \
                flag_cartesiank)
        
        @property
        def k_shift(self):
            """
            Element k_shift ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 114
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_kpoints_py__array__k_shift(self._handle)
            if array_handle in self._arrays:
                k_shift = self._arrays[array_handle]
            else:
                k_shift = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_kpoints_py__array__k_shift)
                self._arrays[array_handle] = k_shift
            return k_shift
        
        @k_shift.setter
        def k_shift(self, k_shift):
            self.k_shift[...] = k_shift
        
        @property
        def mysystem(self):
            """
            Element mysystem ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 115
            
            """
            return _tbfitpy_mod.f90wrap_kpoints_py__get__mysystem(self._handle)
        
        @mysystem.setter
        def mysystem(self, mysystem):
            _tbfitpy_mod.f90wrap_kpoints_py__set__mysystem(self._handle, mysystem)
        
        @property
        def nkpoint(self):
            """
            Element nkpoint ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 116
            
            """
            return _tbfitpy_mod.f90wrap_kpoints_py__get__nkpoint(self._handle)
        
        @nkpoint.setter
        def nkpoint(self, nkpoint):
            _tbfitpy_mod.f90wrap_kpoints_py__set__nkpoint(self._handle, nkpoint)
        
        @property
        def nline(self):
            """
            Element nline ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 116
            
            """
            return _tbfitpy_mod.f90wrap_kpoints_py__get__nline(self._handle)
        
        @nline.setter
        def nline(self, nline):
            _tbfitpy_mod.f90wrap_kpoints_py__set__nline(self._handle, nline)
        
        @property
        def n_ndiv(self):
            """
            Element n_ndiv ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 117
            
            """
            return _tbfitpy_mod.f90wrap_kpoints_py__get__n_ndiv(self._handle)
        
        @n_ndiv.setter
        def n_ndiv(self, n_ndiv):
            _tbfitpy_mod.f90wrap_kpoints_py__set__n_ndiv(self._handle, n_ndiv)
        
        @property
        def idiv_mode(self):
            """
            Element idiv_mode ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 118
            
            """
            return _tbfitpy_mod.f90wrap_kpoints_py__get__idiv_mode(self._handle)
        
        @idiv_mode.setter
        def idiv_mode(self, idiv_mode):
            _tbfitpy_mod.f90wrap_kpoints_py__set__idiv_mode(self._handle, idiv_mode)
        
        @property
        def kreduce(self):
            """
            Element kreduce ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 119
            
            """
            return _tbfitpy_mod.f90wrap_kpoints_py__get__kreduce(self._handle)
        
        @kreduce.setter
        def kreduce(self, kreduce):
            _tbfitpy_mod.f90wrap_kpoints_py__set__kreduce(self._handle, kreduce)
        
        @property
        def kfilenm(self):
            """
            Element kfilenm ftype=character(len=132) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 120
            
            """
            return _tbfitpy_mod.f90wrap_kpoints_py__get__kfilenm(self._handle)
        
        @kfilenm.setter
        def kfilenm(self, kfilenm):
            _tbfitpy_mod.f90wrap_kpoints_py__set__kfilenm(self._handle, kfilenm)
        
        @property
        def ribbon_kfilenm(self):
            """
            Element ribbon_kfilenm ftype=character(len=132) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 121
            
            """
            return _tbfitpy_mod.f90wrap_kpoints_py__get__ribbon_kfilenm(self._handle)
        
        @ribbon_kfilenm.setter
        def ribbon_kfilenm(self, ribbon_kfilenm):
            _tbfitpy_mod.f90wrap_kpoints_py__set__ribbon_kfilenm(self._handle, \
                ribbon_kfilenm)
        
        @property
        def kline_type(self):
            """
            Element kline_type ftype=character(len=8) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 122
            
            """
            return _tbfitpy_mod.f90wrap_kpoints_py__get__kline_type(self._handle)
        
        @kline_type.setter
        def kline_type(self, kline_type):
            _tbfitpy_mod.f90wrap_kpoints_py__set__kline_type(self._handle, kline_type)
        
        @property
        def kunit(self):
            """
            Element kunit ftype=character(len=1) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 123
            
            """
            return _tbfitpy_mod.f90wrap_kpoints_py__get__kunit(self._handle)
        
        @kunit.setter
        def kunit(self, kunit):
            _tbfitpy_mod.f90wrap_kpoints_py__set__kunit(self._handle, kunit)
        
        @property
        def kpoint(self):
            """
            Element kpoint ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 124
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_kpoints_py__array__kpoint(self._handle)
            if array_handle in self._arrays:
                kpoint = self._arrays[array_handle]
            else:
                kpoint = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_kpoints_py__array__kpoint)
                self._arrays[array_handle] = kpoint
            return kpoint
        
        @kpoint.setter
        def kpoint(self, kpoint):
            self.kpoint[...] = kpoint
        
        @property
        def kline(self):
            """
            Element kline ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 125
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_kpoints_py__array__kline(self._handle)
            if array_handle in self._arrays:
                kline = self._arrays[array_handle]
            else:
                kline = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_kpoints_py__array__kline)
                self._arrays[array_handle] = kline
            return kline
        
        @kline.setter
        def kline(self, kline):
            self.kline[...] = kline
        
        @property
        def kpoint_reci(self):
            """
            Element kpoint_reci ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 126
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_kpoints_py__array__kpoint_reci(self._handle)
            if array_handle in self._arrays:
                kpoint_reci = self._arrays[array_handle]
            else:
                kpoint_reci = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_kpoints_py__array__kpoint_reci)
                self._arrays[array_handle] = kpoint_reci
            return kpoint_reci
        
        @kpoint_reci.setter
        def kpoint_reci(self, kpoint_reci):
            self.kpoint_reci[...] = kpoint_reci
        
        @property
        def ndiv(self):
            """
            Element ndiv ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 127
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_kpoints_py__array__ndiv(self._handle)
            if array_handle in self._arrays:
                ndiv = self._arrays[array_handle]
            else:
                ndiv = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_kpoints_py__array__ndiv)
                self._arrays[array_handle] = ndiv
            return ndiv
        
        @ndiv.setter
        def ndiv(self, ndiv):
            self.ndiv[...] = ndiv
        
        @property
        def k_name(self):
            """
            Element k_name ftype=character(len=8) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 128
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_kpoints_py__array__k_name(self._handle)
            if array_handle in self._arrays:
                k_name = self._arrays[array_handle]
            else:
                k_name = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_kpoints_py__array__k_name)
                self._arrays[array_handle] = k_name
            return k_name
        
        @k_name.setter
        def k_name(self, k_name):
            self.k_name[...] = k_name
        
        @property
        def k_name2(self):
            """
            Element k_name2 ftype=character(len=8) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 129
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_kpoints_py__array__k_name2(self._handle)
            if array_handle in self._arrays:
                k_name2 = self._arrays[array_handle]
            else:
                k_name2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_kpoints_py__array__k_name2)
                self._arrays[array_handle] = k_name2
            return k_name2
        
        @k_name2.setter
        def k_name2(self, k_name2):
            self.k_name2[...] = k_name2
        
        @property
        def k_name_index(self):
            """
            Element k_name_index ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 130
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_kpoints_py__array__k_name_index(self._handle)
            if array_handle in self._arrays:
                k_name_index = self._arrays[array_handle]
            else:
                k_name_index = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_kpoints_py__array__k_name_index)
                self._arrays[array_handle] = k_name_index
            return k_name_index
        
        @k_name_index.setter
        def k_name_index(self, k_name_index):
            self.k_name_index[...] = k_name_index
        
        @property
        def kdist(self):
            """
            Element kdist ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 131
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_kpoints_py__array__kdist(self._handle)
            if array_handle in self._arrays:
                kdist = self._arrays[array_handle]
            else:
                kdist = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_kpoints_py__array__kdist)
                self._arrays[array_handle] = kdist
            return kdist
        
        @kdist.setter
        def kdist(self, kdist):
            self.kdist[...] = kdist
        
        def __str__(self):
            ret = ['<kpoints_py>{\n']
            ret.append('    flag_klinemode : ')
            ret.append(repr(self.flag_klinemode))
            ret.append(',\n    flag_kgridmode : ')
            ret.append(repr(self.flag_kgridmode))
            ret.append(',\n    flag_gamma : ')
            ret.append(repr(self.flag_gamma))
            ret.append(',\n    flag_reciprocal : ')
            ret.append(repr(self.flag_reciprocal))
            ret.append(',\n    flag_cartesiank : ')
            ret.append(repr(self.flag_cartesiank))
            ret.append(',\n    k_shift : ')
            ret.append(repr(self.k_shift))
            ret.append(',\n    mysystem : ')
            ret.append(repr(self.mysystem))
            ret.append(',\n    nkpoint : ')
            ret.append(repr(self.nkpoint))
            ret.append(',\n    nline : ')
            ret.append(repr(self.nline))
            ret.append(',\n    n_ndiv : ')
            ret.append(repr(self.n_ndiv))
            ret.append(',\n    idiv_mode : ')
            ret.append(repr(self.idiv_mode))
            ret.append(',\n    kreduce : ')
            ret.append(repr(self.kreduce))
            ret.append(',\n    kfilenm : ')
            ret.append(repr(self.kfilenm))
            ret.append(',\n    ribbon_kfilenm : ')
            ret.append(repr(self.ribbon_kfilenm))
            ret.append(',\n    kline_type : ')
            ret.append(repr(self.kline_type))
            ret.append(',\n    kunit : ')
            ret.append(repr(self.kunit))
            ret.append(',\n    kpoint : ')
            ret.append(repr(self.kpoint))
            ret.append(',\n    kline : ')
            ret.append(repr(self.kline))
            ret.append(',\n    kpoint_reci : ')
            ret.append(repr(self.kpoint_reci))
            ret.append(',\n    ndiv : ')
            ret.append(repr(self.ndiv))
            ret.append(',\n    k_name : ')
            ret.append(repr(self.k_name))
            ret.append(',\n    k_name2 : ')
            ret.append(repr(self.k_name2))
            ret.append(',\n    k_name_index : ')
            ret.append(repr(self.k_name_index))
            ret.append(',\n    kdist : ')
            ret.append(repr(self.kdist))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("tbfitpy_mod.weight_py")
    class weight_py(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=weight_py)
        
        
        Defined at tbfitpy_mod.f90 lines 133-145
        
        """
        def __init__(self, handle=None):
            """
            self = Weight_Py()
            
            
            Defined at tbfitpy_mod.f90 lines 133-145
            
            
            Returns
            -------
            this : Weight_Py
            	Object to be constructed
            
            
            Automatically generated constructor for weight_py
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _tbfitpy_mod.f90wrap_weight_py_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Weight_Py
            
            
            Defined at tbfitpy_mod.f90 lines 133-145
            
            Parameters
            ----------
            this : Weight_Py
            	Object to be destructed
            
            
            Automatically generated destructor for weight_py
            """
            if self._alloc:
                _tbfitpy_mod.f90wrap_weight_py_finalise(this=self._handle)
        
        @property
        def flag_weight_default(self):
            """
            Element flag_weight_default ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 134
            
            """
            return _tbfitpy_mod.f90wrap_weight_py__get__flag_weight_default(self._handle)
        
        @flag_weight_default.setter
        def flag_weight_default(self, flag_weight_default):
            _tbfitpy_mod.f90wrap_weight_py__set__flag_weight_default(self._handle, \
                flag_weight_default)
        
        @property
        def flag_weight_orb(self):
            """
            Element flag_weight_orb ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 135
            
            """
            return _tbfitpy_mod.f90wrap_weight_py__get__flag_weight_orb(self._handle)
        
        @flag_weight_orb.setter
        def flag_weight_orb(self, flag_weight_orb):
            _tbfitpy_mod.f90wrap_weight_py__set__flag_weight_orb(self._handle, \
                flag_weight_orb)
        
        @property
        def efile_ef(self):
            """
            Element efile_ef ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 136
            
            """
            return _tbfitpy_mod.f90wrap_weight_py__get__efile_ef(self._handle)
        
        @efile_ef.setter
        def efile_ef(self, efile_ef):
            _tbfitpy_mod.f90wrap_weight_py__set__efile_ef(self._handle, efile_ef)
        
        @property
        def mysystem(self):
            """
            Element mysystem ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 137
            
            """
            return _tbfitpy_mod.f90wrap_weight_py__get__mysystem(self._handle)
        
        @mysystem.setter
        def mysystem(self, mysystem):
            _tbfitpy_mod.f90wrap_weight_py__set__mysystem(self._handle, mysystem)
        
        @property
        def itarget_e_start(self):
            """
            Element itarget_e_start ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 138
            
            """
            return _tbfitpy_mod.f90wrap_weight_py__get__itarget_e_start(self._handle)
        
        @itarget_e_start.setter
        def itarget_e_start(self, itarget_e_start):
            _tbfitpy_mod.f90wrap_weight_py__set__itarget_e_start(self._handle, \
                itarget_e_start)
        
        @property
        def read_energy_column_index(self):
            """
            Element read_energy_column_index ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 139
            
            """
            return \
                _tbfitpy_mod.f90wrap_weight_py__get__read_energy_column_index(self._handle)
        
        @read_energy_column_index.setter
        def read_energy_column_index(self, read_energy_column_index):
            _tbfitpy_mod.f90wrap_weight_py__set__read_energy_column_index(self._handle, \
                read_energy_column_index)
        
        @property
        def read_energy_column_index_dn(self):
            """
            Element read_energy_column_index_dn ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 140
            
            """
            return \
                _tbfitpy_mod.f90wrap_weight_py__get__read_energy_column_index_dn(self._handle)
        
        @read_energy_column_index_dn.setter
        def read_energy_column_index_dn(self, read_energy_column_index_dn):
            _tbfitpy_mod.f90wrap_weight_py__set__read_energy_column_index_dn(self._handle, \
                read_energy_column_index_dn)
        
        @property
        def ie_cutoff(self):
            """
            Element ie_cutoff ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 141
            
            """
            return _tbfitpy_mod.f90wrap_weight_py__get__ie_cutoff(self._handle)
        
        @ie_cutoff.setter
        def ie_cutoff(self, ie_cutoff):
            _tbfitpy_mod.f90wrap_weight_py__set__ie_cutoff(self._handle, ie_cutoff)
        
        @property
        def efilenmu(self):
            """
            Element efilenmu ftype=character(len=132) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 142
            
            """
            return _tbfitpy_mod.f90wrap_weight_py__get__efilenmu(self._handle)
        
        @efilenmu.setter
        def efilenmu(self, efilenmu):
            _tbfitpy_mod.f90wrap_weight_py__set__efilenmu(self._handle, efilenmu)
        
        @property
        def efilenmd(self):
            """
            Element efilenmd ftype=character(len=132) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 143
            
            """
            return _tbfitpy_mod.f90wrap_weight_py__get__efilenmd(self._handle)
        
        @efilenmd.setter
        def efilenmd(self, efilenmd):
            _tbfitpy_mod.f90wrap_weight_py__set__efilenmd(self._handle, efilenmd)
        
        @property
        def efile_type(self):
            """
            Element efile_type ftype=character(len=16) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 144
            
            """
            return _tbfitpy_mod.f90wrap_weight_py__get__efile_type(self._handle)
        
        @efile_type.setter
        def efile_type(self, efile_type):
            _tbfitpy_mod.f90wrap_weight_py__set__efile_type(self._handle, efile_type)
        
        @property
        def wt(self):
            """
            Element wt ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 145
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_weight_py__array__wt(self._handle)
            if array_handle in self._arrays:
                wt = self._arrays[array_handle]
            else:
                wt = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_weight_py__array__wt)
                self._arrays[array_handle] = wt
            return wt
        
        @wt.setter
        def wt(self, wt):
            self.wt[...] = wt
        
        def __str__(self):
            ret = ['<weight_py>{\n']
            ret.append('    flag_weight_default : ')
            ret.append(repr(self.flag_weight_default))
            ret.append(',\n    flag_weight_orb : ')
            ret.append(repr(self.flag_weight_orb))
            ret.append(',\n    efile_ef : ')
            ret.append(repr(self.efile_ef))
            ret.append(',\n    mysystem : ')
            ret.append(repr(self.mysystem))
            ret.append(',\n    itarget_e_start : ')
            ret.append(repr(self.itarget_e_start))
            ret.append(',\n    read_energy_column_index : ')
            ret.append(repr(self.read_energy_column_index))
            ret.append(',\n    read_energy_column_index_dn : ')
            ret.append(repr(self.read_energy_column_index_dn))
            ret.append(',\n    ie_cutoff : ')
            ret.append(repr(self.ie_cutoff))
            ret.append(',\n    efilenmu : ')
            ret.append(repr(self.efilenmu))
            ret.append(',\n    efilenmd : ')
            ret.append(repr(self.efilenmd))
            ret.append(',\n    efile_type : ')
            ret.append(repr(self.efile_type))
            ret.append(',\n    wt : ')
            ret.append(repr(self.wt))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("tbfitpy_mod.poscar_py")
    class poscar_py(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=poscar_py)
        
        
        Defined at tbfitpy_mod.f90 lines 147-171
        
        """
        def __init__(self, handle=None):
            """
            self = Poscar_Py()
            
            
            Defined at tbfitpy_mod.f90 lines 147-171
            
            
            Returns
            -------
            this : Poscar_Py
            	Object to be constructed
            
            
            Automatically generated constructor for poscar_py
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _tbfitpy_mod.f90wrap_poscar_py_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Poscar_Py
            
            
            Defined at tbfitpy_mod.f90 lines 147-171
            
            Parameters
            ----------
            this : Poscar_Py
            	Object to be destructed
            
            
            Automatically generated destructor for poscar_py
            """
            if self._alloc:
                _tbfitpy_mod.f90wrap_poscar_py_finalise(this=self._handle)
        
        @property
        def flag_selective(self):
            """
            Element flag_selective ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 148
            
            """
            return _tbfitpy_mod.f90wrap_poscar_py__get__flag_selective(self._handle)
        
        @flag_selective.setter
        def flag_selective(self, flag_selective):
            _tbfitpy_mod.f90wrap_poscar_py__set__flag_selective(self._handle, \
                flag_selective)
        
        @property
        def flag_direct(self):
            """
            Element flag_direct ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 149
            
            """
            return _tbfitpy_mod.f90wrap_poscar_py__get__flag_direct(self._handle)
        
        @flag_direct.setter
        def flag_direct(self, flag_direct):
            _tbfitpy_mod.f90wrap_poscar_py__set__flag_direct(self._handle, flag_direct)
        
        @property
        def flag_cartesian(self):
            """
            Element flag_cartesian ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 150
            
            """
            return _tbfitpy_mod.f90wrap_poscar_py__get__flag_cartesian(self._handle)
        
        @flag_cartesian.setter
        def flag_cartesian(self, flag_cartesian):
            _tbfitpy_mod.f90wrap_poscar_py__set__flag_cartesian(self._handle, \
                flag_cartesian)
        
        @property
        def mysystem(self):
            """
            Element mysystem ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 151
            
            """
            return _tbfitpy_mod.f90wrap_poscar_py__get__mysystem(self._handle)
        
        @mysystem.setter
        def mysystem(self, mysystem):
            _tbfitpy_mod.f90wrap_poscar_py__set__mysystem(self._handle, mysystem)
        
        @property
        def neig(self):
            """
            Element neig ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 152
            
            """
            return _tbfitpy_mod.f90wrap_poscar_py__get__neig(self._handle)
        
        @neig.setter
        def neig(self, neig):
            _tbfitpy_mod.f90wrap_poscar_py__set__neig(self._handle, neig)
        
        @property
        def neig_total(self):
            """
            Element neig_total ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 153
            
            """
            return _tbfitpy_mod.f90wrap_poscar_py__get__neig_total(self._handle)
        
        @neig_total.setter
        def neig_total(self, neig_total):
            _tbfitpy_mod.f90wrap_poscar_py__set__neig_total(self._handle, neig_total)
        
        @property
        def neig_target(self):
            """
            Element neig_target ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 154
            
            """
            return _tbfitpy_mod.f90wrap_poscar_py__get__neig_target(self._handle)
        
        @neig_target.setter
        def neig_target(self, neig_target):
            _tbfitpy_mod.f90wrap_poscar_py__set__neig_target(self._handle, neig_target)
        
        @property
        def nbasis(self):
            """
            Element nbasis ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 155
            
            """
            return _tbfitpy_mod.f90wrap_poscar_py__get__nbasis(self._handle)
        
        @nbasis.setter
        def nbasis(self, nbasis):
            _tbfitpy_mod.f90wrap_poscar_py__set__nbasis(self._handle, nbasis)
        
        @property
        def neig_eff(self):
            """
            Element neig_eff ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 156
            
            """
            return _tbfitpy_mod.f90wrap_poscar_py__get__neig_eff(self._handle)
        
        @neig_eff.setter
        def neig_eff(self, neig_eff):
            _tbfitpy_mod.f90wrap_poscar_py__set__neig_eff(self._handle, neig_eff)
        
        @property
        def init_erange(self):
            """
            Element init_erange ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 157
            
            """
            return _tbfitpy_mod.f90wrap_poscar_py__get__init_erange(self._handle)
        
        @init_erange.setter
        def init_erange(self, init_erange):
            _tbfitpy_mod.f90wrap_poscar_py__set__init_erange(self._handle, init_erange)
        
        @property
        def fina_erange(self):
            """
            Element fina_erange ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 158
            
            """
            return _tbfitpy_mod.f90wrap_poscar_py__get__fina_erange(self._handle)
        
        @fina_erange.setter
        def fina_erange(self, fina_erange):
            _tbfitpy_mod.f90wrap_poscar_py__set__fina_erange(self._handle, fina_erange)
        
        @property
        def nband(self):
            """
            Element nband ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 159
            
            """
            return _tbfitpy_mod.f90wrap_poscar_py__get__nband(self._handle)
        
        @nband.setter
        def nband(self, nband):
            _tbfitpy_mod.f90wrap_poscar_py__set__nband(self._handle, nband)
        
        @property
        def n_spec(self):
            """
            Element n_spec ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 160
            
            """
            return _tbfitpy_mod.f90wrap_poscar_py__get__n_spec(self._handle)
        
        @n_spec.setter
        def n_spec(self, n_spec):
            _tbfitpy_mod.f90wrap_poscar_py__set__n_spec(self._handle, n_spec)
        
        @property
        def n_atom(self):
            """
            Element n_atom ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 161
            
            """
            return _tbfitpy_mod.f90wrap_poscar_py__get__n_atom(self._handle)
        
        @n_atom.setter
        def n_atom(self, n_atom):
            _tbfitpy_mod.f90wrap_poscar_py__set__n_atom(self._handle, n_atom)
        
        @property
        def max_orb(self):
            """
            Element max_orb ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 162
            
            """
            return _tbfitpy_mod.f90wrap_poscar_py__get__max_orb(self._handle)
        
        @max_orb.setter
        def max_orb(self, max_orb):
            _tbfitpy_mod.f90wrap_poscar_py__set__max_orb(self._handle, max_orb)
        
        @property
        def title(self):
            """
            Element title ftype=character(len=132) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 163
            
            """
            return _tbfitpy_mod.f90wrap_poscar_py__get__title(self._handle)
        
        @title.setter
        def title(self, title):
            _tbfitpy_mod.f90wrap_poscar_py__set__title(self._handle, title)
        
        @property
        def gfilenm(self):
            """
            Element gfilenm ftype=character(len=132) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 164
            
            """
            return _tbfitpy_mod.f90wrap_poscar_py__get__gfilenm(self._handle)
        
        @gfilenm.setter
        def gfilenm(self, gfilenm):
            _tbfitpy_mod.f90wrap_poscar_py__set__gfilenm(self._handle, gfilenm)
        
        @property
        def system_name(self):
            """
            Element system_name ftype=character(len=40) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 165
            
            """
            return _tbfitpy_mod.f90wrap_poscar_py__get__system_name(self._handle)
        
        @system_name.setter
        def system_name(self, system_name):
            _tbfitpy_mod.f90wrap_poscar_py__set__system_name(self._handle, system_name)
        
        @property
        def a_scale(self):
            """
            Element a_scale ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 166
            
            """
            return _tbfitpy_mod.f90wrap_poscar_py__get__a_scale(self._handle)
        
        @a_scale.setter
        def a_scale(self, a_scale):
            _tbfitpy_mod.f90wrap_poscar_py__set__a_scale(self._handle, a_scale)
        
        @property
        def a_latt(self):
            """
            Element a_latt ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 167
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_poscar_py__array__a_latt(self._handle)
            if array_handle in self._arrays:
                a_latt = self._arrays[array_handle]
            else:
                a_latt = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_poscar_py__array__a_latt)
                self._arrays[array_handle] = a_latt
            return a_latt
        
        @a_latt.setter
        def a_latt(self, a_latt):
            self.a_latt[...] = a_latt
        
        @property
        def b_latt(self):
            """
            Element b_latt ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 168
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_poscar_py__array__b_latt(self._handle)
            if array_handle in self._arrays:
                b_latt = self._arrays[array_handle]
            else:
                b_latt = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_poscar_py__array__b_latt)
                self._arrays[array_handle] = b_latt
            return b_latt
        
        @b_latt.setter
        def b_latt(self, b_latt):
            self.b_latt[...] = b_latt
        
        @property
        def nelect(self):
            """
            Element nelect ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 169
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_poscar_py__array__nelect(self._handle)
            if array_handle in self._arrays:
                nelect = self._arrays[array_handle]
            else:
                nelect = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_poscar_py__array__nelect)
                self._arrays[array_handle] = nelect
            return nelect
        
        @nelect.setter
        def nelect(self, nelect):
            self.nelect[...] = nelect
        
        @property
        def n_orbital(self):
            """
            Element n_orbital ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 170
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_poscar_py__array__n_orbital(self._handle)
            if array_handle in self._arrays:
                n_orbital = self._arrays[array_handle]
            else:
                n_orbital = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_poscar_py__array__n_orbital)
                self._arrays[array_handle] = n_orbital
            return n_orbital
        
        @n_orbital.setter
        def n_orbital(self, n_orbital):
            self.n_orbital[...] = n_orbital
        
        @property
        def orb_index(self):
            """
            Element orb_index ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 171
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_poscar_py__array__orb_index(self._handle)
            if array_handle in self._arrays:
                orb_index = self._arrays[array_handle]
            else:
                orb_index = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_poscar_py__array__orb_index)
                self._arrays[array_handle] = orb_index
            return orb_index
        
        @orb_index.setter
        def orb_index(self, orb_index):
            self.orb_index[...] = orb_index
        
        def __str__(self):
            ret = ['<poscar_py>{\n']
            ret.append('    flag_selective : ')
            ret.append(repr(self.flag_selective))
            ret.append(',\n    flag_direct : ')
            ret.append(repr(self.flag_direct))
            ret.append(',\n    flag_cartesian : ')
            ret.append(repr(self.flag_cartesian))
            ret.append(',\n    mysystem : ')
            ret.append(repr(self.mysystem))
            ret.append(',\n    neig : ')
            ret.append(repr(self.neig))
            ret.append(',\n    neig_total : ')
            ret.append(repr(self.neig_total))
            ret.append(',\n    neig_target : ')
            ret.append(repr(self.neig_target))
            ret.append(',\n    nbasis : ')
            ret.append(repr(self.nbasis))
            ret.append(',\n    neig_eff : ')
            ret.append(repr(self.neig_eff))
            ret.append(',\n    init_erange : ')
            ret.append(repr(self.init_erange))
            ret.append(',\n    fina_erange : ')
            ret.append(repr(self.fina_erange))
            ret.append(',\n    nband : ')
            ret.append(repr(self.nband))
            ret.append(',\n    n_spec : ')
            ret.append(repr(self.n_spec))
            ret.append(',\n    n_atom : ')
            ret.append(repr(self.n_atom))
            ret.append(',\n    max_orb : ')
            ret.append(repr(self.max_orb))
            ret.append(',\n    title : ')
            ret.append(repr(self.title))
            ret.append(',\n    gfilenm : ')
            ret.append(repr(self.gfilenm))
            ret.append(',\n    system_name : ')
            ret.append(repr(self.system_name))
            ret.append(',\n    a_scale : ')
            ret.append(repr(self.a_scale))
            ret.append(',\n    a_latt : ')
            ret.append(repr(self.a_latt))
            ret.append(',\n    b_latt : ')
            ret.append(repr(self.b_latt))
            ret.append(',\n    nelect : ')
            ret.append(repr(self.nelect))
            ret.append(',\n    n_orbital : ')
            ret.append(repr(self.n_orbital))
            ret.append(',\n    orb_index : ')
            ret.append(repr(self.orb_index))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("tbfitpy_mod.hopping_py")
    class hopping_py(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=hopping_py)
        
        
        Defined at tbfitpy_mod.f90 lines 173-217
        
        """
        def __init__(self, handle=None):
            """
            self = Hopping_Py()
            
            
            Defined at tbfitpy_mod.f90 lines 173-217
            
            
            Returns
            -------
            this : Hopping_Py
            	Object to be constructed
            
            
            Automatically generated constructor for hopping_py
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _tbfitpy_mod.f90wrap_hopping_py_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Hopping_Py
            
            
            Defined at tbfitpy_mod.f90 lines 173-217
            
            Parameters
            ----------
            this : Hopping_Py
            	Object to be destructed
            
            
            Automatically generated destructor for hopping_py
            """
            if self._alloc:
                _tbfitpy_mod.f90wrap_hopping_py_finalise(this=self._handle)
        
        @property
        def flag_efield(self):
            """
            Element flag_efield ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 174
            
            """
            return _tbfitpy_mod.f90wrap_hopping_py__get__flag_efield(self._handle)
        
        @flag_efield.setter
        def flag_efield(self, flag_efield):
            _tbfitpy_mod.f90wrap_hopping_py__set__flag_efield(self._handle, flag_efield)
        
        @property
        def flag_efield_frac(self):
            """
            Element flag_efield_frac ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 175
            
            """
            return _tbfitpy_mod.f90wrap_hopping_py__get__flag_efield_frac(self._handle)
        
        @flag_efield_frac.setter
        def flag_efield_frac(self, flag_efield_frac):
            _tbfitpy_mod.f90wrap_hopping_py__set__flag_efield_frac(self._handle, \
                flag_efield_frac)
        
        @property
        def flag_efield_cart(self):
            """
            Element flag_efield_cart ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 176
            
            """
            return _tbfitpy_mod.f90wrap_hopping_py__get__flag_efield_cart(self._handle)
        
        @flag_efield_cart.setter
        def flag_efield_cart(self, flag_efield_cart):
            _tbfitpy_mod.f90wrap_hopping_py__set__flag_efield_cart(self._handle, \
                flag_efield_cart)
        
        @property
        def onsite_tolerance(self):
            """
            Element onsite_tolerance ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 177
            
            """
            return _tbfitpy_mod.f90wrap_hopping_py__get__onsite_tolerance(self._handle)
        
        @onsite_tolerance.setter
        def onsite_tolerance(self, onsite_tolerance):
            _tbfitpy_mod.f90wrap_hopping_py__set__onsite_tolerance(self._handle, \
                onsite_tolerance)
        
        @property
        def efield(self):
            """
            Element efield ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 178
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__efield(self._handle)
            if array_handle in self._arrays:
                efield = self._arrays[array_handle]
            else:
                efield = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__efield)
                self._arrays[array_handle] = efield
            return efield
        
        @efield.setter
        def efield(self, efield):
            self.efield[...] = efield
        
        @property
        def efield_origin(self):
            """
            Element efield_origin ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 179
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__efield_origin(self._handle)
            if array_handle in self._arrays:
                efield_origin = self._arrays[array_handle]
            else:
                efield_origin = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__efield_origin)
                self._arrays[array_handle] = efield_origin
            return efield_origin
        
        @efield_origin.setter
        def efield_origin(self, efield_origin):
            self.efield_origin[...] = efield_origin
        
        @property
        def efield_origin_cart(self):
            """
            Element efield_origin_cart ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 180
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__efield_origin_cart(self._handle)
            if array_handle in self._arrays:
                efield_origin_cart = self._arrays[array_handle]
            else:
                efield_origin_cart = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__efield_origin_cart)
                self._arrays[array_handle] = efield_origin_cart
            return efield_origin_cart
        
        @efield_origin_cart.setter
        def efield_origin_cart(self, efield_origin_cart):
            self.efield_origin_cart[...] = efield_origin_cart
        
        @property
        def mysystem(self):
            """
            Element mysystem ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 181
            
            """
            return _tbfitpy_mod.f90wrap_hopping_py__get__mysystem(self._handle)
        
        @mysystem.setter
        def mysystem(self, mysystem):
            _tbfitpy_mod.f90wrap_hopping_py__set__mysystem(self._handle, mysystem)
        
        @property
        def n_neighbor(self):
            """
            Element n_neighbor ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 182
            
            """
            return _tbfitpy_mod.f90wrap_hopping_py__get__n_neighbor(self._handle)
        
        @n_neighbor.setter
        def n_neighbor(self, n_neighbor):
            _tbfitpy_mod.f90wrap_hopping_py__set__n_neighbor(self._handle, n_neighbor)
        
        @property
        def max_nn_pair(self):
            """
            Element max_nn_pair ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 183
            
            """
            return _tbfitpy_mod.f90wrap_hopping_py__get__max_nn_pair(self._handle)
        
        @max_nn_pair.setter
        def max_nn_pair(self, max_nn_pair):
            _tbfitpy_mod.f90wrap_hopping_py__set__max_nn_pair(self._handle, max_nn_pair)
        
        @property
        def nnfilenm(self):
            """
            Element nnfilenm ftype=character(len=132) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 184
            
            """
            return _tbfitpy_mod.f90wrap_hopping_py__get__nnfilenm(self._handle)
        
        @nnfilenm.setter
        def nnfilenm(self, nnfilenm):
            _tbfitpy_mod.f90wrap_hopping_py__set__nnfilenm(self._handle, nnfilenm)
        
        @property
        def nnfilenmo(self):
            """
            Element nnfilenmo ftype=character(len=132) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 185
            
            """
            return _tbfitpy_mod.f90wrap_hopping_py__get__nnfilenmo(self._handle)
        
        @nnfilenmo.setter
        def nnfilenmo(self, nnfilenmo):
            _tbfitpy_mod.f90wrap_hopping_py__set__nnfilenmo(self._handle, nnfilenmo)
        
        @property
        def flag_site_cindex(self):
            """
            Element flag_site_cindex ftype=logical pytype=bool
            
            
            Defined at tbfitpy_mod.f90 line 186
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__flag_site_cindex(self._handle)
            if array_handle in self._arrays:
                flag_site_cindex = self._arrays[array_handle]
            else:
                flag_site_cindex = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__flag_site_cindex)
                self._arrays[array_handle] = flag_site_cindex
            return flag_site_cindex
        
        @flag_site_cindex.setter
        def flag_site_cindex(self, flag_site_cindex):
            self.flag_site_cindex[...] = flag_site_cindex
        
        @property
        def i_coord(self):
            """
            Element i_coord ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 187
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__i_coord(self._handle)
            if array_handle in self._arrays:
                i_coord = self._arrays[array_handle]
            else:
                i_coord = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__i_coord)
                self._arrays[array_handle] = i_coord
            return i_coord
        
        @i_coord.setter
        def i_coord(self, i_coord):
            self.i_coord[...] = i_coord
        
        @property
        def j_coord(self):
            """
            Element j_coord ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 188
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__j_coord(self._handle)
            if array_handle in self._arrays:
                j_coord = self._arrays[array_handle]
            else:
                j_coord = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__j_coord)
                self._arrays[array_handle] = j_coord
            return j_coord
        
        @j_coord.setter
        def j_coord(self, j_coord):
            self.j_coord[...] = j_coord
        
        @property
        def rij(self):
            """
            Element rij ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 189
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__rij(self._handle)
            if array_handle in self._arrays:
                rij = self._arrays[array_handle]
            else:
                rij = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__rij)
                self._arrays[array_handle] = rij
            return rij
        
        @rij.setter
        def rij(self, rij):
            self.rij[...] = rij
        
        @property
        def r(self):
            """
            Element r ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 190
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__r(self._handle)
            if array_handle in self._arrays:
                r = self._arrays[array_handle]
            else:
                r = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__r)
                self._arrays[array_handle] = r
            return r
        
        @r.setter
        def r(self, r):
            self.r[...] = r
        
        @property
        def dij(self):
            """
            Element dij ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 191
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__dij(self._handle)
            if array_handle in self._arrays:
                dij = self._arrays[array_handle]
            else:
                dij = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__dij)
                self._arrays[array_handle] = dij
            return dij
        
        @dij.setter
        def dij(self, dij):
            self.dij[...] = dij
        
        @property
        def dij0(self):
            """
            Element dij0 ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 192
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__dij0(self._handle)
            if array_handle in self._arrays:
                dij0 = self._arrays[array_handle]
            else:
                dij0 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__dij0)
                self._arrays[array_handle] = dij0
            return dij0
        
        @dij0.setter
        def dij0(self, dij0):
            self.dij0[...] = dij0
        
        @property
        def i_sign(self):
            """
            Element i_sign ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 193
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__i_sign(self._handle)
            if array_handle in self._arrays:
                i_sign = self._arrays[array_handle]
            else:
                i_sign = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__i_sign)
                self._arrays[array_handle] = i_sign
            return i_sign
        
        @i_sign.setter
        def i_sign(self, i_sign):
            self.i_sign[...] = i_sign
        
        @property
        def j_sign(self):
            """
            Element j_sign ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 194
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__j_sign(self._handle)
            if array_handle in self._arrays:
                j_sign = self._arrays[array_handle]
            else:
                j_sign = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__j_sign)
                self._arrays[array_handle] = j_sign
            return j_sign
        
        @j_sign.setter
        def j_sign(self, j_sign):
            self.j_sign[...] = j_sign
        
        @property
        def r_nn(self):
            """
            Element r_nn ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 195
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__r_nn(self._handle)
            if array_handle in self._arrays:
                r_nn = self._arrays[array_handle]
            else:
                r_nn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__r_nn)
                self._arrays[array_handle] = r_nn
            return r_nn
        
        @r_nn.setter
        def r_nn(self, r_nn):
            self.r_nn[...] = r_nn
        
        @property
        def r0_nn(self):
            """
            Element r0_nn ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 196
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__r0_nn(self._handle)
            if array_handle in self._arrays:
                r0_nn = self._arrays[array_handle]
            else:
                r0_nn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__r0_nn)
                self._arrays[array_handle] = r0_nn
            return r0_nn
        
        @r0_nn.setter
        def r0_nn(self, r0_nn):
            self.r0_nn[...] = r0_nn
        
        @property
        def local_charge(self):
            """
            Element local_charge ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 197
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__local_charge(self._handle)
            if array_handle in self._arrays:
                local_charge = self._arrays[array_handle]
            else:
                local_charge = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__local_charge)
                self._arrays[array_handle] = local_charge
            return local_charge
        
        @local_charge.setter
        def local_charge(self, local_charge):
            self.local_charge[...] = local_charge
        
        @property
        def local_moment(self):
            """
            Element local_moment ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 198
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__local_moment(self._handle)
            if array_handle in self._arrays:
                local_moment = self._arrays[array_handle]
            else:
                local_moment = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__local_moment)
                self._arrays[array_handle] = local_moment
            return local_moment
        
        @local_moment.setter
        def local_moment(self, local_moment):
            self.local_moment[...] = local_moment
        
        @property
        def local_moment_rot(self):
            """
            Element local_moment_rot ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 199
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__local_moment_rot(self._handle)
            if array_handle in self._arrays:
                local_moment_rot = self._arrays[array_handle]
            else:
                local_moment_rot = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__local_moment_rot)
                self._arrays[array_handle] = local_moment_rot
            return local_moment_rot
        
        @local_moment_rot.setter
        def local_moment_rot(self, local_moment_rot):
            self.local_moment_rot[...] = local_moment_rot
        
        @property
        def i_atom(self):
            """
            Element i_atom ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 200
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__i_atom(self._handle)
            if array_handle in self._arrays:
                i_atom = self._arrays[array_handle]
            else:
                i_atom = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__i_atom)
                self._arrays[array_handle] = i_atom
            return i_atom
        
        @i_atom.setter
        def i_atom(self, i_atom):
            self.i_atom[...] = i_atom
        
        @property
        def j_atom(self):
            """
            Element j_atom ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 201
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__j_atom(self._handle)
            if array_handle in self._arrays:
                j_atom = self._arrays[array_handle]
            else:
                j_atom = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__j_atom)
                self._arrays[array_handle] = j_atom
            return j_atom
        
        @j_atom.setter
        def j_atom(self, j_atom):
            self.j_atom[...] = j_atom
        
        @property
        def i_matrix(self):
            """
            Element i_matrix ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 202
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__i_matrix(self._handle)
            if array_handle in self._arrays:
                i_matrix = self._arrays[array_handle]
            else:
                i_matrix = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__i_matrix)
                self._arrays[array_handle] = i_matrix
            return i_matrix
        
        @i_matrix.setter
        def i_matrix(self, i_matrix):
            self.i_matrix[...] = i_matrix
        
        @property
        def j_matrix(self):
            """
            Element j_matrix ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 203
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__j_matrix(self._handle)
            if array_handle in self._arrays:
                j_matrix = self._arrays[array_handle]
            else:
                j_matrix = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__j_matrix)
                self._arrays[array_handle] = j_matrix
            return j_matrix
        
        @j_matrix.setter
        def j_matrix(self, j_matrix):
            self.j_matrix[...] = j_matrix
        
        @property
        def n_class(self):
            """
            Element n_class ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 204
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__n_class(self._handle)
            if array_handle in self._arrays:
                n_class = self._arrays[array_handle]
            else:
                n_class = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__n_class)
                self._arrays[array_handle] = n_class
            return n_class
        
        @n_class.setter
        def n_class(self, n_class):
            self.n_class[...] = n_class
        
        @property
        def sk_index_set(self):
            """
            Element sk_index_set ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 205
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__sk_index_set(self._handle)
            if array_handle in self._arrays:
                sk_index_set = self._arrays[array_handle]
            else:
                sk_index_set = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__sk_index_set)
                self._arrays[array_handle] = sk_index_set
            return sk_index_set
        
        @sk_index_set.setter
        def sk_index_set(self, sk_index_set):
            self.sk_index_set[...] = sk_index_set
        
        @property
        def cc_index_set(self):
            """
            Element cc_index_set ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 206
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__cc_index_set(self._handle)
            if array_handle in self._arrays:
                cc_index_set = self._arrays[array_handle]
            else:
                cc_index_set = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__cc_index_set)
                self._arrays[array_handle] = cc_index_set
            return cc_index_set
        
        @cc_index_set.setter
        def cc_index_set(self, cc_index_set):
            self.cc_index_set[...] = cc_index_set
        
        @property
        def n_nn(self):
            """
            Element n_nn ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 207
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__n_nn(self._handle)
            if array_handle in self._arrays:
                n_nn = self._arrays[array_handle]
            else:
                n_nn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__n_nn)
                self._arrays[array_handle] = n_nn
            return n_nn
        
        @n_nn.setter
        def n_nn(self, n_nn):
            self.n_nn[...] = n_nn
        
        @property
        def j_nn(self):
            """
            Element j_nn ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 208
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__j_nn(self._handle)
            if array_handle in self._arrays:
                j_nn = self._arrays[array_handle]
            else:
                j_nn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__j_nn)
                self._arrays[array_handle] = j_nn
            return j_nn
        
        @j_nn.setter
        def j_nn(self, j_nn):
            self.j_nn[...] = j_nn
        
        @property
        def l_onsite_param_index(self):
            """
            Element l_onsite_param_index ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 209
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__l_onsite_param_index(self._handle)
            if array_handle in self._arrays:
                l_onsite_param_index = self._arrays[array_handle]
            else:
                l_onsite_param_index = \
                    f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__l_onsite_param_index)
                self._arrays[array_handle] = l_onsite_param_index
            return l_onsite_param_index
        
        @l_onsite_param_index.setter
        def l_onsite_param_index(self, l_onsite_param_index):
            self.l_onsite_param_index[...] = l_onsite_param_index
        
        @property
        def stoner_i_param_index(self):
            """
            Element stoner_i_param_index ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 210
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__stoner_i_param_index(self._handle)
            if array_handle in self._arrays:
                stoner_i_param_index = self._arrays[array_handle]
            else:
                stoner_i_param_index = \
                    f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__stoner_i_param_index)
                self._arrays[array_handle] = stoner_i_param_index
            return stoner_i_param_index
        
        @stoner_i_param_index.setter
        def stoner_i_param_index(self, stoner_i_param_index):
            self.stoner_i_param_index[...] = stoner_i_param_index
        
        @property
        def local_u_param_index(self):
            """
            Element local_u_param_index ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 211
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__local_u_param_index(self._handle)
            if array_handle in self._arrays:
                local_u_param_index = self._arrays[array_handle]
            else:
                local_u_param_index = \
                    f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__local_u_param_index)
                self._arrays[array_handle] = local_u_param_index
            return local_u_param_index
        
        @local_u_param_index.setter
        def local_u_param_index(self, local_u_param_index):
            self.local_u_param_index[...] = local_u_param_index
        
        @property
        def plus_u_param_index(self):
            """
            Element plus_u_param_index ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 212
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__plus_u_param_index(self._handle)
            if array_handle in self._arrays:
                plus_u_param_index = self._arrays[array_handle]
            else:
                plus_u_param_index = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__plus_u_param_index)
                self._arrays[array_handle] = plus_u_param_index
            return plus_u_param_index
        
        @plus_u_param_index.setter
        def plus_u_param_index(self, plus_u_param_index):
            self.plus_u_param_index[...] = plus_u_param_index
        
        @property
        def soc_param_index(self):
            """
            Element soc_param_index ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 213
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__soc_param_index(self._handle)
            if array_handle in self._arrays:
                soc_param_index = self._arrays[array_handle]
            else:
                soc_param_index = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__soc_param_index)
                self._arrays[array_handle] = soc_param_index
            return soc_param_index
        
        @soc_param_index.setter
        def soc_param_index(self, soc_param_index):
            self.soc_param_index[...] = soc_param_index
        
        @property
        def ci_orb(self):
            """
            Element ci_orb ftype=character(len=8) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 214
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__ci_orb(self._handle)
            if array_handle in self._arrays:
                ci_orb = self._arrays[array_handle]
            else:
                ci_orb = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__ci_orb)
                self._arrays[array_handle] = ci_orb
            return ci_orb
        
        @ci_orb.setter
        def ci_orb(self, ci_orb):
            self.ci_orb[...] = ci_orb
        
        @property
        def cj_orb(self):
            """
            Element cj_orb ftype=character(len=8) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 215
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__cj_orb(self._handle)
            if array_handle in self._arrays:
                cj_orb = self._arrays[array_handle]
            else:
                cj_orb = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__cj_orb)
                self._arrays[array_handle] = cj_orb
            return cj_orb
        
        @cj_orb.setter
        def cj_orb(self, cj_orb):
            self.cj_orb[...] = cj_orb
        
        @property
        def p_class(self):
            """
            Element p_class ftype=character(len=2) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 216
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__p_class(self._handle)
            if array_handle in self._arrays:
                p_class = self._arrays[array_handle]
            else:
                p_class = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__p_class)
                self._arrays[array_handle] = p_class
            return p_class
        
        @p_class.setter
        def p_class(self, p_class):
            self.p_class[...] = p_class
        
        @property
        def site_cindex(self):
            """
            Element site_cindex ftype=character(len=20) pytype=str
            
            
            Defined at tbfitpy_mod.f90 line 217
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_hopping_py__array__site_cindex(self._handle)
            if array_handle in self._arrays:
                site_cindex = self._arrays[array_handle]
            else:
                site_cindex = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_hopping_py__array__site_cindex)
                self._arrays[array_handle] = site_cindex
            return site_cindex
        
        @site_cindex.setter
        def site_cindex(self, site_cindex):
            self.site_cindex[...] = site_cindex
        
        def __str__(self):
            ret = ['<hopping_py>{\n']
            ret.append('    flag_efield : ')
            ret.append(repr(self.flag_efield))
            ret.append(',\n    flag_efield_frac : ')
            ret.append(repr(self.flag_efield_frac))
            ret.append(',\n    flag_efield_cart : ')
            ret.append(repr(self.flag_efield_cart))
            ret.append(',\n    onsite_tolerance : ')
            ret.append(repr(self.onsite_tolerance))
            ret.append(',\n    efield : ')
            ret.append(repr(self.efield))
            ret.append(',\n    efield_origin : ')
            ret.append(repr(self.efield_origin))
            ret.append(',\n    efield_origin_cart : ')
            ret.append(repr(self.efield_origin_cart))
            ret.append(',\n    mysystem : ')
            ret.append(repr(self.mysystem))
            ret.append(',\n    n_neighbor : ')
            ret.append(repr(self.n_neighbor))
            ret.append(',\n    max_nn_pair : ')
            ret.append(repr(self.max_nn_pair))
            ret.append(',\n    nnfilenm : ')
            ret.append(repr(self.nnfilenm))
            ret.append(',\n    nnfilenmo : ')
            ret.append(repr(self.nnfilenmo))
            ret.append(',\n    flag_site_cindex : ')
            ret.append(repr(self.flag_site_cindex))
            ret.append(',\n    i_coord : ')
            ret.append(repr(self.i_coord))
            ret.append(',\n    j_coord : ')
            ret.append(repr(self.j_coord))
            ret.append(',\n    rij : ')
            ret.append(repr(self.rij))
            ret.append(',\n    r : ')
            ret.append(repr(self.r))
            ret.append(',\n    dij : ')
            ret.append(repr(self.dij))
            ret.append(',\n    dij0 : ')
            ret.append(repr(self.dij0))
            ret.append(',\n    i_sign : ')
            ret.append(repr(self.i_sign))
            ret.append(',\n    j_sign : ')
            ret.append(repr(self.j_sign))
            ret.append(',\n    r_nn : ')
            ret.append(repr(self.r_nn))
            ret.append(',\n    r0_nn : ')
            ret.append(repr(self.r0_nn))
            ret.append(',\n    local_charge : ')
            ret.append(repr(self.local_charge))
            ret.append(',\n    local_moment : ')
            ret.append(repr(self.local_moment))
            ret.append(',\n    local_moment_rot : ')
            ret.append(repr(self.local_moment_rot))
            ret.append(',\n    i_atom : ')
            ret.append(repr(self.i_atom))
            ret.append(',\n    j_atom : ')
            ret.append(repr(self.j_atom))
            ret.append(',\n    i_matrix : ')
            ret.append(repr(self.i_matrix))
            ret.append(',\n    j_matrix : ')
            ret.append(repr(self.j_matrix))
            ret.append(',\n    n_class : ')
            ret.append(repr(self.n_class))
            ret.append(',\n    sk_index_set : ')
            ret.append(repr(self.sk_index_set))
            ret.append(',\n    cc_index_set : ')
            ret.append(repr(self.cc_index_set))
            ret.append(',\n    n_nn : ')
            ret.append(repr(self.n_nn))
            ret.append(',\n    j_nn : ')
            ret.append(repr(self.j_nn))
            ret.append(',\n    l_onsite_param_index : ')
            ret.append(repr(self.l_onsite_param_index))
            ret.append(',\n    stoner_i_param_index : ')
            ret.append(repr(self.stoner_i_param_index))
            ret.append(',\n    local_u_param_index : ')
            ret.append(repr(self.local_u_param_index))
            ret.append(',\n    plus_u_param_index : ')
            ret.append(repr(self.plus_u_param_index))
            ret.append(',\n    soc_param_index : ')
            ret.append(repr(self.soc_param_index))
            ret.append(',\n    ci_orb : ')
            ret.append(repr(self.ci_orb))
            ret.append(',\n    cj_orb : ')
            ret.append(repr(self.cj_orb))
            ret.append(',\n    p_class : ')
            ret.append(repr(self.p_class))
            ret.append(',\n    site_cindex : ')
            ret.append(repr(self.site_cindex))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("tbfitpy_mod.energy_py")
    class energy_py(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=energy_py)
        
        
        Defined at tbfitpy_mod.f90 lines 219-230
        
        """
        def __init__(self, handle=None):
            """
            self = Energy_Py()
            
            
            Defined at tbfitpy_mod.f90 lines 219-230
            
            
            Returns
            -------
            this : Energy_Py
            	Object to be constructed
            
            
            Automatically generated constructor for energy_py
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _tbfitpy_mod.f90wrap_energy_py_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Energy_Py
            
            
            Defined at tbfitpy_mod.f90 lines 219-230
            
            Parameters
            ----------
            this : Energy_Py
            	Object to be destructed
            
            
            Automatically generated destructor for energy_py
            """
            if self._alloc:
                _tbfitpy_mod.f90wrap_energy_py_finalise(this=self._handle)
        
        @property
        def mysystem(self):
            """
            Element mysystem ftype=integer(kind=sp) pytype=int
            
            
            Defined at tbfitpy_mod.f90 line 220
            
            """
            return _tbfitpy_mod.f90wrap_energy_py__get__mysystem(self._handle)
        
        @mysystem.setter
        def mysystem(self, mysystem):
            _tbfitpy_mod.f90wrap_energy_py__set__mysystem(self._handle, mysystem)
        
        @property
        def e(self):
            """
            Element e ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 221
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_energy_py__array__e(self._handle)
            if array_handle in self._arrays:
                e = self._arrays[array_handle]
            else:
                e = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_energy_py__array__e)
                self._arrays[array_handle] = e
            return e
        
        @e.setter
        def e(self, e):
            self.e[...] = e
        
        @property
        def de(self):
            """
            Element de ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 222
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_energy_py__array__de(self._handle)
            if array_handle in self._arrays:
                de = self._arrays[array_handle]
            else:
                de = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_energy_py__array__de)
                self._arrays[array_handle] = de
            return de
        
        @de.setter
        def de(self, de):
            self.de[...] = de
        
        @property
        def dorb(self):
            """
            Element dorb ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 223
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_energy_py__array__dorb(self._handle)
            if array_handle in self._arrays:
                dorb = self._arrays[array_handle]
            else:
                dorb = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_energy_py__array__dorb)
                self._arrays[array_handle] = dorb
            return dorb
        
        @dorb.setter
        def dorb(self, dorb):
            self.dorb[...] = dorb
        
        @property
        def orb(self):
            """
            Element orb ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 224
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_energy_py__array__orb(self._handle)
            if array_handle in self._arrays:
                orb = self._arrays[array_handle]
            else:
                orb = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_energy_py__array__orb)
                self._arrays[array_handle] = orb
            return orb
        
        @orb.setter
        def orb(self, orb):
            self.orb[...] = orb
        
        @property
        def v(self):
            """
            Element v ftype=complex(kind=dp) pytype=complex
            
            
            Defined at tbfitpy_mod.f90 line 225
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_energy_py__array__v(self._handle)
            if array_handle in self._arrays:
                v = self._arrays[array_handle]
            else:
                v = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_energy_py__array__v)
                self._arrays[array_handle] = v
            return v
        
        @v.setter
        def v(self, v):
            self.v[...] = v
        
        @property
        def sv(self):
            """
            Element sv ftype=complex(kind=dp) pytype=complex
            
            
            Defined at tbfitpy_mod.f90 line 226
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_energy_py__array__sv(self._handle)
            if array_handle in self._arrays:
                sv = self._arrays[array_handle]
            else:
                sv = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_energy_py__array__sv)
                self._arrays[array_handle] = sv
            return sv
        
        @sv.setter
        def sv(self, sv):
            self.sv[...] = sv
        
        @property
        def e_band(self):
            """
            Element e_band ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 227
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_energy_py__array__e_band(self._handle)
            if array_handle in self._arrays:
                e_band = self._arrays[array_handle]
            else:
                e_band = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_energy_py__array__e_band)
                self._arrays[array_handle] = e_band
            return e_band
        
        @e_band.setter
        def e_band(self, e_band):
            self.e_band[...] = e_band
        
        @property
        def e_tot(self):
            """
            Element e_tot ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 228
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_energy_py__array__e_tot(self._handle)
            if array_handle in self._arrays:
                e_tot = self._arrays[array_handle]
            else:
                e_tot = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_energy_py__array__e_tot)
                self._arrays[array_handle] = e_tot
            return e_tot
        
        @e_tot.setter
        def e_tot(self, e_tot):
            self.e_tot[...] = e_tot
        
        @property
        def f_occ(self):
            """
            Element f_occ ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 229
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _tbfitpy_mod.f90wrap_energy_py__array__f_occ(self._handle)
            if array_handle in self._arrays:
                f_occ = self._arrays[array_handle]
            else:
                f_occ = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _tbfitpy_mod.f90wrap_energy_py__array__f_occ)
                self._arrays[array_handle] = f_occ
            return f_occ
        
        @f_occ.setter
        def f_occ(self, f_occ):
            self.f_occ[...] = f_occ
        
        @property
        def e_f(self):
            """
            Element e_f ftype=real(kind=dp) pytype=float
            
            
            Defined at tbfitpy_mod.f90 line 230
            
            """
            return _tbfitpy_mod.f90wrap_energy_py__get__e_f(self._handle)
        
        @e_f.setter
        def e_f(self, e_f):
            _tbfitpy_mod.f90wrap_energy_py__set__e_f(self._handle, e_f)
        
        def __str__(self):
            ret = ['<energy_py>{\n']
            ret.append('    mysystem : ')
            ret.append(repr(self.mysystem))
            ret.append(',\n    e : ')
            ret.append(repr(self.e))
            ret.append(',\n    de : ')
            ret.append(repr(self.de))
            ret.append(',\n    dorb : ')
            ret.append(repr(self.dorb))
            ret.append(',\n    orb : ')
            ret.append(repr(self.orb))
            ret.append(',\n    v : ')
            ret.append(repr(self.v))
            ret.append(',\n    sv : ')
            ret.append(repr(self.sv))
            ret.append(',\n    e_band : ')
            ret.append(repr(self.e_band))
            ret.append(',\n    e_tot : ')
            ret.append(repr(self.e_tot))
            ret.append(',\n    f_occ : ')
            ret.append(repr(self.f_occ))
            ret.append(',\n    e_f : ')
            ret.append(repr(self.e_f))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def init(comm, pinpt_py, ppram_py, pkpts_py, pwght_py, pgeom_py, nn_table_py, \
        edft_py, etba_py):
        """
        init(comm, pinpt_py, ppram_py, pkpts_py, pwght_py, pgeom_py, nn_table_py, \
            edft_py, etba_py)
        
        
        Defined at tbfitpy_mod.f90 lines 234-299
        
        Parameters
        ----------
        comm : int
        pinpt_py : Incar_Py
        ppram_py : Params_Py
        pkpts_py : Kpoints_Py
        pwght_py : Weight_Py
        pgeom_py : Poscar_Py
        nn_table_py : Hopping_Py
        edft_py : Energy_Py
        etba_py : Energy_Py
        
        """
        _tbfitpy_mod.f90wrap_init(comm=comm, pinpt_py=pinpt_py._handle, \
            ppram_py=ppram_py._handle, pkpts_py=pkpts_py._handle, \
            pwght_py=pwght_py._handle, pgeom_py=pgeom_py._handle, \
            nn_table_py=nn_table_py._handle, edft_py=edft_py._handle, \
            etba_py=etba_py._handle)
    
    @staticmethod
    def pso(comm, pinpt_py, ppram_py, pkpts_py, pwght_py, pgeom_py, nn_table_py, \
        edft_py, etba_py, iseed, pso_miter):
        """
        pso(comm, pinpt_py, ppram_py, pkpts_py, pwght_py, pgeom_py, nn_table_py, \
            edft_py, etba_py, iseed, pso_miter)
        
        
        Defined at tbfitpy_mod.f90 lines 301-345
        
        Parameters
        ----------
        comm : int
        pinpt_py : Incar_Py
        ppram_py : Params_Py
        pkpts_py : Kpoints_Py
        pwght_py : Weight_Py
        pgeom_py : Poscar_Py
        nn_table_py : Hopping_Py
        edft_py : Energy_Py
        etba_py : Energy_Py
        iseed : int
        pso_miter : int
        
        """
        _tbfitpy_mod.f90wrap_pso(comm=comm, pinpt_py=pinpt_py._handle, \
            ppram_py=ppram_py._handle, pkpts_py=pkpts_py._handle, \
            pwght_py=pwght_py._handle, pgeom_py=pgeom_py._handle, \
            nn_table_py=nn_table_py._handle, edft_py=edft_py._handle, \
            etba_py=etba_py._handle, iseed=iseed, pso_miter=pso_miter)
    
    @staticmethod
    def eig(comm, pinpt_py, ppram_py, pkpts_py, pgeom_py, nn_table_py, etba_py):
        """
        eig(comm, pinpt_py, ppram_py, pkpts_py, pgeom_py, nn_table_py, etba_py)
        
        
        Defined at tbfitpy_mod.f90 lines 347-378
        
        Parameters
        ----------
        comm : int
        pinpt_py : Incar_Py
        ppram_py : Params_Py
        pkpts_py : Kpoints_Py
        pgeom_py : Poscar_Py
        nn_table_py : Hopping_Py
        etba_py : Energy_Py
        
        """
        _tbfitpy_mod.f90wrap_eig(comm=comm, pinpt_py=pinpt_py._handle, \
            ppram_py=ppram_py._handle, pkpts_py=pkpts_py._handle, \
            pgeom_py=pgeom_py._handle, nn_table_py=nn_table_py._handle, \
            etba_py=etba_py._handle)
    
    @staticmethod
    def toten(comm, pinpt_py, pkpts_py, pgeom_py, etba_py):
        """
        toten(comm, pinpt_py, pkpts_py, pgeom_py, etba_py)
        
        
        Defined at tbfitpy_mod.f90 lines 380-436
        
        Parameters
        ----------
        comm : int
        pinpt_py : Incar_Py
        pkpts_py : Kpoints_Py
        pgeom_py : Poscar_Py
        etba_py : Energy_Py
        
        """
        _tbfitpy_mod.f90wrap_toten(comm=comm, pinpt_py=pinpt_py._handle, \
            pkpts_py=pkpts_py._handle, pgeom_py=pgeom_py._handle, \
            etba_py=etba_py._handle)
    
    @staticmethod
    def fit(comm, pinpt_py, ppram_py, pkpts_py, pwght_py, pgeom_py, nn_table_py, \
        edft_py, etba_py):
        """
        fit(comm, pinpt_py, ppram_py, pkpts_py, pwght_py, pgeom_py, nn_table_py, \
            edft_py, etba_py)
        
        
        Defined at tbfitpy_mod.f90 lines 438-519
        
        Parameters
        ----------
        comm : int
        pinpt_py : Incar_Py
        ppram_py : Params_Py
        pkpts_py : Kpoints_Py
        pwght_py : Weight_Py
        pgeom_py : Poscar_Py
        nn_table_py : Hopping_Py
        edft_py : Energy_Py
        etba_py : Energy_Py
        
        """
        _tbfitpy_mod.f90wrap_fit(comm=comm, pinpt_py=pinpt_py._handle, \
            ppram_py=ppram_py._handle, pkpts_py=pkpts_py._handle, \
            pwght_py=pwght_py._handle, pgeom_py=pgeom_py._handle, \
            nn_table_py=nn_table_py._handle, edft_py=edft_py._handle, \
            etba_py=etba_py._handle)
    
    @staticmethod
    def init_incar_py(ifilenm, nsystem):
        """
        pinpt_py = init_incar_py(ifilenm, nsystem)
        
        
        Defined at tbfitpy_mod.f90 lines 521-598
        
        Parameters
        ----------
        ifilenm : str array
        nsystem : int
        
        Returns
        -------
        pinpt_py : Incar_Py
        
        """
        pinpt_py = _tbfitpy_mod.f90wrap_init_incar_py(ifilenm=ifilenm, nsystem=nsystem)
        pinpt_py = \
            f90wrap.runtime.lookup_class("tbfitpy_mod.incar_py").from_handle(pinpt_py, \
            alloc=True)
        return pinpt_py
    
    @staticmethod
    def init_params_py():
        """
        ppram_py = init_params_py()
        
        
        Defined at tbfitpy_mod.f90 lines 749-787
        
        
        Returns
        -------
        ppram_py : Params_Py
        
        """
        ppram_py = _tbfitpy_mod.f90wrap_init_params_py()
        ppram_py = \
            f90wrap.runtime.lookup_class("tbfitpy_mod.params_py").from_handle(ppram_py, \
            alloc=True)
        return ppram_py
    
    @staticmethod
    def copy_params_best(self, imode):
        """
        copy_params_best(self, imode)
        
        
        Defined at tbfitpy_mod.f90 lines 789-797
        
        Parameters
        ----------
        ppram_py : Params_Py
        imode : int
        
        """
        _tbfitpy_mod.f90wrap_copy_params_best(ppram_py=self._handle, imode=imode)
    
    @staticmethod
    def init_kpoints_py():
        """
        pkpts_py = init_kpoints_py()
        
        
        Defined at tbfitpy_mod.f90 lines 1026-1042
        
        
        Returns
        -------
        pkpts_py : Kpoints_Py
        
        """
        pkpts_py = _tbfitpy_mod.f90wrap_init_kpoints_py()
        pkpts_py = \
            f90wrap.runtime.lookup_class("tbfitpy_mod.kpoints_py").from_handle(pkpts_py, \
            alloc=True)
        return pkpts_py
    
    @staticmethod
    def init_weight_py():
        """
        pwght_py = init_weight_py()
        
        
        Defined at tbfitpy_mod.f90 lines 1169-1182
        
        
        Returns
        -------
        pwght_py : Weight_Py
        
        """
        pwght_py = _tbfitpy_mod.f90wrap_init_weight_py()
        pwght_py = \
            f90wrap.runtime.lookup_class("tbfitpy_mod.weight_py").from_handle(pwght_py, \
            alloc=True)
        return pwght_py
    
    @staticmethod
    def init_poscar_py():
        """
        pgeom_py = init_poscar_py()
        
        
        Defined at tbfitpy_mod.f90 lines 1226-1233
        
        
        Returns
        -------
        pgeom_py : Poscar_Py
        
        """
        pgeom_py = _tbfitpy_mod.f90wrap_init_poscar_py()
        pgeom_py = \
            f90wrap.runtime.lookup_class("tbfitpy_mod.poscar_py").from_handle(pgeom_py, \
            alloc=True)
        return pgeom_py
    
    @staticmethod
    def init_hopping_py():
        """
        nn_table_py = init_hopping_py()
        
        
        Defined at tbfitpy_mod.f90 lines 1322-1360
        
        
        Returns
        -------
        nn_table_py : Hopping_Py
        
        """
        nn_table_py = _tbfitpy_mod.f90wrap_init_hopping_py()
        nn_table_py = \
            f90wrap.runtime.lookup_class("tbfitpy_mod.hopping_py").from_handle(nn_table_py, \
            alloc=True)
        return nn_table_py
    
    @staticmethod
    def init_energy_py():
        """
        e_py = init_energy_py()
        
        
        Defined at tbfitpy_mod.f90 lines 1715-1727
        
        
        Returns
        -------
        e_py : Energy_Py
        
        """
        e_py = _tbfitpy_mod.f90wrap_init_energy_py()
        e_py = f90wrap.runtime.lookup_class("tbfitpy_mod.energy_py").from_handle(e_py, \
            alloc=True)
        return e_py
    
    @staticmethod
    def print_param_py(self, ppram_py, pfileoutnm):
        """
        print_param_py(self, ppram_py, pfileoutnm)
        
        
        Defined at tbfitpy_mod.f90 lines 1832-1882
        
        Parameters
        ----------
        pinpt_py : Incar_Py
        ppram_py : Params_Py
        pfileoutnm : str
        
        """
        _tbfitpy_mod.f90wrap_print_param_py(pinpt_py=self._handle, \
            ppram_py=ppram_py._handle, pfileoutnm=pfileoutnm)
    
    @staticmethod
    def print_weight(self, wfileoutnm):
        """
        print_weight(self, wfileoutnm)
        
        
        Defined at tbfitpy_mod.f90 lines 1884-1904
        
        Parameters
        ----------
        pwght_py : Weight_Py
        wfileoutnm : str
        
        """
        _tbfitpy_mod.f90wrap_print_weight(pwght_py=self._handle, wfileoutnm=wfileoutnm)
    
    @staticmethod
    def print_target(self, pkpts_py, edft_py, pwght_py, pgeom_py, tfileoutnm):
        """
        print_target(self, pkpts_py, edft_py, pwght_py, pgeom_py, tfileoutnm)
        
        
        Defined at tbfitpy_mod.f90 lines 1906-1924
        
        Parameters
        ----------
        pinpt_py : Incar_Py
        pkpts_py : Kpoints_Py
        edft_py : Energy_Py
        pwght_py : Weight_Py
        pgeom_py : Poscar_Py
        tfileoutnm : str
        
        """
        _tbfitpy_mod.f90wrap_print_target(pinpt_py=self._handle, \
            pkpts_py=pkpts_py._handle, edft_py=edft_py._handle, \
            pwght_py=pwght_py._handle, pgeom_py=pgeom_py._handle, tfileoutnm=tfileoutnm)
    
    @staticmethod
    def print_etba(self, pkpts_py, etba_py, edft_py, pwght_py, pgeom_py, suffix, \
        flag_use_overlap):
        """
        print_etba(self, pkpts_py, etba_py, edft_py, pwght_py, pgeom_py, suffix, \
            flag_use_overlap)
        
        
        Defined at tbfitpy_mod.f90 lines 1926-1952
        
        Parameters
        ----------
        pinpt_py : Incar_Py
        pkpts_py : Kpoints_Py
        etba_py : Energy_Py
        edft_py : Energy_Py
        pwght_py : Weight_Py
        pgeom_py : Poscar_Py
        suffix : str
        flag_use_overlap : bool
        
        """
        _tbfitpy_mod.f90wrap_print_etba(pinpt_py=self._handle, \
            pkpts_py=pkpts_py._handle, etba_py=etba_py._handle, edft_py=edft_py._handle, \
            pwght_py=pwght_py._handle, pgeom_py=pgeom_py._handle, suffix=suffix, \
            flag_use_overlap=flag_use_overlap)
    
    @staticmethod
    def set_verbose(iverbose_mod):
        """
        set_verbose(iverbose_mod)
        
        
        Defined at tbfitpy_mod.f90 lines 1954-1964
        
        Parameters
        ----------
        iverbose_mod : int
        
        """
        _tbfitpy_mod.f90wrap_set_verbose(iverbose_mod=iverbose_mod)
    
    _dt_array_initialisers = []
    

pyfit = Pyfit()

