import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
#from random import random
import random
import warnings
from tqdm import tqdm 
import time
import sys
import os
import gc
import torch 
from mpi4py import MPI
from tbfitpy_mod_mpi import pyfit
#from tbfitpy_mod import pyfit
warnings.filterwarnings("ignore")

# last update: 14.04.2022 HJ Kim
# last update: 21.06.2022 HJ Kim
# last update: 16.07.2022 HJ Kim - add nsystem=4,5
# last update: 22.07.2022 HJ Kim - add CSA class
# last update: 31.07.2022 HJ Kim - add parameter constraint function
# last update: 09.08.2022 HJ Kim - fixed parameters are to be excluded (CSA class)
# last update: 10.08.2022 HJ Kim - update parameter constraint function

# IMPORT NOTE:
# if you want to run tbfitpy_mod_mpi with MPI implementation, 
# uncomment "from mpi4py import MPI " and, "from tbfitpy_mod_mpi ..."
# comment out "from tbfitpy_mod import pyfit".
# If you want to run serial tbfitpy_mod, 
# comment out "from mpi4py ..." and "from tbfitpy_mod_mpi ..."
class mycolor:
    orb_index = {'s'            :0,
                 'px'           :0,
                 'py'           :1,
                 'pz'           :2,
                 'dxy'          :0,
                 'dyz'          :1,
                 'dxz'          :2,
                 'dx2-y2'       :3,
                 'dx2'          :3,
                 'dz2'          :4,
                 'fxz2'         :0,
                 'fyz2'         :1,
                 'fz3'          :2,
                 'fxyz'         :3,
                 'fz(x2-y2)'    :4,
                 'fz'           :4,
                 'fx(x2-3y2)'   :5,
                 'fx'           :5,
                 'fy(3x2-y2)'   :6,
                 'fy'           :6,
                 }
    marker_index = {'s'            :'o',
                    'px'           :'>',
                    'py'           :'^',
                    'pz'           :'d',
                    'dxy'          :'x',
                    'dyz'          :'|',
                    'dxz'          :'_',
                    'dx2-y2'       :'+',
                    'dx2'          :'+',
                    'dz2'          :'8',
                    'fxz2'         :'H',
                    'fyz2'         :'h',
                    'fz3'          :'D',
                    'fxyz'         :',',
                    'fz(x2-y2)'    :'X',
                    'fz'           :'X',
                    'fx(x2-3y2)'   :'*',
                    'fx'           :'*',
                    'fy(3x2-y2)'   :'p',
                    'fy'           :'p',
                 }

    color_index =  {'s'            :['black','darkgrey'],
                    'px'           :['navy','lightcoral'],
                    'py'           :['darkgreen','firebrik'],
                    'pz'           :['darkcyan','darkorange'],
                    'dxy'          :['mediumpurple','darkgoldenrod'],
                    'dyz'          :['deepskyblue','yellow'],
                    'dxz'          :['limegreen','hotpink'],
                    'dx2-y2'       :['indigo','peachpuff'],
                    'dx2'          :['indigo','peachpuff'],
                    'dz2'          :['b','r'],
                    'fxz2'         :['olivedrab','tan'],
                    'fyz2'         :['darkslateblue','rosybrown'],
                    'fz3'          :['steelblue','sandybrown'],
                    'fxyz'         :['darkorchid','khaki'],
                    'fz(x2-y2)'    :['lightgreen','fuchsia'],
                    'fz'           :['lightgreen','fuchsia'],
                    'fx(x2-3y2)'   :['lightblue','gold'],
                    'fx'           :['lightblue','gold'],
                    'fy(3x2-y2)'   :['aqua','lightgoldenrodyellow'],
                    'fy'           :['aqua','lightgoldenrodyellow'],
                 }
class myfont:
    font = {'family': 'sans-serif',
            'color':  'black',
            'weight': 'normal',
            'size': 16,
            }


class pytbfit:
    def __init__(self, mpicomm =None, filenm = 'INCAR-TB'):

        self.pfile_parse = None
        self.reduce_overlap = None
        self.reduce_hopping = None
    
       #self.reduce_overlap = 1.0
       #self.reduce_hopping = 1.0
        if len(sys.argv) >=2:
            for i in range( len(sys.argv) ):
                if str(sys.argv[i]) == '-p' :
                    self.pfile_parse = str(sys.argv[i+1])
                elif str(sys.argv[i]) == '-red_ovl':
                    self.reduce_overlap = float(sys.argv[i+1])
                elif str(sys.argv[i]) == '-red_hop':
                    self.reduce_hopping = float(sys.argv[i+1])

        if mpicomm is not None :
            self.comm  = mpicomm
            self.fcomm = self.comm.py2f()
        else:
            self.comm  = None
            self.fcomm = 0

        if type(filenm) is tuple:
            pass
        elif type(filenm) is str:
            filenm=[filenm]
        elif type(filenm) is set:
            pass
        elif type(filenm) is dict:
            filenm=list(filenm.values())
        elif type(filenm) is list:
            filenm=list(filenm)

        flist=list(filenm)
        nsystem = len(flist) ; self.nsystem = nsystem
        self.filenm=np.empty((nsystem,132),dtype='c')
        for i in range(nsystem):
            f=flist[i]
            self.filenm[i]=f+' '*(132-len(f))
        self.pinpt = pyfit.init_incar_py( self.filenm.T , nsystem=nsystem )

        self.ppram = pyfit.init_params_py()    
        self.pkpts = pyfit.init_kpoints_py()
        self.pwght = pyfit.init_weight_py()
        self.pgeom = pyfit.init_poscar_py()
        self.hopping = pyfit.init_hopping_py()
        self.edft  = pyfit.init_energy_py()
        self.etba  = pyfit.init_energy_py()
        self.orbfit= False

        if nsystem >= 2:
            self.pkpts2= pyfit.init_kpoints_py()
            self.pwght2= pyfit.init_weight_py()
            self.pgeom2= pyfit.init_poscar_py()
            self.hopping2= pyfit.init_hopping_py()
            self.edft2 = pyfit.init_energy_py()
            self.etba2 = pyfit.init_energy_py()

        if nsystem >= 3:
            self.pkpts3= pyfit.init_kpoints_py()
            self.pwght3= pyfit.init_weight_py()
            self.pgeom3= pyfit.init_poscar_py()
            self.hopping3= pyfit.init_hopping_py()
            self.edft3 = pyfit.init_energy_py()
            self.etba3 = pyfit.init_energy_py()

        if nsystem >= 4:
            self.pkpts4= pyfit.init_kpoints_py()
            self.pwght4= pyfit.init_weight_py()
            self.pgeom4= pyfit.init_poscar_py()
            self.hopping4= pyfit.init_hopping_py()
            self.edft4 = pyfit.init_energy_py()
            self.etba4 = pyfit.init_energy_py()

        if nsystem >= 5:
            self.pkpts5= pyfit.init_kpoints_py()
            self.pwght5= pyfit.init_weight_py()
            self.pgeom5= pyfit.init_poscar_py()
            self.hopping5= pyfit.init_hopping_py()
            self.edft5 = pyfit.init_energy_py()
            self.etba5 = pyfit.init_energy_py()


    def init(self, verbose=False, orbfit=False, myid=0, pfilenm=None, red_ovl=None, red_hop=None):
        # if pfilenm is specified, use it. 
        # if pfilenm is not specified but remotely specified by -p option, use this (2nd option)
        # if pfilenm is not specified and -p is also not specified, PFILE in INCAR-TB is used (3rd option)
        if (self.pfile_parse != None) and (pfilenm is None) :
            pfilenm = self.pfile_parse

        if (self.reduce_overlap != None) and (red_ovl is None):
            red_ovl = self.reduce_overlap

        if (self.reduce_hopping != None) and (red_hop is None):
            red_hop = self.reduce_hopping

        if verbose is True:
            self.pinpt.iverbose = 1
        else :
            self.pinpt.iverbose = 2

        if orbfit is True:
            self.pinpt.flag_fit_orbital_parse = True
        else:
            self.pinpt.flag_fit_orbital_parse = False
        self.orbfit = orbfit

        pfilenm_parse = np.empty((132),dtype='c')
        if pfilenm != None:
            pfilenm_parse = pfilenm + ' '*(132-len(pfilenm))
            self.pinpt.flag_pfile_parse = True
        else:
            self.pinpt.flag_pfile_parse = False

        if red_ovl !=None:
            self.pinpt.flag_reduce_overlap_parse = True
            self.pinpt.reduce_overlap_parse = red_ovl
        else:
            self.pinpt.flag_reduce_overlap_parse = False       

        if red_hop !=None:
            self.pinpt.flag_reduce_hopping_parse = True
            self.pinpt.reduce_hopping_parse = red_hop
        else:
            self.pinpt.flag_reduce_hopping_parse = False       

        # initialize
        if self.nsystem == 1:
            pyfit.init(self.fcomm, self.pinpt,self.ppram,self.pkpts,self.pwght,self.pgeom,self.hopping,self.edft,self.etba, 
                       pfilenm_parse)
            self.energy_target, self.orb_target = self.load_band(self.pwght, self.pgeom)
            self.i_specs , self.i_atoms = self.orb_index(self.pgeom)
        elif self.nsystem == 2:
            pyfit.init2(self.fcomm, self.pinpt,self.ppram,self.pkpts,self.pwght,self.pgeom,self.hopping,self.edft,self.etba,
                                                          self.pkpts2,self.pwght2,self.pgeom2,self.hopping2,self.edft2,self.etba2,
                                                          pfilenm_parse)
            self.energy_target , self.orb_target  = self.load_band(self.pwght , self.pgeom )
            self.energy_target2, self.orb_target2 = self.load_band(self.pwght2, self.pgeom2)
            self.i_specs , self.i_atoms = self.orb_index(self.pgeom)
            self.i_specs2, self.i_atoms2= self.orb_index(self.pgeom2)
        elif self.nsystem == 3:
            pyfit.init3(self.fcomm, self.pinpt,self.ppram,self.pkpts,self.pwght,self.pgeom,self.hopping,self.edft,self.etba,
                                                          self.pkpts2,self.pwght2,self.pgeom2,self.hopping2,self.edft2,self.etba2,
                                                          self.pkpts3,self.pwght3,self.pgeom3,self.hopping3,self.edft3,self.etba3, 
                                                          pfilenm_parse)
            self.energy_target , self.orb_target  = self.load_band(self.pwght , self.pgeom )
            self.energy_target2, self.orb_target2 = self.load_band(self.pwght2, self.pgeom2)
            self.energy_target3, self.orb_target3 = self.load_band(self.pwght3, self.pgeom3)
            self.i_specs , self.i_atoms = self.orb_index(self.pgeom)
            self.i_specs2, self.i_atoms2= self.orb_index(self.pgeom2)
            self.i_specs3, self.i_atoms3= self.orb_index(self.pgeom3)
        elif self.nsystem == 4:
            pyfit.init4(self.fcomm, self.pinpt,self.ppram,self.pkpts,self.pwght,self.pgeom,self.hopping,self.edft,self.etba,
                                                          self.pkpts2,self.pwght2,self.pgeom2,self.hopping2,self.edft2,self.etba2,
                                                          self.pkpts3,self.pwght3,self.pgeom3,self.hopping3,self.edft3,self.etba3,
                                                          self.pkpts4,self.pwght4,self.pgeom4,self.hopping4,self.edft4,self.etba4,
                                                          pfilenm_parse)
            self.energy_target , self.orb_target  = self.load_band(self.pwght , self.pgeom )
            self.energy_target2, self.orb_target2 = self.load_band(self.pwght2, self.pgeom2)
            self.energy_target3, self.orb_target3 = self.load_band(self.pwght3, self.pgeom3)
            self.energy_target4, self.orb_target4 = self.load_band(self.pwght4, self.pgeom4)
            self.i_specs , self.i_atoms = self.orb_index(self.pgeom)
            self.i_specs2, self.i_atoms2= self.orb_index(self.pgeom2)
            self.i_specs3, self.i_atoms3= self.orb_index(self.pgeom3)
            self.i_specs4, self.i_atoms4= self.orb_index(self.pgeom4)

        elif self.nsystem == 5:
            pyfit.init5(self.fcomm, self.pinpt,self.ppram,self.pkpts,self.pwght,self.pgeom,self.hopping,self.edft,self.etba,
                                                          self.pkpts2,self.pwght2,self.pgeom2,self.hopping2,self.edft2,self.etba2,
                                                          self.pkpts3,self.pwght3,self.pgeom3,self.hopping3,self.edft3,self.etba3,
                                                          self.pkpts4,self.pwght4,self.pgeom4,self.hopping4,self.edft4,self.etba4,
                                                          self.pkpts5,self.pwght5,self.pgeom5,self.hopping5,self.edft5,self.etba5,
                                                          pfilenm_parse)
            self.energy_target , self.orb_target  = self.load_band(self.pwght , self.pgeom )
            self.energy_target2, self.orb_target2 = self.load_band(self.pwght2, self.pgeom2)
            self.energy_target3, self.orb_target3 = self.load_band(self.pwght3, self.pgeom3)
            self.energy_target4, self.orb_target4 = self.load_band(self.pwght4, self.pgeom4)
            self.energy_target5, self.orb_target5 = self.load_band(self.pwght5, self.pgeom5)
            self.i_specs , self.i_atoms = self.orb_index(self.pgeom)
            self.i_specs2, self.i_atoms2= self.orb_index(self.pgeom2)
            self.i_specs3, self.i_atoms3= self.orb_index(self.pgeom3)
            self.i_specs4, self.i_atoms4= self.orb_index(self.pgeom4)
            self.i_specs5, self.i_atoms5= self.orb_index(self.pgeom5)

        self.cost_history = []
        self.niter = 0
        self.cost = 0.0
        self.cost_orb = 0.0
        self.cost_total = self.cost + self.cost_orb
        self.myid = myid

    def orb_index(self,pgeom):
        ii = 0
       #c_orbital = np.array([]) # initialize list of orbital
        i_specs   = {} # initialize dictionary 
        i_atoms   = {}
        for iatom in range(pgeom.n_atom):
            spec_id = pgeom.spec[iatom] - 1 # -1 is due to the pgeom.spec is from fortran 
            spec_key= str(pgeom.c_spec.T[spec_id].view('S8'),'utf-8').strip()
            atom_key= str(iatom)
            for iorb in range(pgeom.n_orbital[iatom]):
                c_orb = str(pgeom.c_orbital[:,iorb,iatom].view('S8'),'utf-8').strip()
               #c_orbital = np.append(c_orbital,c_orb)
                key = spec_key+':'+c_orb
                key_atom = 'AT'+atom_key+':'+c_orb
                i_atoms[key_atom]=[ii]
                try:
                    i_specs[key].append(ii)
                    ii = ii + 1
                except KeyError:
                    i_specs[key]=[ii]
                    ii = ii + 1
                
        return i_specs, i_atoms

    def load_band(self,pwght, pgeom):
        fname_target = pwght.efilenmu.decode('utf-8').strip()
        ntarget = pgeom.neig_target
        ef = pwght.efile_ef
        
        data = np.loadtxt(fname_target)
        energy_target = data[:,1].reshape(ntarget,-1) - ef
        if data.shape[1] >= 15:
            orb_target = data[:,2:].reshape(ntarget,-1,data.shape[1]-2)
        else:
            orb_target = None
        return energy_target, orb_target 

    def constraint_param(self,iparam_type=None, iparam=None, fix=None, ub=None, lb=None):
        # Note: iparam -> param index, ppram.param_const[:,iparam], 1 ~ self.ppram.nparam
        # note: ppram.param_const[:, :].shape = (5,nparam)
        #       i=0 : 'is same as'
        #       i=1 : 'is lower than' (.le.) : maximum/upper bound  ! <=  20
        #       i=2 : 'is greater than' (.ge.) : minimum/lower bound  ! >= -20 or >= 0.001 (if scale factor)
        #       i=3 : 'is fixed' : fixed during fitting (1 if fixed, 0 if not fixed)
        #       i=4 : if 'fixed', original value will be copied to PINPT%param_const(i=5,:)
        if iparam_type is None:
            pass
        else:
            # Note: iparam_type = 1 -> onsite        (e_, ...) 
            #       iparam_type = 2 -> hopping       (sss_, sps_, pps_, ..., etc)
            #       iparam_type = 3 -> hopping_scale (s_sss_, s_pps_, s_pps_, ..., etc)
            #       iparam_type = 4 -> overlap       (o_sss_, o_pps_, o_pps_, ..., etc)
            #       iparam_type = 5 -> overlap_scale (os_sss_, os_pps_, os_sps_, ..., etc)
            if isinstance(iparam_type, int):
                valid_param_idx=np.ma.masked_equal(self.ppram.iparam_type,iparam_type).mask            
                self.set_param_fix(fix=fix,param_idx=valid_param_idx)
                self.set_param_bound(bnds=[lb,ub],param_idx=valid_param_idx)
            elif isinstance(iparam_type, list) or isinstance(iparam_type, tuple):
                for itype in iparam_type:
                    valid_param_idx=np.ma.masked_equal(self.ppram.iparam_type,itype).mask
                    self.set_param_fix(fix=fix,param_idx=valid_param_idx)
                    self.set_param_bound(bnds=[lb,ub],param_idx=valid_param_idx)

        if iparam is not None:
            self.set_param_fix(fix=fix,param_idx=iparam)

    def set_param_bound(self, bnds=None, param_idx=None):
        if bnds is not None:
            lb = bnds[0] # lower bound
            ub = bnds[1] # upper bound
            if ub is not None:
                self.ppram.param_const[1,param_idx] = ub
            if lb is not None:
                self.ppram.param_const[2,param_idx] = lb

    def set_param_fix(self, fix=None, param_idx=None): 
        if fix is True:
            self.ppram.param_const[3,param_idx] = 1.0
            self.ppram.param_const[4,param_idx] = self.ppram.param[param_idx]
        elif fix is False:
            self.ppram.param_const[3,param_idx] = 0.0

    def fit(self, verbose=False, miter=None, tol=None, pso_miter=None, 
            method='lmdif', pso_options=None, n_particles=None, iseed=None, sigma=400, sigma_orb=4000):
        '''
        Fit the parameters.
        Following methods are supported:
         - Minimization function via LMFIT method
           *lmdif: Levenberg-Marquardt method with MINPACK subroutine modified by H.-J. Kim (TBFIT)
              (see details in: https://github.com/Infant83/TBFIT)

         - Particle Swarm Optimization (PSO) scheme
           *mypso: A Global-best Particle Swarm Optimization algorithm
           *mypso.lmdif: same as pso.leastsq but much faster

           Note on all PSO based altorhtim: one should provide n_particle and pso_options when call
                n_particles: number of particles (the swarm size), integer
                pso_options: velocity and position update policy, dictionary, 
                             ex: pso_options = {'c1': 1.0, 'c2': 1.0, 'w': 2.0'}
                             for details, see: J. Kennedy and R. Eberhart, Particle Swarm Optimization,
                                               IEEE, Piscataway, NJ, 1995, p. 1942.
        NONE: For global parameter fitting, mypso is suggested. 
              To utilize lmdif method with mypso, mypso.lmdif can be used. 

        - sigma : smearing function for cost function. Default = 400

        '''

        if verbose is not None:
            if verbose is True:
                self.pinpt.iverbose = 1
            else:
                self.pinpt.iverbose = 2

        flag_tbfit = self.pinpt.flag_tbfit
        self.pinpt.flag_tbfit = -1

        t0 = time.time()

        if tol is not None: self.pinpt.ptol = tol
        if tol is not None: self.pinpt.ftol = tol

        if miter is not None: self.pinpt.miter = miter
        self.n_particles = n_particles if n_particles  is not None else 30
        if pso_options is None: pso_options = {'c1': 1.0, 'c2': 2.0, 'w':2.0} # set these values as default

        if method == 'lmdif' :
            self.cost_history = np.full( (self.pinpt.miter), 0.0)
            self.pinpt.flag_tbfit = True
            if self.nsystem == 1:
                pyfit.fit(self.fcomm, self.pinpt, self.ppram, self.pkpts, self.pwght, 
                                      self.pgeom, self.hopping, self.edft, self.etba)
                max_wt = np.max(self.pwght.wt)
                sum_cost     = np.sum(abs(self.etba.de   * (self.pwght.wt/max_wt) ))
                if self.orbfit is True:
                    sum_cost_orb = np.sum(abs(self.etba.dorb * (self.pwght.wt/max_wt) ))
            elif self.nsystem == 2:
                pyfit.fit2(self.fcomm, self.pinpt, self.ppram, self.pkpts, self.pwght, 
                                       self.pgeom, self.hopping, self.edft, self.etba,
                                       self.pkpts2, self.pwght2,self.pgeom2, self.hopping2, self.edft2, self.etba2)
                max_wt = np.max(self.pwght.wt)
                max_wt2= np.max(self.pwght2.wt)
                sum_cost = np.sum(abs(self.etba.de * (self.pwght.wt/max_wt) )) + np.sum(abs(self.etba2.de * (self.pwght2.wt/max_wt2)))
                if self.orbfit is True:
                    sum_cost_orb = np.sum(abs(self.etba.dorb * (self.pwght.wt/max_wt) )) + np.sum(abs(self.etba2.dorb * (self.pwght2.wt/max_wt2)))
            elif self.nsystem == 3:
                pyfit.fit3(self.fcomm, self.pinpt, self.ppram, self.pkpts, self.pwght, 
                                       self.pgeom, self.hopping, self.edft, self.etba,
                                       self.pkpts2, self.pwght2,self.pgeom2, self.hopping2, self.edft2, self.etba2,
                                       self.pkpts3, self.pwght3,self.pgeom3, self.hopping3, self.edft3, self.etba3)
                max_wt = np.max(self.pwght.wt)
                max_wt2= np.max(self.pwght2.wt)
                max_wt3= np.max(self.pwght3.wt)
                sum_cost = np.sum(abs(self.etba.de  * (self.pwght.wt/max_wt) )) +\
                           np.sum(abs(self.etba2.de * (self.pwght2.wt/max_wt2) )) +\
                           np.sum(abs(self.etba3.de * (self.pwght3.wt/max_wt3) ))
                if self.orbfit is True:
                    sum_cost_orb = np.sum(abs(self.etba.dorb  * (self.pwght.wt/max_wt)   )) +\
                                   np.sum(abs(self.etba2.dorb * (self.pwght2.wt/max_wt2) )) +\
                                   np.sum(abs(self.etba3.dorb * (self.pwght3.wt/max_wt3) ))
            elif self.nsystem == 4:
                pyfit.fit4(self.fcomm, self.pinpt, self.ppram, self.pkpts, self.pwght,
                                       self.pgeom, self.hopping, self.edft, self.etba,
                                       self.pkpts2, self.pwght2,self.pgeom2, self.hopping2, self.edft2, self.etba2,
                                       self.pkpts3, self.pwght3,self.pgeom3, self.hopping3, self.edft3, self.etba3,
                                       self.pkpts4, self.pwght4,self.pgeom4, self.hopping4, self.edft4, self.etba4)
                max_wt = np.max(self.pwght.wt)
                max_wt2= np.max(self.pwght2.wt)
                max_wt3= np.max(self.pwght3.wt)
                max_wt4= np.max(self.pwght4.wt)
                sum_cost = np.sum(abs(self.etba.de  * (self.pwght.wt/max_wt) )) +\
                           np.sum(abs(self.etba2.de * (self.pwght2.wt/max_wt2) )) +\
                           np.sum(abs(self.etba3.de * (self.pwght3.wt/max_wt3) )) +\
                           np.sum(abs(self.etba4.de * (self.pwght4.wt/max_wt4) ))
                if self.orbfit is True:
                    sum_cost_orb = np.sum(abs(self.etba.dorb  * (self.pwght.wt/max_wt)   )) +\
                                   np.sum(abs(self.etba2.dorb * (self.pwght2.wt/max_wt2) )) +\
                                   np.sum(abs(self.etba3.dorb * (self.pwght3.wt/max_wt3) )) +\
                                   np.sum(abs(self.etba4.dorb * (self.pwght4.wt/max_wt4) ))
            elif self.nsystem == 5:
                pyfit.fit5(self.fcomm, self.pinpt, self.ppram, self.pkpts, self.pwght,
                                       self.pgeom, self.hopping, self.edft, self.etba,
                                       self.pkpts2, self.pwght2,self.pgeom2, self.hopping2, self.edft2, self.etba2,
                                       self.pkpts3, self.pwght3,self.pgeom3, self.hopping3, self.edft3, self.etba3,
                                       self.pkpts4, self.pwght4,self.pgeom4, self.hopping4, self.edft4, self.etba4,
                                       self.pkpts5, self.pwght5,self.pgeom5, self.hopping5, self.edft5, self.etba5)
                max_wt = np.max(self.pwght.wt)
                max_wt2= np.max(self.pwght2.wt)
                max_wt3= np.max(self.pwght3.wt)
                max_wt4= np.max(self.pwght4.wt)
                max_wt5= np.max(self.pwght5.wt)
                sum_cost = np.sum(abs(self.etba.de  * (self.pwght.wt/max_wt)   )) +\
                           np.sum(abs(self.etba2.de * (self.pwght2.wt/max_wt2) )) +\
                           np.sum(abs(self.etba3.de * (self.pwght3.wt/max_wt3) )) +\
                           np.sum(abs(self.etba4.de * (self.pwght4.wt/max_wt4) )) +\
                           np.sum(abs(self.etba5.de * (self.pwght5.wt/max_wt5) ))
                if self.orbfit is True:
                    sum_cost_orb = np.sum(abs(self.etba.dorb  * (self.pwght.wt/max_wt)   )) +\
                                   np.sum(abs(self.etba2.dorb * (self.pwght2.wt/max_wt2) )) +\
                                   np.sum(abs(self.etba3.dorb * (self.pwght3.wt/max_wt3) )) +\
                                   np.sum(abs(self.etba4.dorb * (self.pwght4.wt/max_wt4) )) +\
                                   np.sum(abs(self.etba5.dorb * (self.pwght5.wt/max_wt5) ))

            self.cost_history = self.ppram.cost_history
            self.cost     = 100.0 - np.exp( -( sum_cost / sigma)**2)*100.0
            if self.orbfit is True: self.cost_orb = 100.0 - np.exp( -( sum_cost_orb / sigma_orb)**2)*100.0
            self.cost_total = self.cost + self.cost_orb

        elif method == 'mypso' or method == 'mypso.lmdif':
            self.pinpt.flag_tbfit = True
            if iseed is None: iseed = 123
            self.pinpt.flag_pso_with_lmdif = True if method == 'mypso.lmdif' else False
            if pso_miter is None: 
              pso_miter = self.ppram.pso_miter
            self.ppram.pso_c1 = pso_options['c1']
            self.ppram.pso_c2 = pso_options['c2']
            self.ppram.pso_w  = pso_options['w' ]
            self.ppram.pso_nparticles = self.n_particles
            self.cost_history   = np.full( (self.ppram.pso_miter), 0.0)
            self.cost_history_particle = np.full( (self.ppram.pso_miter, self.n_particles), 0.0)
            
            if self.nsystem == 1:
                pyfit.pso(self.fcomm, self.pinpt, self.ppram, self.pkpts, self.pwght,
                                      self.pgeom, self.hopping, self.edft, self.etba, iseed, pso_miter)
            elif self.nsystem == 2:
                pyfit.pso2(self.fcomm, self.pinpt, self.ppram, self.pkpts, self.pwght,
                                       self.pgeom, self.hopping, self.edft, self.etba, 
                                       self.pkpts2,self.pwght2,self.pgeom2,self.hopping2,self.edft2,self.etba2,iseed,pso_miter)
            elif self.nsystem == 3:
                pyfit.pso3(self.fcomm, self.pinpt, self.ppram, self.pkpts, self.pwght,
                                       self.pgeom, self.hopping, self.edft, self.etba, 
                                       self.pkpts2,self.pwght2,self.pgeom2,self.hopping2,self.edft2,self.etba2,
                                       self.pkpts3,self.pwght3,self.pgeom3,self.hopping3,self.edft3,self.etba3,iseed,pso_miter)
            elif self.nsystem == 4:
                pyfit.pso4(self.fcomm, self.pinpt, self.ppram, self.pkpts, self.pwght,
                                       self.pgeom, self.hopping, self.edft, self.etba,
                                       self.pkpts2,self.pwght2,self.pgeom2,self.hopping2,self.edft2,self.etba2,
                                       self.pkpts3,self.pwght3,self.pgeom3,self.hopping3,self.edft3,self.etba3,
                                       self.pkpts4,self.pwght4,self.pgeom4,self.hopping4,self.edft4,self.etba4,iseed,pso_miter)
            elif self.nsystem == 5:
                pyfit.pso5(self.fcomm, self.pinpt, self.ppram, self.pkpts, self.pwght,
                                       self.pgeom, self.hopping, self.edft, self.etba,
                                       self.pkpts2,self.pwght2,self.pgeom2,self.hopping2,self.edft2,self.etba2,
                                       self.pkpts3,self.pwght3,self.pgeom3,self.hopping3,self.edft3,self.etba3,
                                       self.pkpts4,self.pwght4,self.pgeom4,self.hopping4,self.edft4,self.etba4,
                                       self.pkpts5,self.pwght5,self.pgeom5,self.hopping5,self.edft5,self.etba5,iseed,pso_miter)

            self.cost_history   = self.ppram.pso_cost_history
            self.cost_history_particle = self.ppram.pso_cost_history_i

        if self.myid == 0 and verbose is True:
            print('Time elapsed for ',method,' method :', time.time() - t0, ' sec')
            sys.stdout.flush()

        self.pinpt.flag_tbfit = flag_tbfit

        #return self.ppram.param

    def get_eig(self, verbose=None, sys=None):
        if verbose is not None:
            if verbose is True:
                self.pinpt.iverbose = 1
            else:
                self.pinpt.iverbose = 2

        flag_tbfit = self.pinpt.flag_tbfit
        self.pinpt.flag_tbfit = 0
        
        if sys is None:
            if self.nsystem == 1:
                pyfit.eig(self.fcomm, self.pinpt, self.ppram, self.pkpts,self.pgeom,self.hopping,self.etba)
            elif self.nsystem == 2:
                pyfit.eig2(self.fcomm, self.pinpt, self.ppram, self.pkpts ,self.pgeom ,self.hopping ,self.etba ,
                                                              self.pkpts2,self.pgeom2,self.hopping2,self.etba2)
            elif self.nsystem == 3:
                pyfit.eig3(self.fcomm, self.pinpt, self.ppram, self.pkpts ,self.pgeom ,self.hopping ,self.etba ,
                                                              self.pkpts2,self.pgeom2,self.hopping2,self.etba2,
                                                              self.pkpts3,self.pgeom3,self.hopping3,self.etba3)
            elif self.nsystem == 4:
                pyfit.eig4(self.fcomm, self.pinpt, self.ppram, self.pkpts ,self.pgeom ,self.hopping ,self.etba ,
                                                              self.pkpts2,self.pgeom2,self.hopping2,self.etba2,
                                                              self.pkpts3,self.pgeom3,self.hopping3,self.etba3,
                                                              self.pkpts4,self.pgeom4,self.hopping4,self.etba4)
            elif self.nsystem == 5:
                pyfit.eig5(self.fcomm, self.pinpt, self.ppram, self.pkpts ,self.pgeom ,self.hopping ,self.etba ,
                                                              self.pkpts2,self.pgeom2,self.hopping2,self.etba2,
                                                              self.pkpts3,self.pgeom3,self.hopping3,self.etba3,
                                                              self.pkpts4,self.pgeom4,self.hopping4,self.etba4,
                                                              self.pkpts5,self.pgeom5,self.hopping5,self.etba5)

        elif sys is not None:
            if sys == 1:
                pyfit.eig(self.fcomm, self.pinpt, self.ppram, self.pkpts,self.pgeom,self.hopping,self.etba)
            elif sys == 2:
                pyfit.eig(self.fcomm, self.pinpt, self.ppram, self.pkpts2,self.pgeom2,self.hopping2,self.etba2)
            elif sys == 3:
                pyfit.eig(self.fcomm, self.pinpt, self.ppram, self.pkpts3,self.pgeom3,self.hopping3,self.etba3)
            elif sys == 4:
                pyfit.eig(self.fcomm, self.pinpt, self.ppram, self.pkpts4,self.pgeom4,self.hopping4,self.etba4)
            elif sys == 5:
                pyfit.eig(self.fcomm, self.pinpt, self.ppram, self.pkpts5,self.pgeom5,self.hopping5,self.etba5)

        self.pinpt.flag_tbfit = flag_tbfit

    def generate_TBdata(self, filenm='tb_data.pt', ndata=1000, N_fit=50, tol=1e-5, method='lmdif', myid=0):
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        for i in range(ndata):  # assumes that each loop, init loads "random" parameters. -> check PARAM_FIT.dat with 'r' option
            self.init(myid=myid) # it assumes that fitting is done
            self.fit(miter = N_fit, tol=tol, method=method)
            band_fitted = torch.Tensor(self.etba.e).view(1,1,len(self.etba.e),-1).to(device)
            if i == 0:
                band_target = torch.Tensor(self.edft.e).view(1,1,len(self.edft.e), -1).to(device)
                band = band_fitted
                targ = band_target
            else:
                band = torch.cat((band_fitted, band), 0)
            if myid == 0:
                print(f'band:{i}, {band_all.size()}, {targ.size()}')
                sys.stdout.flush()
        if myid == 0:
            tb_data = {'TBA': band, 'DFT': targ}
            torch.save(tb_data, filenm)
            
    def copy_param_best(self, imode=1):
        pyfit.copy_params_best(self.ppram, imode)
       #NOTE:
       #if(imode .eq. 1) then
       #  PPRAM_PY%param_best = PPRAM_PY%param 
       #elseif(imode .eq. 2) then
       #  PPRAM_PY%param = PPRAM_PY%param_best
       #endif 

    def toten(self, eltemp=None, nelect=None):
        # rule for nelect: counting number of states upto certain level of energy
        # -> nelect (no soc) * 2 = nelect (soc)

        if eltemp is not None:
            self.pinpt.electronic_temperature = eltemp
        
        if nelect is not None:
            self.pgeom.nelect=nelect
        else:
            if self.pgeom.nelect[1] <= -1:
                print("The total nergy calculation is requested but number of electons are not explicitly specified.")
                print("Please set 'nelect'. current nelect=",nelect)
                print("Exit...",dir())
                exit()
        pyfit.toten(self.fcomm, self.pinpt, self.pkpts, self.pgeom, self.etba)

    def get_kpath(self, pkpts) :
        KNAME=[]
        KPTS=[]
        kdist = pkpts.kdist
        for i in range( pkpts.nline + 1):
            kname=str( pkpts.k_name2.T[i].view('S8'), 'utf-8').strip()
            KNAME.append(r"${}$".format(kname))
           #KNAME.append( str( pkpts.k_name2.T[i].view('S8'), 'utf-8').strip() )
            KPTS.append(pkpts.kdist[pkpts.k_name_index[i]-1])
        return kdist, KNAME, KPTS

    def plot_line(self,kdist,en, ax=None, ls='-',lw=.6, lc='black', label=None, args=None):
        if ax is None: ax = plt.gca()
        idx = np.where( (en >= args['yin']) & (en <= args['yen']))[0]
        font=args['font'] if 'font' in args else None

        for ie in range(idx.min(), idx.max()+1):
            mylabel = label if ie == idx.min() else None
            ax.plot(kdist, en[ie], linestyle=ls, linewidth=lw, color=lc, label=mylabel)
            
        if ('yin' in args) and ('yen' in args) :
            ax.set_ylim(args['yin'],args['yen'])
            ystep = args['ystep'] if 'ystep' in args else 5.0
            ax.set_yticks(np.arange(args['yin'], args['yen']+.000001, step=ystep))
            ax.tick_params(axis='y',labelsize=font['size'])
        ax.set_xlim(0, kdist[-1])
        if 'kpts'   in args : ax.set_xticks(args['kpts'])
        if 'kname'  in args : ax.set_xticklabels(args['kname'],fontdict=font)
        if 'xlabel' in args : ax.set_xlabel(args['xlabel'],fontdict=font)
        if 'ylabel' in args : ax.set_ylabel(args['ylabel'],fontdict=font)
        if 'title'  in args : ax.set_title(args['title'],fontdict=font)
        return ax

    def plot_dots(self, kdist, en, ax=None, wt=None, lc='black', marker='o', ms=1, lw=.5, label=None, args=None):
        if ax is None: ax = plt.gca()
        if (type(wt) == type(en)) and (wt.size == en.size):
            flag_wt_array = True
            maxwt = np.max(wt)
           #print("MAXWT : ",maxwt)
        else:
            flag_wt_array = False
        font=args['font'] if 'font' in args else None
       #marker_size = args['ms'] if 'ms' in args else 1.0
        marker_size = ms
        yin=args['yin'] if 'yin' in args else -20.0
        yen=args['yen'] if 'yen' in args else  10.0
        idx = np.where( (en >= yin) & (en <= yen))[0]
        for ie in range(idx.min(), idx.max()+1):
            if ie == idx.min():
                mylabel = r'${}$'.format(label) 
            else:
                mylabel = None
            if flag_wt_array is True:
                ax.scatter(kdist, en[ie], s=wt[ie,:]/maxwt*marker_size, color=lc, 
                           label=mylabel, marker=marker, linewidth=lw)
            else:
                ax.scatter(kdist, en[ie], s=wt, color=lc, 
                           label=mylabel, marker=marker, linewidth=lw)

        ax.set_ylim(yin,yen)
        ystep = args['ystep'] if 'ystep' in args else 5.0
        ax.set_yticks(np.arange(yin, yen+.000001, step=ystep))
        ax.tick_params(axis='y',labelsize=font['size'])
        ax.set_xlim(0, kdist[-1])
        if 'kpts'   in args : ax.set_xticks(args['kpts'])
        if 'kname'  in args : ax.set_xticklabels(args['kname'],fontdict=font)
        if 'xlabel' in args : ax.set_xlabel(args['xlabel'],fontdict=font)
        if 'ylabel' in args : ax.set_ylabel(args['ylabel'],fontdict=font)
        if 'title'  in args : ax.set_title(args['title'],fontdict=font)
        return ax

    def get_proj_weight(self, v2=None, i_atoms=None, i_specs=None, proj=None):
        atoms=None # atom info (integer list), atom index starts from 1, not from 0
        specs=None # specis info (character list)
        corbs=None # orbital info (character list)
        wt   = np.zeros(v2.shape[1:])

        # note: if atoms is not None, specs is ignored

        try:
            atoms=proj['atoms']
        except KeyError:
            try:
                specs=proj['specs']
            except KeyError:
                pass
        finally:
            try:
                corbs = proj['orb']
            except KeyError:    
                pass
              
        if (atoms is None) and (specs is None) and (corbs is None):
            print(" No proj target is provided or inappropriatly specified 'proj' dictionary!")
            print(" PROVIDED proj element: ", proj)
            print(" Ex):  projs = {'proj1':{'atoms':[1,2,5,6,...], 'orb':['px','pz',...]},")
            print("                'proj2':{'atoms':[4,8,9,10,..], 'orb':['px','pz',...]},")
            print("                'proj3':{'specs':['In1','In2'], 'orb':[     'pz'    ]}}")
            print(" Important note: atoms index should start from 1 not 0.")
           #exit()
            return wt

        if (atoms is not None) and (specs is None):
            for iatom in atoms:
                if corbs is not None:
                    for corb in corbs:
                        key  = 'AT'+str(iatom-1)+':'+corb
                        ikey = i_atoms[key]
                        wt   = wt + sum(v2[ikey])
                else:
                    key   = 'AT'+str(iatom)+':'
                    ikey_ = [val for keys, val in i_atoms.items() if key in keys]
                    ikey  = [] ; [ikey.extend(i) for i in ikey_]
                    wt    = wt + sum(v2[ikey])

        elif (specs is not None) and (atoms is None):
            for spec in specs:
                if corbs is not None:
                    for corb in corbs:
                        key  = spec+':'+corb
                        ikey = i_specs[key]
                        wt   = wt + sum(v2[ikey])
                else:
                    key   = spec+':'
                    ikey_ = [val for keys, val in i_specs.items() if key in keys]
                    ikey  = [] ; [ikey.extend(i) for i in ikey_] 
                    wt    = wt + sum(v2[ikey])

        elif (atoms is None) and (specs is None):
            for corb in corbs:
                key   = ':'+corb
                ikey_ = [val for keys, val in i_specs.items() if key in keys]
                ikey  = [] ; [ikey.extend(i) for i in ikey_]
                wt    = wt + sum(v2[ikey])

        return wt

    def get_system(self,isystem=1):
        if isystem == 1:
            etarget=self.energy_target ; etba = self.etba.e ; edft = self.edft.e ; wt = self.pwght.wt
            title   = self.pgeom.title.decode('utf-8').strip()
            kdist, KNAME, KPTS = self.get_kpath(self.pkpts)
        elif isystem == 2:
            etarget=self.energy_target2 ; etba = self.etba2.e ; edft = self.edft2.e ; wt = self.pwght2.wt
            title   = self.pgeom2.title.decode('utf-8').strip()
            kdist, KNAME, KPTS = self.get_kpath(self.pkpts2)
        elif isystem == 3:
            etarget=self.energy_target3 ; etba = self.etba3.e ; edft = self.edft3.e ; wt = self.pwght3.wt
            title   = self.pgeom3.title.decode('utf-8').strip()
            kdist, KNAME, KPTS = self.get_kpath(self.pkpts3)
        elif isystem == 4:
            etarget=self.energy_target4 ; etba = self.etba4.e ; edft = self.edft4.e ; wt = self.pwght4.wt
            title   = self.pgeom4.title.decode('utf-8').strip()
            kdist, KNAME, KPTS = self.get_kpath(self.pkpts4)
        elif isystem == 5:
            etarget=self.energy_target5 ; etba = self.etba5.e ; edft = self.edft5.e ; wt = self.pwght5.wt
            title   = self.pgeom5.title.decode('utf-8').strip()
            kdist, KNAME, KPTS = self.get_kpath(self.pkpts5)

        return etarget, etba, edft, wt, title, kdist, KNAME, KPTS

    def get_proj_ldos(self, isystem=1):
        if isystem == 1:
            proj_ldos = self.etba.v2
        elif isystem == 2:
            proj_ldos = self.etba2.v2
        elif isystem == 3:
            proj_ldos = self.etba3.v2
        elif isystem == 4:
            proj_ldos = self.etba4.v2
        elif isystem == 5:
            proj_ldos = self.etba5.v2

        return proj_ldos

    def set_proj_c(self, projs):
        proj_c=[]
        proj_m=[]
        i = 1
        orbs = list(mycolor.color_index)
        for proj in projs:
            proj_c.append(mycolor.color_index[orbs[i]][0])
            proj_m.append(mycolor.marker_index[orbs[i]])
            i = i + 1
            
        return proj_c, proj_m

    def set_proj_m(self,projs):
        proj_m=[]
        i = 1
        orbs = list(mycolor.color_index)
        for proj in projs:
            proj_m.append(mycolor.marker_index[orbs[i]])
            i = i + 1
        return proj_m

    def get_proj_index(self,isystem=1):
        if isystem == 1:
            i_specs = self.i_specs  ; i_atoms = self.i_atoms
        if isystem == 2:
            i_specs = self.i_specs2 ; i_atoms = self.i_atoms2
        if isystem == 3:
            i_specs = self.i_specs3 ; i_atoms = self.i_atoms3
        if isystem == 4:
            i_specs = self.i_specs4 ; i_atoms = self.i_atoms4
        if isystem == 5:
            i_specs = self.i_specs5 ; i_atoms = self.i_atoms5

        return i_specs, i_atoms

    def get_grid_ratio(self,nsystem=1,sys=None):
        grid_ratio = []
        for isystem in range(nsystem):
            if sys is not None:
                if isystem+1 is not sys:
                    continue
            if isystem+1 == 1:
                grid_ratio.append(self.pkpts.kdist[-1])
            elif isystem+1 == 2:
                grid_ratio.append(self.pkpts2.kdist[-1])
            elif isystem+1 == 3:
                grid_ratio.append(self.pkpts3.kdist[-1])
            elif isystem+1 == 4:
                grid_ratio.append(self.pkpts4.kdist[-1])
            elif isystem+1 == 5:
                grid_ratio.append(self.pkpts5.kdist[-1])

        return grid_ratio

    def plot_band(self, figsize=(5,6), ef=0.0, yin=-20., yen=10, ystep=5., plot_target=True,
                 title=' ', ylabel='Energy (eV)', xlabel=' ', fout='band.pdf', ms=1,
                 lw=0.6, marker_lw=0.5, lc='red', projs=None, proj_c=None, sys=None,
                 font_size=16, legend_loc='lower right', legend_font_size = 10, legend_marker_size=30):

        font=myfont.font ; font['size']=font_size
       #args={'ylabel':ylabel,'xlabel':xlabel,'yen':yen,'yin':yin,'ystep':ystep,'ms':ms,'font':font}
        args={'ylabel':ylabel,'xlabel':xlabel,'yen':yen,'yin':yin,'ystep':ystep,'ms':ms,'font':font}
        gr = self.get_grid_ratio(nsystem=self.nsystem,sys=sys)
        if sys is not None or self.nsystem == 1:
            fig, axes = plt.subplots(1,1, figsize=figsize, gridspec_kw={'width_ratios': gr})
        else:
            fig, axes = plt.subplots(1,self.nsystem, figsize=figsize, gridspec_kw={'width_ratios': gr})
        if (projs is not None): 
            if proj_c is None:
                proj_c, proj_m=self.set_proj_c(projs)
            else:
                proj_m = self.set_proj_m(projs)

        for isystem in range(self.nsystem):
            if sys is not None:
                if isystem+1 is not sys:
                    continue
                elif isystem+1 is sys:
                    myax = axes
                    isys = 0
            else:
                if self.nsystem > 1:
                    myax = axes[isystem]
                else:
                    myax = axes
                isys = isystem
            if self.pinpt.iverbose == 1:
                print(" --> PLOT system: ", isystem+1,isys)
            etarget,etba,edft,wt,mysys,kdist,args['kname'],args['kpts']=self.get_system(isystem=isystem+1)
            args['title']=title.strip()+':'+mysys
            if plot_target is True and bool(self.pinpt.flag_tbfit):
                ax = self.plot_line(kdist, etarget, myax, ls='-', lw=lw,lc='lightgray',label='DFT',args=args)
                sc = self.plot_dots(kdist, edft, myax, wt=wt, lc='black', marker='o', ms=ms, label="weight",args=args)
            ax = self.plot_line(kdist,etba, myax, ls='-',lw=lw,lc=lc,label='TBA',args=args)
            if (projs is not None) and bool(self.pinpt.flag_get_orbital):
                # NOTE: proj = {'proj_set1':{'atoms':[1,2,5,6,...], 'orb':['px','pz',...]},
                #               'proj_set2':{'atoms':[4,8,9,10,..], 'orb':['px','pz',...]},
                #               'proj_set3':{'specs':['In1','In2'], 'orb':[     'pz'    ]}}
                i_specs, i_atoms = self.get_proj_index(isystem=isystem+1)
                proj_c_it = iter(proj_c)
                proj_m_it = iter(proj_m)
                proj_ldos = self.get_proj_ldos(isystem=isystem+1)
                for proj in projs: 
                    wt_proj = self.get_proj_weight(v2=proj_ldos, i_atoms=i_atoms, i_specs=i_specs, proj=projs[proj])
                    sc = self.plot_dots(kdist, etba, myax, wt=wt_proj, label=str(proj),
                                        lw=marker_lw, marker=next(proj_m_it), ms=ms*2, lc=next(proj_c_it), args=args)

        fig.tight_layout(pad=.05)
#       plt.legend()
        if (projs is not None) and bool(self.pinpt.flag_get_orbital):
            lgnd = plt.legend(loc=legend_loc, scatterpoints=1, fontsize=legend_font_size,
                              prop={'family':'sans-serif', 'size':legend_font_size})
            ii = 0
           #for _ in projs:
            for _ in lgnd.legendHandles:
                lgnd.legendHandles[ii]._sizes = [legend_marker_size] # adjusting marker legend size for s
                ii = ii + 1

        plt.savefig(fout,bbox_inches='tight',pad_inches=0) 
        if self.pinpt.iverbose == 1:
            print("\n #- Plotting band structure: ",fout, " is written")
        plt.clf()

    def plot_pso_cost_history(self, figsize=(5,6), yin=-20,yen=20, ystep=5.,
                                    title=' ', ylabel='Best cost function', xlabel='PSO iteration', 
                                    fout='COST_HISTORY.pdf', bestn=5):
        yin  = np.min(self.cost_history_particle)
        yen  = np.max(self.cost_history_particle)
        delta = yen  - yin 
        yin  = yin - delta*0.05
        yen  = yen + delta*0.05
        plt.figure(figsize=figsize)
        plt.ylim(yin, yen)
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.title(title.strip())
        best_particle_index = self.cost_history_particle[-1,:].argmin()
        worst_particle_index = self.cost_history_particle[-1,:].argmax()
        best_nparticle_index = self.cost_history_particle[-1,:].argsort()[:bestn]
        worst_nparticle_index = self.cost_history_particle[-1,:].argsort()[::-1][:bestn]
        cmap = mpl.cm.GnBu
        cmap_best = mpl.cm.get_cmap('autumn_r')
        cmap_worst= mpl.cm.get_cmap('Greys_r')

        for i in range( self.n_particles ):
            if i == worst_particle_index :
                lw = 1
                zorder = 555
                color  = 'black'
                mylabel= None
            elif i in best_nparticle_index:
                lw = 1
                order, = np.where(best_nparticle_index == i) ; order = int(order)
                zorder = 999-order
                color = cmap_best(            order/float(bestn))
                if order < 5 :
                    mylabel = 'best:'+str(order)+'-th,ID:'+str(i)
                    plt.text(-self.ppram.niter*0.05, self.cost_history_particle[0,i], str(i), fontsize=5, zorder=zorder)
                    if order == 0:
                        lw = 2
                else:
                    mylabel = None
            else:
                lw = 1
                zorder = 1
                color  = cmap(i/ float(self.n_particles))
                mylabel= None

            plt.plot( range( self.ppram.niter ), 
                      self.cost_history_particle[:,i],'-',linewidth=.6*lw,color= color,zorder=zorder, label=mylabel)
        plt.legend()
        plt.savefig(fout,bbox_inches='tight',pad_inches=0)
        plt.clf()

    def plot_pso_pbest(self, figsize=(5,6), yin=-20, yen=20, ystep=5.,
                             title=' ', ylabel='PSO Iteration', xlabel='Parameter index', 
                              fout='PARAM_particles.pdf', bestn=5):
        zin  = np.min(self.ppram.param)
        zen  = np.max(self.ppram.param)
        delta = zen  - zin
        zin  = zin - delta*0.05
        zen  = zen + delta*0.05

        best_particle_index = self.cost_history_particle[-1,:].argmin()
        worst_particle_index = self.cost_history_particle[-1,:].argmax()
        best_nparticle_index = self.cost_history_particle[-1,:].argsort()[:bestn]
        worst_nparticle_index = self.cost_history_particle[-1,:].argsort()[::-1][:bestn]
        cmap = mpl.cm.GnBu
        cmap_best = mpl.cm.get_cmap('autumn_r')
        cmap_worst= mpl.cm.get_cmap('Greys_r')

        plt.figure(figsize=figsize)
        plt.ylim(zin, zen)
        plt.ylabel('Parameter values')
        plt.xlabel(xlabel)
        plt.title(title.strip())
        for i in range( self.n_particles ):
            if i in best_nparticle_index:
                lw = 1
                order, = np.where(best_nparticle_index == i) ; order = int(order)
                zorder = 999-order
                color = cmap_best( order/float(bestn))
                if order < 5:
                    mylabel = 'best:'+str(order)+'-th,ID:'+str(i)
                    plt.text((-self.ppram.nparam_free)*0.05, self.ppram.pso_pbest_history[-1,i,0], str(i), fontsize=5, zorder=zorder)
                    if order == 0:
                        lw = 2
                else:
                    mylabel = None
            elif i in worst_nparticle_index:
                lw = 1
                order, = np.where(worst_nparticle_index == i) ; order = int(order)
                zorder = 555-order
                color = cmap_worst( order/float(bestn))
                mylabel= None
                if order == 0:
                    lw = 2
            else:
                lw = 1
                zorder = 1
                color  = cmap(i/ float(self.n_particles))
                mylabel= None

            plt.plot( range( self.ppram.nparam_free ), 
                self.ppram.pso_pbest_history[-1,i,:],'-',linewidth=.6*lw,color= color ,zorder=zorder, label=mylabel)
        plt.legend()
        plt.savefig(fout,bbox_inches='tight',pad_inches=0)
        plt.clf()

    def print_param(self, param_out=None):
        if param_out is not None:
            pfileoutnm = param_out+' '*(132-len(param_out))
        else:
            pfileoutnm = self.ppram.pfileoutnm

        pyfit.print_param_py(self.pinpt, self.ppram, pfileoutnm)
        if self.pinpt.iverbose == 1:
            print("\n #- Writing parameter file: ",pfileoutnm.strip(), " is written")

    def print_weight(self, weight_out=None):
        if weight_out is not None:
           wfileoutnm = weight_out+' '*(132-len(weight_out))
        else:
           wfileoutnm = 'WEIGHT.dat'+' '*(132-len('WEIGHT.dat'))

        pyfit.print_weight(self.pwght, wfileoutnm)
        if self.pinpt.iverbose == 1:
            print("\n #- Writing weight information file: ",wfileoutnm.strip(), " is written")

    def load_weight(self, weight_in=None):
        if weight_in is not None:
            wfileinnm = weight_in 
        else:
            wfileinnm = 'WEIGHT.dat'

        self.pwght.wt = np.loadtxt(wfileinnm).T

    def select_param(self, param_set, i):
        k = round(self.ppram.param_const[0,i])
        if k == 0 :
            ifix = round(self.ppram.param_const[3, i])
            if ifix == 1 :
                params =  self.ppram.param_const[4,i]
            else:
                params =  param_set[i]
        elif k >= 1 :
            ifix = round(self.ppram.param_const[3, i])
            if ifix == 1:
                params = self.ppram.param_const[4,i]
            else:
                params = param_set[ k-1 ]
        return params

    def print_fit(self, suffix='None'):
        # print target band structure
        if suffix is not None:
            tfileoutnm = 'band_structure_DFT.'+suffix+'.dat' + ' '*(132-len('band_structure_DFT.'+suffix+'.dat'))
        else:
            tfileoutnm = 'band_structure_DFT.dat' + ' '*(132-len('band_structure_DFT.dat'))
        
        pyfit.print_target(self.pinpt, self.pkpts, self.edft, self.pwght, self.pgeom, tfileoutnm)

        # print fitted band structure
        if suffix is not None:
            suffixx = '.'+suffix + ' '*(132-len(suffix)+1)
        else:
            suffixx = ' '*132

        use_overlap = self.ppram.flag_use_overlap
        pyfit.print_etba(self.pinpt, self.pkpts, self.etba, self.edft, self.pwght, self.pgeom, 
                         suffixx, use_overlap)
    
    def save(self, title=' ', plot_band=None, plot_target=False, target=None, band=None, param=None, weight=None, 
                   cost_history=None, cost_history_particle=None , pso_param_history = None, yin=-20., yen=10, ystep=5.,
                   figsize=(5,6),lw=.6, marker_lw=.5, marker_size=1, lc='red', projs=None, proj_c=None, sys=None,
                   font_size=16, legend_loc='lower right', legend_font_size = 10, legend_marker_size=30):
        if self.myid != 0:
            return

        if plot_band is True:
            if title is not None:
                fout = 'BAND.'+title.strip()+'.pdf'
            else:
                fout = 'BAND.pdf'
            self.plot_band(title=title.strip(), figsize=figsize, fout=fout, yin=yin, yen=yen, ystep=ystep, 
                          plot_target=plot_target, lw=lw, marker_lw=marker_lw, ms=marker_size, 
                          lc=lc, projs=projs, proj_c=proj_c, sys=sys,
                          font_size=font_size, legend_loc=legend_loc, legend_font_size = legend_font_size, 
                          legend_marker_size=legend_marker_size)

        if target is True:
            if title is not None:
                tfileoutnm = 'band_structure_DFT.'+title.strip()+'.dat' + ' '*(132-len('band_structure_DFT.'+title.strip()+'.dat'))
            else:
                tfileoutnm = 'band_structure_DFT.dat' + ' '*(132-len('band_structure_DFT.dat'))
            pyfit.print_target(self.pinpt, self.pkpts, self.edft, self.pwght, self.pgeom, tfileoutnm)

        if band is True:
            if title is not None:
                suffix = '.'+title.strip() + ' '*(132-len(title.strip())+1)
            else:
                suffix = ' '*132
            use_overlap = self.ppram.flag_use_overlap
            pyfit.print_etba(self.pinpt, self.pkpts, self.etba, self.edft, self.pwght, self.pgeom,
                             suffix, use_overlap)

        if param is True:
            if (title is not None) and (title.strip() != '') :
                param_out = 'PARAM_FIT.new.'+title.strip()+'.dat'
            else:
                param_out = 'PARAM_FIT.new.dat'
            self.print_param(param_out=param_out)
    
        if weight is True:
            if title is not None:
                weight_out = 'WEIGHT.'+title.strip()+'.dat'
            else:
                weight_out = 'WEIGHT.dat'
            self.print_weight(weight_out = weight_out)

        if cost_history is True:
            if title is not None:
                cost_out_dat = 'COST_HISTORY.'+title.strip()+'.dat'
              # cost_out_pdf = 'COST_HISTORY.'+title.strip()+'.pdf'
            else:
                cost_out_dat = 'COST_HISTORY.dat'
              # cost_out_pdf = 'COST_HISTORY.pdf'
            if self.ppram.niter == 0:
                niter = self.ppram.pso_miter
            else:
                niter = self.ppram.niter

           # deactivated due to the error in Infant-fzj: HJKim (16.09.2021)
           #plot_cost_history(self.cost_history[:niter], title=title.strip())
           #plt.savefig(cost_out_pdf, bbox_inches='tight', pad_inches=0)

            cost_out = open(cost_out_dat, 'w+')
            for ii in range(niter):
                cost_out.write(" %d    %.9f \n" %(ii, self.cost_history[ii]))
            cost_out.close()

        if cost_history_particle is True:
            if title is not None:
                cost_out_dat = 'COST_HISTORY_particles.'+title.strip()+'.dat'
                cost_out_pdf = 'COST_HISTORY_particles.'+title.strip()+'.pdf'
            else:
                cost_out_dat = 'COST_HISTORY_particles.dat'
                cost_out_pdf = 'COST_HISTORY_particles.pdf'
            if self.ppram.niter == 0:
                niter = self.ppram.pso_miter
            else:
                niter = self.ppram.niter
            
            bestn = int(self.n_particles * 0.2)
            if(bestn > 5 ): bestn = 5 
            self.plot_pso_cost_history(title=title.strip(), fout = cost_out_pdf, bestn=bestn)

            cost_out = open(cost_out_dat, 'w+')
            for jj in range(self.n_particles):
                cost_out.write("# Particle : %d  \n" %(jj))
                for ii in range(niter):
                    cost_out.write("%d  %.9f \n"%(ii, self.cost_history_particle[ii,jj]))
                cost_out.write("  \n")
            cost_out.close()

        if pso_param_history is True:
            if title is not None:
                param_out_dat = 'PARAM_best_particles.'+title.strip()+'.dat'
                param_out_pdf = 'PARAM_best_particles.'+title.strip()+'.pdf'
            else:
                param_out_dat = 'PARAM_best_particles.dat'
                param_out_pdf = 'PARAM_best_particles.pdf'
            if self.ppram.niter == 0:
                niter = self.ppram.pso_miter
            else:
                niter = self.ppram.niter

            bestn_ = int(self.n_particles * 0.2)
            if(bestn_ > 5): # and self.n_particles > 5 : 
                bestn = 5 
            else:
                bestn = bestn_

            self.plot_pso_pbest(title=title.strip(), fout = param_out_pdf, bestn=bestn)

            best_nparticle_index = self.cost_history_particle[-1,:].argsort()[:bestn_]
            param_out = open(param_out_dat, 'w+')
            for i in range( bestn_):
                ii = best_nparticle_index[i]
                cost = self.cost_history_particle[-1,ii]
                param_out.write("# Best Particle : rank=%d , ID=%d, COST=%.9f  \n" %(i, ii, cost))
                param_set = self.ppram.param 
                for ip in range(self.ppram.nparam_free):
                    param_set[self.ppram.iparam_free[ip]-1] = self.ppram.pso_pbest_history[-1,ii,ip]
                for j in range( self.ppram.nparam):
                    pj = self.select_param(param_set, j)
                    param_out.write("%5d  %20.9f   # %5d %20.9f\n" %(j+1, pj, ii, cost))
                param_out.write("  \n")
            param_out.close()

    
    def plot_param(self, fname=None, nplot=None, iplot=0, nparam = None, diffmin = 10.0,
                   fout='PARAM_best_particles.pdf',
                   title='Best particles',
                   figsize=(10,6), 
                   plot = False) :
        self.pset_id = []
        f = open(fname)
        data = np.loadtxt(f) ; data = data.reshape(data.shape)
        zin  = np.min(data[:,1])
        zen  = np.max(data[:,1])
        if nparam is None:
            nparam = int(np.max(data[:,0]))

        if nplot is None:
            nplot = int((data.shape[0]/nparam) * 0.1)

        #cmap = mpl.cm.BuGn
        cmap = mpl.cm.tab20c
        plt.style.use('dark_background')
        plt.figure(figsize=figsize)
        plt.ylim(zin, zen)
        plt.ylabel('Parameter values')
        plt.xlabel('Parameter index')
        plt.title(title.strip())

        ii = 0 ; jj = 0
        for iset in range(iplot, iplot+nplot):
            params= data[iset*nparam:(iset+1)*nparam].reshape(nparam,2)[:,1]
            if iset == iplot:
                params_old = params
                ii += 1
            diff  = np.linalg.norm(params_old - params) / np.linalg.norm(params_old) * 100.0
            if diff >= diffmin:
                ii += 1
            params_old = params

        self.pset = np.ndarray((ii, nparam))
        
        for iset in range(iplot, iplot+nplot):
            lw = 1
            zorder = nplot - iset
            mylabel = None
            index = range(nparam)
            params= data[iset*nparam:(iset+1)*nparam].reshape(nparam,2)[:,1]
            
            if iset == 0:
                params_old = params
                self.pset_id.append(iset)

            diff  = np.linalg.norm(params_old - params) / np.linalg.norm(params_old) * 100.0
            if diff >= diffmin:
                mylabel = 'PID: '+str(iset)
                self.pset_id.append(iset)
                color = cmap( float(jj)/float(ii)      )
                plt.plot(index, params, '-', linewidth=.6*lw, color = color, zorder=zorder, label=mylabel)
                self.pset[jj,:] = params
                jj += 1
            elif iset == iplot:
                mylabel = 'PID: 0'
                color = cmap(  0                       )
                plt.plot(index, params, '-', linewidth=.6*lw, color = color, zorder=zorder, label=mylabel)
                self.pset[jj,:] = params
                jj += 1
            else:
                color = color
        
            params_old = params


        plt.legend(ncol=6)
        plt.savefig(fout,bbox_inches='tight', pad_inches=0)
        plt.clf()
        return jj

class csa_tools:
    def __init__(self,Info):
        self.npop = Info.npop
        self.ndim = Info.ndim

        try:
            self.apath = Info.apath.rstrip()
        except AttributeError:
            self.apath = os.getcwd()
        if self.apath[-1] == '/' :
            self.apath=self.apath[:-1]

        self.rmin = -20.0
        self.rmax =  20.0

        try:
            self.fname_dump = Info.fname_dump
        except AttributeError:
            self.fname_dump = 'dump.txt'

        try:
            self.ntotal = max(Info.nmax,self.npop + self.npop*self.ndim)
        except AttributeError:
            self.ntotal = self.npop + 500 * self.npop * self.ndim

        try:
            self.csa_soldier_command=Info.soldier_command.rstrip()
            if self.csa_soldier_command[-1] != '&' :
                self.csa_soldier_command = self.csa_soldier_command + ' &'
        except AttributeError:
            self.csa_soldier_command='srun --nodes=1  --ntasks=1 --exclusive python CSA_SOLDIER.py &'
        # Note: 
        #       for iffslurm with sinteractive:
        #         srun --mpi=pmi2 --ntasks=1 --exclusive python CSA_SOLDIER.py &
        #       for iffslurm with sbatch :     
        #         srun --nodes=1  --ntasks=1 --exclusive python CSA_SOLDIER.py &
        
    def init(self,ntotal = None):
        cpath = self.apath+'/bnds.txt'
        self.bnds = self.get_bnds(cpath)
        npop  = self.npop
        ndim  = self.ndim
        self.static=np.zeros((npop,ndim)) 
        self.ostatic=np.zeros(npop)
        self.xvector=np.zeros(ndim)
        self.dynamic=np.zeros((npop,ndim))
        self.odynamic=np.zeros(npop) 
        self.yvector=np.zeros(ndim)
        self.staticr=np.zeros((npop,ndim))
        print(npop,' npop',flush=True)
        print(self.apath,' apath',flush=True)
        print(ndim,' ndim',flush=True)

    def append_multiple_lines(self, file_name, lines_to_append):
        with open(file_name, "a+") as file_object:
            appendEOL = False
            file_object.seek(0)
            data = file_object.read(100)
            if len(data) > 0:
                appendEOL = True
            for line in lines_to_append:
                if appendEOL == True:
                    file_object.write("\n")
                else:
                    appendEOL = True
                file_object.write(line)

    def repulsion(self,npop,ndim,static):
        tmp=0.
        for i in range(npop):
             for j in range(npop):
                 if i > j :
                     tmq=0.
                     for k in range(ndim):
                         tmq=tmq+(static[i,k]-static[j,k])**2
                     tmq=np.sqrt(tmq)+1e-13
                     tmp=tmp+1./tmq
        return tmp

    def coulomb(self, old):
        npop = self.npop
        ndim = self.ndim
        bnds = self.bnds
        work=np.zeros((npop,ndim)) ; new=np.zeros((npop,ndim))
        new[:,:]=old[:,:]
        tmp0=self.repulsion(npop,ndim,old)
        for kter in range(npop+10):
            for i in range(1,npop):
                for k in range(ndim):
                    rmin,rmax=bnds[k]
                    work[i,k]=new[i,k]+(rmax-rmin)*np.random.random()/4.
                    if work[i,k] < rmin or work[i,k] > rmax:
                        work[i,k]=rmin+(rmax-rmin)*np.random.random()
            tmp=self.repulsion(npop,ndim,work)
            if tmp < tmp0 :
               tmp0=tmp ; new[:,:]=work[:,:]
        return new

    def gen_directories(self):
        npop = self.npop
        apath = self.apath
        for i in range(npop):
             astring='mkdir '+apath+'/'+str(i).zfill(4) ; os.system(astring)
        for i in range(npop):
             astring='cp '+apath+'/CSA_SOLDIER.py'+' '+apath+'/'+str(i).zfill(4)
             os.system(astring)
             if os.path.isfile(apath+'/bnds.txt') :
                 astring='cp '+apath+'/bnds.txt'+' '+apath+'/'+str(i).zfill(4)
                 os.system(astring)

    def del_directories(self):
        npop = self.npop
        apath = self.apath
        for i in range(npop):
            astring='rm '+apath+'/'+str(i).zfill(4) ; os.system(astring)

    def partial_shuffle(self,l,nitera=5):
        n21=len(l)-1
        for _ in range(nitera):
            a, b = random.randint(0, n21), random.randint(0, n21)
            l[b], l[a] = l[a], l[b]
        return l

    def load_solutions(self,fname='dump.txt'):
        afile=open(fname,'r')
        jline=0
        for line in afile:
            if jline == 0 :
                npop=int(line.split()[0]) ; ndim=int(line.split()[1])
                static=np.zeros((npop,ndim)) ; ostatic=np.zeros(npop)
                dynamic=np.zeros((npop,ndim)) ; odynamic=np.zeros(npop)
                i=0
            if jline > 0 and jline  <= npop :
                if i < npop and i > -1 :
                    for j in range(ndim):
                        static[i,j]=float(line.split()[j])
                i=i+1
                if i == npop :
                    i=0
            if jline > npop and jline  <= npop*2 :
                if i < npop  and i > -1 :
                    for j in range(ndim):
                        dynamic[i,j]=float(line.split()[j])
                i=i+1
            if jline == npop*2+1 :
                for i in range(npop):
                    ostatic[i]=float(line.split()[i])
            if jline == npop*2+2 :
                for i in range(npop):
                    odynamic[i]=float(line.split()[i])
            jline=jline+1
        afile.close()
        return npop,ndim,ostatic,odynamic,static,dynamic

    def gen_trial_solution(self,fname,ncal):
        ndim = self.ndim
        bnds = self.bnds
        alist=[ np.random.random() for j in range(ndim)]
        for j in range(ndim):
            rmin,rmax=bnds[j]
            alist[j]=rmin+(rmax-rmin)*np.random.random()
        fname100=self.fname_dump #'dump.txt'
        if os.path.isfile(fname100) :
            print(fname100,' is present, perturbation',flush=True)
            npop1,ndim1,ostatic1,odynamic1,static1,dynamic1=self.load_solutions(fname100)
            if ndim1 == ndim :
                if ncal >= 0 and ncal <= npop1-1:
                    alist=[ dynamic1[ncal,j] for j in range(ndim)]
                    for j in range(ndim):
                        rmin,rmax=bnds[j]
                        alist[j]=alist[j]+np.random.random()*(rmax-rmin)/4.
                        if alist[j] < rmin or alist[j] > rmax:
                             alist[j]=np.random.random()*(rmax-rmin)+rmin
            else :
               print('check ndim,ndim0 ')
               sys.exit()
        lines_to_append=[]
        lines_to_append.append(str(ndim))
        for j in range(ndim):
            lines_to_append.append(str(alist[j]))
        lines_to_append.append(str(ncal)+'\n')
        self.append_multiple_lines(fname, lines_to_append)

    def write_trial_solution(self,ndim0,xvector,iidd,ncal):
        apath=self.apath
        bnds =self.bnds
        fname=apath+'/'+str(iidd).zfill(4)+'/input.txt'
        gname=apath+'/'+str(iidd).zfill(4)+'/STATUS'
        if os.path.isfile(fname) :
            os.remove(fname)
        isign=1
        if ndim0 < 0:
           isign=-1 
           ndim=-ndim0
           self.gen_trial_solution(fname,ncal)
        if isign == 1:
            ndim=ndim0
            lines_to_append=[]
            lines_to_append.append(str(ndim))
            for j in range(ndim):
                rmin,rmax=bnds[j]
                if xvector[j] < rmin or xvector[j] > rmax:
                    xvector[j]=rmin+(rmax-rmin)*np.random.random()
            for j in range(ndim):
                lines_to_append.append(str(xvector[j]))
            lines_to_append.append(str(ncal)+'\n')
            self.append_multiple_lines(fname, lines_to_append)
        astring='cd '+apath+'/'+str(iidd).zfill(4)+' ; ' + self.csa_soldier_command
        os.system(astring)
        astring='echo "ING" >> '+gname
        os.system(astring)

    def get_solution(self):
        ncal0=0
        npop = self.npop
        ndim = self.ndim
        apath = self.apath
        bnds = self.bnds
        while True:
            idirectory0=np.random.randint(npop)
            fname=apath+'/'+str(idirectory0).zfill(4)+'/STOP'
            gname=apath+'/'+str(idirectory0).zfill(4)+'/STATUS'
            hname=apath+'/'+str(idirectory0).zfill(4)+'/output.txt'
            sol0=np.zeros(ndim)
            for j in range(ndim):
                rmin,rmax=bnds[j]
                sol0[j]=rmin+(rmax-rmin)*np.random.random()
            objective0=9e99
            os.system('sleep 0.1')
            if os.path.isfile(fname) :
                os.system('sleep 0.1')
                if os.path.isfile(gname):
                    astring=' '
                    afile=open(gname,'r')
                    for line in afile:
                        if len(line.split()) > 0:
                            astring=line.split()[0]
                    afile.close()
                    astring=astring.lower()
                    if astring == 'done' :
                        os.remove(fname) ; os.system('sleep 0.1')
                        if os.path.isfile(hname) :
                            afile=open(hname,'r')
                            jline=0
                            for line in afile:
                                if jline == 0:
                                   kkk=int(line.split()[0])
                                if jline <= ndim:
                                   sol0[jline-1]=float(line.split()[0])
                                if jline == ndim+1:
                                   objective0=float(line.split()[0])
                                if jline == ndim+2:
                                   ncal0=int(line.split()[0])
                                jline=jline+1
                            afile.close()
                            os.system('sleep 0.1')
                            os.remove(hname)
                        else:
                            os.system('sleep 0.1')
                        print(objective0,idirectory0,' loss',flush=True)
                        return objective0,sol0,idirectory0,ncal0

    def csadistance(self,x,y):
        tmp=0.
        for j in range(len(x)):
            tmp=tmp+np.abs(x[j]-y[j])
        return tmp

    def mutation(self,xvector):
        ndim = self.ndim
        bnds = self.bnds
        if np.random.random() < 0.10:
            j=np.random.randint(5)+1
            xvector=self.partial_shuffle(xvector, nitera=j)
        if np.random.random() < 0.05:
            xvector=np.flipud(xvector)
        tmq=np.random.random()*0.15+0.25
        for j in range(ndim):
            rmin,rmax=bnds[j]
            if np.random.random() < tmq:
                xvector[j]=xvector[j]+(np.random.random()-0.5)*2.*(rmax-rmin)/2.
        for j in range(ndim):
            rmin,rmax=bnds[j]
            if xvector[j] < rmin or xvector[j] > rmax :
                xvector[j]=rmin+np.random.random()*(rmax-rmin)
        return xvector

    def crossover(self,xvector,yvector):
        ndim = self.ndim
        bnds = self.bnds
        if np.random.random() < 0.05:
            xvector=np.flipud(xvector)
        tmq=np.random.random()*0.25+0.25
        for j in range(ndim):
            if np.random.random() < tmq:
                xvector[j]=yvector[j]
        for j in range(ndim):
            rmin,rmax=bnds[j]
            if xvector[j] < rmin or xvector[j] > rmax :
                xvector[j]=rmin+np.random.random()*(rmax-rmin)
        return xvector

    def selection(self):
        npop0=min(self.npop,50)
        k0=0 ; k1=0
        while k0 == k1 :
            i=np.random.randint(npop0) ;  j=np.random.randint(npop0)  ; k0=j
            if i < j :
                k0=i
            i=np.random.randint(npop0) ;  j=np.random.randint(npop0)  ; k1=j
            if i < j :
                k1=i
        return k0,k1

    def gen_parents(self):
        k0,k1=self.selection()
        return k0,k1

    def save_solutions(self):
        fname=self.fname_dump #'dump.txt'
        npop = self.npop
        ndim = self.ndim
        if os.path.isfile(fname) :
            os.remove(fname)
        lines_to_append=[]
        astring=str(npop)+' ' +str(ndim)
        lines_to_append.append(astring)
        for i in range(npop):
            astring=''
            for j in range(ndim):
                astring=astring+' '+str(self.static[i,j])
            lines_to_append.append(astring)
        for i in range(npop):
            astring=''
            for j in range(ndim):
                astring=astring+' '+str(self.dynamic[i,j])
            lines_to_append.append(astring)
        astring=''
        for i in range(npop):
                astring=astring+' '+str(self.ostatic[i])
        lines_to_append.append(astring)
        astring=''
        for i in range(npop):
                astring=astring+' '+str(self.odynamic[i])
        lines_to_append.append(astring)
        self.append_multiple_lines(fname, lines_to_append)

    def get_bnds(self,cpath):
        if os.path.isfile(cpath) :
            afile=open(cpath,'r')
            j=0 ; i0=0
            for line in afile:
                if j == 0:
                   ndim=int(line.split()[0])
                   if ndim != self.ndim :
                      print('check bnds.txt or ndim')
                      sys.exit()
                   self.ndim=ndim
                   bnds=[ (-1.,1.) for _ in range(ndim)]
                else :
                   if i0 < ndim:
                       bnds[i0]=(float(line.split()[0]),float(line.split()[1]))
                   i0=i0+1
                j=j+1
            afile.close()
        else:
            bnds=[(self.rmin,self.rmax) for _ in range(ndim)]
            print('no bnds.txt file is assumed. MIN and MAX = ',self.rmin,self.rmax)
        bnds=np.array(bnds)
        return bnds

    def csa_initial_step(self):
        nin=0 ; nout=0
        ndim = self.ndim
        npop = self.npop
        ntotal = self.ntotal

        for iidd in range(npop):
            if nin < ntotal :
                self.write_trial_solution(-ndim,self.xvector,iidd,nin) ; nin=nin+1
                for j in range(ndim):
                    rmin,rmax=self.bnds[j]
                    self.xvector[j]=rmin+(rmax-rmin)*np.random.random()
                for j in range(ndim):
                    self.dynamic[iidd,j]=self.xvector[j] ; self.static[iidd,j]=self.xvector[j]
                    self.staticr[iidd,j]=self.xvector[j]

        while nout < npop:
            if nout <  ntotal:
                trial,self.xvector,idledir,ncal0=self.get_solution() ; nout=nout+1
                self.static[nout-1,:]=self.xvector[:] ; self.ostatic[nout-1]=trial
                self.dynamic[nout-1,:]=self.xvector[:] ; self.odynamic[nout-1]=trial
                print(nout,trial,' Loss',flush=True)
    

        if os.path.isfile(self.fname_dump) :
            print(self.fname_dump,' is present, merge',flush=True)
            npop1,ndim1,ostatic1,odynamic1,static1,dynamic1=self.load_solutions(self.fname_dump)
            if ndim == ndim1:
                print('ostatic old')
                print(ostatic1)
                print('ostatic try')
                print(self.ostatic)
                print('odynamic old')
                print(odynamic1)
                print('odynamic try')
                print(self.odynamic)
                i0=npop+npop1
                static2=np.zeros((i0,ndim)) ; ostatic2=np.zeros(i0)
                dynamic2=np.zeros((i0,ndim)) ; odynamic2=np.zeros(i0)
                static2[0:npop,:]=self.static[0:npop,:] ; static2[npop:i0,:]=static1[0:npop1,:]
                dynamic2[0:npop,:]=self.dynamic[0:npop,:] ; dynamic2[npop:i0,:]=dynamic1[0:npop1,:]
                ostatic2[0:npop]=self.ostatic[0:npop] ; ostatic2[npop:i0]=ostatic1[0:npop1]
                odynamic2[0:npop]=self.odynamic[0:npop] ; odynamic2[npop:i0]=odynamic1[0:npop1]
                ind=ostatic2.argsort() ; ostatic2=ostatic2[ind]  ; static2=static2[ind,:]
                ind=odynamic2.argsort() ; odynamic2=odynamic2[ind]  ; dynamic2=dynamic2[ind,:]
                self.ostatic[0:npop]=ostatic2[0:npop] ; self.odynamic[0:npop]=odynamic2[0:npop]
                self.static[0:npop,:]=static2[0:npop,:] ; self.dynamic[0:npop,:]=dynamic2[0:npop,:]
                del ostatic1,static1,ostatic2,static2
                del odynamic1,dynamic1,odynamic2,dynamic2
                del npop1,ndim1
                gc.collect()
        else :
                ind=self.ostatic.argsort() ; self.ostatic=self.ostatic[ind] ; self.static=self.static[ind,:]
                self.dynamic[:,:]=self.static[:,:] ; self.odynamic[:]=self.ostatic[:]

        if nout == npop:
            print('nout,npop ',nout,npop,flush=True)
            ind=self.ostatic.argsort() ; self.ostatic=self.ostatic[ind] ; self.static=self.static[ind,:]
            self.staticr[:,:]=self.static[:,:]
            best=self.odynamic[0] ; davg=0. ; i0=0
            for i in range(npop):
                for j in range(npop):
                    if i > j:
                        davg=davg+self.csadistance(self.static[i,:],self.static[j,:]) ; i0=i0+1
            davg=davg/i0 ; dcut=davg/2.
            print('dcut ',dcut,flush=True)
            self.save_solutions()
            print(self.ostatic,flush=True)
            for iidd in range(npop):
                k0,k1=self.gen_parents() ; self.xvector[:]=self.dynamic[k0,:]
                if np.random.random() < 0.1:
                    self.xvector[:]=self.static[k0,:]
                if np.random.random() < 0.5:
                    self.xvector=self.mutation(self.xvector)
                else:
                    self.yvector[:]=self.dynamic[k1,:]
                    if np.random.random() < 0.1:
                        self.yvector[:]=self.static[k1,:]
                    self.xvector=self.crossover(self.xvector,self.yvector)
                if nin < ntotal :
                    self.write_trial_solution(ndim,self.xvector,iidd,nin) ; nin=nin+1

        self.nin = nin
        self.nout = nout
        self.dcut = dcut
        self.best = best
        self.davg = davg

    def csa_main_step(self):
        ntotal = self.ntotal
        ndim   = self.ndim
        npop   = self.npop
        nin    = self.nin
        nout   = self.nout
        dcut   = self.dcut
        best   = self.best
        davg   = self.davg

        while True:
            if nout <  ntotal:
                trial,self.xvector,idledir,ncal0=self.get_solution() ; nout=nout+1
                print(nout,trial,' Loss',flush=True)
                lupdate=False ; dsmall=9e99 ; i0=npop-1
                for i in range(npop):
                    tmp=self.csadistance(self.xvector,self.dynamic[i,:])
                    if dsmall > tmp:
                        i0=i ; dsmall=tmp
                if dcut > dsmall :
                    if self.odynamic[i0] > trial:
                        tmq=self.odynamic[i0] ; self.odynamic[i0]=trial  ;  self.dynamic[i0,:]=self.xvector[:]
                        lupdate=True
                        print(tmq,trial,tmq-trial,' old type',nin,nout,flush=True)
                else:
                    i0=npop-1
                    if self.odynamic[i0] > trial:
                        tmq=self.odynamic[i0] ; self.odynamic[i0]=trial ;  self.dynamic[i0,:]=self.xvector[:]
                        lupdate=True
                        print(tmq,trial,tmq-trial,' new type',nin,nout,flush=True)
                ind = self.odynamic.argsort() ; self.odynamic=self.odynamic[ind] ; self.dynamic=self.dynamic[ind,:]
                if best > self.odynamic[0] :
                    tmq=best ; best=self.odynamic[0] ; self.xvector[:]=self.dynamic[0,:]
                    print('best ',tmq,best,tmq-best,flush=True)
                if lupdate :
                    self.save_solutions()
                if nout > npop*2 and np.mod(nout,npop) == 0:
                    self.staticr[0,:]=self.dynamic[0,:]
                    self.staticr=self.coulomb(self.staticr)
                    dcut=dcut*0.99
                if dcut < davg/5. :
                    dcut=davg/5.
                if nout == ntotal:
                    return
                if nin < ntotal :
                    k0,k1=self.gen_parents() ; self.xvector[:]=self.dynamic[k0,:]
                    if np.random.random() < 0.1:
                        self.xvector[:]=self.static[k0,:]
                        if np.random.random() < 0.1:
                            self.xvector[:]=self.staticr[k0,:]
                    if np.random.random() < 0.5:
                        self.xvector=self.mutation(self.xvector)
                    else:
                        self.yvector[:]=self.dynamic[k1,:]
                        if np.random.random() < 0.1:
                            self.yvector[:]=self.static[k1,:]
                            if np.random.random() < 0.1:
                                self.yvector[:]=self.staticr[k1,:]
                        self.xvector=self.crossover(self.xvector,self.yvector)
                    self.write_trial_solution(ndim,self.xvector,idledir,nin) ; nin=nin+1
            else :
                return

    def csa_run(self):
        self.csa_initial_step()
        self.csa_main_step()


    # for CSA_SOLDIER
class csa_soldier_tools:

    def __init__(self,param_type, param_const):
        # iparam_type : array of parameter type: mytb.ppram.iparam_type
        # check only valid parameters: onsite and hopping. 
        #                              default value of overlap and scale parameters will be used.
        # Be aware that csa_soldier_tools.nparam is not always equal to nparam from pytbfit.ppram.nparam (depend on the constraint such as 'fix')
        #self.valid_param_idx = np.ma.masked_less_equal(param_type,2).mask
        param_fix = param_const[3,:]

        # parameters with ( param_type <=2  .AND. parameters not fixed )
        self.valid_param_idx = np.logical_and(np.ma.masked_less_equal(param_type,2).mask, [np.invert(np.bool(x)) for x in np.int64(param_fix)])
        self.nparam = np.count_nonzero(self.valid_param_idx)

    def read_csa_input(self,fname='input.txt'):
        if os.path.isfile(fname):
            ncal = 0
            afile = open('input.txt','r')
            jline = 0; i0 = 0
            for line in afile:
                if jline == 0 :
                    ndim = int(line.split()[0])
                    xvector = np.zeros(ndim)
                else:
                    if i0 < ndim:
                        xvector[i0]=float(line.split()[0])
                    i0 = i0 + 1
                if jline == ndim+1:
                    if len(line.split()) > 0:
                        ncal = int(line.split()[0])
                jline = jline + 1
            afile.close()
        else:
            print('input.txt is not preset',flush=True)
            ndim = self.nparam
            xvector = np.zeros(ndim)
            for i in range(ndim):
                xvector[i] = (np.random.random()-0.5)*2.*5.
            ncal = 0

        self.ndim = ndim
        self.xvector = xvector
        self.ncal = ncal
       #return ndim, xvector, ncal
        return xvector

    def write_csa_output(self, obj = None):
        file_name = 'output.txt'
        ndim = self.ndim
        xvector = self.xvector
        ncal = self.ncal

        if os.path.isfile(file_name):
            os.remove(file_name)
        lines_to_append=[]
        lines_to_append.append(str(ndim))
        for j in range(ndim):
            lines_to_append.append(str(xvector[j]))
        lines_to_append.append(str(obj))
        lines_to_append.append(str(ncal))
        self.append_multiple_lines(file_name, lines_to_append)


    def check_obj_rank(self,obj,dump_path):
        if os.path.isfile(dump_path):
            with open(dump_path, "r") as file:
                first_line = file.readline()
                for last_line in file:
                    pass
            npop = int(first_line.split()[0])
            obj_last = float(last_line.split()[npop-1])
            if obj<obj_last:
                return True
            else:
                return False
        else:
            return True

    def append_multiple_lines(self, file_name, lines_to_append):
        with open(file_name, "a+") as file_object:
            appendEOL = False
            file_object.seek(0)
            data = file_object.read(100)
            if len(data) > 0:
                appendEOL = True
            for line in lines_to_append:
                if appendEOL == True:
                    file_object.write("\n")
                else:
                    appendEOL = True
                file_object.write(line)

    def update_csa_result(self, obj, dump_path):
        if self.check_obj_rank(obj,dump_path) is True:
            j=np.random.randint(90000)
            file_name='solution_'+str(self.ncal).zfill(8)+'_'+str(obj)+'.txt'
            if os.path.isfile(file_name):
                os.remove(file_name)
            lines_to_append=[]
            for j in range(self.ndim):
                lines_to_append.append(str(self.xvector[j]))
            lines_to_append.append(str(obj))
            lines_to_append.append(str(self.ncal))
            self.append_multiple_lines(file_name, lines_to_append)
            flag_save_result = True
            return flag_save_result
        else:
            flag_save_result = False
            return flag_save_result

    def update_status(self,sleep=1.0):
        time.sleep(sleep)
        os.system('touch STOP')
        os.system('echo "DONE" >> STATUS')

