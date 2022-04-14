import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from random import random
import warnings
from tqdm import tqdm 
import time
import sys
import torch 
#from mpi4py import MPI
#from tbfitpy_mod_mpi import pyfit
from tbfitpy_mod import pyfit
warnings.filterwarnings("ignore")

# last update: 14.04.2022 HJ Kim

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

        self.cost_history = []
        self.niter = 0
        self.cost = 0.0
        self.cost_orb = 0.0
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

#       if verbose is True:
#           self.pinpt.iverbose = 1
#       else :
#           self.pinpt.iverbose = 2

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
                sum_cost = np.sum(abs(self.etba.de))
                sum_cost_orb = np.sum(abs(self.etba.dorb))
                self.cost     = 100.0 - np.exp( -( sum_cost / sigma)**2)*100.0
                if self.orbfit is True:
                    self.cost_orb = 100.0 - np.exp( -(np.sum(abs(self.etba.dorb)) / sigma_orb)**2)*100.0
            elif self.nsystem == 2:
                pyfit.fit2(self.fcomm, self.pinpt, self.ppram, self.pkpts, self.pwght, 
                                       self.pgeom, self.hopping, self.edft, self.etba,
                                       self.pkpts2, self.pwght2,self.pgeom2, self.hopping2, self.edft2, self.etba2)
                sum_cost = np.sum(abs(self.etba.de)) + np.sum(abs(self.etba2.de))
                sum_cost_orb = np.sum(abs(self.etba.dorb)) + np.sum(abs(self.etba2.dorb))
                self.cost     = 100.0 - np.exp( -( sum_cost / sigma)**2)*100.0
            elif self.nsystem == 3:
                pyfit.fit3(self.fcomm, self.pinpt, self.ppram, self.pkpts, self.pwght, 
                                       self.pgeom, self.hopping, self.edft, self.etba,
                                       self.pkpts2, self.pwght2,self.pgeom2, self.hopping2, self.edft2, self.etba2,
                                       self.pkpts3, self.pwght3,self.pgeom3, self.hopping3, self.edft3, self.etba3)
                sum_cost = np.sum(abs(self.etba.de)) + np.sum(abs(self.etba2.de)) + np.sum(abs(self.etba3.de))
                sum_cost_orb = np.sum(abs(self.etba.dorb)) + np.sum(abs(self.etba2.dorb)) + np.sum(abs(self.etba3.dorb))

            self.cost_history = self.ppram.cost_history
            self.cost     = 100.0 - np.exp( -( sum_cost / sigma)**2)*100.0
            if self.orbfit is True: self.cost_orb = 100.0 - np.exp( -( sum_cost_orb / sigma_orb)**2)*100.0

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

            self.cost_history   = self.ppram.pso_cost_history
            self.cost_history_particle = self.ppram.pso_cost_history_i

        if self.myid == 0 and verbose is True:
            print('Time elapsed for ',method,' method :', time.time() - t0, ' sec')
            sys.stdout.flush()

        self.pinpt.flag_tbfit = flag_tbfit

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
        elif sys is not None:
            if sys == 1:
                pyfit.eig(self.fcomm, self.pinpt, self.ppram, self.pkpts,self.pgeom,self.hopping,self.etba)
            elif sys == 2:
                pyfit.eig(self.fcomm, self.pinpt, self.ppram, self.pkpts2,self.pgeom2,self.hopping2,self.etba2)
            elif sys == 3:
                pyfit.eig(self.fcomm, self.pinpt, self.ppram, self.pkpts3,self.pgeom3,self.hopping3,self.etba3)

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
        return etarget, etba, edft, wt, title, kdist, KNAME, KPTS

    def get_proj_ldos(self, isystem=1):
        if isystem == 1:
            proj_ldos = self.etba.v2
        elif isystem == 2:
            proj_ldos = self.etba2.v2
        elif isystem == 3:
            proj_ldos = self.etba3.v2
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
