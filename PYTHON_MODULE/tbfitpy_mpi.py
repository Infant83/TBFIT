import numpy as np
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters
import pyswarms as ps
from pyswarms.utils.plotters import plot_cost_history
from random import random
import warnings
from tqdm import tqdm 
import time
import sys
from mpi4py import MPI
from tbfitpy_mod_mpi import pyfit
#from tbfitpy_mod import pyfit
warnings.filterwarnings("ignore")

# IMPORT NOTE:
# if you want to run tbfitpy_mod_mpi with MPI implementation, 
# uncomment "from mpi4py import MPI " and, "from tbfitpy_mod_mpi ..."
# comment out "from tbfitpy_mod import pyfit".
# If you want to run serial tbfitpy_mod, 
# comment out "from mpi4py ..." and "from tbfitpy_mod_mpi ..."

class pytbfit:
    def __init__(self, mpicomm =None, filenm = 'INCAR-TB'):
        self.filenm = filenm+' '*(132-len(filenm))

        self.pinpt = pyfit.init_incar_py( self.filenm , nsystem=1 )
        self.ppram = pyfit.init_params_py()    
        self.pkpts = pyfit.init_kpoints_py()
        self.pwght = pyfit.init_weight_py()
        self.pgeom = pyfit.init_poscar_py()
        self.hopping = pyfit.init_hopping_py()
        self.edft  = pyfit.init_energy_py()
        self.etba  = pyfit.init_energy_py()
        if mpicomm is not None :
           #self.fcomm = mpicomm
            self.comm  = mpicomm
            self.fcomm = self.comm.py2f()
        else:
            self.comm  = None
            self.fcomm = 0

    def init(self, verbose=False, myid=0):

        if verbose is True:
            self.pinpt.iverbose = 1
        else :
            self.pinpt.iverbose = 2

        # initialize
        pyfit.init(self.fcomm, self.pinpt, self.ppram, self.pkpts, self.pwght, 
                               self.pgeom, self.hopping, self.edft, self.etba)
        self.cost_history = []
        self.niter = 0

        self.myid = myid

    def cost(self):
        self.get_eig() # get eig
        cost = ((self.pwght.wt*(self.etba.e - self.edft.e))**2).ravel()
        return cost
    
    def residual(self, params):
        for iparam in range(self.ppram.nparam_free):
            pname = list(params.valuesdict())[iparam]
            param = params[pname].value
            self.ppram.param[self.ppram.iparam_free[iparam]-1] = param
        return  self.cost()
    
    def iter_cb(self, params, iiter, resid):
        self.cost_history.append ( sum(abs(resid)) )
        self.niter += 1
        return None

    def cost_pso(self, params):
        cost=[]
        for i in range(self.n_particles):
            for iparam in range(self.ppram.nparam_free):
                self.ppram.param[self.ppram.iparam_free[iparam]-1] = params[i, iparam]
            cost.append( np.sqrt(sum( (self.cost() ) )) )
        return cost

    def cost_psoleastsq(self, params):
        cost=[]
        for i in range(self.n_particles):
            for iparam in range(self.ppram.nparam_free):
                self.ppram.param[self.ppram.iparam_free[iparam]-1] = params[i, iparam]
            
            myparam = self.get_free_param()
            out_temp = minimize(self.residual, myparam, method='leastsq')
            for iparam in range(self.ppram.nparam_free):
                pname = list(out_temp.params.valuesdict())[iparam]
                self.ppram.param[self.ppram.iparam_free[iparam]-1] = out_temp.params[pname].value

           #cost.append( sum(abs(self.cost())) )
            cost.append( np.sqrt(sum( (self.cost() ) )) )

        return cost

    def cost_pso_single(self, params):
        for iparam in range(self.ppram.nparam_free):
            self.ppram.param[self.ppram.iparam_free[iparam]-1] = params[iparam]
       #return sum(abs(self.cost()))
        return np.sqrt(sum( (self.cost() ) ))

    def cost_pso_leastsq(self, params):
        cost=[]
        for i in range(self.n_particles):
            for iparam in range(self.ppram.nparam_free):
                self.ppram.param[self.ppram.iparam_free[iparam]-1] = params[i, iparam]
            myparam = self.get_free_param()
            self.fit_out = minimize(self.residual, myparam, method='leastsq')
            for iparam in range(self.ppram.nparam_free):
                pname = list(self.fit_out.params.valuesdict())[iparam]
                self.ppram.param[self.ppram.iparam_free[iparam]-1] = self.fit_out.params[pname].value
            cost.append( self.fit_out.chisqr )
        return cost

    def update_pso_pos(self, pos, vel, bounds):
        next_pos = pos + vel
       #for i in range(self.n_particles):
       #    for iparam in range(self.ppram.nparam_free):
       #        if next_pos[i,iparam] < bounds[0][iparam] :
       #            next_pos[i,iparam] = bounds[0][iparam] + bounds[0][iparam] - next_pos[i,iparam]

       #        if next_pos[i,iparam] > bounds[1][iparam]:
       #            next_pos[i,iparam] = bounds[1][iparam] + bounds[1][iparam] - next_pos[i,iparam]

        return next_pos

    def update_pso_vel(self, pos, vel, p_best, g_best, c1=0.5, c2=1.0, w=0.75):
        for i in range(self.n_particles):
            for j in range(self.ppram.nparam_free):
                r1 = np.random.uniform() ; r2 = np.random.uniform()
                if self.comm is not None:
                    r1 = self.comm.bcast(r1, root=0)
                    r2 = self.comm.bcast(r2, root=0)
                vel[i,j] = w * vel[i,j] + c1 * r1 * (p_best[i,j] - pos[i,j]) + c2 * r2 * (g_best[j] - pos[i,j])
       #for i in range(self.ppram.nparam_free):
       #    r1 = np.random.uniform() ; r2 = np.random.uniform()
       #    if self.comm is not None:
       #        r1 = self.comm.bcast(r1, root=0)
       #        r2 = self.comm.bcast(r2, root=0)
       #    vel[i] = w * vel[i] + c1 * r1 * (p_best[i] - pos[i] ) + c2 * r2 * ( g_best[i] - pos[i] )
        return vel

    def fit(self, verbose=False, miter=None, method='leastsq', pso_options=None, n_particles=None, iseed=None):
        '''
        Fit the parameters.
        Following methods are supported:
           *lmdif: Levenberg-Marquardt method with MINPACK subroutine modified by H.-J. Kim (TBFIT)
              (see details in: https://github.com/Infant83/TBFIT)

         - Minimization function via LMFIT module 
              (see details in: https://lmfit.github.io/lmfit-py/intro.html)
           *leastsq: Levenberg-Marquardt (default)
           *least_squares: Least-Squares minimization, using Trust Region Reflective method
           *bfgs: BFGS

         - Particle Swarm Optimization (PSO) method
           *pso: A Global-best Particle Swarm Optimization algorithm
           *pso.leastsq: PSO with leastsq method (better performance and higher computational load)
           *mypso: same as pso, but much faster.
           *mypso.lmdif: same as pso.leastsq but much faster

         - Particle Swarm Optimization (PSO) scheme via Pyswars module
              (see details in: https://pyswarms.readthedocs.io/en/latest/index.html)
           *gbest.pso: A Global-best Particle Swarm Optimization algorithm
           *lbest.pso: A Local-best Particle Swarm Optimization algorithm
            Note: To use this method, one might pass 'option' variable as well, which defines
                  cognitive, social, and inertia parameter used in PSO algorithm.
                  In addition, one also need 'n_particles' defining size of swarm of parameters.
                  The details can be found in https://pyswarms.readthedocs.io/en/latest/index.html
    
           Note on all PSO based altorhtim: one should provide n_particle and pso_options when call
                n_particles: number of particles (the swarm size), integer
                pso_options: velocity and position update policy, dictionary, 
                             ex: pso_options = {'c1': 1.0, 'c2': 1.0, 'w': 2.0'}
                             for details, see: J. Kennedy and R. Eberhart, Particle Swarm Optimization,
                                              􏰀IEEE, Piscataway, NJ, 1995􏰁, p. 1942.
        NONE: I suggest to use lmdif method for Levenberg-Marquardt method for local fitting 
              than leastsq and least_squares. lmdif is much faster that the others.
              For global parameter fitting, mypso is suggested. 
              To utilize lmdif method with mypso, mypso.lmdif can be used. 
              gbest.pso, lbest.pso does not support MPI, and pso is slower than mypso.
        '''

        if verbose is True:
            self.pinpt.iverbose = 1
        else :
            self.pinpt.iverbose = 2

        t0 = time.time()

        if miter is not None: self.pinpt.miter = miter
        self.n_particles = n_particles if n_particles  is not None else 30
        if pso_options is None: pso_options = {'c1': 1.0, 'c2': 2.0, 'w':2.0} # set these values as default

        if method == 'lmdif' :
            self.cost_history = np.full( (self.pinpt.miter), 0.0)
            pyfit.fit(self.fcomm, self.pinpt, self.ppram, self.pkpts, self.pwght, 
                                  self.pgeom, self.hopping, self.edft, self.etba)

            self.cost_history = self.ppram.cost_history

        elif method == 'mypso' or method == 'mypso.lmdif':
            # iseed : random seed
            if iseed is None: iseed = 123
            ilmdif = 1 if method == 'mypso.lmdif' else 0

            self.cost_history = np.full( (self.pinpt.miter), 0.0)
            self.ppram.pso_c1 = pso_options['c1']
            self.ppram.pso_c2 = pso_options['c2']
            self.ppram.pso_w  = pso_options['w' ]
            self.ppram.pso_nparticles = self.n_particles
            
            pyfit.pso(self.fcomm, self.pinpt, self.ppram, self.pkpts, self.pwght,
                                  self.pgeom, self.hopping, self.edft, self.etba, iseed, ilmdif)

            self.cost_history = self.ppram.pso_cost_history

        elif method == 'pso' or method == 'pso.leastsq':
           # initialize
            self.pinpt.iverbose = 2

            # set bound
            max_bounds = np.zeros(self.ppram.nparam_free)
            min_bounds = np.zeros(self.ppram.nparam_free)
            bounds     = (np.zeros(self.ppram.nparam_free),np.zeros(self.ppram.nparam_free))
            for iparam in range(self.ppram.nparam_free):
                max_bounds[iparam]=self.ppram.param_const[1, self.ppram.iparam_free[iparam]-1]
                min_bounds[iparam]=self.ppram.param_const[2, self.ppram.iparam_free[iparam]-1]
            bounds=(min_bounds, max_bounds)
        
            # allocate and initialize pos and vel with random noise
            pos = np.ndarray((self.n_particles, self.ppram.nparam_free))
            vel = np.ndarray((self.n_particles, self.ppram.nparam_free))
            for iparam in range(self.ppram.nparam_free):
                for i in range(self.n_particles):
                    cond = False
                    while cond is False:
                        r = (random()*2. - 1.) * 5.
                        if self.comm is not None: r = self.comm.bcast(r, root=0)

                        pos[i,iparam] = self.ppram.param[self.ppram.iparam_free[iparam]-1] + r
                        vel[i,iparam] = r
                        if pos[i,iparam] > max_bounds[iparam] or pos[i,iparam] < min_bounds[iparam]:
                            cond = False
                        else:
                            cond = True

            cost = self.cost_pso(pos)
            g_best = pos[ cost.index(min(cost)) ]
            cg_best = min(cost)
            p_best  = pos
            cp_best = cost

            # optimize with PSO algorithm
            if verbose is True:
                pbar = tqdm(range(self.pinpt.miter)) if self.myid == 0 else range(self.pinpt.miter)
            else:
                pbar = range(self.pinpt.miter)

            for _ in pbar:
                if self.myid == 0 and verbose is True:
                    pbar.set_description("PSO: BEST COST = %.8f " % cg_best)
                vel  = self.update_pso_vel(pos, vel, p_best, g_best, 
                                           pso_options['c1'], pso_options['c2'], pso_options['w'])
                pos  = self.update_pso_pos(pos, vel, bounds)
                
                if method == 'pso.leastsq':
                    cost = self.cost_psoleastsq(pos)
                elif method == 'pso':
                    cost = self.cost_pso(pos)

                if min(cost) < cg_best : 
                    cg_best = min(cost)
                    g_best = pos[cost.index(min(cost))]
                for i in range(self.n_particles):
                    if cost[i] < cp_best[i]:
                        p_best[i] = pos[i]
                        cp_best[i] = cost[i]

                self.cost_history.append(cg_best)
   
            # save best post and mininum cost 
            for iparam in range(self.ppram.nparam_free):
                self.ppram.param[self.ppram.iparam_free[iparam]-1] = g_best[iparam]
           #self.fit(method='leastsq')
            if method == 'pso.leastsq':
                self.fit(method='leastsq')
            elif method == 'pso':
                self.get_eig() # save eigenvalue with global best parameters (pos)
            self.etba.de = self.etba.e - self.edft.e

        elif method == 'gbest.pso' or method == 'lbest.pso' or method == 'gbest.pso.leastsq':
            # NOTE: This methods should be run with non-parallel tbfitpy module, tbfitpy_mod
            max_bounds=np.zeros(self.ppram.nparam_free)
            min_bounds=np.zeros(self.ppram.nparam_free)
            bounds=(np.zeros(self.ppram.nparam_free),np.zeros(self.ppram.nparam_free))
            for iparam in range(self.ppram.nparam_free):
                max_bounds[iparam]=self.ppram.param_const[1, self.ppram.iparam_free[iparam]-1]
                min_bounds[iparam]=self.ppram.param_const[2, self.ppram.iparam_free[iparam]-1]
            bounds=(min_bounds, max_bounds)
    
            init_pos=np.ndarray((self.n_particles, self.ppram.nparam_free))
            for iparam in range(self.ppram.nparam_free):
                for i in range(self.n_particles):
                    cond = True 
                    while cond is True:
                        r = (random()*2. - 1.) * 5.
                        if self.comm is not None:
                            r = self.comm.bcast(r, root=0)

                        init_pos[i,iparam] = self.ppram.param[self.ppram.iparam_free[iparam]-1] + r
                        if init_pos[i,iparam] < max_bounds[iparam] and init_pos[i,iparam] > min_bounds[iparam]:
                            cond = False

            if method == 'gbest.pso' or method == 'gbest.pso.leastsq':
                optimizer = ps.single.GlobalBestPSO(n_particles=self.n_particles, 
                                                    dimensions=self.ppram.nparam_free, 
                                                    options=pso_options, 
                                                    bounds=bounds, 
                                                    init_pos=init_pos)
            elif method == 'lbest.pso':
                optimizer = ps.single.LocalBestPSO(n_particles=self.n_particles, 
                                                    dimensions=self.ppram.nparam_free, 
                                                    options=pso_options, 
                                                    bounds=bounds, 
                                                    init_pos=init_pos)
            
            if method == 'gbest.pso' or method == 'lbest.pso' :
                cost, pos = optimizer.optimize(self.cost_pso, iters=self.pinpt.miter, verbose=True)
            elif method == 'gbest.pso.leastsq':
                cost, pos = optimizer.optimize(self.cost_pso_leastsq, iters=self.pinpt.miter,verbose=True)
                
                params = self.get_free_param()
                self.fit_out = minimize(self.residual, params, method=method)
                for iparam in range(self.ppram.nparam_free):
                    pname = list(self.fit_out.params.valuesdict())[iparam]
                    self.ppram.param[self.ppram.iparam_free[iparam]-1] = self.fit_out.params[pname].value
                
            for iparam in range(self.ppram.nparam_free):
                self.ppram.param[self.ppram.iparam_free[iparam]-1] = pos[iparam]

            self.get_eig()
            
            self.cost_history =  optimizer.cost_history

        elif method == 'leastsq' : 
            params  = self.get_free_param()

            self.fit_out = minimize(self.residual, params, method=method, iter_cb=self.iter_cb)

            for iparam in range(self.ppram.nparam_free):
                pname = list(self.fit_out.params.valuesdict())[iparam]
                self.ppram.param[self.ppram.iparam_free[iparam]-1] = self.fit_out.params[pname].value
            self.get_eig()

            self.etba.de = self.etba.e - self.edft.e

        else :
            params = self.get_free_param()
            self.fit_out = minimize(self.residual, params, method=method)
            for iparam in range(self.ppram.nparam_free):
                pname = list(self.fit_out.params.valuesdict())[iparam]
                self.ppram.param[self.ppram.iparam_free[iparam]-1] = self.fit_out.params[pname].value

            self.etba.de = self.etba.e - self.edft.e


        if self.myid == 0 and verbose is True:
            print('Time elapsed for ',method,' method :', time.time() - t0, ' sec')
            sys.stdout.flush()

    def get_eig(self):
        pyfit.eig(self.fcomm, self.pinpt, self.ppram, self.pkpts, 
                              self.pgeom, self.hopping, self.etba)

    def get_free_param(self):
        params = Parameters()
        for iparam in range(self.ppram.nparam_free):
            pname      = str(self.ppram.param_name.T[ self.ppram.iparam_free[iparam]-1 ].view('S40'),'utf-8').strip()
            param      = self.ppram.param[self.ppram.iparam_free[iparam]-1]
            params.add(name=pname, value=param)

            # NONE: With "leastsq" and "least_squares" method, setting parameter with bounds results in error.
            #       It is not clear what is the origin for that but for the safety, just leave it default (-inf,inf)
            #max_bound  = self.ppram.param_const[1, self.ppram.iparam_free[iparam]-1]
            #min_bound  = self.ppram.param_const[2, self.ppram.iparam_free[iparam]-1]
            #params.add(name=pname, value=param, min=min_bound, max=max_bound)

        return params

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


    def plot_fit(self, figsize=(5,6), ef=0.0, yin=-20., yen=10, ystep=5.,
                 title=' ', ylabel='Energy (eV)', xlabel=' ', fout='band.pdf'):
        # get k-path name and position
        self.KNAME=[] ; self.KPTS=[]
        for i in range( self.pkpts.nline + 1 ):
            #print(self.KNAME)
            self.KNAME.append( str( self.pkpts.k_name2.T[i].view('S8'), 'utf-8').strip() )
            self.KPTS.append(self.pkpts.kdist[self.pkpts.k_name_index[i]-1])

        maxwt = np.max(self.pwght.wt)
        plt.figure(figsize=figsize)
        plt.ylim(yin, yen)
        plt.xlim(0, self.pkpts.kdist[-1])
        plt.xticks(self.KPTS, self.KNAME)
        plt.yticks( np.arange(yin, yen+.000001, step=ystep) )
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.title(title.strip())

        for i in range( self.etba.e.shape[0] ) :
            mylabel = 'DFT' if i == 0 else None
            plt.plot(self.pkpts.kdist,self.edft.e[i,:],'-',linewidth=.6,color='black', label=mylabel)
            plt.scatter(self.pkpts.kdist,self.edft.e[i,:], s=self.pwght.wt[i,:]/maxwt*10.0, color='black')
        for i in range( self.etba.e.shape[0] ) :
            mylabel = 'TBA' if i == 0 else None
            plt.plot(self.pkpts.kdist,self.etba.e[i,:],'-',linewidth=.6,color='red' , label=mylabel)

        plt.legend()
        #plt.savefig(fout,bbox_inches='tight',transparent=False,pad_inches=0)
        plt.savefig(fout,bbox_inches='tight',pad_inches=0)

    def print_param(self, param_out=None):
        if param_out is not None:
            pfileoutnm = param_out+' '*(132-len(param_out))
        else:
            pfileoutnm = self.ppram.pfileoutnm

        pyfit.print_param_py(self.pinpt, self.ppram, pfileoutnm)

    def print_weight(self, weight_out=None):
        if weight_out is not None:
           wfileoutnm = weight_out+' '*(132-len(weight_out))
        else:
           wfileoutnm = 'WEIGHT.dat'+' '*(132-len('WEIGHT.dat'))

        pyfit.print_weight(self.pwght, wfileoutnm)

    def load_weight(self, weight_in=None):
        if weight_in is not None:
            wfileinnm = weight_in 
        else:
            wfileinnm = 'WEIGHT.dat'

        self.pwght.wt = np.loadtxt(wfileinnm).T

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
    
    def save(self, title=None, plot_fit=None, target=None, band=None, param=None, weight=None, cost_history=None ):
        if self.myid != 0:
            return

        if plot_fit is True:
            if title is not None:
                fout = 'BAND.'+title.strip()+'.pdf'
            else:
                fout = 'BAND.pdf'
            self.plot_fit(title=title.strip(), fout=fout)

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
            if title is not None:
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
                cost_out_pdf = 'COST_HISTORY.'+title.strip()+'.pdf'
            else:
                cost_out_dat = 'COST_HISTORY.dat'
                cost_out_pdf = 'COST_HISTORY.pdf'
            if self.ppram.niter == 0:
                niter = self.pinpt.miter
            else:
                niter = self.ppram.niter

            plot_cost_history(self.cost_history[:niter], title=title.strip())
            plt.savefig(cost_out_pdf, bbox_inches='tight', pad_inches=0)

            cost_out = open(cost_out_dat, 'w+')
            for ii in range(niter):
                cost_out.write(" %d    %.9f \n" %(ii, self.cost_history[ii]))
            cost_out.close()
