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

# last update: 06.01.2022 HJ Kim

# IMPORT NOTE:
# if you want to run tbfitpy_mod_mpi with MPI implementation, 
# uncomment "from mpi4py import MPI " and, "from tbfitpy_mod_mpi ..."
# comment out "from tbfitpy_mod import pyfit".
# If you want to run serial tbfitpy_mod, 
# comment out "from mpi4py ..." and "from tbfitpy_mod_mpi ..."

class pytbfit:
    def __init__(self, mpicomm =None, filenm = 'INCAR-TB'):
        if mpicomm is not None :
            self.comm  = mpicomm
            self.fcomm = self.comm.py2f()
        else:
            self.comm  = None
            self.fcomm = 0

        self.filenm = filenm+' '*(132-len(filenm))
        self.pinpt = pyfit.init_incar_py( self.filenm , nsystem=1 )
        self.ppram = pyfit.init_params_py()    
        self.pkpts = pyfit.init_kpoints_py()
        self.pwght = pyfit.init_weight_py()
        self.pgeom = pyfit.init_poscar_py()
        self.hopping = pyfit.init_hopping_py()
        self.edft  = pyfit.init_energy_py()
        self.etba  = pyfit.init_energy_py()
        self.orbfit= False

    def init(self, verbose=False, orbfit=False, myid=0):

        if verbose is True:
            self.pinpt.iverbose = 1
        else :
            self.pinpt.iverbose = 2

        if orbfit is True:
            self.pinpt.flag_fit_orbital_parse = True
        else:
            self.pinpt.flag_fit_orbital_parse = False
        self.orbfit = orbfit

        # initialize
        pyfit.init(self.fcomm, self.pinpt, self.ppram, self.pkpts, self.pwght, 
                               self.pgeom, self.hopping, self.edft, self.etba)
        self.cost_history = []
        self.niter = 0

        self.cost = 0.0
        self.cost_orb = 0.0

        self.myid = myid

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

        if verbose is True:
            self.pinpt.iverbose = 1
        else :
            self.pinpt.iverbose = 2

        t0 = time.time()

        if tol is not None: self.pinpt.ptol = tol
        if tol is not None: self.pinpt.ftol = tol

        if miter is not None: self.pinpt.miter = miter
        self.n_particles = n_particles if n_particles  is not None else 30
        if pso_options is None: pso_options = {'c1': 1.0, 'c2': 2.0, 'w':2.0} # set these values as default

        if method == 'lmdif' :
            self.cost_history = np.full( (self.pinpt.miter), 0.0)
            self.pinpt.flag_tbfit = True
            pyfit.fit(self.fcomm, self.pinpt, self.ppram, self.pkpts, self.pwght, 
                                  self.pgeom, self.hopping, self.edft, self.etba)

            self.cost_history = self.ppram.cost_history
            self.cost     = np.exp( -(np.sum(abs(self.etba.de)) / sigma)**2)*100
            if self.orbfit is True:
                self.cost_orb = np.exp( -(np.sum(abs(self.etba.dorb)) / sigma_orb)**2)*100
    
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
            
            pyfit.pso(self.fcomm, self.pinpt, self.ppram, self.pkpts, self.pwght,
                                  self.pgeom, self.hopping, self.edft, self.etba, iseed, pso_miter)

            self.cost_history   = self.ppram.pso_cost_history
            self.cost_history_particle = self.ppram.pso_cost_history_i

        if self.myid == 0 and verbose is True:
            print('Time elapsed for ',method,' method :', time.time() - t0, ' sec')
            sys.stdout.flush()

    def get_eig(self, verbose=None):
        if verbose is True:
            self.pinpt.iverbose = 1
        else :
            self.pinpt.iverbose = 2

        pyfit.eig(self.fcomm, self.pinpt, self.ppram, self.pkpts, 
                              self.pgeom, self.hopping, self.etba)

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


    def plot_fit(self, figsize=(5,6), ef=0.0, yin=-20., yen=10, ystep=5., plot_model=None,
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
        if plot_model is not None:
            mylabel = 'MODEL' if i == 0 else None
            for i in range(self.etba.e.shape[0]):
                plt.plot(self.pkpts.kdist, self.model.e[i,:],'--', linewidth=.6, color='green', label=mylabel)

        plt.legend()
        #plt.savefig(fout,bbox_inches='tight',transparent=False,pad_inches=0)
        plt.savefig(fout,bbox_inches='tight',pad_inches=0)

        plt.clf()

    def plot_pso_cost_history(self, figsize=(5,6), yin=-20,yen=20, ystep=5.,
                                    title=' ', ylabel='Best cost function', xlabel='PSO iteration', fout='COST_HISTORY.pdf', bestn=5):
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
    
    def save(self, title=None, plot_fit=None, plot_model=None, target=None, band=None, param=None, weight=None, 
                   cost_history=None, cost_history_particle=None , pso_param_history = None):
        if self.myid != 0:
            return

        if plot_fit is True:
            if title is not None:
                fout = 'BAND.'+title.strip()+'.pdf'
            else:
                fout = 'BAND.pdf'
            self.plot_fit(title=title.strip(), fout=fout, plot_model=None)

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
