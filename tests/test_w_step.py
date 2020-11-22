from __future__ import print_function
from orphics import maps,io,cosmology,mpi,stats
from enlib import enmap
import numpy as np
import os,sys


ells = np.arange(2,6000)
low_acc = False
    
def dcl(step=0.1):

    
    params = cosmology.defaultCosmology.copy()
    params['w0'] = -1 + step/2.
    cc = cosmology.Cosmology(params,low_acc=low_acc,pickling=False,skipPower=True,skip_growth=True,lmax=7000)
    ucltt_up = cc.theory.uCl('TT',ells)
    lcltt_up = cc.theory.lCl('TT',ells)
    uclte_up = cc.theory.uCl('TE',ells)
    lclte_up = cc.theory.lCl('TE',ells)
    uclee_up = cc.theory.uCl('EE',ells)
    lclee_up = cc.theory.lCl('EE',ells)
    
    params = cosmology.defaultCosmology.copy()
    params['w0'] = -1 - step/2.
    cc = cosmology.Cosmology(params,low_acc=low_acc,pickling=False,skipPower=True,skip_growth=True,lmax=7000)
    ucltt_dn = cc.theory.uCl('TT',ells)
    lcltt_dn = cc.theory.lCl('TT',ells)
    uclte_dn = cc.theory.uCl('TE',ells)
    lclte_dn = cc.theory.lCl('TE',ells)
    uclee_dn = cc.theory.uCl('EE',ells)
    lclee_dn = cc.theory.lCl('EE',ells)

    dcltt = (ucltt_up-ucltt_dn)/step
    dlcltt = (lcltt_up-lcltt_dn)/step
    dclte = (uclte_up-uclte_dn)/step
    dclee = (uclee_up-uclee_dn)/step

    return dcltt,dclte,dclee#,dlcltt



#wsteps = [0.005,0.01,0.05,0.1,0.3]
wsteps = [0.000001,0.000005,0.00001,0.00005,0.0001]

comm = mpi.MPI.COMM_WORLD
rank = comm.Get_rank()
numcores = comm.Get_size()
Njobs = len(wsteps)
num_each,each_tasks = mpi.mpi_distribute(Njobs,numcores)
if rank==0: print ("At most ", max(num_each) , " tasks...")
my_tasks = each_tasks[rank]

s = stats.Stats(comm)
if True:
    for task in my_tasks:
        #dcltt,dlcltt = dcl(wsteps[task])
        dcltt,dclte,dclee = dcl(wsteps[task])
        vec = np.append(np.array((wsteps[task],)),dcltt)
        vec2 = np.append(np.array((wsteps[task],)),dclte)
        vec3 = np.append(np.array((wsteps[task],)),dclee)
        s.add_to_stats("vec",vec)
        s.add_to_stats("vec2",vec2)
        s.add_to_stats("vec3",vec3)
        if rank==0: print ("Rank 0 done with task ", task+1, " / " , len(my_tasks))

    s.get_stats()

if rank==0:

    
    #np.savetxt("dcls.txt",s.vectors['vec'])
    vecs = s.vectors['vec']
    vecs2 = s.vectors['vec2']
    vecs3 = s.vectors['vec3']
    #vecs = np.loadtxt("dcls.txt")

    print(vecs)
    print(vecs.shape)

    pl = io.Plotter()
    for row in vecs:
        step = row[0]
        dcl = row[1:]
        pl.add(ells,dcl*ells**2.,label="%.3f"%step)
    pl.legend(loc='upper right')
    pl.done(io.dout_dir+"wsteptt.png")


    pl = io.Plotter()
    for row in vecs2:
        step = row[0]
        dcl = row[1:]
        pl.add(ells,dcl*ells**2.,label="%.3f"%step)
    pl.legend(loc='upper right')
    pl.done(io.dout_dir+"wstepte.png")


    pl = io.Plotter()
    for row in vecs3:
        step = row[0]
        dcl = row[1:]
        pl.add(ells,dcl*ells**2.,label="%.3f"%step)
    pl.legend(loc='upper right')
    pl.done(io.dout_dir+"wstepee.png")
    
