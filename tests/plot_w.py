from __future__ import print_function
from orphics import maps,io,cosmology,mpi,stats
from enlib import enmap
import numpy as np
import os,sys


ells = np.arange(2,3000)
low_acc = True
    
def dcl(w):

    
    params = cosmology.defaultCosmology.copy()
    params['w0'] = w
    cc = cosmology.Cosmology(params,low_acc=low_acc,pickling=False,skipPower=True,skip_growth=True)
    ucltt = cc.theory.uCl('TT',ells)
    lcltt = cc.theory.lCl('TT',ells)

    uclte = cc.theory.uCl('TE',ells)
    lclte = cc.theory.lCl('TE',ells)

    uclee = cc.theory.uCl('EE',ells)
    lclee = cc.theory.lCl('EE',ells)
    
    return ucltt,lcltt,uclte,lclte,uclee,lclee



ws = np.linspace(-1.2,-0.8,6)

comm = mpi.MPI.COMM_WORLD
rank = comm.Get_rank()
numcores = comm.Get_size()
Njobs = len(ws)
num_each,each_tasks = mpi.mpi_distribute(Njobs,numcores)
if rank==0: print ("At most ", max(num_each) , " tasks...")
my_tasks = each_tasks[rank]

s = stats.Stats(comm)
if True:
    for task in my_tasks:
        ucltt,lcltt,uclte,lclte,uclee,lclee = dcl(ws[task])
        vec = np.append(np.array((ws[task],)),ucltt)
        vec2 = np.append(np.array((ws[task],)),lcltt)
        vecte = np.append(np.array((ws[task],)),uclte)
        vecte2 = np.append(np.array((ws[task],)),lclte)
        vecee = np.append(np.array((ws[task],)),uclee)
        vecee2 = np.append(np.array((ws[task],)),lclee)
        s.add_to_stats("vec",vec)
        s.add_to_stats("vec2",vec2)
        s.add_to_stats("vecte",vecte)
        s.add_to_stats("vecte2",vecte2)
        s.add_to_stats("vecee",vecee)
        s.add_to_stats("vecee2",vecee2)
        if rank==0: print ("Rank 0 done with task ", task+1, " / " , len(my_tasks))

    s.get_stats()

if rank==0:

    
    vecs = s.vectors['vec']
    vecs2 = s.vectors['vec2']
    vecste = s.vectors['vecte']
    vecste2 = s.vectors['vecte2']
    vecsee = s.vectors['vecee']
    vecsee2 = s.vectors['vecee2']
    #vecs = np.loadtxt("dcls.txt")

    print(vecs)
    print(vecs.shape)

    pl = io.Plotter()
    for row,row2 in zip(vecs,vecs2):
        step = row[0]
        dcl = row[1:]
        dcl2 = row2[1:]
        pl.add(ells,dcl*ells**2.,label="%.3f"%step)
        #pl.add(ells,dcl2*ells**2.,label="%.3f"%step,ls="--")
    pl.legend(loc='upper right')
    pl.done(io.dout_dir+"wtt.png")


    pl = io.Plotter()
    for row,row2 in zip(vecsee,vecsee2):
        step = row[0]
        dcl = row[1:]
        dcl2 = row2[1:]
        pl.add(ells,dcl*ells**2.,label="%.3f"%step)
        #pl.add(ells,dcl2*ells**2.,label="%.3f"%step,ls="--")
    pl.legend(loc='upper right')
    pl.done(io.dout_dir+"wee.png")
    

    pl = io.Plotter()
    for row,row2 in zip(vecste,vecste2):
        step = row[0]
        dcl = row[1:]
        dcl2 = row2[1:]
        pl.add(ells,dcl*ells**2.,label="%.3f"%step)
        #pl.add(ells,dcl2*ells**2.,label="%.3f"%step,ls="--")
    pl.legend(loc='upper right')
    pl.done(io.dout_dir+"wte.png")
    
