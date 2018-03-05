#!~/usr/bin/python
import formatxfer as fx
import numpy as np
import multiprocessing
from functools import partial

def get_surf_atoms(limits):
    aa=fx.info('example/dump.000001677.cfg','cfg',2)
    #aa=fx.info(a,'cfg',tot_num)
    cell=aa.cell
    data=aa.data     
    n_type=aa.atom_type_num
    tot_num=int(len(data))
    surface_atoms=[]
    l2 = [data[i] for i in range(tot_num) if data[i][1] > limits[0] and  data[i][1] < limits[1] and  data[i][2] > limits[2] and  data[i][2] < limits[3] and  data[i][0]==1.07868200e+02 ]
    l2 = sorted(l2,key=lambda x:x[3],reverse=True)
    if len(l2)==0: pass
    else:
        return(l2[0])

def get_surf_atoms_process(data,tot_num,limits):
    surface_atoms=[]
    l2 = [data[i] for i in range(tot_num) if data[i][1] > limits[0] and data[i][1] < limits[1] and data[i][2] > limits[2] and data[i][2] < limits[3] and data[i][0]==1.07868200e+02 ]
    l2 = sorted(l2,key=lambda x:x[3],reverse=True)
    if len(l2)==0: pass
    else:
        surface_atoms.append(l2[0])
    return(np.asarray(surface_atoms))

def surface_atoms_filter(a,N_type,N_space):
    aa=fx.info(a,'cfg',2)
    cell=aa.cell
    data=aa.data     
    n_type=aa.atom_type_num
    tot_num=int(len(data))
    xx, yy = np.meshgrid(np.linspace(0.0, 1.0, N_space-1),np.linspace(0.0, 1.0, N_space-1))
    M = np.array([[ x1, x1+1.0/float(N_space), x2, x2+1.0/float(N_space)] for x1, x2 in zip(np.ravel(xx), np.ravel(yy))])
    pool = multiprocessing.Pool(processes=8)
    surface_atoms = pool.map(get_surf_atoms,M)
    surface_atoms=[x for x in surface_atoms if x is not None]
    surface_atoms=np.asarray(surface_atoms)

    with open("surf.dat",'w') as fout:
        for i in range(3):
            fout.write("%12.8f   %12.8f   %12.8f\n"%(cell[i][0],cell[i][1],cell[i][2]))
        for i in range(len(surface_atoms)):
            fout.write("%12.8f   %12.8f   %12.8f   %12.8f\n"%(surface_atoms[i][0],surface_atoms[i][1],surface_atoms[i][2],surface_atoms[i][3]))

    #print np.asarray(surface_atoms) 
    #print(len(surface_atoms))
    aa.data=np.asarray(surface_atoms)
    aa.tot_num=int(len(aa.data))
    aa.atom_type_num=[n_type[1]]
    aa.get_cfg_file_one('surface.cfg')
    h_bar=np.mean(surface_atoms[:,3])
    return

def surface_statistics(a):
    with open(a,'r') as fin:
        mylines=fin.readlines()
        cell=np.zeros((3,3))
        data=np.zeros((len(mylines)-3,4))
        for i in range(3):
            cell[i][0]=float(mylines[i].split()[0])
            cell[i][1]=float(mylines[i].split()[1])
            cell[i][2]=float(mylines[i].split()[2])
        for i in range(len(mylines)-3):
            for j in range(4):
                data[i][j]=float(mylines[3+i].split()[j])
    print data 
    return



surface_atoms_filter('example/dump.000001677.cfg',2,20)
surface_statistics('surf.dat')
