#!~/usr/bin/python
import formatxfer as fx
import numpy as np

def surface_roughness(a,N_type,N_space):
    aa=fx.info(a,'cfg',2)
    cell=aa.cell
    data=aa.data     
    n_type=aa.atom_type_num
    tot_num=int(sum(n_type))
    #data=np.delete(data,np.arange(n_type[0]))
    surface_atoms=[]
    for ii in range(N_space):
        xmin=ii/float(N_space)
        xmax=(ii+1)/float(N_space)
        for jj in range(N_space):
            ymin=jj/float(N_space)
            ymax=(jj+1)/float(N_space)
            l2 = [data[i] for i in range(tot_num) if data[i][1] > xmin and data[i][1] < xmax and data[i][2] > ymin and data[i][2] < ymax and data[i][0]==1.07868200e+02 ]
            l2 = sorted(l2,key=lambda x:x[3],reverse=True)
            if len(l2)==0: pass
            else:
                surface_atoms.append(l2[0])
    aa.data=np.asarray(surface_atoms)
    aa.tot_num=int(len(aa.data))
    aa.atom_type_num=[n_type[1]]
    aa.get_cfg_file_one('surface.cfg')
    return

surface_roughness('example/dump.000001677.cfg',2,100)
