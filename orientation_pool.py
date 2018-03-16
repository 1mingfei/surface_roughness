#!~/usr/bin/python
import formatxfer as fx
import numpy as np
import multiprocessing

def all_NNs(a,N_type,cutoff):
    aa=fx.info(a,'cfg',N_type)
    cell=aa.cell
    data1=aa.data     
    data1=[x for x in data1 if (x[0]==1.07868200e+02)]
    tot_num=int(len(data1))
    tmp=np.arange(len(data1))
    tmp=tmp.reshape(len(data1),1)
    data = np.concatenate((data1,tmp),axis=1)
    print data

    pool = multiprocessing.Pool(processes=8)


    return



all_NNs('example/dump.000001677.cfg',2,4.0)

