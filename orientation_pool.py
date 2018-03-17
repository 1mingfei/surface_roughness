#!~/usr/bin/python
import formatxfer as fx
import get_NN
import numpy as np
import multiprocessing
		
def all_NNs(a,N_type,cutoff):
    aa=fx.info(a,'cfg',N_type)
    cell=aa.cell
    data1=aa.data     
    #data1=[x for x in data1 if (x[0]==1.07868200e+02)]
    tot_num=int(len(data1))
    tmp=np.arange(len(data1))
    tmp=tmp.reshape(len(data1),1)
    data = np.concatenate((data1,tmp),axis=1)

    pool = multiprocessing.Pool(processes=8)
    query_list=[5561]
    NN_list = [pool.apply_async(get_NN.get_NN, args=(cell,data,i,cutoff)).get() for i in query_list]
    print NN_list 

    return



all_NNs('example/dump.000001677.cfg',2,4.0)
