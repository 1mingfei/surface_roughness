#!~/usr/bin/python
import formatxfer as fx
import get_NN
import numpy as np
import multiprocessing

def id_orientation(nn_list):
    data=[x for x in nn_list if (x[2]>=nn_list[0][2])]
    data=sorted(data,key=lambda x:x[2],reverse=True)
    #data=np.asarray(data[:5])
    data=np.asarray(data)
    return data


		
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
    query_list=[0]
    NN_list = [pool.apply_async(get_NN.get_NN, args=(cell,data,i,cutoff)).get() for i in query_list]
    print NN_list[0]
    print id_orientation(NN_list[0][1])

    return



#all_NNs('example/dump.000001677.cfg',2,4.0)
all_NNs('example/Ag_111.cfg',1,4.0)

