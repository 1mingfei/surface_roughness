#!~/usr/bin/python
import formatxfer as fx
import get_NN
import numpy as np
import multiprocessing

'''
100 010 001 12 NNs
[[ 0.         -0.70710678  0.70710678]
 [ 0.          0.70710678  0.70710678]
 [ 0.70710678  0.          0.70710678]
 [-0.70710678  0.          0.70710678]
 [ 0.70710678 -0.70710678  0.        ]
 [-0.70710678 -0.70710678  0.        ]
 [ 0.70710678  0.70710678  0.        ]
 [-0.70710678  0.70710678  0.        ]
 [ 0.         -0.70710678 -0.70710678]
 [ 0.          0.70710678 -0.70710678]
 [ 0.70710678  0.         -0.70710678]
 [-0.70710678  0.         -0.70710678]]
'''

def id_orientation(nn_list):
    data=[x for x in nn_list if (x[2]>=nn_list[0][2])]
    data=sorted(data,key=lambda x:x[2],reverse=True)
    #data=np.asarray(data[:5])
    data=np.asarray(data)
    return data

def plane_normal(p1,p2,p3):
    v1=p3-p1
    v2=p2-p1
    cross_prod=np.cross(v1,v2)
    norm=np.linalg.norm(cross_prod)
    return cross_prod/norm

def norm_12NNs(nn_list):
    if len(nn_list)<>13:
        print "error while normailizing 12 NN vectors"
    else:
        data=[]
        for i in range(1,13):
            a=(nn_list[i]-nn_list[0])/np.linalg.norm(nn_list[i]-nn_list[0])
            data.append(a)
        data=sorted(data,key=lambda x:x[2],reverse=True)
        data=np.asarray(data)
    return data

def rotation_M_ref_001(data):
    ref=np.asarray(
        [[ 0.        ,-0.70710678, 0.70710678],
         [ 0.        , 0.70710678, 0.70710678],
         [ 0.70710678, 0.        , 0.70710678],
         [-0.70710678, 0.        , 0.70710678],
         [ 0.70710678,-0.70710678, 0.        ],
         [-0.70710678,-0.70710678, 0.        ],
         [ 0.70710678, 0.70710678, 0.        ],
         [-0.70710678, 0.70710678, 0.        ],
         [ 0.        ,-0.70710678,-0.70710678],
         [ 0.        , 0.70710678,-0.70710678],
         [ 0.70710678, 0.        ,-0.70710678],
         [-0.70710678, 0.        ,-0.70710678]]
        )
    print(ref)
    print(data)
    #data[0],data[1],data[2]=data[1],data[2],data[0]
    #a=np.linalg.lstsq(ref[:3],data[:3])
    a=np.linalg.lstsq(ref[:3],data[:3])
    return a[0]

def rotation_M_ref_111(data):
    ref=np.asarray(
        [[ 0.        , 0.70710678, 0.70710678],
         [-0.80178373,-0.26726126, 0.53452247],
         [ 0.80178373,-0.26726126, 0.53452247],
         [ 0.70710678,-0.70710678, 0.        ],
         [-1.        , 0.        , 0.        ],
         [ 1.        , 0.        , 0.        ],
         [ 0.70710678, 0.70710678, 0.        ],
         [-0.70710678, 0.70710678, 0.        ],
         [-0.70710678,-0.70710678, 0.        ],
         [-0.80178373, 0.26726126,-0.53452247],
         [ 0.80178373, 0.26726126,-0.53452247],
         [ 0.        ,-0.70710678,-0.70710678]]
        )
    a=np.linalg.lstsq(ref,data)
    return a[0]

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
    #print NN_list[0][1]
    #a=id_orientation(NN_list[0][1])
    a=norm_12NNs(NN_list[0][1])
    b0=rotation_M_ref_001(a)  #rotation matrix ref to [100][010][001]
    #b1=rotation_M_ref_111(a)  #rotation matrix ref to [0-11][-211][111]
    print b0
    print np.linalg.det(b0)

    return



#all_NNs('example/dump.000001677.cfg',2,4.0)
all_NNs('example/Ag_111.cfg',1,4.0)
#all_NNs('example/Ag_001.cfg',1,4.0)

