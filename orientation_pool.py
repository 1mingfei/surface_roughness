#!~/usr/bin/python
import formatxfer as fx
import get_NN
import numpy as np
import multiprocessing
from itertools import permutations

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
    data=np.asarray(data)
    return data

def plane_normal(p1,p2,p3):
    v1=p3-p1
    v2=p2-p1
    cross_prod=np.cross(v1,v2)
    norm=np.linalg.norm(cross_prod)
    return cross_prod/norm

def norm_NN_vectors(nn_list):
    if len(nn_list)==13:
        data=[]
        for i in range(1,13):
            a=(nn_list[i]-nn_list[0])/np.linalg.norm(nn_list[i]-nn_list[0])
            data.append(a)
        data=sorted(data,key=lambda x:x[2],reverse=True)
        data=np.asarray(data)
    elif len(nn_list)==7: 
        data=[]
        for i in range(1,7):
            a=(nn_list[i]-nn_list[0])/np.linalg.norm(nn_list[i]-nn_list[0])
            data.append(a)
        data=sorted(data,key=lambda x:x[2],reverse=True)
        data=np.asarray(data)
    else:
        print "error while normailizing NN vectors (#of 1NN <>12 or #of 2NN <>6)"
        print len(nn_list)
    return data

def orthogonal_2nn(nn_list):
    a=norm_NN_vectors(nn_list)[:3]
    #print a
    b=[[],[],[]]
    l = list(permutations(range(0, 3)))
    sco=[]
    for i in l:
        b[0]=a[i[0]]
        b[1]=a[i[1]]
        b[2]=a[i[2]]
        sco.append(np.linalg.det(b))
    print sco
    i= np.argmax(sco)
    c=l[i]
    b[0]=a[c[0]]
    b[1]=a[c[1]]
    b[2]=a[c[2]]
    b=np.asarray(b)
    print "det of V_2nn:"
    print np.linalg.det(b) 
    return b

def rotation_M_2NN_ref_001(data):
    ref=np.asarray(
        [[ 1.        , 0.        , 0.        ],
         [ 0.        , 1.        , 0.        ],
         [ 0.        , 0.        , 1.        ]]
        )
    a=np.linalg.lstsq(data,ref)
    return a[0]

def decompose_R(a):
    theta_x=np.rad2deg(np.arctan2(a[2][1],a[2][2]))
    theta_y=np.rad2deg(np.arctan2(-a[2][0],np.sqrt(a[2][1]**2.0+a[2][2]**2.0)))
    theta_z=np.rad2deg(np.arctan2(a[1][0],a[0][0]))
    if theta_x <0 : theta_x+=180.0
    if theta_y <0 : theta_y+=180.0
    if theta_z <0 : theta_z+=180.0
    return [theta_x,theta_y,theta_z]
    


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
    #print(ref)
    #print(data)
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
    query_list=[79951,71056]
    NN_list = [pool.apply_async(get_NN.get_NN, args=(cell,data,i,cutoff)).get() for i in query_list]
    #print(NN_list[0][0])
    #print(NN_list[0][1])
    #print(NN_list[0][2])
    #print(NN_list[0][3])

    #a=id_orientation(NN_list[0][1])
    #a=norm_NN_vectors(NN_list[0][3])
    b=orthogonal_2nn(NN_list[0][3])
    b0=rotation_M_2NN_ref_001(b)  #rotation matrix ref to [100][010][001]
    b1=decompose_R(b0)
    print b1

    return



all_NNs('example/dump.000150000.cfg',2,[4.0,5.0])
all_NNs('example/Ag_111.cfg',1,[4.0,5.0])
#all_NNs('example/Ag_001.cfg',1,[4.0,5.0])

