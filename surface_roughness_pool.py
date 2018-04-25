#!~/usr/bin/python
import formatxfer as fx
import numpy as np
import multiprocessing
from numpy import fft
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

def get_surf_atoms(cell,data,n_type,limits):
    tot_num=int(len(data))
    l2 = [data[i] for i in range(tot_num) if data[i][1] > limits[0] and  data[i][1] < limits[1] and  data[i][2] > limits[2] and  data[i][2] < limits[3] and  data[i][0]==1.07868200e+02 ]
    l2 = sorted(l2,key=lambda x:x[3],reverse=True)
    if len(l2)==0: pass
    else:
        return(l2[0])

def cutoff_surf_atoms(data,cutoffs):
    tot_num=int(len(data))
    l2 = [data[i] for i in range(tot_num) if data[i][3] > cutoffs[0] and  data[i][3] < cutoffs[1] and  data[i][0]==1.07868200e+02 ]
    l2 = sorted(l2,key=lambda x:x[3],reverse=True)
    if len(l2)==0: pass
    else:
        return(l2)

def surface_atoms_filter(a,N_type,N_space,estimate):
    aa=fx.info(a,'cfg',N_type)
    cell=aa.cell
    data=aa.data     
    n_type=aa.atom_type_num
    tot_num=int(len(data))

    data_new=cutoff_surf_atoms(data,estimate)
    aa.data=np.asarray(data_new)
    aa.tot_num=int(len(aa.data))
    aa.atom_type_num=[n_type[1]]
    aa.get_cfg_file_one('cutoff.cfg')


    xx, yy = np.meshgrid(np.linspace(0.0, 1.0, N_space-1),np.linspace(0.0, 1.0, N_space-1))
    M = np.array([[ x1, x1+1.0/float(N_space), x2, x2+1.0/float(N_space)] for x1, x2 in zip(np.ravel(xx), np.ravel(yy))])
    pool = multiprocessing.Pool(processes=8)
    surface_atoms = [pool.apply_async(get_surf_atoms, args=(cell,data_new,n_type,M[i])).get() for i in range(len(M))]
    #surface_atoms = [pool.apply_async(get_surf_atoms, args=(cell,data,n_type,M[i])).get() for i in range(len(M))]
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
    return

def surface_statistics(a,N_space):
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

    lowest_h=np.amin(data, axis=0)[3] 
    h_bar=np.mean(data[:,3])-lowest_h
    w=np.std(data[:,3])
    print "mean height:",h_bar*cell[2][2]
    print "RMS ave of height devieation Rq(std):",w*cell[2][2]
    tmp=np.zeros((len(data),1))
    for i in range(len(data)):
        tmp[i]=abs(data[i][3]-h_bar)
        if tmp[i]<0.0 : tmp[i]+=1.0
    print "Arithmetic average of the absolute values",np.mean(tmp)*cell[2][2]
    for i in range(len(data)):
        data[i][3]-=lowest_h
        if data[i][3]<0.0 : data[i][3]+=1.0
    samples=[]
    M=np.zeros((N_space,N_space))
    for ii in range(N_space):
        xmin=ii/float(N_space)
        xmax=(ii+1)/float(N_space)
        for jj in range(N_space):
            ymin=jj/float(N_space)
            ymax=(jj+1)/float(N_space)
            flag1=0
            for i in range(len(data)):
                if (data[i][1]>xmin and data[i][1]<xmax and data[i][2] > ymin and data[i][2] < ymax):
                    samples.append(np.asarray(data[i]))
                    M[ii][jj]=float(data[i][3])
                    flag1=1
                    break
            if flag1==0:
                samples.append(np.array([1.07868200e+02,(xmin+xmax)/2.0,(ymin+ymax)/2.0,0.0]))
                M[ii][jj]=0.0
    samples=np.asarray(samples)
    M=np.asarray(M)
    #----------generate cfg--------------
    aa=fx.info('surface.cfg','cfg',1)
    aa.data=np.asarray(samples)
    aa.tot_num=int(len(aa.data))
    aa.atom_type_num=[aa.tot_num]
    aa.get_cfg_file_one('sample.cfg')
    #-----done generate cfg--------------
    #Fk = fft.fft2(M)/N_space
    #Fk = fft.fft2(M)/N_space # Fourier coefficients (divided by n) nu = fft.fftfreq(n,dx) # Natural frequencies
    #Fk = fft.fftshift(Fk) # Shift zero freq to center
    #freqs = fft.fftfreq(N_space)
    #freqs = fft.fftshift(freqs)
    ##nu = fft.fftshift(nu) # Shift zero freq to center
    #psd2D = np.abs( Fk )**2
    #X, Y = np.meshgrid(freqs, freqs)
    #fig = plt.figure()
    #ax = plt.axes(projection='3d')
    #ax.contour3D(X, Y, psd2D, 50, cmap='binary')
    #fig.savefig('1.png')
    return
'''
sampling the surface and statistics
'''
N=50
surface_atoms_filter('example/dump.000001677.cfg',2,N,[0.6,1.0])
surface_statistics('surf.dat',N)
