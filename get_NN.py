import numpy as np
#define a function to calculate dist btw two atoms
def calc_dist(lst,n1,n2,H):
    if  lst[n2][1]-lst[n1][1]>=0.5:
        a=(lst[n1][1]-lst[n2][1]+1)
    elif  lst[n2][1]-lst[n1][1]<(-0.5):
        a=(lst[n1][1]-lst[n2][1]-1)
    else:
        a=((lst[n1][1]-lst[n2][1]))
    if  lst[n2][2]-lst[n1][2]>=0.5:
        b=((lst[n1][2]-lst[n2][2]+1))
    elif  lst[n2][2]-lst[n1][2]<(-0.5):
        b=((lst[n1][2]-lst[n2][2]-1))
    else:
        b=((lst[n1][2]-lst[n2][2]))
    if  lst[n2][3]-lst[n1][3]>=0.5:
        c=((lst[n1][3]-lst[n2][3]+1))
    elif  lst[n2][3]-lst[n1][3]<(-0.5):
        c=((lst[n1][3]-lst[n2][3]-1))
    else:
        c=((lst[n1][3]-lst[n2][3]))
    tmp=[a,b,c]
    tmp_cart=np.dot(tmp,H)
    aa,bb,cc=tmp_cart[0],tmp_cart[1],tmp_cart[2]
    return np.sqrt(aa**2+bb**2+cc**2)

def get_NN(cell,data,i,cutoff):
    lst=[]
    for j in range(len(data)):
        if (calc_dist(data,i,j,cell) <= cutoff) and (i<>j):
            lst.append(j)
    return len(lst),lst
