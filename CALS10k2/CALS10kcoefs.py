#/usr/bin/env python3
#--*--coding:utf-8--*--

import pandas as pd
import numpy as np

'''
this is rewrite from CALS10kcoef.f from Monika Konte
'''
class CALS_coesf():
    def __init__(self,fAdres=None,time=None):
        lm,nm,nspl,tknts,gt = self.dataLoad(fAdres=fAdres)
        nleft = self.interv(tknts,time,nspl)
        spl = self.bspline(tknts,time,nspl,nleft)
        self.coefData = self.coefs(spl,gt,nleft)

    def dataLoad(self, fAdres=None):
        fo = open(fAdres)
        file = fo.readlines()
        line1 = np.array(file[1].split()).astype(np.float64)
        lm = line1[0]
        nm = line1[1]
        nspl = line1[2]
        tknts = line1[3::]
        gt = np.array(file[2].split()).astype(np.float64)
        gt.shape = (303,120)
        gt = gt.T
        fo.close()
        return lm,nm,nspl,tknts,gt

    def interv(self,tknts,time,nspl):
        if time >tknts[3] or time<tknts[-4]:
            for i in np.arange(4,len(tknts)-3):
                if time<tknts[i]:
                    nleft=i-1
                    break
        else:
            print('time not fit range')
            pass
        return nleft

    def bspline(self,tknts,t,nspl,nleft):
        spl = [1]*5
        #print(spl)
        deltal,deltar = [0]*4,[0]*4
        jorder = 4
        for j in np.arange(1,jorder):
            #print(tknts[nleft-j+1])
            deltar[j] = tknts[nleft+j]-t
            deltal[j] = t-tknts[nleft-j+1]
            saved = 0
            #print(j)
            for i in np.arange(1,j+1):
                #print(deltar[i],deltal[j+1-i])
                #print(deltal[j+1-i])
                #print(spl[i])
                term = spl[i]/(deltar[i]+deltal[j+1-i])
                spl[i] = saved+deltar[i]*term
                saved = deltal[j+1-i]*term
                #print(deltar[i])
            spl[j+1] = saved
        return spl[1::]

    def coefs(self,spl,gt,nleft):
        gg=[]
        for k in np.arange(1,121):
            g=0
            for j in np.arange(0,4):
                #print(k-1,j+nleft-3)
                #print(spl[j])
                g = g + spl[j]*gt[k-1,j+nleft-3]
            gg.append(g)

        #print(gg)

        h0 = 0
        coefData = []
        for l in np.arange(1,11):
            m=0
            k=l**2
            coefData.append([l,m,gg[k-1],h0])
            for m in np.arange(1,l+1):
                k=l**2+2*m-1
                coefData.append([l,m,gg[k-1],gg[k]])
        coefData = np.array(coefData)
        coefData = pd.DataFrame({'degree':coefData[:,0],
                                'order':coefData[:,1],
                                'g':coefData[:,2],
                                 'h':coefData[:,3]})
        return coefData

def main():
    fAdres = './CALS10k.2'
    time = 1900
    coefs = CALS_coesf(fAdres=fAdres, time=time)
    print(coefs.coefData)


if __name__ == '__main__':
    main()
