#/usr/bin/env python3
#--*--coding:utf-8--*--

import math as mt
import numpy as np
from ggfcoefs import CALS_coesf

'''
this is rewrite from CALS10kfield.f from Monika Konte
'''

class cField():
    def __init__(self,fAdres=None,time=None,lat=None,lon=None):
        nMax = 10
        df = CALS_coesf(fAdres='./GGF100k',time=time).coefData
        self.field = self.magnetF(lat,lon,nMax,df)

    def coords(self, h,lat):
        #theta = mt.radians(90-lat)
        b1=40680925
        b2=40408585
        #pi=3.14159265
        theta = mt.pi/2 - mt.radians(90-lat)
        clat = mt.cos(theta)
        slat = mt.sin(theta)
        one = b1*clat**2
        two = b2*slat**2
        three = one + two
        four = mt.sqrt(three)
        r = mt.sqrt(h*(h+2*four)+(b1*one+b2*two)/three)
        cd = (h+four)/r
        sd = (b1-b2)/four*slat*clat/r
        sinth = slat*cd - clat*sd
        costh = clat*cd + slat*sd
        theta = mt.pi/2 - mt.atan2(sinth,costh)

        return cd,sd,theta,r

    def magnetF(self,lat,lon,nMax,df):
        '''
        Jeremy Davis, 2004, Mathmatical Modelling of Earth's Magnetic Field
        '''

        if abs(abs(lat)-90)==0:
            lat = lat - lat/abs(lat)*10**-5
            #print(lat)

        cd,sd,theta,r = self.coords(h=0,lat=lat)
        phi = mt.radians(lon)
        a=6371.2

        B = []

        z,x,y = 0,0,0
        P11 = 1
        P10 = P11
        dP11 = 0
        dP10 = dP11

        for m in np.arange(0,nMax,step=1):
            for n in np.arange(1,nMax,step=1):
                if m<=n:
                    if n==m:
                        P2 = mt.sin(theta)*P11
                        dP2 = mt.sin(theta)*dP11 + mt.cos(theta)*P11
                        P11=P2; P10=P11; P20=0;
                        dP11=dP2; dP10=dP11; dP20=0;
                    elif n==1:
                        P2 = mt.cos(theta)*P10
                        dP2 = mt.cos(theta)*dP10 - mt.sin(theta)*P10
                        P20=P10; P10=P2;
                        dP20=dP10; dP10=dP2;
                    else:
                        K = ((n-1)**2-m**2)/((2*n-1)*(2*n-3))
                        P2 = mt.cos(theta)*P10 - K*P20
                        dP2 = mt.cos(theta)*dP10 - mt.sin(theta)*P10 - K*dP20
                        P20=P10; P10=P2;
                        dP20=dP10; dP10=dP2;

                    df_t = df[(df['degree']==n) & (df['order']==m)]
                    g = df_t['g'].values[0]
                    h = df_t['h'].values[0]

                    z = z + (a/r)**(n+2)*(-n-1)*((g*mt.cos(m*phi)+h*mt.sin(m*phi))*P2)
                    x = x + (a/r)**(n+2)*((g*mt.cos(m*phi)+h*mt.sin(m*phi))*dP2)
                    y = y + (a/r)**(n+2)*(m*(g*mt.sin(m*phi)-h*mt.cos(m*phi))*P2)

        y = y/mt.sin(theta)
        xs = x
        x = x*cd + z*sd
        z = z*cd - xs*sd
        #ys = -y/mt.sin(theta)

        h = mt.sqrt(x**2+y**2)
        f = mt.sqrt(h**2+z**2)
        i = mt.degrees(mt.asin(z/f))
        d = mt.degrees(mt.atan2(y,x))
        return [d,i,f/1000]


def main():
    fAdres = './CALS10k.2'
    time = 1000
    lon = 30
    lat = 40

    print(cField(fAdres,time,lat,lon).field)

if __name__ == '__main__':
    main()
