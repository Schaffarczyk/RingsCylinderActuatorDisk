# Analysis of the FWRV-AD runs 011/012/013 (geometric panelling) and 014 (cosine), v3.4.1, and of run 010 (cosine, v3.4).
# Run in this directory: python3 fwrv-edge-analysis.py  -> decade tables on stdout, runs.npy for fwrv-edge-blocks.py
import numpy as np, sys
def load(f):
    d=np.loadtxt(f,comments='#',skiprows=5)
    return d
runs={'010':'schaffarczyk010quadWake.dat','011':'schaffarczyk011quadWake.dat','012':'schaffarczyk012quadWake.dat','013':'schaffarczyk013quadWake.dat','014':'schaffarczyk014quadWake.dat'}
out={}
for k,f in runs.items():
    d=load(f); z,r,g,L,c,s,th,vz,vr,vm=d[:,:10].T
    rho=np.hypot(z,r-1)
    ang=np.degrees(np.unwrap(np.arctan2(s,c)[::-1])[::-1])
    Vin=vm-np.abs(g)/2; Vout=vm+np.abs(g)/2
    print('===== run',k,' N=%d  L1=%.3e  rho_min=%.3e  z1=%.3e  n(z<0)=%d  Rinf=%.5f'%(len(z),L[0],rho[0],z[0],(z<0).sum(),r[-1]))
    print('   decade          n   theta   dth/dec  beta_geo  gamma    Vin    Vout  Vo/Vi  beta_vel   Vtm*|g|-4/9')
    rows=[]
    for e in range(-12,0):
        lo,hi=10.0**e,10.0**(e+1); m=(rho>=lo)&(rho<hi)
        if m.sum()>=3:
            p=np.polyfit(np.log10(rho[m]),ang[m],1); bg=-p[0]/np.degrees(np.log(10))
            bv=np.log(Vout[m].mean()/Vin[m].mean())/(2*np.pi)
            rows.append([e,e+1,m.sum(),ang[m].mean(),p[0],bg,g[m].mean(),Vin[m].mean(),Vout[m].mean(),Vout[m].mean()/Vin[m].mean(),bv,(vm[m]*np.abs(g[m])-4/9).max()])
            print('   [1e%3d,1e%3d) %5d %7.1f %7.1f  %.4f  %.4f  %.4f  %.4f  %.3f  %.4f   %.1e'%tuple(rows[-1]))
    for t in (90,180):
        i=np.where(ang>=t)[0]
        if len(i): print('   theta>=%d first at rho=%.2e (z=%.2e)'%(t,rho[i.max()],z[i.max()]))
        else: print('   theta max %.1f at rho=%.2e'%(ang.max(),rho[np.argmax(ang)]))
    out[k]=dict(rho=rho,ang=ang,g=g,Vin=Vin,Vout=Vout,z=z,r=r,L=L,vm=vm)
np.save('runs.npy',out,allow_pickle=True)
