# Rows of Table 'fwrvedge' and Table 'adiabatic' and the pgfplots coordinates of Fig. 'fwrvedge' from run 013 (G10+).
# Run after fwrv-edge-analysis.py: python3 fwrv-edge-blocks.py -> fwrv-edge-blocks.txt
import numpy as np
o=np.load('runs.npy',allow_pickle=True).item()
d=o['013']; rho,ang,g,Vin,Vout,z,r=d['rho'],d['ang'],d['g'],d['Vin'],d['Vout'],d['z'],d['r']
rows=[]
for e in range(-12,-1):
    lo,hi=10.0**e,10.0**(e+1); m=(rho>=lo)&(rho<hi)
    p=np.polyfit(np.log10(rho[m]),ang[m],1); bg=-p[0]/np.degrees(np.log(10))
    bv=np.log(Vout[m].mean()/Vin[m].mean())/(2*np.pi)
    rows.append(dict(e=e,n=m.sum(),th=ang[m].mean(),sl=-p[0],bg=bg,vi=Vin[m].mean(),vo=Vout[m].mean(),bv=bv))
# centred differences of ln G over decade means (as in tab:adiabatic)
G=np.array([np.sqrt(R['vi']*R['vo']) for R in rows]); lg=np.array([R['e']+0.5 for R in rows])
dG=-np.gradient(np.log(G),lg)
# ODE
f=lambda b:-(b**2/np.pi)*np.tanh(2*np.pi*b)
def integ(x0,b0,xs,h=-1e-3):
    x,b,res=x0,b0,{}
    for xt in xs:
        while x>xt:
            k1=f(b);k2=f(b+0.5*h*k1);k3=f(b+0.5*h*k2);k4=f(b+h*k3); b+=h/6*(k1+2*k2+2*k3+k4); x+=h
        res[xt]=b
    return res
i0=[R['e'] for R in rows].index(-4); b0=rows[i0]['bv']
pts=[np.log(10**(e+0.5)) for e in range(-4,-13,-1)]; sol=integ(np.log(3.16e-4),b0,pts)
print('b0=',b0)
out=[]
out.append('%% ---- tab:fwrvedge rows (run G10+ = 013, v3.4.1) ----')
for R in rows[::-1]:
    out.append('$10^{%d}$--$10^{%d}$ & %d & $%.1f^\\circ$ & $%.1f^\\circ$ & %.3f & %.3f & %.3f & %.2f & %.3f\\\\'%(R['e'],R['e']+1,R['n'],R['th'],R['sl'],R['bg'],R['vi'],R['vo'],R['vo']/R['vi'],R['bv']))
out.append('%% ---- tab:adiabatic rows ----')
for k,R in enumerate(rows[::-1]):
    if R['e']>-4: continue
    bo=sol[np.log(10**(R['e']+0.5))]
    out.append('$10^{%d}$--$10^{%d}$ & %.3f & %.3f & %.3f & %.3f & %.3f\\\\'%(R['e'],R['e']+1,R['bg'],R['bv'],bo,dG[len(rows)-1-k],R['bg']**2*np.log(10)))
# figure coordinates: theta(rho) thinned log-uniformly from 013
def thin(x,y,n=90):
    lx=np.log10(x); grid=np.linspace(lx.min(),lx.max(),n); idx=[np.argmin(abs(lx-gq)) for gq in grid]; idx=sorted(set(idx))
    return ' '.join('(%.3e,%.2f)'%(x[i],y[i]) for i in idx)
m=rho<=0.5
out.append('%% ---- fig left: theta_s(rho), run 013 ----')
out.append(thin(rho[m],ang[m]))
out.append('%% ---- fig right: beta_geo, beta_vel decade midpoints ----')
out.append(' '.join('(%.1e,%.4f)'%(10**(R['e']+0.5),R['bg']) for R in rows))
out.append(' '.join('(%.1e,%.4f)'%(10**(R['e']+0.5),R['bv']) for R in rows))
out.append('%% ---- fig right: adiabatic ODE curve from 3.16e-4 inward ----')
xs=[np.log(10**x) for x in np.arange(-3.5,-12.51,-0.25)]; s2=integ(np.log(3.16e-4),b0,xs)
out.append(' '.join('(%.2e,%.4f)'%(np.exp(x),s2[x]) for x in xs))
# crossing points
for t in (90,180):
    i=np.where(ang>=t)[0]; print('theta>=%d at rho=%.2e z=%.2e r-1=%.2e'%(t,rho[i.max()],z[i.max()],r[i.max()]-1))
print('n z<0:',(z<0).sum(),' n r<1:',(r<1).sum(),' min r-1 = %.2e at rho=%.2e'%((r-1).min(),rho[np.argmin(r-1)]))
print('ODE point values: 1e-10 %.3f 1e-11 %.3f 1e-12 %.3f 1e-13 %.3f 1e-14 %.3f'%tuple(integ(np.log(3.16e-4),b0,[np.log(10.0**k) for k in (-10,-11,-12,-13,-14)]).values()))
# decade-point values of dz where sheet crosses z=0
i=np.where(z<0)[0].max(); print('z<0 starts at rho=%.2e'%rho[i])
# fig 654 data: r-1 vs rho, z>0, -z for 013 thinned
out.append('%% ---- fig spiral right panel: (rho, r-1) ----'); out.append(thin(rho[m],r[m]-1,70).replace('%.2f','%.3e'))
open('fwrv-edge-blocks.txt','w').write('\n'.join(out))
print('\n'.join(out[:30]))
