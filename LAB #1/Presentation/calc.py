import numpy as np, json
from scipy import signal
fp=1e6; BW=3e6; fs=fp+BW; N=4; Rp=0.5; Rs=50
Ca=np.array([15e-12,15e-12]); Cb=np.array([15e-12,15e-12]); lowQfirst=True
b,a=signal.ellip(N,Rp,Rs,2*np.pi*fp,btype='low',analog=True)
z,p,k=signal.tf2zpk(b,a)
Hdc=abs(np.polyval(b,0)/np.polyval(a,0))
pp=p[p.imag>0]; zz=z[z.imag>0]; nst=len(pp)
w0_all=abs(pp); Q_all=w0_all/(-2*pp.real)
ix=np.argsort(Q_all); Q_all=Q_all[ix]; w0_all=w0_all[ix]
wz_all=np.sort(abs(zz)); wz_pair=wz_all[::-1]
ordr=np.arange(nst) if lowQfirst else np.arange(nst)[::-1]
w0=w0_all[ordr]; Q=Q_all[ordr]; wz=wz_pair[ordr]
g=np.ones(nst); g[-1]=Hdc; k0=g; k2=g/wz**2
R1=Q/(w0*Ca); R3=1/(w0*Ca); R4=1/(w0*Cb*k0); R5=-1/(w0*Cb); C4=-(w0**2)*k2*Ca
f=np.logspace(3,8,5000); s=1j*2*np.pi*f
Hst=[]
for i in range(nst):
    num=-(C4[i]/Ca[i])*s**2+1/(R3[i]*R4[i]*Ca[i]*Cb[i])
    den=s**2+s/(R1[i]*Ca[i])-1/(R3[i]*R5[i]*Ca[i]*Cb[i])
    Hst.append(num/den)
Hst=np.array(Hst); Hch=Hst.prod(axis=0)
_,Href=signal.freqs(b,a,2*np.pi*f)
err=np.max(abs(20*np.log10(abs(Hch))-20*np.log10(abs(Href))))
kB=1.380649e-23; T=300
fi=np.logspace(0,10,200000); si=1j*2*np.pi*fi
r=[abs(R1),abs(R3),abs(R4),abs(R5)]
Pres=np.zeros((nst,4)); Pcum_tot=[]; Pcum_all=np.zeros_like(fi)
def cumtrapz(y,x):
    return np.concatenate([[0],np.cumsum(0.5*(y[1:]+y[:-1])*np.diff(x))])
for i in range(nst):
    D=si**2/w0[i]**2+si/(w0[i]*Q[i])+1
    Gn=np.ones_like(fi,dtype=complex)
    for j in range(i+1,nst):
        Gn=Gn*(k0[j]+k2[j]*si**2)/(si**2/w0[j]**2+si/(w0[j]*Q[j])+1)
    Hr=np.array([(si/(w0[i]*Q[i]))/D,(si/w0[i]*Cb[i]/Ca[i])/D,k0[i]*np.ones_like(fi)/D,np.ones_like(fi)/D])*Gn
    Sv=4*kB*T*np.array([r[0][i],r[1][i],r[2][i],r[3][i]])[:,None]
    Sout=abs(Hr)**2*Sv*2
    Pres[i]=np.trapezoid(Sout,fi,axis=1)
    Pcum_all+=sum(cumtrapz(Sout[m],fi) for m in range(4))
Ptot=Pres.sum()
out=dict(Hdc=Hdc,err=err,f0=(w0/2/np.pi).tolist(),Q=Q.tolist(),fz=(wz/2/np.pi).tolist(),k0=k0.tolist(),
 R1=r[0].tolist(),R3=r[1].tolist(),R4=r[2].tolist(),R5=r[3].tolist(),C4=abs(C4).tolist(),
 Pres=Pres.tolist(),Ptot=Ptot,rms_uV=np.sqrt(Ptot)*1e6)
json.dump(out,open('calc.json','w'),indent=1)
np.savez('curves.npz',f=f,Href=Href,Hch=Hch,Hst=Hst,fi=fi,Pcum=Pcum_all)
print(json.dumps(out,indent=1))
