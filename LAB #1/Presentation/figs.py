import numpy as np, json, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'cmss10','mathtext.fontset':'cm','axes.unicode_minus':False,
 'font.size':8,'axes.linewidth':0.6,'lines.linewidth':1.2,'axes.grid':True,'grid.color':'#d0d0d0','grid.linewidth':0.4,
 'legend.fontsize':6.5,'legend.framealpha':0.9})
d=np.load('curves.npz'); c=json.load(open('calc.json'))
f=d['f']; Hr=20*np.log10(abs(d['Href'])); Hc=20*np.log10(abs(d['Hch'])); Hs=20*np.log10(abs(d['Hst']))
RED='#B3071B'; TEAL='#1b7f8e'; ORG='#d9822b'; GR='#555555'; fp=1e6; fs=4e6
def mask(ax,full=True):
    ax.axhline(-0.5,color=GR,ls='--',lw=0.8); 
    if full: ax.axhline(-50,color=GR,ls='--',lw=0.8)
    ax.axvline(fp,color=GR,ls=':',lw=0.8); ax.axvline(fs,color=GR,ls=':',lw=0.8)
# 1 full response
fig,ax=plt.subplots(figsize=(4.0,2.45))
ax.semilogx(f,Hr,color='k',lw=1.8,label='Ideal elliptic')
ax.semilogx(f,Hc,color=RED,ls='--',lw=1.1,label='From components')
ax.semilogx(f,Hs[0],color=TEAL,lw=0.9,label='Biquad 1 (low $Q$)')
ax.semilogx(f,Hs[1],color=ORG,lw=0.9,label='Biquad 2 (high $Q$)')
mask(ax); ax.set_ylim(-100,20); ax.set_xlim(1e3,1e8)
ax.text(1.15e3,-47,r'$-50$ dB',fontsize=6.5,color=GR); ax.text(1.15e3,-3.5,r'$-0.5$ dB',fontsize=6.5,color=GR)
ax.set_xlabel('Frequency [Hz]'); ax.set_ylabel('Magnitude [dB]'); ax.legend(loc='lower left')
fig.tight_layout(pad=0.3); fig.savefig('figs/fig_resp.pdf'); plt.close(fig)
# 2 zooms
fig,(a1,a2)=plt.subplots(1,2,figsize=(5.8,2.2))
a1.semilogx(f,Hr,color='k',lw=1.8); a1.semilogx(f,Hc,color=RED,ls='--',lw=1.1)
a1.axhline(0,color=GR,ls='--',lw=0.8); a1.axhline(-0.5,color=GR,ls='--',lw=0.8)
a1.set_xlim(1e3,fp); a1.set_ylim(-0.6,0.1); a1.set_title('Pass band (ripple 0 to $-0.5$ dB)',fontsize=8)
a2.semilogx(f,Hr,color='k',lw=1.8,label='Ideal elliptic'); a2.semilogx(f,Hc,color=RED,ls='--',lw=1.1,label='From components')
a2.axhline(-50,color=GR,ls='--',lw=0.8); a2.set_xlim(fs,1e8); a2.set_ylim(-100,-40); a2.set_title('Stop band (below $-50$ dB)',fontsize=8)
a2.legend(loc='lower right')
for a in (a1,a2): a.set_xlabel('Frequency [Hz]'); a.set_ylabel('Magnitude [dB]')
fig.tight_layout(pad=0.3); fig.savefig('figs/fig_zoom.pdf'); plt.close(fig)
# 3 cumulative noise (rms)
fi=d['fi']; P=d['Pcum']; m=(fi>=1e4)&(fi<=1e10)
fig,ax=plt.subplots(figsize=(2.9,2.45))
ax.semilogx(fi[m],np.sqrt(P[m])*1e6,color=RED,lw=1.5,label='Resistors (cumulative)')
ax.axhline(120,color='k',ls='--',lw=0.9); ax.text(1.3e4,124,r'spec: $120\ \mu\mathrm{V_{rms}}$',fontsize=6.5)
ax.axhline(c['rms_uV'],color=GR,ls=':',lw=0.9); ax.text(1.3e4,c['rms_uV']-14,r'%.1f $\mu\mathrm{V_{rms}}$'%c['rms_uV'],fontsize=6.5,color=GR)
ax.set_ylim(0,140); ax.set_xlabel('Upper integration limit [Hz]'); ax.set_ylabel(r'Output noise [$\mu\mathrm{V_{rms}}$]')
fig.tight_layout(pad=0.3); fig.savefig('figs/fig_noise.pdf'); plt.close(fig)
# tables
P=np.array(c['Pres']); tot=P.sum()
def kohm(x): return f'{x/1e3:.2f}'
t=r'''\begin{tabular}{lcc}
\toprule
 & Biquad 1 & Biquad 2 \\
\midrule
$f_0$ [MHz] & %.3f & %.3f \\
$Q$ & %.3f & %.3f \\
$f_z$ [MHz] & %.3f & %.3f \\
$k_0$ & %.3f & %.3f \\
\bottomrule
\end{tabular}'''%(c['f0'][0]/1e6,c['f0'][1]/1e6,c['Q'][0],c['Q'][1],c['fz'][0]/1e6,c['fz'][1]/1e6,c['k0'][0],c['k0'][1])
open('tab_biquads.tex','w').write(t)
t=r'''\begin{tabular}{lcc}
\toprule
 & Biquad 1 & Biquad 2 \\
\midrule
$C_a=C_b$ [pF] & 15 & 15 \\
$R_1$ [k$\Omega$] & %s & %s \\
$R_3$ [k$\Omega$] & %s & %s \\
$R_4$ [k$\Omega$] & %s & %s \\
$|R_5|$ [k$\Omega$] & %s & %s \\
$|C_4|$ [pF] & %.3f & %.3f \\
\bottomrule
\end{tabular}'''%(kohm(c['R1'][0]),kohm(c['R1'][1]),kohm(c['R3'][0]),kohm(c['R3'][1]),kohm(c['R4'][0]),kohm(c['R4'][1]),kohm(c['R5'][0]),kohm(c['R5'][1]),c['C4'][0]*1e12,c['C4'][1]*1e12)
open('tab_components.tex','w').write(t)
rows=''
for i in range(2):
    rows+=r'Biquad %d & '%(i+1)+' & '.join(f'{np.sqrt(P[i,m])*1e6:.1f}' for m in range(4))+f' & {np.sqrt(P[i].sum())*1e6:.1f} \\\\\n'
t=r'''\begin{tabular}{lccccc}
\toprule
{[$\mu$V$_{\rm rms}$ at filter output]} & $R_1$ & $R_3$ & $R_4$ & $R_5$ & Biquad total \\
\midrule
'''+rows+r'''\midrule
Filter total & \multicolumn{5}{c}{%.1f} \\
\bottomrule
\end{tabular}'''%c['rms_uV']
open('tab_noise.tex','w').write(t); print(t)
