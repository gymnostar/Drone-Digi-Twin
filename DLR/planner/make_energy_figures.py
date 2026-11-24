"""Energy-model figures (repo root):  .venv/bin/python DLR/planner/make_energy_figures.py"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from battery_analysis import *  # noqa: F401,F403
FIG=str(Path(__file__).resolve().parent / 'fig') + '/'
import os; os.makedirs(FIG,exist_ok=True)
C=dict(blade="#0072B2",induced="#D55E00",parasite="#009E73",total="#222222")
plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,"axes.grid":True,"grid.alpha":.25})
V=np.linspace(0,30,301)

# 1 power breakdown
fig,ax=plt.subplots(figsize=(7.5,4.6))
for d,ls,lab in ((E.ZENG2019,"-","Zeng 2019 (2 kg)"),(E.AIRSIM_GENERIC_QUAD,"--","AirSim generic quad (1 kg)")):
    b,i,p=E.power_components(V,d)
    for y,k,n in ((b,"blade","blade profile"),(i,"induced","induced"),(p,"parasite","parasite")):
        ax.plot(V,y,ls,color=C[k],lw=1.6,label=f"{n}" if d is E.ZENG2019 else None)
    ax.plot(V,b+i+p,ls,color=C["total"],lw=2.2,label=f"total, {lab}")
ax.axvline(E.max_endurance_speed(),color="grey",ls=":",lw=1); ax.axvline(E.max_range_speed(),color="grey",ls=":",lw=1)
ax.text(E.max_endurance_speed()+.2,325,"max endurance\n%.1f m/s"%E.max_endurance_speed(),fontsize=8,color="grey")
ax.text(E.max_range_speed()+.2,325,"max range\n%.1f m/s"%E.max_range_speed(),fontsize=8,color="grey")
ax.set_xlabel("airspeed V (m/s)"); ax.set_ylabel("propulsion power (W)"); ax.set_ylim(0,370)
ax.set_title("Level-flight power vs airspeed (solid Zeng 2019, dashed AirSim set)",fontsize=9.5)
ax.legend(fontsize=8,loc="upper center",ncol=2,frameon=False,bbox_to_anchor=(.5,.88))
fig.tight_layout(); fig.savefig(FIG+"energy_power_breakdown.png",dpi=150); plt.close(fig)

# 2 J/m vs ground speed
fig,ax=plt.subplots(figsize=(7.5,4.6)); vg=np.linspace(3,26,461)
for w,lab,col in ((-5,"headwind 5 m/s","#D55E00"),(0,"still air","#222222"),(5,"tailwind 5 m/s","#0072B2")):
    e=E.power_level(np.abs(vg-w))/vg; ax.plot(vg,e,color=col,lw=2,label=lab)
    i=np.argmin(e); ax.plot(vg[i],e[i],"o",color=col,ms=7)
    ax.annotate(f"{vg[i]:.1f} m/s\n{e[i]:.2f} J/m",(vg[i],e[i]),xytext=(vg[i],2.5),fontsize=8,color=col,ha="center",arrowprops=dict(arrowstyle="-",color=col,lw=.6))
ax.set_ylim(0,30); ax.set_xlabel("ground speed (m/s)"); ax.set_ylabel("propulsion energy per metre (J/m)")
ax.set_title("Zeng 2019: energy per metre of track; dots = optimal ground speed",fontsize=10); ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(FIG+"energy_per_metre_wind.png",dpi=150); plt.close(fig)

# 3 battery losses
sp=np.arange(4,28.01,0.1)
fig,(a1,a2)=plt.subplots(1,2,figsize=(11,4.6))
res={}
for w,lab,col in ((0,"still air","#222222"),(-5,"headwind 5 m/s","#D55E00"),(5,"tailwind 5 m/s","#0072B2")):
    a=sweep(w,speeds=sp); res[w]=a
    ip,ic=np.argmin(a[:,1]),np.argmin(a[:,3])
    a1.plot(a[:,0],a[:,1]/1e3,"--",color=col,lw=1.3); a1.plot(a[:,0],a[:,3]/1e3,"-",color=col,lw=2,label=lab)
    a1.plot(a[ip,0],a[ip,1]/1e3,"o",mfc="white",color=col,ms=6); a1.plot(a[ic,0],a[ic,3]/1e3,"o",color=col,ms=6)
    print(lab,"prop-opt %.2f cell-opt %.2f shift %.2f m/s, E penalty %.3f%%, loss at opt %.2f%%"%(a[ip,0],a[ic,0],a[ic,0]-a[ip,0],100*(a[ip,3]/a[ic,3]-1),100*(a[ic,3]-a[ic,2])/a[ic,3]))
a1.plot([],[],"k--",lw=1.3,label="dashed: propulsion energy"); a1.plot([],[],"k-",lw=2,label="solid: drawn from cells")
a1.set_ylim(0,40); a1.set_xlabel("ground speed (m/s)"); a1.set_ylabel("energy for a 1 km trip (kJ)")
a1.set_title("1 km trip, 6S 5 Ah pack, motor+ESC eff. 0.8\nopen dot: optimum without battery, filled dot: with battery",fontsize=9.5)
a1.legend(frameon=False,fontsize=8,loc="upper right")
a=res[0]; ip,ic=np.argmin(a[:,1]),np.argmin(a[:,3])
a1.annotate("optimum shifts by\n%.2f m/s (still air)\n%.3f%% energy"%(a[ic,0]-a[ip,0],100*(a[ip,3]/a[ic,3]-1)),(a[ic,0],a[ic,3]/1e3),xytext=(8.5,24),fontsize=8,arrowprops=dict(arrowstyle="->",color="grey"))
# panel 2: loss fraction
for lab,pack,s0,col,ls in (("fresh, warm, SOC 100%",B.PACK_6S_5000,1.0,"#222222","-"),("cold/aged (R x3), SOC 100%",B.PACK_COLD,1.0,"#E69F00","-"),("cold/aged (R x3), SOC 20%",B.PACK_COLD,0.2,"#CC79A7","-")):
    a=sweep(0,pack,s0,speeds=sp); a2.plot(a[:,0],100*(a[:,3]-a[:,2])/a[:,3],ls,color=col,lw=2,label=lab)
    ic=np.argmin(a[:,3]); a2.plot(a[ic,0],100*(a[ic,3]-a[ic,2])/a[ic,3],"o",color=col)
a2.set_xlabel("ground speed (m/s), still air"); a2.set_ylabel("battery loss, I$^2$R share of cell energy (%)")
a2.set_title("Battery loss fraction vs speed (dots: energy-optimal speed)",fontsize=9.5); a2.legend(frameon=False,fontsize=8,loc="upper left"); a2.set_ylim(0,None)
fig.tight_layout(); fig.savefig(FIG+"energy_battery_losses.png",dpi=150); plt.close(fig)

# 4 momentum check
fig,(a1,a2)=plt.subplots(2,1,figsize=(7.5,6),sharex=True,gridspec_kw=dict(height_ratios=[3,1.3]))
V=np.linspace(0,30,301); vec=np.stack([V,0*V],-1)
for d,col,lab in ((E.ZENG2019,"#222222","Zeng 2019"),(E.AIRSIM_GENERIC_QUAD,"#0072B2","AirSim set")):
    z=E.power_level(V,d); m=E.power_momentum(vec,d); m3=E.power_momentum(vec,d,mu_coeff=3.0)
    a1.plot(V,z,color=col,lw=2,label=lab+": Zeng closed form"); a1.plot(V,m,"--",color="#D55E00" if d is E.ZENG2019 else "#E69F00",lw=2,label=lab+": momentum + Glauert")
    a2.plot(V,100*(m/z-1),color=col,lw=1.8,label=lab+" (profile coeff 4.65)"); a2.plot(V,100*(m3/z-1),":",color=col,lw=1.5,label=lab+" (profile coeff 3)")
a1.set_ylabel("propulsion power (W)"); a1.legend(frameon=False,fontsize=8); a1.set_title("Zeng closed form vs first-principles momentum theory (level flight)",fontsize=10)
a2.axhline(0,color="grey",lw=.8); a2.set_ylabel("momentum vs Zeng (%)"); a2.set_xlabel("airspeed V (m/s)"); a2.legend(frameon=False,fontsize=7,ncol=2)
fig.tight_layout(); fig.savefig(FIG+"energy_momentum_check.png",dpi=150); plt.close(fig)
