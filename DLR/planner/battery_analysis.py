"""Numbers behind ENERGY_MODEL.md section 4 (battery vs propulsion energy). Imported by make_energy_figures.py."""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, energy as E, battery as B
L=1000.0
def trip(vg, w, pack=B.PACK_6S_5000, soc0=1.0, d=E.ZENG2019, dt=2.0):
    P=float(E.power_level(abs(vg-w),d)); T=L/vg
    r=B.energy_constant_power(P,T,pack,soc0,dt)
    return P*T, r
def sweep(w, pack=B.PACK_6S_5000, soc0=1.0, speeds=np.arange(4,24.01,0.25)):
    out=[]
    for v in speeds:
        ep,r=trip(v,w,pack,soc0); out.append((v,ep,r.e_elec,r.e_cell,r.soc_used))
    return np.array(out)
if __name__=="__main__":
    for name,w in (("still",0),("head5",-5),("tail5",5)):
        # wind convention: v_air = v_g - w ; headwind -> w=-5 along flight direction
        a=sweep(w)
        ip=np.argmin(a[:,1]); ic=np.argmin(a[:,3])
        print(name,"opt prop %.2f m/s %.1f kJ | opt cell %.2f m/s %.1f kJ | cell at prop-opt %.1f kJ  penalty %.3f%%"%(a[ip,0],a[ip,1]/1e3,a[ic,0],a[ic,3]/1e3,a[ip,3]/1e3,100*(a[ip,3]/a[ic,3]-1)),
              "| loss%% at opt %.2f"%(100*(a[ic,3]-a[ic,2])/a[ic,3]), "SOC %.2f%%"%(100*a[ic,4]))
    # integer speeds table, still
    for w in (0,-5,5):
        print("wind",w)
        for v in (4,8,12,16,20,24):
            ep,r=trip(v,w); print("  v=%2d prop %6.1f kJ elec %6.1f cell %6.1f socused %5.2f%% I %.1f A"%(v,ep/1e3,r.e_elec/1e3,r.e_cell/1e3,100*r.soc_used,r.current[0]))
    # harsh conditions
    for lab,pack,soc0 in (("cold R x3, SOC 100",B.PACK_COLD,1.0),("cold R x3, SOC 20",B.PACK_COLD,0.2),("R x3 + 4S",B.Pack(n_series=4,r_scale=3.0,r_cell=0.010),1.0),("RC pack",B.PACK_6S_5000_RC,1.0)):
        a=sweep(0,pack,soc0); ip=np.argmin(a[:,1]); ic=np.argmin(a[:,3])
        print(lab,"opt cell %.2f vs prop %.2f; penalty %.3f%%; loss at opt %.2f%%"%(a[ic,0],a[ip,0],100*(a[ip,3]/a[ic,3]-1),100*(a[ic,3]-a[ic,2])/a[ic,3]))
