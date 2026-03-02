import matplotlib.pyplot as plt 
import numpy as np 

import os 

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.collections import LineCollection
from matplotlib import cm
from matplotlib import gridspec

from mpl_toolkits.axes_grid1.inset_locator import inset_axes


plt.rcParams.update({'font.size': 8})

time = []
rs = []
rb = []
thickness = []
ecc = []

dir_name = "/nfsy5/kihoulou/IcyBreathing/thickness_40km/data_Europa_ecc_0.04_30Myr"

infile = open(dir_name+"/solution_statistics.dat", "r") 
lines = infile.readlines() 

header = True
for line in lines:
    sline = line.split("\t")

    if (header==True):
        header = False
        continue

    time.append(float(sline[0]))
    rs.append(float(sline[1])/1e3 + 0.1) 
    rb.append(float(sline[2])/1e3 - 0.15)
    thickness.append((float(sline[1])- float(sline[2]))/1e3)
    ecc.append(float(sline[6]))

# fig, (ax0, ax1) = plt.subplots(2, 1, layout='constrained', figsize=(10, 6)) 
# fig.set_constrained_layout_pads(w_pad = 0.2, h_pad = 0.2, hspace=1, wspace=0)

# fig = plt.figure(figsize=(7.25,7.25/2.0))
ratio = (7.25/7.08333)*(7.250/7.243)*(7.250/7.249)
fig = plt.figure(figsize=(8*ratio, 4*ratio))
# fig = plt.figure(figsize=(3.55,2))

gs = gridspec.GridSpec(11, 11)
ax0 = fig.add_subplot(gs[0:5, 0:10])
ax1 = fig.add_subplot(gs[6:11, 0:10])
ax2 = fig.add_subplot(gs[6:11, 10:11])


# ax0.text(45, 0.07, "A", fontsize=9, weight = "bold")
# ax1.text(45, 1561, "B", fontsize=9, weight = "bold")

ax2.axis('off')


ax01 = ax0.twinx() 
# ax1a = ax1.twinx() 
ax0.set_xlim(0, 180)
ax0.set_xticks([])

ax1.set_xlim(0, 180)
ax1.set_xticks([])


ax0.tick_params(axis='y',  colors='black')
ax0.set_ylabel("Eccentricity (-)", labelpad = 5)
ax0.set_ylim(0.009, 0.1)
ax0.set_yticks(np.linspace(0.01, 0.1, num=4))
lns1 = ax0.plot(time, ecc, color='red', label="Eccentricity") 


ax01.set_ylabel("Shell thickness (km)", labelpad = 5)
ax01.tick_params(axis='y', colors='black')                                         
ax01.set_ylim(40, 0)
ax01.set_yticks([0,10,20,30,40])
lns2 = ax01.plot(time, thickness, color='black', label="Shell thickness") 

lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax0.legend(lns, labs, ncol=2, loc="upper right", bbox_to_anchor=(0.975, 1.25))

ax1.set_ylabel("Radius (km)", labelpad = 5)
ax1.set_xlim(0, 180)
ax1.set_xlabel("Time (Myr)", labelpad = 5)
ax1.set_xticks(np.linspace(0, 180, num=5))
ax1.set_yticks(np.linspace(1520, 1560, num=3))

bdt_y = []
bdt_x = []

file_bdt = open("plastic_line.dat","w")
for i in range(0,5000):    
    try:
        file = open(dir_name+"/data_steps/solution_"+str(i)+".dat", "r") 
        lines = file.readlines()   

        rr = []
        tt = []
        sigma_t = []
        sigma_II = [] 
        sigma_y = [] 

        header = True
        for line in lines:
            sline = line.split("\t")

            if (header == True):
                header = False
                continue

            rr.append(float(sline[0]))
            tt.append(i/10.0)
            sigma_t.append(float(sline[6])/1e6)
            sigma_II.append(float(sline[7])/1e6)
            sigma_y.append(float(sline[8])/1e6)
        
        if (i%5 == 0):
            r = np.array(rr)
            t = np.array(tt)
            st = np.array(sigma_t)

            points = np.array([t, r]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)

            norm = plt.Normalize(-4, 4)
            lc1 = LineCollection(segments, cmap="coolwarm", norm=norm)
            lc1.set_array(st)
            lc1.set_linewidth(1.5)
            line1 = ax1.add_collection(lc1)

        bdt_x.append(i/10.0)
        bdt_y.append(r[0])

        for j in range(0, len(sigma_y)-1):
            found = False
            for k in range(0, 100):
                r_int = rr[j] - k*(rr[j] - rr[j+1])/100.0
                sigma_II_int = sigma_II[j] - k*(sigma_II[j] - sigma_II[j+1])/100.0
                sigma_y_int = sigma_y[j] - k*(sigma_y[j] - sigma_y[j+1])/100.0

                if (sigma_II_int < 0.99*sigma_y_int):
                    bdt_y[len(bdt_y)-1] = r_int
                    found = True
                    break

            if (found == True):
                break
        
        file_bdt.write((4*"%.7E\t"+"\n")%(bdt_x[len(bdt_y)-1], bdt_y[len(bdt_y)-1], rr[0], abs(bdt_y[len(bdt_y)-1]- rr[0])))    

    except:
        pass

file_bdt.close()

ax1.plot(bdt_x, bdt_y, color='black', linewidth = 1, linestyle='-', dashes=[2,2])
# cb = fig.colorbar(line1, ax=ax2, pad = -0.08, aspect=15, ticks=[-4, 0, 4])
# cb.set_label(r"$\sigma_t\rm\ (MPa)$", labelpad = 5)
# ax1a.set_yticks([])

ax2.set_yticks([])
ax2.set_xticks([])
cmap = plt.cm.coolwarm
norm = plt.Normalize(0,1)
# cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax2, fraction=0.75, aspect = 15, label=r"$\sigma_{\theta\theta},\sigma_{\varphi\varphi}\rm\ (MPa)$", location='right')
cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax2, fraction=0.75, aspect = 15, label="Lateral stress (MPa)", location='right')
cbar.set_ticks(ticks=[0,  0.5, 1], labels=['-4', "0", "4"])
cbar.ax.tick_params(width=1, length=4)

ax1.plot(time, rs, color='black', linewidth=1) 
ax1.plot(time, rb, color='black', linewidth=1)



# plt.savefig("fig5_mpl.pdf", bbox_inches='tight', dpi = 800) 
# plt.savefig("../coupled_code/software_lateral_compression/final_figures_pdf/fig5.pdf", bbox_inches='tight', pad_inches=0, dpi = 1000) 
# plt.savefig("../coupled_code/software_lateral_compression/final_figures_jpg/fig5.jpg", bbox_inches='tight', pad_inches=0,  dpi = 1000)

plt.savefig("fig5_30Myr.png", bbox_inches='tight', pad_inches=0, dpi = 200)
# plt.savefig("fig5_mpl.jpg", bbox_inches='tight', pad_inches=0, dpi = 1000)

# os.system("pdfcrop fig5_mpl.pdf")
# os.system("identify -verbose fig5_mpl.pdf | grep 'Print size'")
# os.system("identify -verbose fig5_mpl-crop.pdf | grep 'Print size'")
