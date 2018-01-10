import numpy as np
import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import tempfile, os

duration = 5
timesteps = np.arange(0, duration, 1)
ncol = 6
nrow = 13
nlay = 2
cmap = plt.cm.gist_rainbow_r
CBlabel = 'values ([-])'
plt_title = 'VALUES'

# create output folder
ws = os.path.join(tempfile.gettempdir(), "__ani")
f = 0
if os.path.exists(ws):
    ws1 = ws
    while os.path.exists(ws1):
        ws1 = '%s%d' % (ws, f)
        f +=1
    os.makedirs(ws1)
    ws = ws1
    del ws1
else:
    os.makedirs(ws)

# create coordinate arrays
x = np.arange(0.5, ncol+1.5, 1)
y = np.arange(0.5, nrow+1.5, 1)
xg,yg = np.meshgrid(x,y)

x = np.arange(1, ncol+1, 1)
y = np.arange(1, nrow+1, 1)
xg1,yg1 = np.meshgrid(x,y)

# create values to be plotted
V = np.zeros([duration,nrow,ncol,nlay], dtype = float)
for d in timesteps:
    for L in range(nlay):
        V[d,:,:,L] = ((np.random.rand(nrow, ncol)+0.5)*xg1+(np.random.rand(nrow, ncol)+0.5)*yg1)/(1.0+2.0*L)
Vmax = np.max(V)
Vmin = np.min(V)
ticks = np.linspace(Vmin,Vmax,5)

# plot initialization
ax= []
fig = plt.figure(num=None, figsize=(11.7, 8.27), dpi=30)
figtitle = fig.suptitle('')
for L in range(nlay):
    ax.append(fig.add_subplot(1,nlay,L+1, axisbg = 'silver'))
    ax[L].xaxis.set_ticks(np.arange(0,ncol+1,1))
    ax[L].yaxis.set_ticks(np.arange(0,nrow+1,1))
    plt.setp(ax[L].get_xticklabels(), fontsize=8)
    plt.setp(ax[L].get_yticklabels(), fontsize=8)
    plt.ylabel('row i', fontsize=10)
    plt.xlabel('col j', fontsize=10)
    ax[L].set_title('layer ' + str(L+1), fontsize = 10)

# plot sequences of grids
ims = []
for i, day in enumerate(timesteps):
    ims.append([])
    figtitle.set_text(plt_title + '\ntime step %s' % (day))
    plt.draw()
    for L in range(nlay):
        Vtmp = V[day,:,:,L]
        ims[i].append(ax[L].pcolormesh(xg, yg, Vtmp, cmap = cmap, vmin = Vmin, vmax = Vmax))
        # plot contours with labels
        #ims[i].append(ax[L].contour(xg1, yg1[::-1], Vtmp[::-1], ticks, colors = 'gray'))
        #ax[L].clabel(ims[i][4*i+L+1], inline=1, fontsize = 6, fmt='%2.2f', colors = 'gray')
        del Vtmp
        # modify axes range
        ax[L].set_ylim(bottom = np.max(yg1), top = np.min(yg1))
        ax[L].axis('scaled')
    # create color bar
    cax = fig.add_axes([0.035, 0.125, 0.025, 0.75])
    CB = fig.colorbar(ims[0][0], extend='both', ticks = ticks, format = '%2.2f', cax = cax,  orientation = 'vertical')
    CB.set_label(CBlabel, fontsize = 12)
    cax.yaxis.set_label_position('left')
    plt.setp(CB.ax.get_yticklabels(), fontsize = 7)
    # save image
    plt_export_fn = os.path.join(ws, '_plt_%s_timestep%05d.png' % (plt_title, day))
    plt.savefig(plt_export_fn)

# save animation
ani = animation.ArtistAnimation(fig, ims, interval=1000*(timesteps[1]-timesteps[0]), repeat_delay=3000, blit=True)
# save(filename, writer=None, fps=None, dpi=None, codec=None, bitrate=None, extra_args=None, metadata=None, extra_anim=None)
# command line
# ffmpeg -r 1 -i _plt_VALUES_timestep%05d.png -s:v 1280x720 -c:v libx264 -profile:v high -crf 23 -pix_fmt yuv420p -r 30 movie.mp4
ani.save(os.path.join(ws, '_plt_%s_mov.mp4' % plt_title), codec = 'libx264') #, fps=1, codec = 'libx264', bitrate = -1, metadata = dict(title='Dynamic array', artist='Matplotlib', comment='FFMpeg'))  # writer=animation.FFMpegWriter,

plt.close('all')

print "Done!\nCheck output in:\n%s" % ws