"""
Matplotlib Animation Example

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
modified by suman bhattacharya

"""
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import MySQLdb as mdb

def findcmfit(tablename, z):
   if z==0.0: strz='0'
   if z==1.0: strz='1'
   if z==2.0: strz='2'
   var='cmfit_'+strz
   ampl,indx=[],[]
   try:
      con = mdb.connect('localhost', 'mysqlu','','testcm');
      cur = con.cursor()
      cur.execute("select amplitude, indx from "+var)
      numrows= int(cur.rowcount)
      for i in range(numrows):
         row= cur.fetchone()
         ampl.append(row[0])
         indx.append(row[1])
         #print "%le %lf" %(row[0], row[1])   
   except mdb.Error, e:
       print "Error %d: %s" % (e.args[0],e.args[1])
       sys.exit(1)
   finally:
      if con:
         con.close()
   return np.array(ampl), np.array(indx)

x0,y0=findcmfit('cmfit_',0.0)
x1,y1=findcmfit('cmfit_',1.0)

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(1, 10), ylim=(0, 0.7))
line, = ax.plot([], [], marker='o', markersize=20,color='r',label="redshift=0")
line1, = ax.plot([], [], marker='o', markersize=20,color='b',label="redshit=1.0")
plt.xlabel(r'amplitude', fontsize=20)
plt.ylabel(r'index', fontsize=20)
ax.legend(loc="upper right", #bbox_to_anchor=[0.7, 1.0],
           ncol=1, shadow=False)
# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    line1.set_data([],[])
    return line,line1

# animation function.  This is called sequentially
def animate(i):
    print "model # = ",i
    line.set_data(x0[i], y0[i])
    line1.set_data(x1[i],y1[i])
    #plt.title("model=%d"%(i))
    return line,line1

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                                frames=37,interval=400, blit=True,repeat=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#anim.save('basic_animation.mp4', fps=30)# """,extra_args=['-vcodec', 'libx264']""")
plt.show()
