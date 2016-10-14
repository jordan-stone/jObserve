import os
import subprocess
from datetime import datetime, timedelta

import numpy as np

from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as mpl
from matplotlib import dates, rcParams
from matplotlib.widgets import RadioButtons, CheckButtons
from matplotlib.image import imread

from astroquery.simbad import Simbad
from ephem.cities import lookup
import ephem
mpl.ion()

########################################################
# PLOTTING STUFF                                       #  
########################################################
this_dir=os.path.dirname(__file__)

def compModern(size=20):
    computerModern = FontProperties(size=size,
                                    fname=this_dir+'/cmunrm.ttf')
    return computerModern

rcParams['text.usetex'] = 'true'
rcParams['legend.numpoints']='1'
rcParams['font.family'] = 'serif'
rcParams['axes.linewidth']= '2'
rcParams['xtick.major.pad'] = '6'
rcParams['xtick.major.size'] = '8'
rcParams['xtick.major.width'] = '3'
rcParams['ytick.major.size'] = '8'
rcParams['ytick.major.width'] = '3'
rcParams['lines.markersize'] = '10'
rcParams['lines.linewidth'] = '2'
rcParams['axes.unicode_minus'] = False

########################################################
# QUERY STUFF                                          #  
########################################################
custom_query=Simbad()
custom_query.add_votable_fields('id(HIP|1)','coo(s)','flux(V)')
custom_query.remove_votable_fields('main_id','coordinates')
custom_query.TIMEOUT = 2

defaultDate=str(datetime.utcnow()).split()[0]

def set_date(newDefaultDate):
    global defaultDate
    defaultDate=newDefaultDate

def get_observatory(observatory='Mount Graham',
                    date=None):
    if observatory == 'Mount Graham':
        #lbt=lookup('Mount Graham International Observatory, Arizona')
        lbt=ephem.Observer()
        lbt.lon='-109:53:51.0'
        lbt.lat='32:41:56.9'
    else:
        print 'LOOKING UP OBSERVATORY'
        print observatory
        lbt=lookup(observatory)
    if date:
        lbt.date=date
    else:
        lbt.date=defaultDate
    return lbt

def getstar(name,harvard=False,observatory='Mount Graham',date=defaultDate):
    if harvard:
        custom_query.SIMBAD_URL=u'http://simbad.cfa.harvard.edu/simbad/sim-script'
    table=custom_query.query_object(name,)
    Name=table['ID_HIP_1'][0]
    RA_s=table['RA_s'][0].replace(' ',':')
    DEC_s=table['DEC_s'][0].replace(' ',':')
    magv=table['FLUX_V'][0]
    readdb_str='%s,f,%s,%s,%f,2000' % (Name, RA_s, DEC_s, magv)
    star=ephem.readdb(readdb_str)
    lbt=get_observatory(observatory,date=date)
    star.compute(lbt)
    return star

def makestar(RA_s,DEC_s,Name='star',magv=999,observatory='Mount Graham'):
    readdb_str='%s,f,%s,%s,%f,2000' % (Name, RA_s, DEC_s, magv)
    star=ephem.readdb(readdb_str)
    lbt=get_observatory(observatory=observatory)
    star.compute(lbt)
    return star

def plot_secz(star,date=None,ax=None,overplot=False,observatory='Mount Graham'):
    lbt=get_observatory(observatory=observatory)
    if date:
        lbt.date=date
    else:
        lbt.date=defaultDate
    sun=ephem.Sun()
    sun.compute(lbt)
    start=sun.set_time
    end=lbt.next_rising(sun)
    times=[]
    seczs=[]
    for tt in np.linspace(start,end,100):
        lbt.date=tt
        star.compute(lbt)
        z=(np.pi/2.)-star.alt
        secz=(1./np.cos(z))
        if secz < 3.5 and secz>0:
            seczs.append(secz)
            times.append(ephem.date(tt).datetime())
    if ax is None:
        f=mpl.figure()
        ax=f.add_subplot(111)
    elif overplot:
        ax=ax.twinx()

    ax.plot(times,seczs,label=star.name)
    ax.set_xlim(start.datetime(),end.datetime())
    ax.set_ylim(0.5,2.2)
    ax.xaxis.set_major_locator(dates.HourLocator())
    ax.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
    ax.figure.autofmt_xdate(bottom=0.2,rotation=30,ha='right')

    ax.legend(loc='best')
    ax.set_xlabel('UTC Time',fontproperties=compModern(18))
    ax.set_ylabel('Airmass',fontproperties=compModern(18))
    ax.set_title(str(lbt.date).split()[0])
    mpl.show()
    return ax

def get_secz(name,date=None,ax=None,harvard=False,observatory='Mount Graham'):
    '''rectilinear'''
    star=getstar(name,harvard=harvard,observatory=observatory)
    if date:
        pass
    else:
        date=defaultDate
    ax=plot_secz(star,date=date,ax=ax,observatory=observatory)
    return ax
    
def plot_parallactic_angle(star,date=None,ax=None,overplot=False,observatory='Mount Graham'):
    lbt=get_observatory(observatory=observatory)
    if date:
        lbt.date=date
    else:
        lbt.date=defaultDate
    sun=ephem.Sun()
    sun.compute(lbt)
    start=sun.set_time
    end=lbt.next_rising(sun)
    times=[]
    qs=[]
    for tt in np.linspace(start,end,100):
        lbt.date=tt
        star.compute(lbt)
        z=(np.pi/2.)-star.alt
        secz=(1./np.cos(z))
        if secz < 3.5 and secz>0:
            qs.append(180./np.pi*star.parallactic_angle())
            times.append(ephem.date(tt).datetime())
    if ax is None:
        f=mpl.figure()
        ax=f.add_subplot(111)
    elif overplot:
        ax=ax.twinx()

    ax.plot(times,qs,linestyle='-',label=star.name)
    ax.set_xlim(start.datetime(),end.datetime())
    ax.xaxis.set_major_locator(dates.HourLocator())
    ax.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
    ax.figure.autofmt_xdate()
    ax.set_ylim(-180,180)
    ax.legend(loc='best')
    ax.set_xlabel('UTC Time',fontproperties=compModern(18))
    ax.set_ylabel('Parallactic Angle [$^{\circ}$]',fontproperties=compModern(18))
    mpl.show()
    return ax

def get_parallactic_angle(name,date=None,ax=None,harvard=False,observatory='Mount Graham'):
    '''rectilinear'''
    star=getstar(name,harvard=harvard,observatory=observatory)
    if date:
        pass
    else:
        date=defaultDate
    ax=plot_parallactic_angle(star,date=date,ax=ax,observatory=observatory)
    return ax

def fix_az(az):
    if az > 0 and az<=180:
        return az
    if az > 180:
        return az-360

def make_az_data(star, date=None, observatory='Mount Graham'):
    lbt = get_observatory(observatory=observatory)
    print lbt
    if date:
        lbt.date=date
    else:
        lbt.date=defaultDate
    sun=ephem.Sun()
    sun.compute(lbt)

    #start on the hour before sunset
    set_time=list(sun.set_time.tuple())
    set_time[-2]=0
    set_time[-1]=0
    start=ephem.Date(tuple(set_time))
    end=lbt.next_rising(sun)
    night_minutes=24*60.*(end-start)
    times=[]
    azs=[]
    alts=[]
    for tt in np.linspace(start,end,night_minutes):
        lbt.date=tt
        star.compute(lbt)
        z=(np.pi/2.)-star.alt
        secz=(1./np.cos(z))
        if 180./np.pi*((np.pi/2.)-star.alt) < 90:
            azs.append(star.az)
            alts.append(180./np.pi*((np.pi/2.)-star.alt))
            times.append(ephem.date(tt).datetime())
    return times, alts, azs

def plot_az(star,date=None,ax=None, observatory='Mount Graham'):
    times, alts, azs = make_az_data(star, date=date, observatory=observatory)
    if ax is None:
        f=mpl.figure()
        ax=f.add_subplot(111,projection='polar')
    ax.plot(azs,alts,linestyle='-',label=star.name,markevery=60,marker='o')
    ax.set_theta_offset(np.pi/2.)
    tind = (alts == np.min(alts)).nonzero()[0][0]
    ax.text(azs[tind], alts[tind], str(times[tind]).split()[-1][:5],color=ax.lines[-1].get_color())
    ax.legend()
    return ax 

def get_az(name,date=None,ax=None,harvard=False, observatory='Mount Graham'):
    '''polar'''#this is used
    star=getstar(name,harvard=harvard,observatory=observatory)
    if date:
        pass
    else:
        date=defaultDate
    ax=plot_az(star,date=date,ax=ax,observatory=observatory)
    return ax

def get_summary(name,date=None,ax=None,harvard=False,observatory='Mount Graham'):
    if date:
        pass
    else:
        date=defaultDate
    star=getstar(name,harvard=harvard,observatory=observatory)
    ax=plot_secz(star,date,ax=ax,observatory=observatory)
    ax=plot_parallactic_angle(star,date,ax=ax,overplot=True,observatory=observatory)


class jObserve:
    def __init__(self,pfunc=get_secz,targs=[],date=defaultDate,harvard=False,observatory='Mount Graham'):
        self.func=pfunc
        self.f=mpl.figure(figsize=(10,10))
        self.a=self.f.add_axes([0.30,0.30,0.65,0.65],projection=pfunc.func_doc)
        self.raxes=[]
        self.raxes.append(self.f.add_axes([0.025, 0.90, 0.10, 0.05]))
        self.last_button_bottom=0.90
        self.radios=[]
        self.nowline=None
        self.radios.append(CheckButtons(self.raxes[-1], 
                                  ('now',), 
                                  actives=[True]))
        self.radios[-1].labels[-1].set_color('k')
        self.radios[-1].on_clicked(self.now)
        self.line_names={}
        self.harvard=harvard
        self.observatory=observatory
        self.date=date
        for t in targs:
            self.add_star(t)
        mpl.draw()

    def now(self,label):
        linestyles={'-':'--','--':'-'}
        if self.nowline is None:
            self.a.plot([datetime.utcnow()]*2,[-200,200],color='gray',linewidth=3)
            self.nowline=self.a.lines[-1]
            self.a.plot([datetime.utcnow()+timedelta(0,3600)]*2,[-200,200],color='gray',linewidth=3)
            self.hourline=self.a.lines[-1]
        else:
            self.nowline.set_xdata([datetime.utcnow()]*len(self.nowline.get_xdata()))
            self.nowline.set_linestyle(linestyles[self.nowline.get_linestyle()])
            self.hourline.set_xdata([datetime.utcnow()+timedelta(0,3600)] *
                                    len(self.hourline.get_xdata()))
            self.hourline.set_linestyle(linestyles[self.hourline.get_linestyle()])
        
    def showfunc(self,label):
        self.a.lines[self.line_names[label]].set_visible(not self.a.lines[self.line_names[label]].get_visible())
        try:
            self.a.texts[self.line_names[label]].set_visible(not self.a.texts[self.line_names[label]].get_visible())
        except:
            pass
        mpl.draw()

    def add_star(self,name):
        print self.observatory
        self.a=self.func(name,ax=self.a,harvard=self.harvard,observatory=self.observatory,date=self.date)
        self.a.legend_.set_visible(False)

        if len(self.raxes)==0:
            self.raxes.append(self.f.add_axes([0.025, 0.90, 0.17, 0.05]))
            self.last_button_bottom=0.90
        else:
            self.raxes.append(self.f.add_axes([0.025,self.last_button_bottom-0.05,0.17,0.05]))
            self.last_button_bottom=self.last_button_bottom-0.05

        self.line_names[name]=len(self.a.lines)-1
        self.radios.append(CheckButtons(self.raxes[-1], 
                                  (name,), 
                                  actives=[False]))
        self.radios[-1].labels[-1].set_color(self.a.lines[-1].get_color())
        self.radios[-1].on_clicked(self.showfunc)

        mpl.draw()

class night_summary:
    def __init__(self,targs=[],date=defaultDate,harvard=False,observatory='Mount Graham'):
        self.seczfig=jObserve(pfunc=get_secz,targs=targs,date=date,harvard=harvard,observatory=observatory)
        self.qfig=jObserve(pfunc=get_parallactic_angle,targs=targs,date=date,harvard=harvard,observatory=observatory)
        self.azfig=jObserve(pfunc=get_az,targs=targs,date=date,harvard=harvard,observatory=observatory)

    def add_star(self,name):
        self.seczfig.add_star(name)
        self.qfig.add_star(name)
        self.azfig.add_star(name)

def get_finderchart(*star):
    getband = 'DSS2 Red'
    if len(star) == 1:
        obj=getstar(star[0],date='2000-01-01')
        radeg=obj.ra*180./np.pi
        decdeg=obj.dec*180./np.pi
    if len(star) == 2:
        radeg=360.*ten([float(ii) for ii in star[0].split(':')])/24.
        decdeg=ten([float(ii) for ii in star[1].split(':')])

    myurl = ('http://irsa.ipac.caltech.edu/applications/finderchart/servlet/api?locstr=' + 
              str(radeg)+ 
              '+'
              +str(decdeg)+
              '&subsetsize=18.0&survey=dss&survey=2mass&orientation=left')
# It does not make JPEG images if you ask for marker=true.

# Use wget to retrieve the resulting XML on STDOUT
    wcom = '/usr/bin/wget -q "' + myurl + '" -O -'  # create the command
    child = os.popen(wcom) # spawn the process
    data = child.read()   # capture the stdout
    err = child.close()
#
# Find multiple occurrences of "<jpgurl>" in data
    jcount =  data.count("<jpgurl>")
    remain = data           # copy the data array
    for x in range(0, jcount):  # loop over occurrences
#print x, " ", remain.count("<jpgurl>"), " ", remain.index("<jpgurl>")
        imageband = remain[remain.index("<band>")+6:remain.index("</band>")]
        imagedate = remain[remain.index("<obsdate>")+9:remain.index("</obsdate>")]
        print x, " ", imageband, " ", imagedate
        imageurl = remain[remain.index("<jpgurl>")+8:remain.index("</jpgurl>")]
# reject url if there is no date? (IRSA issue?)
        if len(imagedate) > 2:
# Convert "&amp;" to "&"
            cleanurl = imageurl.replace("&amp;","&")
            cleanband = imageband
            cleandate = imagedate
            if cleanband == getband:
                break #stop with this one -- the preferred band

        newremain = remain[remain.index("</jpgurl>")+9:] # save what is left
        remain = newremain   # remaining data for next pass
    child = subprocess.call(["wget", "-O", "finderchart.png", cleanurl])
    arr = imread('finderchart.png')
    mpl.imshow(np.sum(arr,axis=-1),cmap='gray')
    return cleanurl

