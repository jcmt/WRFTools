import netCDF4 as n4
import numpy as np


class get:
  def __init__(self, path):
    self.path = path
    self.nc = n4.Dataset(path)

  def close(self):
    self.nc.close()

  def variables(self):
    '''
    Prints the variables in the wrfout file.
    '''

    varname = np.asarray([v for v in self.nc.variables])
    return varname

  def dim(self):
    '''
    Returns an array containing the domain dimensios
    time, levels, latitude, longitude
    '''

    nt = len(self.nc.dimensions['Time'])
    nx = len(self.nc.dimensions['west_east'])
    ny = len(self.nc.dimensions['south_north'])
    nz = len(self.nc.dimensions['bottom_top'])

    return np.array([nt, nz, ny, nx])

  def lat(self):
    '''
    Returns the latitude array
    '''

    LAT = self.nc.variables['XLAT'][0, ]
    return LAT

  def lon(self):
    '''
    Returns the longitude array
    '''

    LON = self.nc.variables['XLONG'][0, ]
    return LON

  def height(self, tstep=None, nlev=':', ny=':', nx=':'):
    '''
    Returns the height of the model levels at a given time

    usage:
        height(tstep)

        tstep is the time instant, if not specified all the written times
        will be used
    '''

    if not tstep:
      Z = self.getvar('PH', tstep=':', nlev=nlev, ny=ny, nx=nx) / 9.81

    else:
      Z = self.getvar('PH', tstep=tstep, nlev=nlev, ny=ny, nx=nx) / 9.81
    return Z

  def time(self, tstep=None):
    '''
    Returns a datetime
    '''

    tstart = self.nc.SIMULATION_START_DATE
    if not tstep:
      t = self.nc.variables['XTIME'][...]
    else:
      t = self.nc.variables['XTIME'][tstep]
    TIME = n4.num2date(t, units='minutes since ' + tstart, calendar='gregorian')
    return TIME

  def getvar(self, var, tstep=':', nlev=':', ny=':', nx=':'):
    '''
    Returns the data from a given variable in the wrfout file

    usage:
       getvar(var, tstep)

       var is a string with the variable name (example 'U10')
       tstep is the time instant, if not specified all the written times will
       be used

    Warning:
       For the variable P, PH and T, their base state will be added
    '''

    if len(self.nc.variables[var].dimensions) == 4:
      SLICE = str(tstep) + ',' + str(nlev) + ',' + str(ny) + ',' + str(nx)
    elif len(self.nc.variables[var].dimensions) == 3:
      SLICE = str(tstep) + ',' + str(ny) + ',' + str(nx)
    elif len(self.nc.variables[var].dimensions) == 2:
      SLICE = str(ny) + ',' + str(nx)

    VAR = eval("self.nc.variables['" + var + "'][" + SLICE + "]")

    if (var == 'P') | (var == 'PH'):
      VAR += eval("self.nc.variables['" + var + "B'][" + SLICE + "]")
    if (var == 'T'):
      VAR += 300

    return VAR

  def SLP(self, tstep=0):
    '''
    Calculation of Sea-level pressure.

    usage:
       WRF_SLP = SLP(tstep)

       tstep is the time instant, if not specified the first written time
       will be used

    From NCL fortran source code wrf_user.f
    '''

    PR = self.getvar('P', tstep=tstep)
    TH = self.getvar('T', tstep=tstep)
    QVAPOR = self.getvar('QVAPOR', tstep=tstep)
    ELEVATION = self.height(tstep=tstep)

    #constants:
    R=287.04
    G=9.81
    GAMMA=0.0065
    TC=273.16+17.05
    PCONST=10000
    c = 2.0/7.0
    #calculate TK:
    TK = TH*np.power(PR*.00001,c)
    #Find least z that is PCONST Pa above the surface
    #Sweep array from bottom to top
    s = np.shape(PR)      #size of the input array
    ss = [s[1],s[2]] # shape of 2d arrays
    WRF_SLP = np.empty(ss,np.float32)
    LEVEL = np.empty(ss,np.int32)
    # Ridiculous MM5 test:
    RIDTEST = np.empty(ss,np.int32)
    PLO = np.empty(ss, np.float32)
    ZLO = np.empty(ss,np.float32)
    TLO = np.empty(ss,np.float32)
    PHI = np.empty(ss,np.float32)
    ZHI = np.empty(ss,np.float32)
    THI = np.empty(ss,np.float32)
    LEVEL[:,:] = -1
    for K in range(s[0]):
      KHI = np.minimum(K+1, s[0]-1)
      LEVNEED = np.logical_and(np.less(LEVEL,0), np.less(PR[K,:,:] , PR[0,:,:] - PCONST))
      LEVEL[LEVNEED]=K
      PLO=np.where(LEVNEED,PR[K,:,:],PLO[:,:])
      TLO=np.where(LEVNEED,TK[K,:,:]*(1.+0.608*QVAPOR[K,:,:]), TLO[:,:])
      ZLO=np.where(LEVNEED,ELEVATION[K,:,:],ZLO[:,:])
      PHI=np.where(LEVNEED,PR[KHI,:,:],PHI[:,:])
      THI=np.where(LEVNEED,TK[KHI,:,:]*(1.+0.608*QVAPOR[KHI,:,:]), THI[:,:])
      ZHI=np.where(LEVNEED,ELEVATION[KHI,:,:],ZHI[:,:])
    P_AT_PCONST = PR[0,:,:]-PCONST
    T_AT_PCONST = THI - (THI-TLO)*np.log(P_AT_PCONST/PHI)*np.log(PLO/PHI)
    Z_AT_PCONST = ZHI - (ZHI-ZLO)*np.log(P_AT_PCONST/PHI)*np.log(PLO/PHI)
    T_SURF = T_AT_PCONST*np.power((PR[0,:,:]/P_AT_PCONST),(GAMMA*R/G))
    T_SEA_LEVEL = T_AT_PCONST + GAMMA*Z_AT_PCONST
    RIDTEST = np.logical_and(T_SURF <= TC, T_SEA_LEVEL >= TC)
    T_SEA_LEVEL = np.where(RIDTEST, TC, TC - .005*(T_SURF -TC)**2)
    Z_HALF_LOWEST=ELEVATION[0,:,:]
    WRF_SLP = 0.01*(PR[0,:,:]*np.exp(2.*G*Z_HALF_LOWEST/(R*(T_SEA_LEVEL+T_SURF))))

    return WRF_SLP

  def ETH(self, tstep=0):
    '''
    Program to calculate equivalent potential temperature.

    usage:
       WRF_ETH = ETH(tsetp)

       tstep is the time instant, if not specified the first written time
       will be used

    From NCL/Fortran source code DEQTHECALC in eqthecalc.f
    '''

    PRESS = self.getvar('P', tstep=tstep)
    TH = self.getvar('T', tstep=tstep)
    QVAPOR = self.getvar('QVAPOR', tstep=tstep)

    c = 2.0/7.0
    EPS = 0.622
    GAMMA = 287.04/1004.0
    GAMMAMD = 0.608 -0.887
    TLCLC1 = 2840.0
    TLCLC2 = 3.5
    TLCLC3 = 4.805
    TLCLC4 = 55.0
    THTECON1 = 3376.0
    THTECON2 = 2.54
    THTECON3 = 0.81
    #calculate Temp. in Kelvin
    PRESS *= 0.01
    TK = TH*np.power(PRESS*.001, c)
    Q = np.maximum(QVAPOR, 1.e-15)
    E = Q*PRESS/(EPS+Q)
    TLCL = TLCLC4+ TLCLC1/(np.log(np.power(TK, TLCLC2)/E)-TLCLC3)
    EXPNT = (THTECON1/TLCL - THTECON2)*Q*(1.0+THTECON3*Q)
    WRF_ETH = TK*np.power(1000.0/PRESS, GAMMA*(1.0+GAMMAMD*Q))*np.exp(EXPNT)

    return WRF_ETH


  def RH(self, tstep=0):
    '''
    Calculation of relative humidity.

    usage:
       WRF_RH = RH(tstep)

       tstep is the time instant, if not specified the first written time
       will be used

    From NCL formula in wrf_user.f
    '''

    PRESS = self.getvar('P', tstep=tstep)
    TH = self.getvar('T', tstep=tstep)
    QVAPOR = self.getvar('QVAPOR', tstep=tstep)

    c = 2.0/7.0
    SVP1 = 0.6112
    SVP2 = 17.67
    SVPT0 = 273.15
    SVP3 = 29.65
    EP_3 = 0.622
    TK = TH * np.power(PRESS * .00001, c)
    ES = 10 * SVP1 * np.exp(SVP2 * (TK-SVPT0) / (TK-SVP3))
    QVS = EP_3 * ES / (0.01*PRESS - (1.-EP_3) * ES)
    WRF_RH = 100.0 * np.maximum(np.minimum(QVAPOR/QVS, 1.0), 0)

    return WRF_RH

  def SHEAR(self, tstep=0, level1=200., level2=850., leveltype='pressure'):
    '''
    Program calculates horizontal wind shear

    usage:
       SHR = SHEAR(tstep, level1, level2)

       tstep is the time instant, if not specified the first written time
       will be used
       level1 is the top level to consider for (200 hPa)
       level2 is the bottom level to consider (850 hPa)

    From NCAR VAPOR python utils
    '''

    if leveltype == 'pressure':
      print(leveltype)
      PR = self.getvar('P', tstep=tstep)
      U = self.getvar('U', tstep=tstep)
      V = self.getvar('V', tstep=tstep)

      PR *= 0.01
      uinterp1 = interp3d(U, PR, level1)
      uinterp2 = interp3d(U, PR, level2)
      vinterp1 = interp3d(V, PR, level1)
      vinterp2 = interp3d(V, PR, level2)
      result = (uinterp1-uinterp2)*(uinterp1-uinterp2)+(vinterp1-vinterp2)*(vinterp1-vinterp2)
      result = np.sqrt(result)

    elif leveltype == 'eta':
      print(leveltype)
      uinterp1 = self.getvar('U', tstep=tstep, nlev=level1, nx=':-1')
      uinterp2 = self.getvar('U', tstep=tstep, nlev=level2, nx=':-1')
      vinterp1 = self.getvar('V', tstep=tstep, nlev=level1, ny=':-1')
      vinterp2 = self.getvar('V', tstep=tstep, nlev=level2, ny=':-1')
      result = (uinterp1-uinterp2)*(uinterp1-uinterp2)+(vinterp1-vinterp2)*(vinterp1-vinterp2)
      result = np.sqrt(result)

    return result

  def TD(self, tstep=0):
    '''
    Calculation of dewpoint temperature based on WRF variables.

    usage:
       WRFTD = TD(tstep)

       tstep is the time instant, if not specified the first written time
       will be used
    '''
    #Let PR = 0.1*(P+PB) (pressure in hPa)
    #and QV = MAX(QVAPOR,0)
    #Where TDC = QV*PR/(0.622+QV)
    # TDC = MAX(TDC,0.001)
    #Formula is (243.5*log(TDC) - 440.8)/(19.48-log(TDC))

    P = self.getvar('P', tstep=tstep)
    QVAPOR = self.getvar('QVAPOR', tstep=tstep)

    QV = np.maximum(QVAPOR, 0.0)
    TDC = 0.01 * QV * P / (0.622 + QV)
    TDC = np.maximum(TDC, 0.001)
    WRF_TD = (243.5 * np.log(TDC) - 440.8)/(19.48 - np.log(TDC))

    return WRF_TD


  def TK(self, tstep=0):
    '''
    Calculation of temperature in degrees kelvin using WRF variables.

    usage:
       TMP = TK(tstep)

       tstep is the time instant, if not specified the first written time
       will be used
    '''

    #Formula is (T+300)*((P+PB)*10**(-5))**c,
    #Where c is 287/(7*287*.5) = 2/7

    P = self.getvar('P', tstep=tstep)
    TH = self.getvar('T', tstep=tstep)

    c = 2.0/7.0
    WRF_TK = TH * np.power(P * .00001, c)

    return WRF_TK

  def BRUNT(self, tstep=0):
    '''
    Calculation of Brunt-Vaisala frequency.

    usage:
       BV = BRUNT(tstep)

       tstep is the time instant, if not specified the first written time
       will be used
    '''

    THETA  =  self.getvar('T', tstep=tstep) * (1 + 0.61 * self.getvar('QVAPOR', tstep=tstep))
    Z = self.height(tstep=tstep)
    nz = self.dim()[1]
    g = 9.81

    BRUNT = np.zeros(shape=self.dim()[1:])
    for i in range(nz-1):
      BRUNT[i, :, :] = (g/THETA[i, :, :]) * ((THETA[i+1, :, :] - THETA[i, :, :]) / (Z[i+1, :, :] - Z[i, :, :]))
 
    return BRUNT

  def RI(self, tstep=0):
    '''
    Calculation of Richardson Number.

    usage:
       ri = RI(tstep)

       tstep is the time instant, if not specified the first written time
       will be used
    '''

    THETA  =  self.getvar('T', tstep=tstep) * (1 + 0.61 * self.getvar('QVAPOR', tstep=tstep))
    Z = self.height(tstep=tstep)
    U  = self.getvar('U', tstep=tstep, nx=':-1')
    V  = self.getvar('V', tstep=tstep, ny=':-1')
    nz = self.dim()[1]
    g = 9.81
    Td = 9.8 / 1000.# The dry adiabatic lapse rate 9.8 K/km

    RI = np.zeros(shape=self.dim()[1:])
    for i in range(nz-1):
      RI[i, :, :] = (g*((THETA[i+1, :, :] - THETA[i, :, :]) + Td * (Z[i+1, :, :] - Z[i, :, :])) * (Z[i+1, :, :] - Z[i, :, :])) / \
           (THETA[i, :, :] * ((U[i+1, :, :] - U[i, :, :])**2 + (V[i+1, :, :] - V[i, :, :])**2))

    return RI

  def pcolor(self, VAR, tstep=None, colorbar=True, level=0, pcolor=False, norm=None, coastcolor='k', **kargs):
    '''
    lat-lon plot on a base map

    usage:
       pcolor(VAR, colormap, colorbar, tstep, level, shading, norm)

       VAR is a wrfout variable (string) or a 2D numpy array
       if VAR is tring a tstep and level must be given to acquire the
       variable. IF NOT the first level and time will be used
       shading can be one of: flat (default), interp (contourf) or None
       (pcolor)
    '''
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    import ticks

    if not tstep:
      return "A time step must be specified..."
    else:
      if isinstance(VAR, str):
        if len(self.nc.variables[VAR].dimensions) == 4:
          VAR = self.getvar(VAR, tstep=tstep, nlev=level, ny=':', nx=':')
        elif len(self.nc.variables[VAR].dimensions) == 3:
          VAR = self.getvar(VAR, tstep=tstep)
        elif len(self.nc.variables[VAR].dimensions) == 2:
          VAR = self.getvar(VAR)

      if self.nc.MAP_PROJ == 1:
        proj = 'lcc'
      elif self.nc.MAP_PROJ == 3:
        proj = 'merc'
      else:
        return('Projection not suported')

      lat_1 = self.nc.TRUELAT1
      lat_2 = self.nc.TRUELAT2
      lon_0 = self.nc.CEN_LON
      lat_0 = self.nc.CEN_LAT
      llcrnrlat = self.lat().min()
      urcrnrlat = self.lat().max()
      llcrnrlon = self.lon().min()
      urcrnrlon = self.lon().max()

      res = 'i'
      if self.nc.DX < 25000:
        res = 'h'

      plt.figure()
      ax = plt.axes()
      m = Basemap(projection=proj, llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, \
                    llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon, lat_1=lat_1, \
                    lat_2=lat_2, lat_0=lat_0, lon_0=lon_0, resolution=res, area_thresh=10000)

      m.drawcoastlines(color=coastcolor, linewidth=2)
      m.drawcountries(color=coastcolor, linewidth=1.5)

      parallels = ticks.loose_label(self.lat().min(),self.lat().max())
      m.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=14)
      meridians = ticks.loose_label(self.lon().min(),self.lon().max())
      m.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=14)

      x, y = m(self.lon(), self.lat())

      if not pcolor:
        if not norm:
          levels = np.linspace(VAR.min(), VAR.max(), 200)
        else:
          levels = np.linspace(norm.min(), norm.max(), 200)
        cs = ax.contourf(x, y, VAR, levels=levels, **kargs)
      else:
        cs = ax.pcolormesh(x, y, VAR, **kargs)

      ax.set_title(self.time(tstep=tstep))

      if colorbar:
        fmt = plt.matplotlib.ticker.FormatStrFormatter("%.1f")
        if not norm:
          clev = np.linspace(np.round(VAR.min()), np.round(VAR.max()), 10, endpoint=True)
        else:
          clev = np.linspace(np.round(norm.min()), np.round(norm.max()), 10, endpoint=True)
        cbar = m.colorbar(cs, location='right', ticks=clev, format=fmt, pad='5%')
        cbar.ax.tick_params(labelsize=12)
      return ax, m

  def CrossPcolor(self, VAR, tstep=1, latitude=None, longitude=None, colorbar=True, \
                  norm=None, ymax=20000, ymin=0, pcolor=True, lev=None, **kargs):
    import matplotlib.pyplot as plt


    '''
    VErtical cross section plot

    usage:
       pcolor(VAR, latitude, longitude, colormap, colorbar, tstep, level, shading, norm,
       lev)

       VAR is a wrfout variable (string) or a 2D numpy array
       if VAR is tring a tstep and level must be given to acquire the
       variable. IF NOT the first level and time will be used
       shading can be one of: flat (default), interp (contourf) or None
       (pcolor)
    '''


    plt.figure()
    ax = plt.axes(axisbg='grey')

    if not latitude and not longitude:
      return('A latitude and longitude range must be chosen...')

    elif not latitude:
      pos_lon = np.argmin(abs(self.lon()[1, :] - longitude))
      pos_lat = slice(0, np.size(self.lat(), axis=0))
      y = self.height(tstep=tstep, nlev=':', ny=pos_lat, nx=pos_lon)
      x = np.tile(self.lat()[:, pos_lon], (self.dim()[1], 1))
      xlabel = 'Latitude ($^\circ$)'

    elif not longitude:
      pos_lon = slice(0, np.size(self.lon(), axis=1))
      pos_lat = np.argmin(abs(self.lat()[:, 1] - latitude))
      y = self.height(tstep=tstep, ny=pos_lat, nx=pos_lon)
      x = np.tile(self.lon()[pos_lat, :], (self.dim()[1], 1))
      xlabel = 'Longitude ($^\circ$)'

    else:
      return('I cant deal with this.. yet!!!')

    if isinstance(VAR, str):
      if len(self.nc.variables[VAR].dimensions) == 4:
        VAR = self.getvar(VAR, tstep=tstep, nlev=':', ny=pos_lat, nx=pos_lon)
      elif len(self.nc.variables[VAR].dimensions) == 3:
        VAR = self.getvar(VAR, tstep=tstep, ny=pos_lat, nx=pos_lon)
      elif len(self.nc.variables[VAR].dimensions) == 2:
        VAR = self.getvar(VAR, ny=pos_lat, nx=pos_lon)

    else:
      VAR = np.squeeze(VAR[:, pos_lat, pos_lon])

    if not pcolor:
      if not lev:
        levels = np.linspace(VAR.min(), VAR.max(), 200)
      else:
        levels = np.linspace(lev[0], lev[1], 100)

      cs = ax.contourf(x[0:VAR.shape[0], :], y[0:VAR.shape[0], :], VAR, norm=norm, levels=levels, **kargs)
    else:
      cs = ax.pcolormesh(x[0:VAR.shape[0], :], y[0:VAR.shape[0], :], VAR, norm=norm, **kargs)

    if colorbar:
      fmt = plt.matplotlib.ticker.FormatStrFormatter("%.1f")
      if not lev:
        clev = np.linspace(np.round(VAR.min()), np.round(VAR.max()), 10)
      else:
        clev = np.linspace(lev[0], lev[1], 10)

      cbar = plt.colorbar(cs, ticks=clev, format=fmt, norm=norm)

    ax.set_title(self.time()[tstep])
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Height (m)')
    ax.set_ylim([ymin, ymax])
    ax.set_xlim(x.min(), x.max())

    return ax

  def stloc(self, latitude, longitude):
    '''
    Returns nearst the grid points to location

    usage:
        nlat, nlon = stloc(latitude, longitude)
    '''

    pos_lat = np.argmin(abs(self.lat()[:, 1] - latitude))
    pos_lon = np.argmin(abs(self.lon()[1, :] - longitude))
    return pos_lat, pos_lon

  def interpvar(self, VAR, ilon, ilat, tstep=None, nlev=None):
    from matplotlib import tri

    '''
    Interpolates variable to latitude longitude (2D)

    usage:
        IVAR = interpvar(VAR, longitude, latitude, tstep, nlev)
    '''

    intp = tri.delaunay.Triangulation(self.lat().flatten(), self.lon().flatten())

    if isinstance(VAR, str):
      if (not nlev) & (tstep):
        data =  self.getvar(VAR, tstep=tstep)[0, ]
      elif (nlev) & (not tstep):
        data =  self.getvar(VAR, nlev=nlev)
      elif (not nlev) & (not tstep):
        pass
    else:
      data = VAR

    IVAR = np.zeros(shape=(data.shape[0], ilat.size, ilon.size))
    for i in range(data.shape[0]):
      IVAR[i, ] = intp.nn_interpolator(data[i, ].flatten())(ilat, ilon)
    return np.squeeze(IVAR)


  def sounding(self, tstep=0, lat=38.28, lon=-28.24):
    import SkewT
    import thermodynamics

 
    '''
    Creates a virtual sounding
    '''  

    ny, nx = self.stloc(lat, lon)

    data={}
    data['pres'] = self.getvar('P', tstep=tstep, nlev=':', ny=ny, nx=nx) * 1e-2
    data['MIXR'] = self.getvar('QVAPOR', tstep=tstep, nlev=':', ny=ny, nx=nx) * 1000.
    data['temp'] = thermodynamics.TempK(self.getvar('T', tstep=tstep, nlev=':', ny=ny, nx=nx), data['pres'] * 1e2) - 273.15
    data['dwpt'] = thermodynamics.MixR2VaporPress(data['MIXR']/1000. ,data['pres'])
    data['sknt'] = self.WSPEED(tstep=tstep, nlev=':', ny=ny, nx=nx) * 1.94384449
    data['drct'] = self.WDIR(tstep=tstep, nlev=':', ny=ny, nx=nx)
    data['hght'] = self.getvar('HGT', tstep=tstep, nlev=':', ny=ny, nx=nx)
    data['RELH'] = self.RH(tstep=tstep)[:, ny, nx]
    data['StationNumber'] = 'Location = ' + str(lat) + ' N ' + str(lon) + ' E'
    data['SoundingDate'] = self.time(tstep=tstep).strftime('%Y/%m/%d %H:%M')
    data['THTV'] =  thermodynamics.ThetaV(data['temp']+273.15, data['pres']*1e2, thermodynamics.VaporPressure(data['dwpt']))
#    data['THTA'] = 
#    data['THTE'] = 

    return SkewT.Sounding(data=data)


  def WSPEED(self, tstep=0, nlev=0, ny=':', nx=':'):

    if nx == ':':
      u = self.getvar('U', tstep=tstep, nlev=nlev, ny=ny, nx=':-1')
    elif ny == ':':
      v = self.getvar('V', tstep=tstep, nlev=nlev, ny=':-1', nx=nx)
    else: 
      u = self.getvar('U', tstep=tstep, nlev=nlev, ny=ny, nx=nx)
      v = self.getvar('V', tstep=tstep, nlev=nlev, ny=ny, nx=nx)

    return np.sqrt(u*u + v*v)


  def WDIR(self, tstep=0, nlev=0, ny=':', nx=':'):

    if nx == ':':
      u = self.getvar('U', tstep=tstep, nlev=nlev, ny=ny, nx=':-1')
    elif ny == ':':
      v = self.getvar('V', tstep=tstep, nlev=nlev, ny=':-1', nx=nx)
    else:
      u = self.getvar('U', tstep=tstep, nlev=nlev, ny=ny, nx=nx)
      v = self.getvar('V', tstep=tstep, nlev=nlev, ny=ny, nx=nx)

    wspeed = self.WSPEED(tstep=tstep, nlev=nlev, ny=ny, nx=nx)

    dir_rad = np.arctan2(u/wspeed, v/wspeed)
    dir_trig = (dir_rad * 180/np.pi) + 180
    dir_cardinal = 90 - dir_trig

    return dir_cardinal


def myround(x, base=5):
  x *= 100
  y = int(base * round(float(x)/base))
  y /= 100.0

  return y


def interp3d(A, PR, val):
  s = np.shape(PR)	#size of the input arrays
  ss = [s[1], s[2]]     #shape of 2d arrays
  interpVal = np.empty(ss, np.float32)
  ratio = np.zeros(ss, np.float32)

  #  the LEVEL value is determine the lowest level where P<=val
  LEVEL = np.empty(ss, np.int32)
  LEVEL[:, :] = -1 #value where PR<=val has not been found
  for K in range(s[0]):
    #LEVNEED is true if this is first time PR<val.
    LEVNEED = np.logical_and(np.less(LEVEL, 0), np.less(PR[K, :, :], val))
    LEVEL[LEVNEED] = K
    ratio[LEVNEED] = (val-PR[K, LEVNEED]) / (PR[K-1, LEVNEED] - PR[K, LEVNEED])
    interpVal[LEVNEED] = ratio[LEVNEED] * A[K, LEVNEED] + (1-ratio[LEVNEED]) * A[K-1, LEVNEED]
    LEVNEED = np.greater(LEVEL, 0)
  # Set unspecified values to value of A at top of data:
  LEVNEED = np.less(LEVEL, 0)
  interpVal[LEVNEED] = A[s[0]-1, LEVNEED]

  return interpVal

def interp_delaunay(x0,y0,v0,x1, y1):
  from matplotlib import tri

  intp = tri.delaunay.Triangulation(x0.flatten(), y0.flatten())
  v1 = intp.nn_interpolator(v0.flatten())(x1, y1)

  return v1

def var_border(v,di=1,dj=1):
  '''
  Border of 2d numpy array
  di,dj is the interval between points along columns and lines
  Corner points are kept even with di and dj not 1
  '''
  j,i=v.shape
  if (di,dj)==(1,1):
    xb=np.arange(2*i+2*j,dtype=v.dtype)
    yb=np.arange(2*i+2*j,dtype=v.dtype)
    xb[0:j] = v[:,0]
    xb[j:j+i] = v[-1,:]
    xb[j+i:j+i+j] = np.flipud(v[:,-1])
    xb[j+i+j:] = np.flipud(v[0,:])
  else:
    # ensure corner points are kept!!
    tmp1 = v[::dj,0]
    tmp2 = v[-1,::di]
    tmp3 = np.flipud(v[:,-1])[::dj]
    tmp4 = np.flipud(v[0,:])[::di]
    xb=np.concatenate((tmp1,tmp2,tmp3,tmp4))
  return xb


def plot_domains(dlist):
  from mpl_toolkits.basemap import Basemap
  import matplotlib.pyplot as plt
  import ticks

  wrf = {}
  for i in range(1, len(dlist)+1):
    wrf['d%02d'%i] = get(dlist[i-1])

  if wrf['d01'].nc.MAP_PROJ == 1:
    proj = 'lcc'
  elif wrf['d01'].nc.MAP_PROJ == 3:
    proj = 'merc'
  else:
    return('Projection not suported')

  lat_1 = wrf['d01'].nc.TRUELAT1
  lat_2 = wrf['d01'].nc.TRUELAT2
  lon_0 = wrf['d01'].nc.CEN_LON
  lat_0 = wrf['d01'].nc.CEN_LAT
  llcrnrlat = wrf['d01'].lat().min() - 5
  urcrnrlat = wrf['d01'].lat().max() + 5
  llcrnrlon = wrf['d01'].lon().min() - 5
  urcrnrlon = wrf['d01'].lon().max() + 5

  plt.figure()
  ax = plt.axes()
  m = Basemap(projection=proj, llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, \
  llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon, lat_1=lat_1, \
  lat_2=lat_2, lat_0=lat_0, lon_0=lon_0, resolution='i')

  #m.drawcoastlines(color='black', linewidth=2)
  #m.drawcountries(linewidth=1.5)
  m.bluemarble()

  parallels = ticks.loose_label(llcrnrlat, urcrnrlat)
  m.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=14)
  meridians = ticks.loose_label(llcrnrlon, urcrnrlon)
  m.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=14)

  for i in range(1,len(dlist)+1):
    xb = var_border(wrf['d%02d'%i].lon())
    yb = var_border(wrf['d%02d'%i].lat())

    x, y = m(xb,yb)
    tx, ty = m(wrf['d%02d'%i].lon()[-1,0], wrf['d%02d'%i].lat()[-1,0]+0.5)
    colors = ['lightblue', 'pink', 'lightgreen', 'lightsalmon', 'silver', 'khaki']

    ax.plot(x,y, lw=2, c=colors[i-1])
    ax.annotate('d%02d'%i, xy=(tx, ty), fontsize=16, color=colors[i-1])


