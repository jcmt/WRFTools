import numpy as np
import datetime
import netCDF4 as n4


def tracer_emi(emi_lat=38.252, emi_lon=-28.53, emi_start=None, emi_stop=None,  domain=1, namelist='namelist.input', wrfinput='wrfinput_d01', emi=1e-6):
  e_we = int(search_namelist(filein=namelist, key='e_we', domain=domain))
  e_sn = int(search_namelist(filein=namelist, key='e_sn', domain=domain))

  start = {}
  end = {}
  for i in ['year', 'month', 'day', 'hour']:
    start[i] = int(search_namelist(filein=namelist, key='start_' + i, domain=domain))
    end[i] = int(search_namelist(filein=namelist, key='end_' + i, domain=domain))

  nt = int((datetime.datetime(end['year'], end['month'], end['day'], end['hour']) - \
       datetime.datetime(start['year'], start['month'], start['day'], start['hour'])).total_seconds() / 3600.)

  start_date = datetime.datetime(start['year'], start['month'], start['day'], start['hour'])
  time = np.array([start_date + datetime.timedelta(hours=x) for x in range(nt)])
  timestr = [time[i].strftime('%Y-%m-%d_%H:%M:%S') for i in range(nt)]

  if emi_start == None:
    emi_start = time[0]
  if emi_stop == None:
    emi_stop = time[-1]

  ini = np.where(time == emi_start)[0]
  fini = np.where(time == emi_stop)[0]

  tracer = np.zeros(shape=[nt, 1, e_sn, e_we])

  lat, lon = get_latlon(wrfinput)

  pos_lat = np.argmin(abs(lat[:, 1] - emi_lat))
  pos_lon = np.argmin(abs(lon[1, :] - emi_lon))

  tracer[ini:fini, :, pos_lat, pos_lon] = emi

  ### SAVE TO NETCDF
  root_grp = n4.Dataset('wrfchemi_d01_' + timestr[0], 'w', format='NETCDF3_64BIT')
  root_grp.DESCRIPTION = 'Tracer emission for WRF-Chem -- Created by J.C.Teixeira'

  root_grp.createDimension('Time', 0)
  root_grp.createDimension('DateStrLen', len(timestr[0]))
  root_grp.createDimension('west_east', e_we)
  root_grp.createDimension('south_north', e_sn)
  root_grp.createDimension('emissions_zdim', tracer.shape[1])

  Times = root_grp.createVariable('Times', 'S1', ('Time', 'DateStrLen'))
  tracer_1 = root_grp.createVariable('E_CO', 'f4', ('Time', 'emissions_zdim', 'south_north', 'west_east'))
  tracer_1.FieldType = np.array(int(104.))
  tracer_1.MemoryOrder = "XYZ"
  tracer_1.description = "TRACER EMISSION RATE" ;
  tracer_1.units = "ug/m3 m/s" ;
  tracer_1.stagger = "" ;
  tracer_1.coordinates = "XLONG XLAT" ;

  for i in range(nt):
    for j in range(len(timestr[0])):
      Times[i,j] = timestr[i][j]
    tracer_1[i, ] = tracer[i, ]

  ### COPY ATTRIBUTES
  nc = n4.Dataset(wrfinput)
  for i in nc.ncattrs():
    setattr(root_grp, i, getattr(nc,i))

  root_grp.close()

  print('frames per outfile = ' + str(nt) + '\n' \
        + 'emission location (nlat, nlon) = '(pos_lat, pos_lon) + '\n' \
        + 'interval  = 1 hour')

def search_namelist(filein='namelist.input', key=None, domain=1):
  import re

  regexp = re.compile(re.escape(key) + r'.*?([0-9.-]+),.*?([0-9.-]+),')
  with open(filein) as f:
    for line in f:
      match = regexp.match(line)
      if match:
        return match.group(domain)

def get_latlon(wrfinput):
  import netCDF4 as n4

  nc = n4.Dataset(wrfinput)
  lon = nc.variables['XLONG'][0, ]
  lat = nc.variables['XLAT'][0, ]

  return lat, lon
