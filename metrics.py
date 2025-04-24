def check_density(density, epsilon=1e-5) :
    """
    args :
        density (xarray) : DataArray (t, depth, lat, lon) with density value for each point of the grid.
    return :
        (float) proportion of points not respecting density decreasing constraint
    """
    density=density.where(density!=0)
    diff = density  - density.shift(depth=-1)
    return (diff > epsilon).mean().data # Proportion of points not respecting decreasing density





def temperature_500m_30NS_metric(temperature,   file_mask):
    """
        Metric Extraction Function :
        Average Temperature at 500m depth between 30N and 30S.
        Unit : °C    
    
    
        Input : 
           -  thetao    : xarray.DataArray
           -  file_mask : xarray.Dataset
        Output : 
           - np.float32 or np.float64 depending on recording precision of simulation files 
    
    """
   
    # Taking Temperature At 500m depth and between 30N and 30S.

    t500_30NS=temperature.sel(depth=500,method='nearest').where(abs(temperature.nav_lat)<30,drop=False)

    # Computing Area Weights from Mask over 30N-30S latitude zone and @500m depth
    e1t=file_mask.e1t.squeeze()
    e2t=file_mask.e2t.squeeze()
    tmask=file_mask.tmask.squeeze()
    area_500m_30NS=e1t*e2t*tmask.sel(depth=500,method='nearest').where(abs(temperature.nav_lat)<30,drop=False)

    #Returning Average Temperature at 500m depth as a numpy scalar
    return ((t500_30NS*area_500m_30NS).sum(dim=["nav_lat","nav_lon"])/area_500m_30NS.sum(dim=["nav_lat","nav_lon"]))



def temperature_BWbox_metric(thetao,   file_mask):
    """
        Metric Extraction Function :
        Average Temperature in a U-shaped "Bottom Water" box corresponding to waters below 3000m or beyond 30 degrees of latitude North and South.
        
        ________________________________________________ _Surface
        | . . . . |__________________________| . . . . |_500m
        | . . . . |                          | . . . . |
        | . . . . |        Deep Water        | . . . . |
        | . . . . |__________________________| . . . . |_3000m
        | . . . . . . . . Bottom Water . . . . . . . . |
        |______________________________________________|_Bottom
        S        30S           Eq.          30N        N
        
        Figure : Schematic Representation of the Bottom Water box used in this metric.

        Unit : °C    
    
        Input : 
           -  thetao    : xarray.DataArray
           -  file_mask : xarray.Dataset               
        Output : 
           - np.float32 or np.float64 depending on recording precision of simulation files 
    
    """

    t_BW=thetao.where(1-(thetao.depth<3000)*(abs(thetao.nav_lat)<30))

    # Computing Area Weights from Mask over Box
    e1t=file_mask.e1t.squeeze()
    e2t=file_mask.e2t.squeeze()
    tmask=file_mask.tmask.squeeze()
    area_BW=e1t*e2t*tmask.where(1-(thetao.depth<3000)*(abs(thetao.nav_lat)<30))

    #Returning Average Temperature on Box
    return ((t_BW*area_BW).sum(dim=["nav_lat","nav_lon","depth"])/area_BW.sum(dim=["nav_lat","nav_lon","depth"]))



def temperature_DWbox_metric(thetao,   file_mask):
    """
        Metric Extraction Function :
        Average Temperature in a "Deep Water" box corresponding to waters between 500m and 3000m depth and 30°N and 30°S.
        
        ________________________________________________ _Surface
        |         |__________________________|         |_500m
        |         | . . . . . . . . . . . . .|         |
        |         | . . . .Deep Water . . . .|         |
        |         |__________________________|         |_3000m
        |                 Bottom Water                 |
        |______________________________________________|_Bottom
        S        30S           Eq.          30N        N
        
        Figure : Schematic Representation of the Deep Water box used in this metric.

        Unit : °C    
    
        Input : 
           -  thetao    : xarray.DataArray
           -  file_mask : xarray.Dataset
        Output : 
           - np.float32 or np.float64 depending on recording precision of simulation files 
    
    """
    e1t=file_mask.e1t.squeeze()
    e2t=file_mask.e2t.squeeze()
    tmask=file_mask.tmask.squeeze()
    t_DW=thetao.where(abs((thetao.depth-1750)<1250)*(abs(thetao.nav_lat)<30))

    # Computing Area Weights from Mask over Box
    area_DW=e1t*e2t*tmask.where(abs((thetao.depth-1750)<1250)*(abs(thetao.nav_lat)<30))


    #Returning Average Temperature on Box
    return ((t_DW*area_DW).sum(dim=["nav_lat","nav_lon","depth"])/area_DW.sum(dim=["nav_lat","nav_lon","depth"]))




## Version 1 of ACC metric : Computes the flux assuming rigid lid (as if ssh didn’t change)

def ACC_Drake_metric(uo, file_mask):
    """
        Metric Extraction Function :
        Antarctic Circumpolar Current Transport at the DINO equivalent of the Drake Passage (x=0).
        Unit : Sv


        Input :
           -    uo        : xarray.DataArray
           -    file_mask : xarray.Dataset
        Output : 
           - np.float32 or np.float64 depending on recording precision of simulation files
    
        
    """
                                    
    umask_Drake=file_mask.umask.isel(nav_lon=0).squeeze()
    e3u=file_mask.e3u_0.squeeze()
    e2u=file_mask.e2u.squeeze()

    # Masking the variables onto the Drake Passage
    
    u_masked=uo.isel(nav_lon=0)*umask_Drake
    e3u_masked=e3u.isel(nav_lon=0)*umask_Drake
    e2u_masked=e2u.isel(nav_lon=0)*umask_Drake

    # Multiplying zonal velocity by the sectional areas (e2u*e3u)
    
    ubar=(u_masked*e3u_masked)
    flux=(e2u_masked*ubar).sum(dim=["nav_lat",'depth'])

    #Returning Total Transport across Drake passage as a numpy scalar (unit : Sv)
    return flux/1e6



### Version 2 of ACC metric : Computes the flux assuming varying ssh, thus needing to recompute e3u variable from e3u_0


def ACC_Drake_metric_2(uo, ssh, file_mask):
    """
        Metric Extraction Function :
        Antarctic Circumpolar Current Transport at the DINO equivalent of the Drake Passage (x=0).
        Unit : Sv


        Input :
           -    uo        : xarray.DataArray
           -    ssh       : xarray.DataArray
           -    file_mask : xarray.Dataset
        Output : 
           - np.float32 or np.float64 depending on recording precision of simulation files
    
        
    """


    e3u_0=file_mask.e3u_0
    e2u=file_mask.e2u
    umask_Drake=file_mask.umask.isel(nav_lon=0)

    # Recomputing e3u, using ssh to refactor the original e3u_0 cell heights)
    
    ssh_u = (ssh + ssh.roll(_nav_lon=-1))/2
    bathy_u   = e3u_0.sum(dim="depth")
    ssumask  = umask_Drake[:,0]
    e3u = e3u_0*(1+ssh_u*ssumask/(bathy_u+1-ssumask))
    
    # Masking the variables onto the Drake Passage
    
    u_masked=uo.isel(nav_lon=0)*umask_Drake
    e3u_masked=e3u.isel(nav_lon=0)*umask_Drake
    e2u_masked=e2u.isel(nav_lon=0)*umask_Drake

    # Multiplying zonal velocity by the sectional areas (e2u*e3u)
    
    ubar=(u_masked*e3u_masked).sum(dim='depth')
    flux=(e2u_masked*ubar).sum()

    #Returning Total Transport across Drake passage as a numpy scalar (unit : Sv)
    return flux.data/1e6



### Intensity of the North-Atlantic SubTropical Gyre (NASTG) computed from the local maximum of the Barotropic Stream Function (BSF)

def NASTG_BSF_max(vo,ssh,file_mask):
    """
        Metric Extraction Function :
        Intensity of the North-Atlantic SubTropical Gyre (NASTG) which contains the Gulf-Stream Current.
        Computed using the Barotropic Stream Function (BSF).
        Unit : Sv
        
 
        Input :
           -    vo        : xarray.DataArray
           -    file_mask : xarray.Dataset
        Output : 
           - np.float32 or np.float64 depending on recording precision of simulation files   
    """

    e3v_0=file_mask.e3v_0.squeeze()
    e1v=file_mask.e1v.squeeze()
    vmask=file_mask.vmask.squeeze()
    # Updating e3v from e3v_0 and SSH
    ssh_v = (ssh + ssh.roll(nav_lat=-1))/2
    bathy_v   = e3v_0.sum(dim="depth")
    ssvmask  = vmask.isel(depth=0)
    e3v = (e3v_0*(1+ssh_v*ssvmask/(bathy_v+1-ssvmask)))

    # Integrating Meridional Transport (e3v*e1v*vo) along depth and X from Western boundary eastward 
    # to get Barotropic Stream Function with the "American continent" as reference point (BSF=0)
    V = (vo * e3v).sum(dim='depth') #  == "Barotropic Velocity" * Bathymetry
    BSF=(V * e1v * ssvmask).cumsum(dim='nav_lon') / 1e6 # Integrating from the West, and converting from m³/s to Sv
    # Selecting 0N-40N window where to search for the maximum, which will correspond to the center of rotation for the gyre
    BSF_NASPG=BSF.where(abs(BSF.nav_lat-20)<20)

    return BSF_NASPG.max(dim=["nav_lat","nav_lon"])


    