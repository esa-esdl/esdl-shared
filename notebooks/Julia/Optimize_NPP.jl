"""
CALCULATE LIGHT USE EFFICIENCY
calc_lue(T1, T2, We, Emax)
 
return: light use efficiency
T1    : temperature scalar
T2    : temperature scalar
We    : moisture scalar
Emax  : maximum light use efficiency
"""
function calc_lue(T1, T2, We, Emax)
    lue=T1*T2*We*Emax
    if T1>0.0 && T2>0.0
        return lue
    else
        return zero(lue)
    end
end

"""
CALCULATE NET PRIMARY PRODUCTIVITY (OR GROSS PRIMARY PRODUCTIVITY, CASE
THE AUTOTROPHIC RESPIRATION IS BEING EXPLICITLY CALCULATED)
calc_npp(FPAR, SLRD, AIRT, We, T1, TOPT, Emax, solar_coversion, ToptA) 

returns           : net primary production
FPAR              : fraction of photosynthetically active radiation 
SLRD              : solar radiation
AIRT              : air temperature
We                : vegetation water stress factor
T1                : temperature stress effect
TOPT              : optimum temperature
Emax              : maximum light use efficiency
solar_conversion  : solar radiation units conversion
ToptA             : parameter responsable for the growing part of the
                    Topt response curve.
ToptB             : parameter responsable for the descending part of the
                    Topt response curve.

"""
function calc_npp(fpar,slrd,tair,we,t1,topt,emax,solar_conversion,topta,toptb)
    
    # CONSTANT
    t2c = 1.1919

    tair<-10.0 && (t1=0.0)
    
    if topta == .2 && toptb == .3    
        a       = 0.2;       
        b       = 0.3;      
        t2      = t2c/(1+exp(a*(topt-10-tair)))/(1+exp(b*(-topt-10+tair)));
    else
        c       = tair < topt ? topta : toptb     
        t2p1    = inv((1+exp(-10*c))*(1+exp(-10*c)))
        t2c1    = inv(t2p1)
        t21     = t2c1/(1+exp(c*(topt-10.0-tair)))/(1+exp(c*(-topt-10+tair)));
    end
    lue     = calc_lue(t1, t2, we, emax);
    ipar    = slrd*fpar*solar_conversion;
    return ipar*lue;
end

function npp_loop(vars,topt)
    npp=zero(calc_npp(vars[1,2],vars[1,3],vars[1,1],1.0,1.0,topt,0.5,0.5,0.2,0.3))
    cou=0
    @inbounds for itime=1:size(vars,1)
        r=calc_npp(vars[itime,2],vars[itime,3],vars[itime,1],1.0,1.0,topt,0.5,0.5,0.2,0.3)
        if !isnan(r) 
            npp+=r
            cou+=1
        end
    end
    cou==0 && (cou=1)
    -npp/cou
end

"""
Predict npp for a whole time series based on fpar, slrd, and tair.

predict_npp(xout,outmask,vars,varmask,topt,topmask)

xout     :: Output vector to store results
outmask  :: Output mask to store resulting mask
vars     :: Matrix with time series of tair, fpar and slrd in its columns
varmask  :: Mask of the input variables
topt     :: Optimal temperature for photosynthesis, 0-dimensional
topmask  :: Mask for optimal photosynthesis temperature
"""
function predict_npp(xout,outmask,vars,varmask,topt,topmask)
    const gvars=vars
    if varmask[1,1] == CABLAB.CubeAPI.OCEAN
        xout[:]=NaN
        outmask[:]=CABLAB.CubeAPI.OCEAN
        return nothing
    end
    @inbounds for itime=1:size(vars,1)
        xout[itime]=calc_npp(vars[itime,2],vars[itime,3],vars[itime,1],1.0,1.0,topt[1],0.5,0.5,0.2,0.3)
        outmask[itime]=isnan(xout[itime]) ? CABLAB.CubeAPI.MISSING : CABLAB.CubeAPI.VALID
    end
end

using Optim
"""
Estimate topt so that integrated NPP is maximised.

optimize_npp(xout,outmask,vars,varmask)

xout     :: Output for optimal temperature (0-dim)
outmask  :: Output mask for optimal temperature (0-dim)
vars     :: Matrix with time series of tair, fpar and slrd in its columns
varmask  :: Mask of the input variables
"""
function optimize_npp(xout,outmask,vars,varmask)
    const gvars=vars
    if (varmask[1,1] & CABLAB.CubeAPI.OCEAN) > 0
        xout[1]=NaN
        outmask[1]=CABLAB.CubeAPI.OCEAN
        return nothing
    end
    getNPPtot(topt)=npp_loop(gvars,topt)
    o=optimize(getNPPtot,0.0,40.0)
    if Optim.converged(o) && Optim.minimizer(o)<39.0 && Optim.minimizer(o)>1.0
        xout[1]=Optim.minimizer(o)
        outmask[1]=CABLAB.CubeAPI.VALID
    else
        xout[1]=NaN
        outmask[1]=CABLAB.CubeAPI.MISSING
    end
end

