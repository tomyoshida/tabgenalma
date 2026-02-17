import casatasks as ct
import casatools as ctools
import numpy as np
from astropy.time import Time

from decimal import Decimal, ROUND_HALF_UP

import os


def generate_alma_obs_table(data, header = True, footer = True):
    
    if header:
        latex_code = [
            r"\begin{deluxetable*}{cccccc}",
            r"\tablecaption{Summary of ALMA Observations. \label{tab:obs}}",
            r"\tablewidth{0pt}",
            r"\tablehead{",
            r"\colhead{Project} & \colhead{P.I.} & \colhead{Date} & \colhead{On-source} & \colhead{Baselines} & \colhead{Frequencies}\\",
            r"\colhead{Code} & \colhead{} & \colhead{} & \colhead{time (min)} & \colhead{(m)} & \colhead{(GHz)\tablenotemark{a}} ",
            r"}",
            r"\startdata"
        ]

    for band_name, eb_list in data.items():
        latex_code.append(r"\hline")
        latex_code.append(rf"\multicolumn{{6}}{{c}}{{{band_name}}} \\")
        latex_code.append(r"\hline")

        last_project = None
        last_pi = None

        for eb in eb_list:
           
            current_project = eb.get('project_code', '')
            current_pi = eb.get('pi', '')
            
            
            disp_project = current_project
            disp_pi = current_pi

            if current_project == last_project:
                disp_project = "~" 
                if current_pi == last_pi:
                    disp_pi = "~"
            
            
            row = (
                f"{disp_project} & "
                f"{disp_pi} & "
                f"{eb.get('date', '')} & "
                f"{eb.get('on_source_time', '')} & "
                f"{eb.get('baselines', '')} & "
                f"{eb.get('frequencies', '')} \\\\"
            )
            latex_code.append(row)
            
            
            last_project = current_project
            last_pi = current_pi
            
    if footer:

        latex_code.extend([
            r"\enddata",
            r"\tablenotetext{a}{Mean frequency of spectral windows.}",
            r"\end{deluxetable*}"
        ])

    return "\n".join(latex_code)

def get_date_tint_spws( summary, obs_id = 0, field_id = 0, array_id = 0 ):
    
    
    n_scans = [k for k in summary[f'observationID={obs_id}'][f'arrayID={array_id}'].keys() if k.startswith('scan=')]
    
    dT = 0.0
    
    spws_all = []
    times_arr = np.array([])
    
    for n_scan in n_scans:
        data = summary[f'observationID={obs_id}'][f'arrayID={array_id}'][f'{n_scan}'][f'fieldID={field_id}']
       
        times = [ data[k]['time'] for k in data.keys() if k.isdigit()]
        
        times_arr = np.append(times_arr, times)
        
        dT += times[-1] - times[0]
        
        spws = data['data description IDs']
        
        spws_all.append(  spws )
        
        
    tstart =  np.min(times_arr)
    
    if np.all(spws_all == spws_all[0]):
        return MJS_to_Date(tstart), dT / 60.0, np.array(spws_all)[0]
    else:
        raise ValueError('Consistency check failed: Multiple SPW patterns detected within a scan.')

    
def MJS_to_Date( MJS ):


    mjd_value = MJS / 3600 / 24

    t = Time(mjd_value, format='mjd')
    
    formatted_date = t.strftime('%Y %b %d')
    
    return formatted_date


def get_baselines( uvw ):
    
    u,v,_ = uvw
    
    uvd = np.sqrt( u**2 + v**2 )
    
    return np.min(uvd), np.max(uvd)

def get_freqs( summary, msmd, spws, obs_id = 0 ):
    
    
    freqs = np.array([])
    for spw in spws:
        freqs = np.append(freqs, msmd.meanfreq( spw )/1e9 )
        
    freq_txt = ''
    i = 0
    for freq in freqs:
        
        if i == 0:
            freq_txt += f'{dround(freq, 0.1)}'
        else:
            freq_txt += f', {dround(freq, 0.1)}'
            
        i+=1
        
    
    return freq_txt
    
def make_one_line(  msmd, summary, obs_id, uvw, observers, project_code = None ):
    
    
    umin, umax = get_baselines( uvw )
    date, dT, spws = get_date_tint_spws( summary, obs_id )
    freqs = get_freqs( summary, msmd, spws, obs_id )
    
    observer = observers[obs_id]
    
    res = {
            "project_code": project_code,
            "pi": observer,
            "date":  date,
            "on_source_time": dround(dT, 0.1),
            "baselines": f"{dround(umin, 0)} -- {dround(umax, 0)}",
            "frequencies": freqs
         }
    
    return res
    
def dround( v, d = 0 ):
    num = Decimal(f"{v}")
    rounded = num.quantize(Decimal(f'{d}'), rounding=ROUND_HALF_UP)
    
    return str(rounded)
    
def make_lines(  vis, project_code = None ):
    
    
    tb = ctools.table()
    msmd = ctools.msmetadata()
    
    msmd.open(vis)
    observers = msmd.observers()
    summary = msmd.summary()

    
    tb.open(vis)
    obs_ids = np.unique(tb.getcol('OBSERVATION_ID'))


    res_list = []
    
    for obs_id in obs_ids:
        
        subtb = tb.query(f'OBSERVATION_ID=={obs_id}')
        
        uvw = subtb.getcol('UVW')
        
        
        res_list.append( make_one_line(  msmd, summary, obs_id, uvw, observers, project_code = None ) )
        
        
    tb.close()
    msmd.close()
        
        
    return res_list


def generate_table( vis_list, bands ):
    
    for i, vis in enumerate(vis_list):
        
        print( f'Processing {vis}...' )
        
        data_dict = {  f'Band {bands[i]}' : make_lines( vis, project_code = None ) }
    
    res = generate_alma_obs_table(data_dict)
    print( res )
        
    return res
