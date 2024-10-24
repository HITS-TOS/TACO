import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.interpolate import interp1d


# This code uses the Coefficient of Variation method proposed by Bell et al. (2019) to detect potential
# solar-like oscillators without fitting the background of the power density spectrum (pds)

# Some abbreviations: 
#   - ind : independent (reference to independent spectrum made out of 43 independent bins without overlap)
#   - os  : oversampled (reference to oversampled spectrum made out of 2000 overlapping bins)


################################################################################
#              FUNCTIONS TO COMPUTE CV & CLASSIFY PULSATIONS                   #
################################################################################


def delta_nu(numax):
    # Values taken From Yu et al. (2018)
    A = numax ** 0.764
    dnu = 0.267 * A
    
    return dnu

def cvmax_solarlike(nu):
    # Determined by Bell et al. (2019) using numax of APOKASC 2 stars.
    cv_max = 2.69 * (nu ** 0.154)
    
    return cv_max

def bins_ind(f):
    # f has to be a pandas DataFrame with columns "frequency" and "power"
    f = f[f['frequency'] >= 1].reset_index().drop(columns = ['index'])
    
    central, edges, widths = [0.0001], [], []
    cnte = 1.029 

    n = 0
    while n < len(f['frequency']):
        if n == 0:
            freq = np.array(f['frequency'])[n]
            dnu_bin = delta_nu(freq)

            width05 = dnu_bin / 2
            
            binned_f = f[(f['frequency'] <= (freq + width05) ** cnte ) & (f['frequency'] >= (freq - width05) ** cnte)]

            edges.append(freq)
            central.append(max(np.array(binned_f['frequency'])))
            
            n += 1
            
        
        elif n == 1:
            prev_center = np.array(central)[-1]
            dnu_bin = delta_nu(prev_center)
            up_edge = (prev_center + (dnu_bin / 2)) ** cnte

            edges.append(up_edge)
            widths.append(up_edge - prev_center)

            n += 1


        else:
            prev_edge = np.array(edges)[-1]
            prev_width = np.array(widths)[-1]

            if prev_edge + prev_width <= 283.4:
                pot_center = prev_edge 
                
                fbinned = f[f['frequency'] >=  pot_center].reset_index().drop(columns = ['index'])
                freqs = np.array(fbinned['frequency'])
                    
                m = 0
                while m < len(freqs):
                    freq = freqs[m]
                    logfreq = np.log10(freq)

                    dnu_bin = delta_nu(freq)
                    width05 = dnu_bin / 2

                    low = freq - width05
                    up = freq + width05

                    alpha = cnte * (np.log10(up) - np.log10(low))
                    
                    new_low = 10 ** (logfreq - (alpha / 2))
                    new_up = 10 ** (logfreq + (alpha / 2))

                    
                    if (new_low >= prev_edge) and (new_low < freq):

                        central.append(freq)
                        edges.append(new_up)
                        widths.append(new_up - freq)
                        
                        current_f = f[f['frequency'] >= new_up]
                        if len(current_f) != 0:
                            n = f[f['frequency'] >= new_up].index[0]
                        
                        else:
                            n += len(f)

                        m = len(freqs)
                    
                    elif new_low > freq:
                        break

        
                    else:
                        m += 1

            else:
                n += len(f)  

    return np.array(central), np.array(edges)

def cv_bins_ind(f, edges):
    cvs, bin_size, centers = [], [], []
    freq_res = np.mean(np.diff(np.array(f['frequency'])))
    
    n = 1
    while n < len(edges):
        
        if n < len(edges) - 1:
            left = edges[n - 1]
            right = edges[n]
            #print(left, right)
            binned_f = f[(f['frequency'] <= right) & (f['frequency'] >= left)].reset_index().drop(columns = ['index'])
            
        else:
            left = edges[n - 1]
            
            binned_f = f[f['frequency'] >= left].reset_index().drop(columns = ['index'])

        
        if len(binned_f) == 0:
            if n < len(edges) - 1:
                binned_f = pd.DataFrame(data = {'frequency': np.arange(left, right, freq_res),
                                                'power': np.array([1]*len(np.arange(left, right, freq_res)))})
            else:
                binned_f = pd.DataFrame(data = {'frequency': np.arange(left, 283.4, freq_res),
                                                'power': np.array([1]*len(np.arange(left, 283.4, freq_res)))})
        
        
        binned_power = np.array(binned_f['power'])
        
        cvs.append(np.std(binned_power) / np.mean(binned_power))
        bin_size.append(len(binned_power))

        center = (max(binned_f['frequency']) + min(binned_f['frequency'])) / 2
        centers.append(center)

        n += 1

    # Determine FAP threshold for each bin according to its size. Requires interpolation!
    # From longer table
    fap_binsize = [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]
    fap_01pct = [1.8451482441347544, 1.950064568761642, 1.8427104174792766, 1.6493357526052286,
                 1.4779871899958588, 1.3410580005395578, 1.2295718191143508, 1.157506593963815,
                 1.107075066468454, 1.0723238218380722, 1.0508277867416875, 1.035366347090339,
                 1.0243920088967091]

    # If FAPs.dat table available:
    #table = pd.DataFrame(pd.read_csv('./FAPs.dat', sep = '\s', names = ['nbin', '50', '80', '90', '95', '99', '99.9'], comment='#'))
    #fap_binsize = np.array(table['nbin'])
    #fap_01pct = np.array(table['99.9'])


    interp_funct = interp1d(x = fap_binsize, y = fap_01pct, fill_value = 'extrapolate')
    interp_01pct = interp_funct(bin_size)
    
    return np.array(cvs), np.array(interp_01pct)

def find_peaks(cvs, fap):
    # Identify potential peaks using CV values and FAP threshold  
    # Returns a list with index position of each potential peak
    peaks = []

    for i, j in enumerate(cvs):
        fap_val = fap[i]

        if j > fap_val:
            # Possible peak
            if i < len(cvs) - 1:
                peaks.append(i) 

    return peaks

def cv_bins_os(f):
    f = f[f['frequency'] >= 1].reset_index().drop(columns = ['index'])

    freq = np.array(f['frequency'])
    power = np.array(f['power'])
    
    # In log-scale
    log_freq = np.log10(freq)
    log_power = np.log10(power)
    f['log_freq'] = log_freq
    f['log_power'] = log_power
    
    
    # In Bell et al. (2019) they use 2000 overlapping bins in log-freq space:
    # Divides the range of log-freq in 2000, and calculate DeltaNu for each of these points
    nbins = 2000
    freq_step = (max(log_freq) -  min(log_freq)) / nbins
    
    # Generate frequencies for each bin
    bin_freqs, bin_size = [], []
    
    f0 = min(log_freq)
    for i in range(nbins):
        bin_freqs.append(f0 + (i * freq_step))

    
    # Determine CV values for each of the 2000 bins
    cvs = []
    n = 0
    while n < len(bin_freqs):
        cfreq = 10 ** bin_freqs[n]
        dnu_cfreq = delta_nu(cfreq)

        width05 = dnu_cfreq / 2
   
        binned_f = f[(f['frequency'] <= cfreq + width05) & (f['frequency'] > cfreq - width05)]
        
        binned_power = np.array(binned_f['power'])
        bin_size.append(len(binned_power))
        cvs.append(np.std(binned_power) / np.mean(binned_power))
        
        n += 1

    bin_freqs = 10 ** np.array(bin_freqs)

    # Determine FAP threshold for each bin according to its size. Requires interpolation!
    # From longer table
    fap_binsize = [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]
    fap_01pct = [1.8451482441347544, 1.950064568761642, 1.8427104174792766, 1.6493357526052286,
                 1.4779871899958588, 1.3410580005395578, 1.2295718191143508, 1.157506593963815,
                 1.107075066468454, 1.0723238218380722, 1.0508277867416875, 1.035366347090339,
                 1.0243920088967091]

    # If FAPs.dat table available:
    #table = pd.DataFrame(pd.read_csv('./FAPs.dat', sep = '\s', names = ['nbin', '50', '80', '90', '95', '99', '99.9'], comment='#'))
    #fap_binsize = np.array(table['nbin'])
    #fap_01pct = np.array(table['99.9'])


    interp_funct = interp1d(x = fap_binsize, y = fap_01pct, fill_value = 'extrapolate')
    interp_01pct = interp_funct(bin_size)
    
    return np.array(bin_freqs), np.array(cvs), np.array(interp_01pct)

def filter_peaks(peaks, peaks_freq, peaks_cvs):
    # peaks is an array of index where cvs is higher than FAP
    # According to Bell et al. (2019), solar-like oscillators have CV values smaller than
    #           cvmax ~ 2.69 nu ** 0.154

    non_solar, solar = [], []
    
    for i in peaks:
        cv_max = cvmax_solarlike(peaks_freq[i])

        if peaks_cvs[i] > cv_max:
            non_solar.append(i)

        else:
            solar.append(i)

    return non_solar, solar

def interpolate_spikes(cvs, cvs_freqs, pds):
    # Makes use of the oversampled spectrum

    ################################################################################
    # A. INTERPOLATION OF CV OVERSAMPLED SPECTRUM TO THEN ESTIMATE NUMAX
    ################################################################################

    cv_max = cvmax_solarlike(cvs_freqs)
    non_solar = []
    masked_cvs = np.copy(cvs)


    # Find potential non-solar-like peaks in oversampled spectrum    
    for i, j in enumerate(cvs):
        if j > cv_max[i]:
            non_solar.append(True)
 
        else:
            non_solar.append(False)

    masked_cvs[non_solar] = np.nan
    
    # Remove these peaks from oversampled spectrum (useful to then estimate numax)
    # Linear interpolation:
    non_nan = ~np.isnan(masked_cvs)
    xp = non_nan.ravel().nonzero()[0]
    fp = masked_cvs[~np.isnan(masked_cvs)]
    x = np.isnan(masked_cvs).ravel().nonzero()[0]

    interp_cvs = masked_cvs
    interp_cvs[np.isnan(masked_cvs)] = np.interp(x, xp, fp)


    ################################################################################
    # B. INTERPOLATION OF PDS TO BETTER ESTIMATE NUMAX AND BACKGROUND FIT
    ################################################################################

    # Find the beginning and end of the spike
    spike_freqs = cvs_freqs[non_solar]
    interp_flag = 0
    
    if len(spike_freqs) != 0:
        nspikes = np.where(np.diff(np.log10(spike_freqs)) > 0.001224 * 2)[0]
        interp_flag = len(nspikes) + 1
        
        copy_power = np.copy(np.array(pds['power']))

        n = 0
        while n < len(nspikes) + 1:

            # 1. Select frequency range of the peak
            if len(nspikes) == 0:
                local_spike_freq = spike_freqs

            else:
                if n == 0:     
                    cut = nspikes[n]
                    local_spike_freq = spike_freqs[0 : cut + 1]
                    

                elif n > 0 and n < len(nspikes):
                    local_spike_freq = spike_freqs[nspikes[n - 1] + 1 : nspikes[n] + 1]

                else:
                    local_spike_freq = spike_freqs[nspikes[n - 1] + 1 : ]

            
            # 2. Select region where power really goes higher than noise
            if len(local_spike_freq) == 1:
                mask_local = (pds['frequency'] > local_spike_freq[0] - 0.05) & (pds['frequency'] < local_spike_freq[0] + 0.05)

            else:
                min_spike = min(local_spike_freq)
                max_spike = max(local_spike_freq)
                
                if len(local_spike_freq) == 2:
                    mask_local = (pds['frequency'] >= min_spike - 0.1) & (pds['frequency'] <= max_spike + 0.1)
            
                else:
                    mask_local = (pds['frequency'] >= min_spike) & (pds['frequency'] <= max_spike)
                

            local_spike_power = np.array(pds['power'][mask_local])
            local_median = np.median(local_spike_power)
            
            mask_spike = local_spike_power > 8 * local_median
            local_spike_power[mask_spike] = np.nan

            
            mean_mask = np.mean(local_spike_power[~np.isnan(local_spike_power)])
            sigma_mask = np.std(local_spike_power[~np.isnan(local_spike_power)])
            size_mask = len(local_spike_power[mask_spike])

            
            # Remove the spike from pds
            # Linear interpolation:
            non_nan = ~np.isnan(local_spike_power)
            
            xp = non_nan.ravel().nonzero()[0]
            fp = local_spike_power[~np.isnan(local_spike_power)]
            x = np.isnan(local_spike_power).ravel().nonzero()[0]

            interp_pds = np.copy(local_spike_power)

            # Add some normally distributed noise
            gaussian_noise = np.random.normal(mean_mask, sigma_mask, size_mask)

            # Interpolated spike region + noise            
            interp_pds[np.isnan(local_spike_power)] = np.interp(x, xp, fp) + gaussian_noise

            # Final pds:
            copy_power[mask_local] = interp_pds

            n += 1


        final_pds = pd.DataFrame(data = {'frequency': np.array(pds['frequency']),
                                        'power': copy_power})

        pds = final_pds
    

    return interp_flag, interp_cvs, cvs_freqs, pds

def find_solarlike(freqs, cvs, fap, peaks, cv_peaks):
    # Freqs and cvs are taken from the oversampled spectrum
    idxs = cvs > 1.0

    osc_freqs = freqs[idxs]
    osc_cvs = cvs[idxs]
    osc_fap = fap[idxs]


    plt.rcParams['figure.figsize'] = (12, 6)
    plt.plot(freqs, cvs, color = 'forestgreen')
    plt.axhline(y = 1, xmin = 0, xmax = 1, color = 'k', linestyle = '--')
    plt.plot(osc_freqs, osc_cvs, color = 'deepskyblue')
    plt.plot(osc_freqs, osc_fap, color = 'red', linestyle = 'dotted')
    
    # If there is more than just one oscillation, let's find where these oscillations divide
    division = np.where(np.diff(np.log10(osc_freqs)) > 0.001224 * 3)[0]
    
    potential_solarlike, solarlike_freqs, solarlike_cvs = [], [], []

    n = 0
    while n < len(division) + 1 :
        if len(division) == 0:
            region = osc_freqs
            region_cvs = osc_cvs
            region_fap = osc_fap
            
        else:
            if n == 0:
                last = division[n] + 1
                region  =  osc_freqs[0 : last]
                region_cvs = osc_cvs[0 : last]
                region_fap = osc_fap[0 : last]
                
            elif n > 0 and n < len(division):
                first = division[n - 1] + 1
                last = division[n] + 1
                region  =  osc_freqs[first : last]
                region_cvs = osc_cvs[first : last]
                region_fap = osc_fap[first : last]


            else:
                first = division[n - 1] + 1
                region  =  osc_freqs[first :]
                region_cvs = osc_cvs[first :]
                region_fap = osc_fap[first :]

        
        over_fap = region_cvs > region_fap

        regions = region[over_fap]
        regions_cvs = region_cvs[over_fap]

        if len(regions) > 0:
            division2 = np.where(np.diff(np.log10(regions)) > 0.001224 * 10)[0]
            
            m = 0
            while m < len(division2) + 1 :
                
                if len(division2) > 0:
                    if m == 0:
                        cut = division2[m]
                        region2  =        regions[0 : cut + 1]
                        region_cvs2 = regions_cvs[0 : cut + 1]

                        
                    elif m > 0 and m < len(division2):
                        region2  =        regions[division2[m - 1] + 1: division2[m] + 1]
                        region_cvs2 = regions_cvs[division2[m - 1] + 1: division2[m] + 1]
                        

                    else:
                        region2  =        regions[division2[m - 1] + 1:]
                        region_cvs2 = regions_cvs[division2[m - 1] + 1:]
                        
                        
                else:
                    region2 = regions
                    region_cvs2 = regions_cvs

                    
                numax_region = region2[int(round(len(region2) / 2, 0))]
                cv_numax = region_cvs2[int(round(len(region_cvs2) / 2, 0))]
                

                # According to Bell+2019, the width of solar-like oscillations changes as a function of numax or delta nu
                if numax_region < 100:
                    # For numax < 100muHz (probed by Mosser+12):
                    #delta_env = 0.66 * (numax_region ** 0.88)
                    # But according to new results with TACO I found:
                    delta_env = 0.9* ((3.083e-6 * (numax_region ** 3)) - (0.001801 * (numax_region ** 2)) + (0.4141 * numax_region) - 0.03592)
        
                else:
                    # For numax > 100muHz they propose
                    #delta_env = 4.2 * delta_nu(numax = numax_region)
                    # My estimation:
                    delta_env =  0.9* ((3.083e-6 * (numax_region ** 3)) - (0.001801 * (numax_region ** 2)) + (0.4141 * numax_region) - 0.03592)

                if (max(region2) - min(region2)) * 1.1 > delta_env: 
                    cv_max = 2.69 * ((numax_region) ** 0.154)

                    if (cv_numax < cv_max):
                        
                        solar = False
                        spikes = 0
                        for i, j in enumerate(peaks):
                            if j >= min(region) and j <= max(region):
                                if cv_peaks[i] < cv_max:
                                    solar = True
                                
                                else:
                                    spikes += 1
                                    #solar = False

                        
                        if solar and spikes < 2:
                            potential_solarlike.append(numax_region)
                            solarlike_freqs.append(region2)
                            solarlike_cvs.append(region_cvs2)

                m += 1

        n += 1 

    
    return np.array(potential_solarlike)


def cv_method(pds):

    ################################################################################
    #                                  DATA                                        #
    ################################################################################
    # File with 43 non-overlapping bins 
    file_43bins = '43_nonoverlapping_bins.csv'
    flag_file = False

    if file_43bins not in os.listdir('./'):
        print('File not found. Must create it!')
        flag_file = True

    
    ################################################################################
    #                                   ANALYSIS                                   #
    ################################################################################
    flag, interpolation, solar_peaks, non_solar_peaks, numax_candidates  = [], [], [], [], []


    # Instead of computing the bins for each star, the code will compute them once and use the same bins for all of them.
    if flag_file:
        print('Computing bins')
        central_ind, edges_ind = bins_ind(f = pds)
    
        bins43 = pd.DataFrame(data = {'central': central_ind, 'edges': edges_ind})
        bins43.to_csv('./' + file_43bins, index = False)

    else:
        print('Loading bins file')
        bins43 = pd.DataFrame(pd.read_csv(file_43bins))
        central_ind = np.array(bins43['central'])
        edges_ind = np.array(bins43['edges'])

    central_ind = central_ind[1:]

    # Following Bell+19 method
    # STEP 1:   
    # Compute CV values for non-overlapping bins (independent spectrum)
    cvs_ind, fap_01pct_ind = cv_bins_ind(f = pds, edges = edges_ind)

    # STEP 2:   
    # Compare CV values tot the 0.1 per cent FAP threshold to find potential peaks:
    peaks_fap_ind = find_peaks(cvs_ind, fap_01pct_ind)
    
    # STEP 3:   
    # Compute CV values for overlapping bins (oversampled spectrum)
    freqs_os, cvs_os, fap_01pct_os = cv_bins_os(pds)

    # STEP 4:   
    # Classify peaks to check for contaminating signals
    non_solar, solar = filter_peaks(peaks = peaks_fap_ind, peaks_freq = central_ind, peaks_cvs = cvs_ind)

    # STEP 5:   
    # Interpolate non solar-like peaks
    interp_flag, interp_cvs_os, interp_cvs_freqs, interp_pds = interpolate_spikes(cvs = cvs_os, cvs_freqs = freqs_os, pds = pds)
    
    # STEP 6:   
    # Find potential solar-like oscillations
    candidates = find_solarlike(freqs = interp_cvs_freqs, cvs = interp_cvs_os, fap = fap_01pct_os,
                                                peaks = central_ind[peaks_fap_ind],
                                                cv_peaks = cvs_ind[peaks_fap_ind])


    # Flag multiple power excess:
    # It is equal to the number of potential solar-like oscillators found
    cv_flag = len(candidates)

    
    # Summary data:
    # A. From filter_peaks
    solar_peaks.append(len(solar))
    non_solar_peaks.append(len(non_solar))

    # B. From interpolation
    interpolation.append(interp_flag)

    # C. From find_solarlike_overlapped_bins
    numax_candidates.append(candidates)
    flag.append(cv_flag)
    
    results = pd.DataFrame(data = {'flag_cv': flag,
                                   'numax_estimate': numax_candidates, 
                                   'number_solar_peaks': solar_peaks, 
                                   'number_non_solar_peaks': non_solar_peaks, 
                                   'interp_spikes_flag': interpolation})
    
    return results, interp_pds

    

    