from pspy import so_map, so_window, sph_tools, so_mcm, pspy_utils, so_spectra, so_cov #flat_tools
from pixell import enmap
import time
import numpy as np, pylab as plt

map0_info = {"name":"split0_IQU.fits", "data_type":"IQU", "id":"split0", "cal" : None}
map1_info = {"name":"split1_IQU.fits", "data_type":"IQU", "id":"split1", "cal" : None}


source_mask = {"name":"source_mask.fits", "apo_type":"C1" , "apo_radius": 0.3}
galactic_mask = "mask_galactic_equatorial_car_halfarcmin.fits"

patch = {"patch_type": "Rectangle", "patch_coordinate": [[-5,-20],[5,20]]}

maps_info_list = [map0_info, map1_info]


def compute_ps(patch,
               maps_info_list,
               ps_method = "master",
               error_method_1d = "master",
               type="Dl",
               beam=None,
               binning_file=None,
               bin_size=None,
               source_mask = None,
               galactic_mask = None,
               apo_radius_survey=1,
               compute_T_only = False,
               lmax = 1000,
               lmax_pad = None,
               master_threshold = None):
               
    """Compute spectra

    Parameters
    ----------
    patch: dict
      a dict containing the patch type and coordinates
      if patch_type is "Rectangle" the coordinate are expected to be the 4 corners
      if patch_type is "Disc" we expect the coordinate of the center and the radius in degree
    maps_info_list: list of dicts describing the data maps
      dictionnary should contain the name, the data type ("IQU" or "I") and optionally a calibration factor to apply to the map
      note that all map in the list should have the same data type
    beam: text file
      file describing the beam of the map, expect bl to be the second column and start at l=0 (standard is : l,bl, ...)
    binning_file: text file
      a binning file with three columns bin low, bin high, bin mean
      note that either binning_file or bin_size should be provided
    bin_size: integer
      the bin size
      note that either binning_file or bin_size should be provided
    type: string
      the type of binning, either bin Cl or bin Dl
    source_mask: dict
      a dict containing an optional source mask and its properties
      the dictionnary should contain the name, the type of apodisation and the radius of apodisation
    galactic_mask: fits file
      an optional galactic mask to apply
    apo_radius_survey: float
      the apodisation radius in degree
    ps_method: string
      the method for the computation of the power spectrum
      can be "master", "pseudo", or "2dflat" for now
    error_method_1d: string
      the method for the computation of error
      can be "master" or "knox" for now
    compute_T_only: boolean
      True to compute only T spectra, should always be true for data_type= "I"
    lmax : integer
      the maximum multipole to consider for the spectra computation
    lmax_pad: integer
      the maximum multipole to consider for the mcm computation (optional)
      lmax_pad should always be greater than lmax
    master_threshold: integer
      for approximating the coupling computation, threshold is the max |l1-l2| considered
      in the calculation
    """

    if binning_file is None:
        pspy_utils.create_binning_file(bin_size=bin_size, n_bins=1000, file_name="binning.dat")
        binning_file = "binning.dat"
        
    bin_lo, bin_hi, bin_c, bin_size = pspy_utils.read_binning_file(binning_file, lmax)
    n_bins = len(bin_hi)

    if patch["patch_type"] == "Rectangle":
        car_box = patch["patch_coordinate"]
        window = so_map.read_map(maps_info_list[0]["name"], car_box = car_box)
        if maps_info_list[0]["data_type"]=="IQU":
            window.data = window.data[0]
            window.ncomp = 1
        window.data[:] = 0
        window.data[1:-1, 1:-1] = 1
        apo_type_survey = "Rectangle"
        
    elif patch["patch_type"] == "Disk":
        apo_type_survey = "C1"
        # ToDO

        
    if galactic_mask is not None:
        gal_mask = so_map.read_map(galactic_mask, car_box = car_box)
        window.data *= gal_mask.data
        del gal_mask
    
    for map_info in maps_info_list:
        split = so_map.read_map(map_info["name"], car_box = car_box)
        if (compute_T_only == True) & (map_info["data_type"] == "IQU"):
            split.data = split.data[0]
            split.ncomp = 1
            
        if split.ncomp == 1:
            window.data[split.data==0] = 0
        else:
            for i in range(split.ncomp):
                window.data[split.data[i]==0] = 0

    
    window =  so_window.create_apodization(window,
                                           apo_type=apo_type_survey,
                                           apo_radius_degree=apo_radius_survey)

    if source_mask is not None:
        ps_mask = so_map.read_map(source_mask["name"], car_box = car_box)
        ps_mask = so_window.create_apodization(ps_mask,
                                               apo_type=source_mask["apo_type"],
                                               apo_radius_degree=source_mask["apo_radius"])
        window.data *= ps_mask.data
        del ps_mask
    
    fsky = enmap.area(window.data.shape, window.data.wcs) / 4. / np.pi
    fsky *= np.mean(window.data)
    
    if beam is not None:
        beam_data = np.loadtxt(beam)
        if compute_T_only == True:
            beam = beam_data[:,1]
        else:
            beam = (beam_data[:,1], beam_data[:,1])
        
    if compute_T_only == True:
        if ps_method == "master":
            mbb_inv, Bbl = so_mcm.mcm_and_bbl_spin0(window,
                                                    binning_file,
                                                    bl1=beam,
                                                    lmax=lmax,
                                                    type=type,
                                                    niter=0,
                                                    lmax_pad=lmax_pad,
                                                    threshold=master_threshold)
                                                    
        elif ps_method == "pseudo":
            mbb_inv = np.identity(n_bins)
            mbb_inv *= 1/fsky
    else:
        window = (window, window)
        if ps_method == "master":
            print("compute master MCM")
            mbb_inv, Bbl = so_mcm.mcm_and_bbl_spin0and2(window,
                                                        binning_file,
                                                        bl1=beam,
                                                        lmax=lmax,
                                                        type=type,
                                                        niter=0,
                                                        lmax_pad=lmax_pad,
                                                        threshold=master_threshold)
        elif ps_method == "pseudo":
            mbb_inv = {}
            spin_list = ["spin0xspin0", "spin0xspin2", "spin2xspin0"]
            for spin in spin_list:
                mbb_inv[spin] = np.identity(n_bins)
                mbb_inv[spin] *= 1/fsky
            mbb_inv["spin2xspin2"] = np.identity(4 * n_bins)
            mbb_inv["spin2xspin2"] *= 1/fsky

    ht_list = []
    name_list = []

    for map_info in maps_info_list:
    
        split = so_map.read_map(map_info["name"], car_box = car_box)
        
        if (compute_T_only == True) & (map_info["data_type"] == "IQU"):
            split.data = split.data[0]
            split.ncomp = 1

        if map_info["cal"] is not None:
            split.data *= map_info["cal"]
        
        if ps_method == "master" or ps_method == "pseudo":
            print ("SPHT of %s in the patch" % map_info["name"])
            alms = sph_tools.get_alms(split, window, niter=0, lmax=lmax+50)
            ht_list += [alms]
            
        elif ps_method == "2dflat":
            print ("FFT of %s in the patch" % map_info["name"])
            fts = flat_tools.fft(windowed_map, lmax=lmax+50)
            ht_list += [fts]

        name_list += [map_info["id"]]
        
    split_num = np.arange(len(maps_info_list))

    if (compute_T_only == True):
        spectra = None
    else:
        spectra = ["TT", "TE", "TB", "ET", "BT", "EE", "EB", "BE", "BB"]

    ps_dict = {}
    spec_name_list = []

    for name1, ht1, c1 in zip(name_list, ht_list, split_num):
        for name2, ht2, c2 in zip(name_list, ht_list, split_num):
            if c1 > c2: continue
            
            spec_name = "%sx%s" % (name1, name2)

            if ps_method == "master" or ps_method == "pseudo":
                l, ps = so_spectra.get_spectra(ht1, ht2, spectra=spectra)
                lb, ps_dict[spec_name] = so_spectra.bin_spectra(l,
                                                                ps,
                                                                binning_file,
                                                                lmax,
                                                                type=type,
                                                                mbb_inv=mbb_inv,
                                                                spectra=spectra)
            
            spec_name_list += [spec_name]
            
            
    if (error_method_1d is not None) & (ps_method is not "2dflat"):
    
        cov_dict = {}
        
        if compute_T_only == True:
            if error_method_1d == "Knox":
                for name in spec_name_list:
                    m1, m2 = name.split("x")
                    prefac = 1 / ((2*lb+1) * fsky * bin_size)
                    cov_dict[name] = np.diag(prefac * (ps_dict["%sx%s"%(m1,m1)] * ps_dict["%sx%s"%(m2,m2)] + ps_dict["%sx%s"%(m1,m2)]**2))
                    
            elif error_method_1d == "master":
                print("compute master error")
                coupling_dict = so_cov.cov_coupling_spin0(window, lmax, niter=0, threshold=master_threshold)
                coupling = so_cov.bin_mat(coupling_dict["TaTcTbTd"], binning_file, lmax)

                for name in spec_name_list:
                    m1, m2 = name.split("x")

                    cov_dict[name] = so_cov.symmetrize(ps_dict["%sx%s"%(m1,m1)]) * so_cov.symmetrize(ps_dict["%sx%s"%(m2,m2)])
                    cov_dict[name] += so_cov.symmetrize(ps_dict["%sx%s"%(m1,m2)]**2)
                    cov_dict[name] *= coupling
                    
                    cov_dict[name] = np.dot(np.dot(mbb_inv, cov_dict[name]), mbb_inv.T)
        else:
        
            if error_method_1d == "Knox":
                for name in spec_name_list:
                    m1, m2 = name.split("x")
                    cov_dict[name] = {}
                    for spec in spectra:
                        X, Y = spec
                        prefac = 1 / ((2 * lb + 1) * fsky * bin_size)
                        cov_dict[name][X+Y] =  np.diag( prefac *  (ps_dict["%sx%s"%(m1,m1)][X+X] * ps_dict["%sx%s"%(m2,m2)][Y+Y] + ps_dict["%sx%s"%(m1,m2)][X+Y]**2))
                        
            elif error_method_1d == "master":
                print("compute master error")
                coupling_dict = so_cov.cov_coupling_spin0(window[0], lmax, niter=0, threshold=master_threshold)
                coupling = so_cov.bin_mat(coupling_dict["TaTcTbTd"], binning_file, lmax)

                for name in spec_name_list:
                    m1, m2 = name.split("x")
                    cov_dict[name] = {}
                    for spec in spectra:
                        X, Y = spec
                        cov_dict[name][X+Y] = so_cov.symmetrize(ps_dict["%sx%s"%(m1,m1)][X+X]) * so_cov.symmetrize(ps_dict["%sx%s"%(m2,m2)][Y+Y])
                        cov_dict[name][X+Y] += so_cov.symmetrize(ps_dict["%sx%s"%(m1,m2)][X+Y]**2)
                        cov_dict[name][X+Y] *= coupling
                        cov_dict[name][X+Y] = np.dot(np.dot(mbb_inv["spin0xspin0"], cov_dict[name][X+Y]), mbb_inv["spin0xspin0"].T)

        return spectra, spec_name_list, lb, ps_dict, cov_dict
        
    else:
    
        return spectra, spec_name_list, lb, ps_dict
   
     
     
     
t = time.time()
spectra, spec_name_list, lb, ps_dict, cov_dict = compute_ps(patch,
                                                            maps_info_list,
                                                            ps_method = "master",
                                                            error_method_1d = "master",
                                                            beam=None,
                                                            binning_file=None,
                                                            bin_size=40,
                                                            type="Dl",
                                                            source_mask = source_mask,
                                                            galactic_mask = galactic_mask,
                                                            apo_radius_survey=1,
                                                            compute_T_only = False,
                                                            lmax = 4000)
print ("time to compute exact Cl: %0.2f"% (time.time()-t))


t = time.time()
spectra, spec_name_list, lb, ps_dict_threshold, cov_dict_threshold = compute_ps(patch,
                                                                              maps_info_list,
                                                                              ps_method = "master",
                                                                              error_method_1d = "master",
                                                                              beam=None,
                                                                              binning_file=None,
                                                                              bin_size=40,
                                                                              type="Dl",
                                                                              source_mask = source_mask,
                                                                              galactic_mask = galactic_mask,
                                                                              apo_radius_survey=1,
                                                                              compute_T_only = False,
                                                                              lmax = 4000,
                                                                              master_threshold= 500)
print ("time to compute threshod Cl: %0.2f"% (time.time()-t))






spec_name_list = ["split0xsplit1"]
if spectra is None:
    for spec_name in spec_name_list:
        plt.semilogy()
        plt.errorbar(lb, ps_dict[spec_name], np.sqrt(np.diag(cov_dict[spec_name])), fmt = ".")
        plt.errorbar(lb - 5, ps_dict_threshold[spec_name], np.sqrt(np.diag(cov_dict_threshold[spec_name])), fmt = ".")

    plt.legend()
    plt.show()
    
else:
    for spec in spectra:
        for spec_name in spec_name_list:
            if spec=="TT":
                plt.semilogy()
            plt.errorbar(lb, ps_dict[spec_name][spec], np.sqrt(np.diag(cov_dict[spec_name][spec])), fmt = ".")
            plt.errorbar(lb - 5, ps_dict_threshold[spec_name][spec], np.sqrt(np.diag(cov_dict_threshold[spec_name][spec])), fmt = ".")

        plt.legend()
        plt.show()

