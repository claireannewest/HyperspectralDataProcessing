import numpy as np
from scipy.optimize import curve_fit
from nptdms import TdmsFile


class HyperspectralAnalysis:

    def __init(self):
        """Defines system paramters
        Keyword arguments:

        """

    def read_tdms(self, tdms_file, datalen):
        ncol = int(tdms_file.properties['strips'])
        bot_pix = int(tdms_file.properties['bottom pixel'])
        top_pix = int(tdms_file.properties['top pixel'])
        nrow = int(bot_pix - top_pix + 1)
        data_cube = np.zeros((int(ncol * nrow), datalen))
        idx = 0
        for group in tdms_file.groups():
            group_name = group.name
            for channel in group.channels():
                # channel_name = channel.name
                # properties = channel.properties
                data = channel[:]
                data_cube[idx, :] = data
                idx = idx + 1
            break
        data_cube = np.reshape(data_cube, ((ncol, nrow, datalen)))
        for group in tdms_file.groups():
            group_name = group.name
            if group_name == 'Spectra':
                continue
            for channel in group.channels():
                # channel_name = channel.name
                # properties = channel.properties
                waves = channel[:]
                break
        pixel_params = [bot_pix, top_pix, ncol, nrow]
        return pixel_params, waves, data_cube

    def process_wcdc(self, dc_raw, bc_raw, pixel_params):
        dc = np.mean(dc_raw, 0)
        bc = np.mean(bc_raw, 0)
        standard = bc - dc  # called stan in matlab
        stan = standard - np.min(standard) + 0.1  # don't think this matters
        bot_pix, top_pix, ncol, nrow = pixel_params
        stanim = np.zeros((ncol, nrow, 670))
        for i in range(ncol):
            if stan.shape[1] == 1340:
                stanim[i, :, :] = stan[top_pix - 1:bot_pix, ::2]
            if stan.shape[1] == 670:
                stanim[i, :, :] = stan[top_pix - 1:bot_pix, :]
        return stanim

    def back_global(self, I_raw, percent):
        bkavg = np.zeros((670,))
        numPixelsBack = int(np.ceil(percent * I_raw.shape[0] * I_raw.shape[1]))
        I_sumwave = np.sum(I_raw, 2)  # sum all wavelengths
        # Normalize I_sumwave such that it ranges [0,1]
        I_norm = ((I_sumwave - np.min(I_sumwave))
                  / np.max(I_sumwave - np.min(I_sumwave)))
        # All unique values of inten_norm, and sort small to large
        I_norm_sorted = np.unique(I_norm)
        # The cutoff value of I_norm_sorted that meets the threshold
        nthsmlst = I_norm_sorted[numPixelsBack]
        # Find all pixels that are less than the cutoff value
        idx_a, idx_b = np.where(I_norm < nthsmlst)
        idx_a_cut = idx_a[:numPixelsBack]
        idx_b_cut = idx_b[:numPixelsBack]
        # Find the average of all the pixels
        bkavg = sum(I_raw[idx_a_cut, idx_b_cut, :]) / numPixelsBack
        bkavg = np.round(bkavg, 5)  # This is to exactly match matlab
        return bkavg, [idx_a_cut, idx_b_cut]

    def back_local_grid(self, specfin, xi, yi, offset):
        back_grid = np.zeros((670, 4))
        back_grid[:, 0] = self.avg3pixels(M=specfin, xi=xi, yi=(yi + offset))
        back_grid[:, 1] = self.avg3pixels(M=specfin, xi=xi, yi=(yi - offset))
        back_grid[:, 2] = self.avg3pixels(M=specfin, xi=(xi + offset), yi=yi)
        back_grid[:, 3] = self.avg3pixels(M=specfin, xi=(xi - offset), yi=yi)
        whichmin = np.argmin(np.sum(back_grid, axis=0))
        back_mingrid = back_grid[:, whichmin]
        return back_grid, whichmin, back_mingrid

    def fit_2Dgauss(self, M, a, x0, c1, y0, c2, d):
        # For fitting local background.
        x, y = M
        return 1E3 * a * np.exp(-((x - x0)**2 / (2 * c1**2) + (y - y0)**2
                                  / (2 * c2**2))) + d * 1E3

    def fitNP_eachwave(self, specfin, waveidx, xi, yi):
        # For fitting local background.
        # Find the background at a given wavelength
        specfin_wind = specfin[xi - 10: xi + 10, yi - 10: yi + 10, waveidx]
        xran = np.linspace(0, specfin_wind.shape[0] - 1, specfin_wind.shape[0])
        yran = np.linspace(0, specfin_wind.shape[1] - 1, specfin_wind.shape[1])
        X, Y = np.meshgrid(xran, yran)
        xdata = np.vstack((X.ravel(), Y.ravel()))
        Z = specfin_wind.ravel()
        p0 = [1., 10., 0.86, 10., 1., 2.]
        bounds = [[0, 0, 0, 0, 0, 0],
                  [1000, 1000, 1000, 1000, 1000, 1000]]
        try:
            popt, pcov = curve_fit(self.fit_2Dgauss,
                                   xdata, Z, p0, bounds=bounds)
            fitdata = self.fit_2Dgauss(xdata, *popt).reshape(len(yran),
                                                             len(xran))
        except (ValueError, RuntimeError):
            fitdata = []
            popt = []
        return specfin_wind, fitdata, popt

    def fit_all_NPs_eachwave(self,
                             inten_raw,
                             positions,
                             background,
                             numPart,
                             numWaves,
                             wavei,
                             wavef):
        # For fitting local background
        # At some wavelengths, we cannot fit a background.
        # In this case, the code just grabs the previously defined background.
        # I'll silence the error messages that arise.
        import warnings
        warnings.simplefilter("ignore")
        for npi in range(numPart):
            countfits = np.zeros(numWaves)
            xi = int(positions[npi, 1])
            yi = int(positions[npi, 2])
            for wi in range(numWaves):
                _, fitdata, popt = self.fitNP_eachwave(specfin=inten_raw,
                                                       waveidx=(wi + wavei),
                                                       xi=xi,
                                                       yi=yi)
                if fitdata == []:
                    bk_global, _ = self.back_global(I_raw=inten_raw,
                                                    percent=0.1)
                    background[npi, wi] = bk_global[wi + wavei]
                else:
                    background[npi, wi] = np.round(popt[-1] * 1E3, 1)
                    countfits[wi] = 1
            countwaves = (str('NP ') + str(npi) + str(': fit ')
                          + str(int(np.sum(countfits) / numWaves * 100))
                          + str('% of wavelengths'))
            print(countwaves)
        return background

    def calc_DFS(self, inten_raw, wc_minus_dc, back):
        specfin = np.round((inten_raw - back) / wc_minus_dc, 5)
        return specfin

    def avg3pixels(self, M, xi, yi):
        sum_M = np.sum(M[xi - 1:xi + 2, yi - 1:yi + 2, :], axis=(0, 1))
        return np.round(sum_M, 5)

    def calc_DFS_localback(self,
                           inten_raw,
                           wc_minus_dc,
                           background,
                           npi,
                           xi,
                           yi,
                           wavei,
                           wavef):
        dfs_avg3 = self.avg3pixels(M=(inten_raw / wc_minus_dc),
                                   xi=xi,
                                   yi=yi,
                                   )
        total = np.round(dfs_avg3[wavei:wavef] - 9 * background[npi, :]
                         / wc_minus_dc[xi, yi, wavei:wavef], 5)
        return total

    def calc_DFS_localgrid(self,
                           specfin,
                           numPart,
                           positions,
                           wavei,
                           wavef,
                           wc_minus_dc,
                           offset,
                           ):
        total = np.zeros((wavef - wavei, numPart))
        for npi in range(numPart):
            xi = int(positions[npi, 1])
            yi = int(positions[npi, 2])
            _, _, back_mingrid = self.back_local_grid(specfin=specfin,
                                                      xi=xi,
                                                      yi=yi,
                                                      offset=offset)
            dfs_avg3 = self.avg3pixels(M=(specfin / wc_minus_dc),
                                       xi=xi,
                                       yi=yi,
                                       )
            total[:, npi] = (dfs_avg3[wavei:wavef] - back_mingrid[wavei:wavef]
                             / wc_minus_dc[xi, yi, wavei:wavef])
        return np.round(total, 5)

    def lorentz(self, wave_eV, A, Gam_eV, wave0_eV):
        return A * 0.5 * Gam_eV / ((wave_eV - wave0_eV)**2 + (0.5 * Gam_eV)**2)

    def fit_spectrum(self, wave, inten):
        A_ran = [1E-3, 5]
        Gam_ran = [0.1, 2]
        wave0_ran = [1, 3]
        energy = 1240 / wave
        popt, pcov = curve_fit(self.lorentz, energy, inten,
                               p0=[0.5, 0.5, 1.7],
                               bounds=[[A_ran[0], Gam_ran[0], wave0_ran[0]],
                                       [A_ran[1], Gam_ran[1], wave0_ran[1]]]
                               )
        if popt[0] in A_ran or popt[1] in Gam_ran or popt[2] in wave0_ran:
            print('Error with Lorentzian fit.')
        fit = self.lorentz(energy, *popt)
        return fit, popt

    def process_data(self, path_int, pathdc, pathbc, path_positions):
        pxl, wave, inten = self.read_tdms(tdms_file=TdmsFile.read(path_int),
                                          datalen=670)
        _, w1, dc_raw = self.read_tdms(tdms_file=TdmsFile.read(pathdc),
                                       datalen=1340)
        _, w2, bc_raw = self.read_tdms(tdms_file=TdmsFile.read(pathbc),
                                       datalen=1340)
        wc_minus_dc = self.process_wcdc(dc_raw=dc_raw,
                                        bc_raw=bc_raw,
                                        pixel_params=pxl)
        positions = np.loadtxt(path_positions)
        numPart = len(positions)
        return numPart, positions, wave, inten, wc_minus_dc

    def calc_all_DFS(self,
                     wave_raw,
                     inten_raw,
                     numPart,
                     positions,
                     wc_minus_dc):
        wavei = 0  # index of starting wavelength
        wavef = 670  # index of end wavelength
        wave = wave_raw[wavei:wavef]
        gam_tot = np.zeros(numPart)
        eres_tot = np.zeros(numPart)
        # fit_tot = np.zeros((wavef - wavei, numPart))
        fit_tot = np.zeros((670, numPart))

        dfs_tot = self.calc_DFS_localgrid(specfin=inten_raw,
                                          numPart=numPart,
                                          positions=positions,
                                          wavei=wavei,
                                          wavef=wavef,
                                          wc_minus_dc=wc_minus_dc,
                                          offset=4)
        for npi in range(0, numPart):
            DFSi = dfs_tot[:, npi]
            fit_tot[:, npi], popt = self.fit_spectrum(wave, DFSi)
            gam_tot[npi] = 1240 * popt[1] / popt[2]**2
            eres_tot[npi] = 1240 / popt[2]
        return wave, dfs_tot, fit_tot, gam_tot, eres_tot
