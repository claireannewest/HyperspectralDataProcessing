import numpy as np


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
        return bot_pix, top_pix, waves, data_cube

    def process_wcdc(self, dc_raw, bc_raw):
        dc = np.mean(dc_raw, 0)
        bc = np.mean(bc_raw, 0)
        standard = bc - dc  # called stan in matlab
        stan = standard - np.min(standard) + 0.1  # don't think this matters
        return stan

    def calc_background_avg(self, I_raw, percent):
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
        return bkavg

    def calc_DFS_onepixel(self, inten_raw, stan, bkavg, top_pix, bot_pix):
        stanim = np.zeros((103, 261, 670))
        for i in range(103):
            stanim[i, :, :] = stan[top_pix - 1:bot_pix, ::2]
        specfin = np.round((inten_raw - bkavg) / stanim, 5)
        return specfin

    def calc_DFS_avg3pixels(self, specfin, xi, yi):
        sum_DFS = np.sum(specfin[xi - 1:xi + 2, yi - 1:yi + 2, :], axis=(0, 1))
        return np.round(sum_DFS, 5)

