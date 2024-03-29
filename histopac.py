#!/usr/bin/env python

from configparser import ConfigParser
from os import path
import sys
from struct import unpack
import logging
import json

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Pango

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk3agg import FigureCanvasGTK3Agg as FigureCanvas
from matplotlib.backends.backend_gtk3 import NavigationToolbar2GTK3 as NavigationToolbar

import click

import gui_params
from ref_spk import Ref_spk


det_num = 4
t_spk_num = 12
tsum_spk_num = t_spk_num // 2

histo_size = 4096 #(chanels)
histo_finisize = 512
histo_fsize = 16384 #(bytes)
histo_en_fnames = ["BUFKA" + str(i) + ".SPK" for i in range(1, 5)]
histo_t_fnames = ["TIME1.SPK", "TIME2.SPK", "TIME3.SPK", "TIME4.SPK",
                  "TIME5.SPK", "TIME6.SPK", "TIME7.SPK", "TIME8.SPK",
                  "TIME9.SPK", "TIME10.SPK", "TIME11.SPK", "TIME12.SPK"]

glob_vukap_ver = 1.0

wd_script = path.dirname(path.realpath(sys.argv[0]))
basic_ini_fname = path.join(wd_script, "ini.json")
refspk_fname = path.join(wd_script, "ref_transitions.json")
ini = None


@click.command()
@click.argument("dir_name", type=click.Path(exists=True))
def main(dir_name):
    """"Load and plot PAC spectra.
    histopac allows to view/set ENergy windows,
    analyze ENergy and Time peaks, analyze Time exponents.
    The pathways for files which are required for the program
    are specified in file {:s}
    Args:\n
        dir_name (str): Directory path with PAC spectra."""
    logging.basicConfig(format="{File %(filename)s, line %(lineno)d} \
%(levelname)s:%(asctime)s:%(message)s",
                        datefmt="%Y/%m/%d %H-%M-%S",
                        level=logging.INFO)

    global ini
    ini = parse_ini(basic_ini_fname)
    global glob_vukap_ver
    if path.exists(path.join(dir_name, "hdrainer.ini")) and \
       path.exists(path.join(dir_name, "dsp.ini")):
        #VUKAP 2.0
        print("The spectra in folder are from VUKAP 2.0")
        glob_vukap_ver = 2.0
        cfg = parse_hdr_ini2(path.join(dir_name, "hdrainer.ini"))
        consts = parse_dsp_ini2(path.join(dir_name, "dsp.ini"))
        ini["hdr_path"] = path.join(dir_name, "hdrainer.ini")
        ini["dsp_path"] = path.join(dir_name, "dsp.ini")
    else:
        #VUKAP 1.0
        print("The spectra in folder are from VUKAP 1.0")
        glob_vukap_ver = 1.0
        try:
            if ini["cfg_path"] == -1:
                cfg = parse_cfg(path.join(dir_name, "cfg.json"))
                #dir_name + "cfg.json" if dir_name[-1] == "/" else dir_name + "/cfg.json")
            else:
                cfg = parse_cfg(ini["cfg_path"])
        except FileNotFoundError:
            cfg = {}
            print("cfg file Not Found in folder with spk")

        try:
            if ini["const_path"] == -1:
                consts = parse_constants(path.join(dir_name, "constants.json"))
                #dir_name + "constants.json" if dir_name[-1] == "/" else dir_name + "/constants.json")
            else:
                consts = parse_constants(ini["const_path"])
        except FileNotFoundError:
            consts = {"T_SCALE": [1.0, 1.0]}
            print("const file Not Found in folder with spk")

    ui = Create_UI(dir_name, consts)

    en_spk, t_spk = get_histos_from_folder(dir_name)

    ui.set_histos(en_spk, t_spk)
    ui.plot_all_histos()

    ui.set_cfg(cfg)

    ui.hide_all_spk_t()

    ui.show_all()
    Gtk.main()


class Analyze_peak():
    def __init__(self, en_spk, x_l, x_r):
        self.en_spk = np.array(en_spk, dtype=np.float64)
        self.x_l = x_l
        self.x_r = x_r

        self.bg = self.calc_bg()

        self.en_spk[self.x_l:self.x_r] -= self.bg
        #make zeros all negative elements in self.en_spk
        self.en_spk[self.en_spk < 0] = 0

        self.integral, self.integral_err = self.calc_integral()
        self.area, self.area_err = self.calc_area()
        self.mean, self.mean_err = self.calc_mean()
        self.fwhm, self.fwhm_err = self.calc_fwhm()
        self.resol = self.calc_resol()

        logging.info("mean = {}, integral = {}, area = {}, fwhm = {}".format(self.mean,
                                                                             self.integral,
                                                                             self.area,
                                                                             self.fwhm))


    def calc_bg(self):
        y_avg_l = sum(self.en_spk[self.x_l:self.x_l + 3]) / 3.0
        y_avg_r = sum(self.en_spk[self.x_r - 3:self.x_r]) / 3.0

        k_bg = (y_avg_r - y_avg_l) / (self.x_r - 1 - self.x_l - 1)
        b_bg = y_avg_l - k_bg * (self.x_l + 1)

        return np.array([k_bg * i + b_bg for i in range(self.x_l, self.x_r)], dtype=np.int64)


    def calc_integral(self):
        integral = np.sum(self.en_spk[self.x_l:self.x_r] + self.bg)
        integral_err = np.sqrt(integral)

        return integral, integral_err


    def calc_area(self):
        area = np.sum(self.en_spk[self.x_l:self.x_r])
        area_err = np.sqrt(area)

        return area, area_err
        

    def calc_mean(self):
        w_sum = np.sum(np.multiply(np.arange(self.x_l, self.x_r), self.en_spk[self.x_l:self.x_r]))
        mean = w_sum / sum(self.en_spk[self.x_l:self.x_r])
        mean_err = np.sqrt(np.sum((self.en_spk[self.x_l:self.x_r] + self.bg) *
                                   (mean * np.ones(self.x_r - self.x_l) - np.arange(self.x_l, self.x_r))**2) /
                             ((self.x_r - self.x_l) * np.sum(self.en_spk[self.x_l:self.x_r])))

        return mean, mean_err


    def calc_fwhm(self):
        arr = self.en_spk[self.x_l:self.x_r]
        mean, _ = self.calc_mean()

        while True:
            max_index = np.argmax(arr)
            #IF BAD PICK exit
            if arr[max_index] == 0:
                return (0, 0)
            #if |MAX - MEAN| too large find another MAX
            if abs(max_index - (mean-self.x_l)) > 10:
                arr[max_index] = 0
            else:
                break

        more_than_max2 = np.where(arr > arr[max_index]/2)[0]

        #find FWHM_left bound and FWHM_right bound
        yp_l = [more_than_max2[0]-1, more_than_max2[0]]
        xp_l = [arr[y] for y in yp_l]
        k_fwhm_l, b_fwhm_l = np.polyfit(xp_l, yp_l, 1)
        fwhm_l = k_fwhm_l * arr[max_index]/2 + b_fwhm_l
        yp_r = [more_than_max2[-1], more_than_max2[-1]+1]
        xp_r = [arr[y] for y in yp_r]
        k_fwhm_r, b_fwhm_r = np.polyfit(xp_r, yp_r, 1)
        fwhm_r = k_fwhm_r * arr[max_index]/2 + b_fwhm_r
        fwhm = fwhm_r-fwhm_l
        self.fwhm_y = arr[max_index]/2.0

        #set vars for plots
        self.fwhm_ch_l = k_fwhm_l * self.fwhm_y + b_fwhm_l#(self.fwhm_y / 2 - b_fwhm_l) / k_fwhm_l
        self.fwhm_ch_r = k_fwhm_r * self.fwhm_y + b_fwhm_r#(self.fwhm_y / 2 - b_fwhm_r) / k_fwhm_r
        print("fwhm_ch_l = {}, fwhm_ch_r = {}".format(self.fwhm_ch_l, self.fwhm_ch_r))
        print(yp_l, xp_l)
        print("k, b = {}, {}".format(k_fwhm_l, b_fwhm_l))
        print("by hand k, b = {}, {}".format((yp_l[0] - yp_l[1]) / (xp_l[0] - xp_l[1]),
                                             yp_l[0] - (yp_l[0] - yp_l[1]) / (xp_l[0] - xp_l[1]) * xp_l[0]))
        print(yp_r, xp_r)
        print("k, b = {}, {}".format(k_fwhm_r, b_fwhm_r))
        print("by hand k, b = {}, {}".format((yp_r[0] - yp_r[1]) / (xp_r[0] - xp_r[1]),
                                             yp_r[0] - (yp_r[0] - yp_r[1]) / (xp_r[0] - xp_r[1]) * xp_r[0]))

        #FWHM Error
        fwhm_err_l, fwhm_err_r = np.zeros(2), np.zeros(2)
        xp_l_err = [x+(abs(x))**0.5 for x in xp_l]
        k, b = np.polyfit(xp_l_err, yp_l, 1)
        fwhm_err_l[0] = k * arr[max_index]/2 + b
        xp_l_err = [x-(abs(x))**0.5 for x in xp_l]
        k, b = np.polyfit(xp_l_err, yp_l, 1)
        fwhm_err_l[1] = k * arr[max_index]/2 + b
        print(fwhm_l, fwhm_err_l[0], fwhm_err_l[1])

        xp_r_err = [x+(x)**0.5 for x in xp_r]
        k, b = np.polyfit(xp_r_err, yp_r, 1)
        fwhm_err_r[0] = k * arr[max_index]/2 + b
        xp_r_err = [x-(x)**0.5 for x in xp_r]
        k, b = np.polyfit(xp_r_err, yp_r, 1)
        fwhm_err_r[1] = k * arr[max_index]/2 + b
        print(fwhm_r, fwhm_err_r[0], fwhm_err_r[1])

        fwhm_err = ((fwhm_err_l[0]-fwhm_err_l[1])**2 + (fwhm_err_r[0]-fwhm_err_r[1])**2)**0.5

        return fwhm, fwhm_err

    def calc_resol(self):
        return self.fwhm / self.mean



class Calibr_en():
    def __init__(self):
        self.ch = [[-1.0, -1.0] for i in range(det_num)]
        self.en = [[-1.0, -1.0] for i in range(det_num)]

        self.set_ch_en_from_file(path.join(wd_script, ini["calibr_en_path"]))

        self.k = det_num * [0]
        self.b = det_num * [0]
        self.calc_k_b()


    def calc_k_b(self):
        for i in range(det_num):
            try:
                self.k[i] = (self.en[i][0] - self.en[i][1]) / (self.ch[i][0] - self.ch[i][1])
            except ZeroDivisionError:
                self.k[i] = 1.0
            self.b[i] = self.en[i][0] - self.k[i] * self.ch[i][0]


    def ch_to_keV(self, det_i, ch):
        return self.k[det_i] * ch + self.b[det_i]


    def keV_to_ch(self, det_i, keV):
        return (keV - self.b[det_i]) / self.k[det_i]


    def set_ch_en_from_file(self, fname):
        with open(fname, "r") as f:
            data = f.read(2048)
            vals = json.loads(data)

            for i in range(det_num):
                for j in range(2):
                    self.ch[i][j] = int(vals["channels"][i][j])
                    self.en[i][j] = float(vals["energies"][i][j])


    def set_ch_en_from_input(self, ch, en):
        self.ch = ch
        self.en = en


    def save_file_ch_en(self, fname):
        with open(fname, "w+") as f:
            d = {"channels": self.ch,
                 "energies": self.en}
            json.dump(d, f)


class Calibr_t():
    def __init__(self, T_SCALE):
        self.t_scale = T_SCALE
        if glob_vukap_ver == 2.0:
            self.adc_ns_per_ch = 4.0
        elif glob_vukap_ver == 1.0:
            self.adc_ns_per_ch = 8.0
        self.ns_per_ch = self.adc_ns_per_ch / ((1.0 + self.t_scale[1]) /
                                               (2.0 * self.t_scale[1] / histo_size) - 0.5 * histo_size)

    def ch_to_ns(self, ch):
        return ch * self.ns_per_ch

    def ns_to_ch(self, ns):
        return ns / self.ns_per_ch
        
            
class Dialog_calibr_en(Gtk.Dialog):
    def __init__(self, parent, calibr_en):
        Gtk.Dialog.__init__(self, "Energy calibration", parent, 0,
                            (Gtk.STOCK_OK, Gtk.ResponseType.OK,
                             Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL))

        self.set_default_size(200, 100)

        label_ch = Gtk.Label("Channels [ch]:")
        label_en = Gtk.Label("Energies [keV]:")
        
        self.entry_ch = [Gtk.Entry() for i in range(det_num)]
        self.entry_en = [Gtk.Entry() for i in range(det_num)]

        for i in range(det_num):
            self.entry_ch[i].set_text(str(calibr_en.ch[i][0]) + ', ' + str(calibr_en.ch[i][1]))
            self.entry_en[i].set_text(str(calibr_en.en[i][0]) + ', ' + str(calibr_en.en[i][1]))

        self.grid = Gtk.Grid()

        self.grid.attach(label_ch, 1, 0, 1, 1)
        self.grid.attach(label_en, 2, 0, 1, 1)
        for i in range(det_num):
            self.grid.attach(Gtk.Label("D{}".format(i+1)), 0, i + 1, 1, 1)
            self.grid.attach(self.entry_ch[i], 1, i + 1, 1, 1)
            self.grid.attach(self.entry_en[i], 2, i + 1, 1, 1)

        box = self.get_content_area()
        box.add(self.grid)

        self.show_all()



class Calc_view_en():
    def __init__(self):
        self.txtview = Gtk.TextView()
        self.txtview.set_editable(False)
        self.txtview.set_size_request(-1, 150)
        self.buf = self.txtview.get_buffer()

        self.bold_tag = self.buf.create_tag("bold", weight=Pango.Weight.BOLD)
        self.large_size_tag = self.buf.create_tag("large_fontsize", size=14*Pango.SCALE)

        self._title_line = 0

        self._ch_analyze_line = 1
        self._integral_line = 2
        self._area_line = 3
        self._mean_line = 4
        self._fwhm_line = 5
        self._resol_line = 6

        self.set_title("Analyze results")
        
        
    def set_text(self, text):
        self.buf.set_text(text)


    def set_analyze_peak(self, analyze, ch_to_phys):
        self.clr_buf()
        self.set_edges(analyze.x_l, analyze.x_r)
        self.set_integral(analyze.integral, analyze.integral_err)
        self.set_area(analyze.area, analyze.area_err)
        self.set_mean(analyze.mean, analyze.mean_err, ch_to_phys)
        self.set_fwhm(analyze.fwhm, analyze.fwhm_err, ch_to_phys)
        self.set_resol(analyze.resol)


    def clr_buf(self):
        start_iter = self.buf.get_start_iter()
        end_inter = self.buf.get_end_iter()
        self.buf.delete(start_iter, end_inter)


    def set_title(self, txt):
        self._insert_txt_at_line(txt, self._title_line)
        self._apply_tag_at_line_offset("bold", self._title_line, len(txt))
        self._apply_tag_at_line_offset("large_fontsize", self._title_line, len(txt))

        
    def set_edges(self, x_l, x_r):
        txt = "Edges: {:d} - {:d} ch\n".format(x_l, x_r)
        self._insert_txt_at_line(txt, self._ch_analyze_line)
        self._apply_tag_at_line_offset("bold", self._ch_analyze_line, len("Edges"))
        self._apply_tag_at_line_offset("large_fontsize", self._ch_analyze_line, len(txt) - 1)


    def set_integral(self, integral_val, integral_err):
        txt = "Integral: {:.1f}\u00b1{:.1f} k\n".format(integral_val / 1000, integral_err / 1000)
        self._insert_txt_at_line(txt, self._integral_line)
        self._apply_tag_at_line_offset("bold", self._integral_line, len("Integral"))
        self._apply_tag_at_line_offset("large_fontsize", self._integral_line, len(txt) - 1)


    def set_area(self, area_val, area_err):
        txt = "Area: {:.1f}\u00b1{:.1f} k\n".format(area_val / 1000, area_err / 1000)
        self._insert_txt_at_line(txt, self._area_line)
        self._apply_tag_at_line_offset("bold", self._area_line, len("Area"))
        self._apply_tag_at_line_offset("large_fontsize", self._area_line, len(txt) - 1)

        
    def set_mean(self, mean_val, mean_err, ch_to_phys):
        txt = "Position: {:.1f}\u00b1{:.1f} ch ({:.1f} keV)\n".format(mean_val, mean_err, ch_to_phys(mean_val))
        self._insert_txt_at_line(txt, self._mean_line)
        self._apply_tag_at_line_offset("bold", self._mean_line, len("Position"))
        self._apply_tag_at_line_offset("large_fontsize", self._mean_line, len(txt) - 1)


    def set_fwhm(self, fwhm_val, fwhm_err, ch_to_phys):
        txt = "FWHM: {:.1f}\u00b1{:.1f} ch ({:.1f} keV)\n".format(fwhm_val, fwhm_err, ch_to_phys(fwhm_val))
        self._insert_txt_at_line(txt, self._fwhm_line)
        self._apply_tag_at_line_offset("bold", self._fwhm_line, len("FWHM"))
        self._apply_tag_at_line_offset("large_fontsize", self._fwhm_line, len(txt) - 1)


    def set_resol(self, resol_val):
        txt = "Resolution (FWHM): {:.1f} %\n".format(100.0*resol_val)
        self._insert_txt_at_line(txt, self._resol_line)
        self._apply_tag_at_line_offset("bold", self._resol_line, len("Resolution (FWHM)"))
        self._apply_tag_at_line_offset("large_fontsize", self._resol_line, len(txt) - 1)


    def _insert_txt_at_line(self, txt, line):
        start_iter = self.buf.get_iter_at_line(line)
        self.buf.insert(start_iter, txt, -1)
        

    def _apply_tag_at_line_offset(self, name_tag, start_line, offset):
        start_iter = self.buf.get_iter_at_line(start_line - 1)
        if offset == -1:
            end_iter = self.buf.get_iter_at_line(start_line)
        else:
            end_iter = self.buf.get_iter_at_line_offset(start_line - 1, offset)
        self.buf.apply_tag_by_name(name_tag, start_iter, end_iter)



class Calc_view_t(Calc_view_en):
    def __init__(self):
        Calc_view_en.__init__(self)

        self._exp_eq_line = 1
        self._exp_A_line = 2
        self._exp_tau_line = 3
        self._exp_B_line = 4

        
    def set_analyze_exp(self, analyze, ch_to_ns):
        super().clr_buf()
        self.set_exp_eq()
        self.set_exp_A(analyze.A)
        self.set_exp_tau(analyze.tau, ch_to_ns)

    def set_exp_eq(self):
        txt = u"y(t) = A*exp(-t/\u03c4) + B\n"
        super()._insert_txt_at_line(txt, self._exp_eq_line)

        super()._apply_tag_at_line_offset("bold", self._exp_eq_line, len(txt) - 1)
        #super()._apply_tag_at_line_offset("large_fontsize", self._exp_eq_line, len(txt) - 1)

    def set_exp_A(self, A):
        txt = "A = {:.0e}\n".format(A)
        super()._insert_txt_at_line(txt, self._exp_A_line)

        self._apply_tag_at_line_offset("bold", self._exp_A_line, len("A"))
        super()._apply_tag_at_line_offset("large_fontsize", self._exp_A_line, len(txt) - 1)


    def set_exp_tau(self, tau, ch_to_ns):
        txt = u"\u03c4 = {:.2f} ch ({:.3f} ns)\n".format(tau, ch_to_ns(tau))
        super()._insert_txt_at_line(txt, self._exp_tau_line)

        self._apply_tag_at_line_offset("bold", self._exp_tau_line, len(u"\u03c4"))
        super()._apply_tag_at_line_offset("large_fontsize", self._exp_tau_line, len(txt) - 1)

    def set_B_line(self, B):
        txt = "B = {:.1f}".format(B)
        super()._insert_txt_at_line(txt, self._exp_B_line)

    def set_mean(self, mean_val, mean_err, ch_to_phys):
        txt = "Position: {:.1f}\u00b1{:.1f} ({:.1f} ns)\n".format(mean_val, mean_err, ch_to_phys(mean_val))
        self._insert_txt_at_line(txt, self._mean_line)
        self._apply_tag_at_line_offset("bold", self._mean_line, len("Position"))
        self._apply_tag_at_line_offset("large_fontsize", self._mean_line, len(txt) - 1)

    def set_fwhm(self, fwhm_val, fwhm_err, ch_to_phys):
        txt = "FWHM: {:.1f}\u00b1{:.1f} ({:.3f} ns)\n".format(fwhm_val, fwhm_err, ch_to_phys(fwhm_val))
        self._insert_txt_at_line(txt, self._fwhm_line)
        self._apply_tag_at_line_offset("bold", self._fwhm_line, len("FWHM"))
        self._apply_tag_at_line_offset("large_fontsize", self._fwhm_line, len(txt) - 1)

        

class Analyze_exp_t():
    def __init__(self, t_spk, x_l, x_r):
        self.t_spk = t_spk
        self.x_l = x_l
        self.x_r = x_r
        self.exp_right = self.is_exp_right()

        self.bg = self.calc_bg()

        self.curve_exp = self.approx_exp()

    def calc_bg(self):
        avg_bins = 8
        if self.exp_right:
            avg = np.sum(self.t_spk[self.x_r-avg_bins//2:self.x_r+avg_bins//2]) / avg_bins
            std = np.std(self.t_spk[self.x_r-avg_bins//2:self.x_r+avg_bins//2]) / avg_bins
        else:
            avg = np.sum(self.t_spk[self.x_l-avg_bins//2:self.x_l+avg_bins//2]) / avg_bins
            std = np.std(self.t_spk[self.x_l-avg_bins//2:self.x_l+avg_bins//2]) / avg_bins

        print("std / avg = {:.6f}".format(std / avg))
        if std / avg > 0.1:
            return np.array((self.x_r - self.x_l) * [avg], dtype=np.int64)
        else:
            return np.array((self.x_r - self.x_l) * [0], dtype=np.int64)


    """
     Approximate exponent with A * exp (-t / tau)
    """
    def approx_exp(self):
        x = np.arange(self.x_l, self.x_r)
        y = np.log(np.clip(self.t_spk[self.x_l:self.x_r] - self.bg,
                           1,
                           np.amax(self.t_spk[self.x_l:self.x_r])))

        p1, p0 = np.polyfit(x, y, 1) # y = p1 * x + p0
        self.tau = -1.0 / p1
        self.A = np.exp(p0)
        self.chi2 = np.sum((y - p1 * x - p0)**2 / (p1 * x + p0))
        print("chi2 of fit = {:.4f}".format(self.chi2))
        
        return np.array([self.A * np.exp(-i / self.tau) for i in np.arange(self.x_l, self.x_r)]) + self.bg


    """
     True if exp is right (i.e. |\)
     False if exp is left (i.e. /|)
    """
    def is_exp_right(self):
        l_avg_bins = np.sum(self.t_spk[self.x_l - 2:self.x_l + 3]) / 5.0
        r_avg_bins = np.sum(self.t_spk[self.x_r - 2:self.x_r + 3]) / 5.0

        return l_avg_bins > r_avg_bins



class Create_UI(Gtk.Window):
    def __init__(self, dir_name="-", consts={}):
        ###Calibration energy spectra
        self.calibr_en = Calibr_en()

        ###Consts (time to ch)
        self.consts = consts

        Gtk.Window.__init__(self)
        self.maximize()
        self.set_title(path.basename(path.abspath(dir_name)) + " - histopac")
        self.connect("delete-event", self.main_quit)
        self.set_icon_from_file(path.join(wd_script, "./icon.svg"))

        box_main = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        self.add(box_main)

        stack = Gtk.Stack()
        stack_switch = Gtk.StackSwitcher()
        stack_switch.set_stack(stack)
        box_main.pack_start(stack_switch, False, False, 0)
        box_main.pack_start(stack, True, True, 0)

        hbox_en = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
        hbox_t = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)

        stack.add_titled(hbox_en, "EN stack", "EN histo")
        stack.add_titled(hbox_t, "T stack", "T histo")

        self.fig_en = Figure(figsize=(5, 4), dpi=100, tight_layout=True)
        self.canvas_en = FigureCanvas(self.fig_en)
        self.canvas_en.mpl_connect("motion_notify_event", self.motion_mpl_en)
        self.canvas_en.mpl_connect("button_press_event", self.press_btn_mpl_en)
        
        self.ax_en = self.fig_en.add_subplot(111)
        self.ax_en.autoscale(False, "both", True)
        
        nav_toolbar_en = NavigationToolbar(self.canvas_en, self)
        
        vbox_mpl_en = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        vbox_mpl_en.pack_start(self.canvas_en, True, True, 0)
        vbox_mpl_en.pack_start(nav_toolbar_en, False, False, 0)
        
        grid_en = Gtk.Grid()
        grid_en.set_row_spacing(10)
        
        ###Check buttons for EN
        vbox_en_spk_chooser = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        grid_en.attach(vbox_en_spk_chooser, 0, 0, 1, 1)

        lbl = Gtk.Label()
        lbl.set_markup("<span font-weight='bold'>ENergy spectra</span>")
        vbox_en_spk_chooser.pack_start(lbl, False, False, 0)
        
        self.check_btn_en = []
        for i in range(0, det_num):
            hbox = Gtk.Box()
            check_btn = Gtk.CheckButton()
            hbox.pack_start(check_btn, False, False, 0)
            lbl = Gtk.Label()
            lbl.set_markup("<span fgcolor='{:s}'>{:s}</span>".
                           format(gui_params.det_colors[i], gui_params.en[i]))
            hbox.pack_start(lbl, False, False, 0)

            self.check_btn_en.append(check_btn)
            self.check_btn_en[-1].set_active(True)
            self.check_btn_en[-1].connect("toggled", self.toggle_check_btn_en, gui_params.en[i])
            vbox_en_spk_chooser.pack_start(hbox, False, False, 0)
    
        vbox_calc_en = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        txtview_width = 270

        self.calc_view_en = Calc_view_en() 
        self.calc_view_en.txtview.set_size_request(txtview_width, -1)
        vbox_calc_en.pack_start(self.calc_view_en.txtview, False, False, 0)
        grid_en.attach(vbox_calc_en, 0, 1, 1, 1)
        
        ##add buttons_cntrl
        grid_btns_cntrl = Gtk.Grid()
        btn_calibr_en = Gtk.Button(label="EN Calibraion")
        btn_calibr_en.connect("clicked", self.click_btn_calibr_en)
        btn_analyze_en = Gtk.Button(label="Analyze")
        btn_analyze_en.connect("clicked", self.click_btn_analyze_en)
        btn_clr_analyze_en = Gtk.Button(label="Clear Analyze")
        btn_clr_analyze_en.connect("clicked", self.click_btn_clr_analyze_en)
        btn_set_lwin_en = Gtk.Button(label="Set Left Win")
        btn_set_lwin_en.connect("clicked", self.click_btn_set_lwin_en)
        btn_set_rwin_en = Gtk.Button(label="Set Right Win")
        btn_set_rwin_en.connect("clicked", self.click_btn_set_rwin_en)
        btn_show_wins_en = Gtk.Button(label="Show Wins")
        btn_show_wins_en.connect("clicked", self.click_btn_en_show_wins)
        btn_log_scale_en = Gtk.Button(label="Log Scale")
        btn_log_scale_en.connect("clicked", self.click_btn_log_scale_en)
        self.combobox_isotopes = Gtk.ComboBoxText()
        ref_spk = Ref_spk()
        ref_spk.set_path(refspk_fname)
        ref_spk.read_spk()
        for k in ref_spk.trans.keys():
            self.combobox_isotopes.append_text(k)
        self.combobox_isotopes.set_active(0)
        btn_show_ref_spk = Gtk.Button(label="Show Ref")
        btn_show_ref_spk.connect("clicked", self.click_btn_show_ref_spk)
        
        grid_btns_cntrl.attach(btn_calibr_en, 0, 0, 1, 1)
        grid_btns_cntrl.attach(btn_analyze_en, 0, 1, 1, 1)
        grid_btns_cntrl.attach(btn_clr_analyze_en, 1, 1, 1, 1)
        grid_btns_cntrl.attach(btn_set_lwin_en, 0, 2, 1, 1)
        grid_btns_cntrl.attach(btn_set_rwin_en, 1, 2, 1, 1)
        grid_btns_cntrl.attach(btn_show_wins_en, 0, 3, 1, 1)
        grid_btns_cntrl.attach(btn_log_scale_en, 1, 3, 1, 1)
        grid_btns_cntrl.attach(self.combobox_isotopes, 0, 4, 1, 1)
        grid_btns_cntrl.attach(btn_show_ref_spk, 1, 4, 1, 1)
        
        grid_en.attach(grid_btns_cntrl, 0, 2, 1, 1)
        
        ##add entry_pointer
        vbox_ptr_en = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        self.entry_ptr_en = Gtk.Entry()
        self.entry_ptr_en.editable = False
        vbox_ptr_en.pack_start(self.entry_ptr_en, False, False, 0)
        grid_en.attach(vbox_ptr_en, 0, 3, 1, 1)
                       
        hbox_en.pack_start(vbox_mpl_en, True, True, 0)
        hbox_en.pack_start(grid_en, False, False, 5)

        ###Figs t spk###
        self.fig_t = Figure(figsize=(5, 4), dpi=100, tight_layout=True)
        self.canvas_t = FigureCanvas(self.fig_t)
        self.canvas_t.mpl_connect("motion_notify_event", self.motion_mpl_t)
        self.canvas_t.mpl_connect("button_press_event", self.press_btn_mpl_t)
        self.canvas_t.mpl_connect("button_release_event", self.resize_mpl_t)

        self.ax_t = self.fig_t.add_subplot(111)
        self.ax_t.autoscale(False, "both", True)
        
        nav_toolbar_t = NavigationToolbar(self.canvas_t, self)
        
        vbox_mp_t = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        vbox_mp_t.pack_start(self.canvas_t, True, True, 0)
        vbox_mp_t.pack_start(nav_toolbar_t, False, False, 0)
        
        grid_t = Gtk.Grid()
        grid_t.set_row_spacing(10)

        #Check buttons for T
        vbox_t_spk_chooser = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        grid_t.attach(vbox_t_spk_chooser, 0, 0, 1, 1)

        lbl = Gtk.Label()
        lbl.set_markup("<span font-weight='bold'>Time spectra</span>")
        vbox_t_spk_chooser.pack_start(lbl, False, False, 0)

        self.check_btn_t = []
        grid_check_btn_t = Gtk.Grid()
        grid_check_btn_t.set_row_spacing(5)
        for i in range(t_spk_num):
            hbox = Gtk.Box()
            check_btn = Gtk.CheckButton()
            hbox.pack_start(check_btn, False, False, 0)
            lbl = Gtk.Label()
            lbl.set_markup("<span fgcolor='{:s}'>{:s}</span> ".
                           format(gui_params.t_spk_colors[i], gui_params.t[i]))
            hbox.pack_start(lbl, False, False, 0)

            self.check_btn_t.append(check_btn)
            self.check_btn_t[-1].set_active(True)
            self.check_btn_t[-1].set_active(True)
            self.check_btn_t[-1].connect("toggled", self.toggle_check_btn_t, gui_params.t[i])
            grid_check_btn_t.attach(hbox, i % 2, i / 2, 1, 1)

        self.check_btn_tsum = []
        for i in range(t_spk_num // 2):
            hbox = Gtk.Box()
            check_btn = Gtk.CheckButton()
            hbox.pack_start(check_btn, False, False, 0)
            lbl = Gtk.Label()
            lbl.set_markup("<span fgcolor='{:s}'>{:s}</span>".
                           format(gui_params.t_spk_colors[2*i], gui_params.tsum[i]))
            hbox.pack_start(lbl, False, False, 0)

            self.check_btn_tsum.append(check_btn)
            self.check_btn_tsum[-1].set_active(False)
            self.check_btn_tsum[-1].connect("toggled", self.toggle_check_btn_t, gui_params.tsum[i])
            grid_check_btn_t.attach(hbox, 2, i, 1, 1)
            
                        
        vbox_t_spk_chooser.pack_start(grid_check_btn_t, False, False, 0)

        vbox_calc_t = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)

        self.calc_view_t = Calc_view_t()
        self.calc_view_t.txtview.set_size_request(txtview_width, -1)
        vbox_calc_t.pack_start(self.calc_view_t.txtview, False, False, 0)
        grid_t.attach(vbox_calc_t, 0, 1, 1, 1)
        
        grid_btns_cntrl = Gtk.Grid()
        btn_analyze_peak_t = Gtk.Button(label="Peak Analyze")
        btn_analyze_peak_t.connect("clicked", self.click_btn_analyze_peak_t)
        btn_analyze_exp_t = Gtk.Button(label="Exp Analyze")
        btn_analyze_exp_t.connect("clicked", self.click_btn_analyze_exp_t)
        btn_clr_analyze_t = Gtk.Button(label="Clear Analyze")
        btn_clr_analyze_t.connect("clicked", self.click_btn_clr_analyze_t)
        btn_log_scale_t = Gtk.Button(label="Log Scale")
        btn_log_scale_t.connect("clicked", self.click_btn_log_scale_t)
        btn_sum_t = Gtk.Button(label="Sum")
        btn_sum_t.connect("clicked", self.click_btn_sum_t)
        
        grid_btns_cntrl.attach(btn_analyze_peak_t, 0, 0, 1, 1)
        grid_btns_cntrl.attach(btn_analyze_exp_t, 1, 0, 1, 1)
        grid_btns_cntrl.attach(btn_clr_analyze_t, 0, 1, 1, 1)
        grid_btns_cntrl.attach(btn_log_scale_t, 0, 2, 1, 1)
        grid_btns_cntrl.attach(btn_sum_t, 1, 2, 1, 1)
                               
        grid_t.attach(grid_btns_cntrl, 0, 2, 1, 1)

        vbox_ptr_t = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        self.entry_ptr_t = Gtk.Entry()
        self.entry_ptr_t.editable = False
        vbox_ptr_t.pack_start(self.entry_ptr_t, False, False, 0)
        grid_t.attach(vbox_ptr_t, 0, 3, 1, 1)

        self.calibr_t = Calibr_t(self.consts["T_SCALE"])
        
        if self.calibr_t.ns_per_ch > 1.0:
            lbl_ns_per_ch = Gtk.Label()
            lbl_ns_per_ch.set_use_markup(True)
            lbl_ns_per_ch.set_markup("<span font='20'>{:.3f} ns/ch</span>".
                                     format(self.calibr_t.ns_per_ch))
        else:
            lbl_ns_per_ch = Gtk.Label()
            lbl_ns_per_ch.set_use_markup(True)
            lbl_ns_per_ch.set_markup("<span font='20'>{:.3f} ps/ch</span>".
                                     format(1e3 * self.calibr_t.ns_per_ch))
        grid_t.attach(lbl_ns_per_ch, 0, 4, 1, 1)
            
        hbox_t.pack_start(vbox_mp_t, True, True, 0)
        hbox_t.pack_start(grid_t, False, False, 5)

        
    def main_quit(self, event, data):
        logging.info("quit | data = {}".format(data))
        #save params in file
        #save check EN and T btns pos
        
        Gtk.main_quit()
        
        
    def count_act_check_btns_en(self):
        num_act_btns = 0
        btn_ind = -1
        for i in range(det_num):
            if self.check_btn_en[i].get_active():
                num_act_btns += 1 
                btn_ind = i

        return num_act_btns, btn_ind

        
    def motion_mpl_en(self, event):
        txt = ""
        num_act_btns, btn_ind = self.count_act_check_btns_en()

        if num_act_btns == 1:
            try:
                x = int(round(event.xdata))
                y = self.en_spk[btn_ind][x]
                txt = "{:d} <{:.0f}> {:d}".format(x, self.calibr_en.ch_to_keV(btn_ind, x), y)
            except TypeError:
                txt = "Out of range"

            self.entry_ptr_en.set_text(txt)


    def press_btn_mpl_en(self, event):
        r_mouse_btn = 3
        num_act_btns, btn_ind = self.count_act_check_btns_en()

        if num_act_btns == 1:
            if event.button == r_mouse_btn:
                try:
                    self.vlines_en

                    if len(self.vlines_en) == 2:
                        self._clr_vlines_en()
                except AttributeError:
                    self.vlines_en = []
                    self.x_vlines_en = []
                    self.txt_vlines_en = []
                     
                x = int(round(event.xdata))
                y = self.en_spk[btn_ind][x]

                self.x_vlines_en.append(x)
                self.vlines_en.append(self.ax_en.vlines(x=x,
                                                        ymin=0,
                                                        ymax=y,
                                                        color="black"))
                #text near vlines

                text_hshift = {2000: 60, 1000: 30, 500: 15, 100: 10, 50: 5, 10: 0}
                x_l, x_r = self.ax_en.get_xlim()
                for k in text_hshift.keys():
                    if k < (x_r - x_l):
                        hshift = text_hshift[k]
                self.txt_vlines_en.append(self.ax_en.text(x = x - hshift,
                                                          y = 0.9*y,
                                                          s = "{:.0f}".format(x),
                                                          rotation = 90,
                                                          fontsize = 16))
                
                self.canvas_en.draw()
        
                                   
    def toggle_check_btn_en(self, btn, name):        
        if name in gui_params.en:
            ind = gui_params.en.index(name)
            btn = self.check_btn_en[ind]
            if btn.get_active():
                self.plot_histo_name(name)
            else:
                self.clr_histo_name(name)
                

    def click_btn_calibr_en(self, btn):
        dialog = Dialog_calibr_en(self, self.calibr_en)
        response = dialog.run()

        if response == Gtk.ResponseType.OK:
            #get data from entry_en and entry_ch
            ch = [[-1.0, -1.0] for i in range(det_num)]
            en = [[-1.0, -1.0] for i in range(det_num)]
            for i in range(det_num):
                for j in range(2):
                    ch[i][j] = float(dialog.entry_ch[i].get_text().split(', ')[j])
                    en[i][j] = float(dialog.entry_en[i].get_text().split(', ')[j])

            logging.info("\nch = {}\nen = {}".format(ch, en))        
            #save data in obj and file
            self.calibr_en.set_ch_en_from_input(ch, en)
            self.calibr_en.save_file_ch_en(ini["calibr_en_path"])
            #calc k, b
            self.calibr_en.calc_k_b()
        elif response == Gtk.ResponseType.CANCEL:
            None

        dialog.destroy()


    def click_btn_analyze_en(self, btn):
        num_act_btns, btn_ind = self.count_act_check_btns_en()

        try:
            self._clr_analyze_en()
        except AttributeError:
            logging.error("AttributeError in click_btn_analyze_en()")
            None
            
        if num_act_btns == 1:
            x_l = min(self.x_vlines_en)
            x_r = max(self.x_vlines_en)

            analyze = Analyze_peak(self.en_spk[btn_ind], x_l, x_r)
            self.analyze_curve_en = self.ax_en.fill_between(range(x_l, x_r+1),
                                                            self.en_spk[btn_ind][x_l:x_r+1],
                                                            color="#42f4ee")
            self.analyze_mean_en = self.ax_en.vlines(analyze.mean,
                                                     0,
                                                     self.en_spk[btn_ind][int(analyze.mean)])
            self.analyze_fwhm_en = self.ax_en.hlines(analyze.fwhm_y + np.mean(analyze.bg),
                                                     x_l + analyze.fwhm_ch_l,
                                                     x_l + analyze.fwhm_ch_r)
            
            self.canvas_en.draw()

            def ch_to_phys(ch):
                return self.calibr_en.ch_to_keV(btn_ind, ch)

            self.calc_view_en.set_analyze_peak(analyze, ch_to_phys)


    def click_btn_clr_analyze_en(self, btn):
        try:
            self._clr_analyze_en()
        except AttributeError:
            logging.error("AttributeError in click_btn_clr_analyze_en()")
            None

        self.canvas_en.draw()
            
                
    def click_btn_set_lwin_en(self, btn):
        if len(self.x_vlines_en) != 2:
            None
        else:
            num_act_btns, btn_ind = self.count_act_check_btns_en()
            
            self.x_vlines_en.sort()
            #cfg lef wins set
            self.cfg["en_range"][btn_ind][0] = self.x_vlines_en[0]
            self.cfg["en_range"][btn_ind][1] = self.x_vlines_en[1]
            #save cfg to file
            if glob_vukap_ver == 1.0:
                save_cfg(self.cfg, ini["cfg_path"])
            elif glob_vukap_ver == 2.0:
                save_hdr_ini2(self.cfg, ini["hdr_path"])

            self._clr_vlines_en()
            self.canvas_en.draw()

            
    def click_btn_set_rwin_en(self, btn):
        if len(self.x_vlines_en) != 2:
            None
        else:
            num_act_btns, btn_ind = self.count_act_check_btns_en()
            
            self.x_vlines_en.sort()
            #cfg lef wins set
            self.cfg["en_range"][btn_ind][2] = self.x_vlines_en[0]
            self.cfg["en_range"][btn_ind][3] = self.x_vlines_en[1]
            #save cfg to file
            if glob_vukap_ver == 1.0:
                save_cfg(self.cfg, ini["cfg_path"])
            elif glob_vukap_ver == 2.0:
                save_hdr_ini2(self.cfg, ini["hdr_path"])
            
            self._clr_vlines_en()
            self.canvas_en.draw()

            
    def click_btn_en_show_wins(self, btn):
        if btn.get_label() == "Show Wins":
            num_act_btns, btn_ind = self.count_act_check_btns_en()

            if num_act_btns == 1:
                btn.set_label("Hide Wins")
                
                lwin_l, lwin_r = self.cfg["en_range"][btn_ind][0], self.cfg["en_range"][btn_ind][1]
                rwin_l, rwin_r = self.cfg["en_range"][btn_ind][2], self.cfg["en_range"][btn_ind][3]
                
                self.lines_lwin, = self.ax_en.bar(lwin_l,
                                                  self.y_en_max0,
                                                  lwin_r - lwin_l,
                                                  color = "#f4be41",
                                                  alpha = 0.4,
                                                  align = "edge")
            
                self.lines_rwin, = self.ax_en.bar(rwin_l,
                                                  self.y_en_max0,
                                                  rwin_r - rwin_l,
                                                  color = "#d3f441",
                                                  alpha = 0.4,
                                                  align = "edge")
                
        else:
            btn.set_label("Show Wins")
                        
            self.lines_lwin.remove()
            self.lines_rwin.remove()

        self.canvas_en.draw()


    def click_btn_log_scale_en(self, btn):
        if btn.get_label() == "Log Scale":
            btn.set_label("Linear Scale")
            self.ax_en.set_yscale("log")
            bottom, top = self.ax_en.get_ylim()
            self.ax_en.set_ylim(bottom, top)
        else:
            btn.set_label("Log Scale")
            self.ax_en.set_yscale("linear")

        self.canvas_en.draw()
        
    def click_btn_show_ref_spk(self, btn):
        num_act_btns, btn_ind = self.count_act_check_btns_en()
        if num_act_btns == 1:
            if btn.get_label() == "Show Ref":
                isotope = self.combobox_isotopes.get_active_text()

                try:
                    ref_spk = Ref_spk()
                    ref_spk.set_path(refspk_fname)
                    ref_spk.read_spk()
                    ref_spk.fill_isotope_trans(isotope)
                    
                    y_min, y_max = self.ax_en.get_ylim()

                    self.vlines_transitions = []
                    for trans in ref_spk.transitions:
                        ch = self.calibr_en.keV_to_ch(btn_ind, trans[0])
                        self.vlines_transitions.append(self.ax_en.vlines(ch,
                                                                         0.02 * y_max,
                                                                         0.98 * trans[1] * y_max,
                                                                         colors="black",
                                                                         linestyles="dotted",
                                                                         linewidths=2))
                    
                    btn.set_label("Hide Ref")
                    self.canvas_en.draw()
                    return
                except KeyError:
                    return

        if btn.get_label() == "Hide Ref":
            for line in self.vlines_transitions:
                line.remove()
                del line

            btn.set_label("Show Ref")        
            self.canvas_en.draw()
                

    def count_act_check_btns_t(self):
        num_act_btns = 0
        btn_ind = -1
        for i in range(t_spk_num):
            if self.check_btn_t[i].get_active():
                num_act_btns += 1 
                btn_ind = i

        return num_act_btns, btn_ind

    def count_act_check_btns_tsum(self):
        num_act_btns = 0
        btn_ind = -1
        for i in range(tsum_spk_num):
            if self.check_btn_tsum[i].get_active():
                num_act_btns += 1 
                btn_ind = i 

        return num_act_btns, btn_ind


    def motion_mpl_t(self, event):
        txt = ""
        num_act_btns, btn_ind = self.count_act_check_btns_t()
        if num_act_btns == 1:
            try:
                x = int(round(event.xdata))
                y = self.t_spk[btn_ind][x]
                txt = "{:d} ({:.3f} ns) {:d}".format(x, x*self.calibr_t.ns_per_ch, y)
            except TypeError:
                txt = "Out of range"
            self.entry_ptr_t.set_text(txt)
            return 1

        num_act_btns, btn_ind = self.count_act_check_btns_tsum()
        if num_act_btns == 1:
            try:
                x = int(round(event.xdata))
                y = self.tsum_spk[btn_ind][x]
                txt = "{:d} ({:.3f} ns) {:d}".format(x, x*self.calibr_t.ns_per_ch, y)
            except TypeError:
                txt = "Out of range"
            self.entry_ptr_t.set_text(txt)
            return 1


    
    def press_btn_mpl_t(self, event):
        r_mouse_btn = 3
        num_act_btns_t, btn_ind_t = self.count_act_check_btns_t()
        num_act_btns_tsum, btn_ind_tsum = self.count_act_check_btns_tsum()
        
        if (num_act_btns_t == 1 and num_act_btns_tsum == 0) or (num_act_btns_t == 0 and num_act_btns_tsum == 1):
            if event.button == r_mouse_btn:
                try:
                    self.vlines_t
                    
                    if len(self.vlines_t) == 2:
                        self._clr_vlines_t()
                except AttributeError:
                    self.vlines_t = []
                    self.x_vlines_t = []
                    self.txt_vlines_t = []
                    
                x = int(round(event.xdata))
                if num_act_btns_t:
                    y = self.t_spk[btn_ind_t][x]
                else:
                    y = self.tsum_spk[btn_ind_tsum][x]

                self.x_vlines_t.append(x)

                self.vlines_t.append(self.ax_t.vlines(x = x,
                                                      ymin = 0,
                                                      ymax = y,
                                                      color = "black"))
                text_hshift = {2000: 60, 1000: 30, 500: 15, 100: 10, 50: 5, 10: 0}
                x_l, x_r = self.ax_t.get_xlim()
                for k in text_hshift.keys():
                    if k < (x_r - x_l):
                        hshift = text_hshift[k]
                self.txt_vlines_t.append(self.ax_t.text(x = x - hshift,
                                                        y = 0.9*y,
                                                        s = "{:.0f}".format(x),
                                                        rotation = 90,
                                                        fontsize = 16))
                self.canvas_t.draw()


    def resize_mpl_t(self, event):
        l_mouse_btn = 1
        if event.button == l_mouse_btn:
            try:
                self.vlines_t
                    
                if len(self.vlines_t) == 2:
                    self._clr_vlines_t()
            except AttributeError:
                self.vlines_t = []
                self.x_vlines_t = []
                self.txt_vlines_t = [] 
        
        
    def click_btn_analyze_peak_t(self, btn):
        logging.info("Click btn analyze peak")
        num_act_btns_t, btn_ind_t = self.count_act_check_btns_t()
        num_act_btns_tsum, btn_ind_tsum = self.count_act_check_btns_tsum()
        
        try:
            self._clr_analyze_t()
        except AttributeError:
            logging.error("AttributeError in click_btn_analyze_exp_t()")
            
        if (num_act_btns_t == 1 and num_act_btns_tsum == 0) or (num_act_btns_t == 0 and num_act_btns_tsum == 1):
            x_l = min(self.x_vlines_t)
            x_r = max(self.x_vlines_t)

            
            if num_act_btns_t:
                analyze = Analyze_peak(self.t_spk[btn_ind_t], x_l, x_r)
                self.analyze_curve_peak_t = self.ax_t.fill_between(np.arange(x_l, x_r+1),
                                                                   self.t_spk[btn_ind_t][x_l:x_r+1],
                                                                   color="#42f4ee")
                self.analyze_mean_t = self.ax_t.vlines(analyze.mean,
                                                         0,
                                                         self.t_spk[btn_ind_t][int(analyze.mean)])
                self.analyze_fwhm_t = self.ax_t.hlines(analyze.fwhm_y + np.mean(analyze.bg),
                                                         x_l + analyze.fwhm_ch_l,
                                                         x_l + analyze.fwhm_ch_r)
            else:
                analyze = Analyze_peak(self.tsum_spk[btn_ind_tsum], x_l, x_r)
                self.analyze_curve_peak_t = self.ax_t.fill_between(np.arange(x_l, x_r+1),
                                                                self.tsum_spk[btn_ind_tsum][x_l:x_r+1],
                                                                color="#42f4ee")
                self.analyze_mean_t = self.ax_t.vlines(analyze.mean,
                                                       0,
                                                       self.tsum_spk[btn_ind_tsum][int(analyze.mean)])
                self.analyze_fwhm_t = self.ax_t.hlines(analyze.fwhm_y + np.mean(analyze.bg),
                                                       x_l + analyze.fwhm_ch_l,
                                                       x_l + analyze.fwhm_ch_r)
                
            self.canvas_t.draw()

            def ch_to_phys(ch):
                try:
                    return self.calibr_t.ch_to_ns(ch)
                except AttributeError:
                    logging.error("AttributeError, ch is {:d}".format(ch))
                    return -1

            self.calc_view_t.set_analyze_peak(analyze, ch_to_phys)
            
        
    def click_btn_analyze_exp_t(self, btn):
        logging.info("Click btn analyze exp")
        num_act_btns_t, btn_ind_t = self.count_act_check_btns_t()
        num_act_btns_tsum, btn_ind_tsum = self.count_act_check_btns_tsum()
        try:
            self._clr_analyze_t()
        except AttributeError:
            logging.error("AttributeError in click_btn_analyze_exp_t()")
            
        if (num_act_btns_t == 1 and num_act_btns_tsum == 0) or (num_act_btns_t == 0 and num_act_btns_tsum == 1):
            x_l = min(self.x_vlines_t)
            x_r = max(self.x_vlines_t)

            if num_act_btns_t:
                analyze = Analyze_exp_t(self.t_spk[btn_ind_t], x_l, x_r)
            else:
                analyze = Analyze_exp_t(self.tsum_spk[btn_ind_tsum], x_l, x_r)
            self.analyze_curve_exp_t = self.ax_t.plot(np.arange(analyze.x_l, analyze.x_r),
                                                      analyze.curve_exp,
                                                      c="#ecff00",
                                                      lw=2)

            self.canvas_t.draw()

            self.calc_view_t.set_analyze_exp(analyze, self.calibr_t.ch_to_ns)
            
            
    def click_btn_clr_analyze_t(self, btn):
        try:
            self._clr_analyze_t()
        except AttributeError:
            logging.error("AttributeError in click_btn_clr_analyze_t()")

        self.canvas_t.draw()


    def click_btn_log_scale_t(self, btn):
        if btn.get_label() == "Log Scale":
            btn.set_label("Linear Scale")
            self.ax_t.set_yscale("log")
            bottom, top = self.ax_t.get_ylim()
            self.ax_t.set_ylim(bottom, top)
        else:
            btn.set_label("Log Scale")
            self.ax_t.set_yscale("linear")

        self.canvas_t.draw()

        
    def click_btn_sum_t(self, btn):
        num_btn, ind_btn = self.count_act_check_btns_t()
        if num_btn == 2:
            ind_btn = []
            for i in range(t_spk_num):
                if self.check_btn_t[i].get_active(): 
                    ind_btn.append(i)
            print("Exactly 2 btns actives")
            print(ind_btn)
            sum_spk_t = np.array(self.t_spk[ind_btn[0]]) + np.array(self.t_spk[ind_btn[1]])
            self.ax_t.plot(sum_spk_t)
            self.canvas_t.draw()

            
    def set_histos(self, en_spk, t_spk):
        self.en_spk = en_spk
        self.t_spk = t_spk
        self.tsum_spk = []
        for i in range(tsum_spk_num):
            self.tsum_spk.append(np.array(t_spk[2*i]) + np.array(t_spk[2*i+1]))

        
    def set_cfg(self, cfg):
        self.cfg = cfg

        
    def set_consts(self, consts):
        self.consts = consts

        
    def plot_histo_name(self, name):
        if name in gui_params.en:
            ind = gui_params.en.index(name)
            self.en_lines[ind].set_ydata(self.en_spk[ind])
            self.canvas_en.draw()
        elif name in gui_params.t:
            ind = gui_params.t.index(name)
            self.t_lines[ind].set_ydata(self.t_spk[ind])
            self.canvas_t.draw()
        elif name in gui_params.tsum:
            ind = gui_params.tsum.index(name)
            tsum = np.array(self.t_spk[2*ind]) + np.array(self.t_spk[2*ind+1])
            self.tsum_lines[ind].set_ydata(tsum)
            self.canvas_t.draw()
            
    def clr_histo_name(self, name):
        y_zeros = [0] * histo_size
        
        for n in gui_params.en:
            ind = gui_params.en.index(n)
            btn = self.check_btn_en[ind]

            if btn.get_active() is False:
                self.en_lines[ind].set_ydata(y_zeros)

        self.canvas_en.draw()

        for n in gui_params.t:
            ind = gui_params.t.index(n)
            btn = self.check_btn_t[ind]

            if btn.get_active() is False:
                self.t_lines[ind].set_ydata(y_zeros)

        for n in gui_params.tsum:
            ind = gui_params.tsum.index(n)
            btn = self.check_btn_tsum[ind]

            if btn.get_active() is False:
                self.tsum_lines[ind].set_ydata(y_zeros)
                
        self.canvas_t.draw()

        
    def plot_all_histos(self):
        self.en_lines = []
        self.t_lines = []
        self.tsum_lines = []
        x = range(histo_size)
        
        for i in range(det_num):
            self.en_lines.append(self.ax_en.plot(x, self.en_spk[i], marker="o",
                                                 ms=3, mew=0,
                                                 color=gui_params.det_mpl_colors[i], lw=1.0)[0])
        self.set_lim_vals_en(0)     
        self.canvas_en.draw()

        for i in range(t_spk_num):
            self.t_lines.append( self.ax_t.plot(x, self.t_spk[i], marker="o",
                                                ms=5, mew=0,
                                                color=gui_params.t_spk_mpl_colors[i], lw=0)[0] )
        for i in range(tsum_spk_num):
            self.tsum_lines.append( self.ax_t.plot(x, self.tsum_spk[i], marker="o",
                                                ms=5, mew=0,
                                                color=gui_params.t_spk_mpl_colors[2*i], lw=0)[0] )

        self.set_lim_vals_t(0)
        
        self.canvas_t.draw()

        
    def set_lim_vals_en(self, flag):
        if flag == 0:
            self.x_en_min0, self.x_en_max0 = 0, histo_size
            self.y_en_min0, self.y_en_max0 = 0, 1.05 * max([max(spk) for spk in self.en_spk]) or 1
            self.ax_en.set_xlim(self.x_en_min0, self.x_en_max0)
            self.ax_en.set_ylim(self.y_en_min0, self.y_en_max0)

            
    def set_lim_vals_t(self, flag):
        if flag == 0:
            self.x_t_min0, self.x_t_max0 = 0, histo_size
            self.y_t_min0, self.y_t_max0 = 0, 1.05 * max([max(spk) for spk in self.t_spk]) or 1
            self.ax_t.set_xlim(self.x_t_min0, self.x_t_max0)
            self.ax_t.set_ylim(self.y_t_min0, self.y_t_max0)
            

    def _clr_analyze_en(self):
        try:
            self.analyze_curve_en.remove()
            self.analyze_mean_en.remove()
            self.analyze_fwhm_en.remove()
        except ValueError:
            logging.error("ValueError in _clr_analyze_en")


    def _clr_analyze_t(self):
        logging.info("_clr_analyze_t()")
        try:
            self.analyze_curve_peak_t.remove()
            self.analyze_mean_t.remove()
            self.analyze_fwhm_t.remove()
        except AttributeError:
            logging.error("AttributeError in _clr_analyze_t")
        except ValueError:
            logging.error("AttributeError in _clr_analyze_t")
            
        try:
            logging.info("before del {}".format(self.analyze_curve_exp_t))
            self.ax_t.lines.remove(self.analyze_curve_exp_t[-1])
            del self.analyze_curve_exp_t[-1]
            logging.info("after del {}".format(self.analyze_curve_exp_t))
        except ValueError:
            logging.error("ValueError in _clr_analyze_t()")
        except IndexError:
            logging.error("IndexError in _clr_analyze_t()")

        self.canvas_t.draw()
            
            
    def _clr_vlines_en(self):
        self.vlines_en[0].remove()
        self.vlines_en[1].remove()

        self.txt_vlines_en[0].remove()
        self.txt_vlines_en[1].remove()
        
        self.vlines_en = []
        self.x_vlines_en = []
        self.txt_vlines_en = []
        

    def _clr_vlines_t(self):
        self.vlines_t[0].remove()
        self.vlines_t[1].remove()

        self.txt_vlines_t[0].remove()
        self.txt_vlines_t[1].remove()
        
        self.vlines_t = []
        self.x_vlines_t = []
        self.txt_vlines_t = []
    

    def toggle_check_btn_t(self, btn, name):
        if name in gui_params.t:
            ind = gui_params.t.index(name)
            btn = self.check_btn_t[ind]
            if btn.get_active():
                self.plot_histo_name(name)
            else:
                self.clr_histo_name(name)
        elif name in gui_params.tsum:
            ind = gui_params.tsum.index(name)
            btn = self.check_btn_tsum[ind]
            if btn.get_active():
                self.plot_histo_name(name)
            else:
                self.clr_histo_name(name)
                

    def hide_all_spk_t(self):
        for i in range(t_spk_num):
            self.check_btn_t[i].set_active(False)
            self.check_btn_t[i].emit("toggled")
        for i in range(tsum_spk_num):
            self.check_btn_tsum[i].set_active(False)
            self.check_btn_tsum[i].emit("toggled")
            
                
def get_histos_from_folder(foldername="./testspk/"):
    foldername = foldername if foldername[-1] == "/" else foldername + "/"
    
    en_spk = []
    for fname in histo_en_fnames:
        with open(foldername + fname, "rb") as fd:
            fd.seek(histo_finisize, 0)
            en_spk.append(fd.read(histo_fsize))
            en_spk[-1] = unpack(str(histo_size) + "i", en_spk[-1])

    t_spk = []
    for fname in histo_t_fnames:
        with open(foldername + fname, "rb") as fd:
            fd.seek(histo_finisize, 0)
            t_spk.append(fd.read(histo_fsize))
            t_spk[-1] = unpack(str(histo_size) + "i", t_spk[-1])
            
    return en_spk, t_spk


def parse_hdr_ini2(path_ini_file):
    config = ConfigParser()
    config.read(path_ini_file)
    cfg = {}
    cfg["histo_folder"] = config["acquisition"]["out-foldername"]
    cfg["en_range"] = []
    cfg["with-signal"] = config["acquisition"]["with-signal"]
    cfg["time-acq"] = config["acquisition"]["time-acq"]
    for i in range(det_num):
        try:
            rang = [int(xs) for xs in config["en-wins-1"]["det-{:d}".format(i+1)].split()]
            cfg["en_range"].append(rang)
        except KeyError:
            rang = [int(xs) for xs in config["acquisition"]["det-{:d}".format(i+1)].split()]
            cfg["en_range"].append(rang)
    return cfg
        
def parse_dsp_ini2(path_ini_file):
    config = ConfigParser()
    config.read(path_ini_file)
    cnst = {}
    cnst["T_SCALE"] = [float(xs) for xs in config["dsp"]["t-scale"].split()]
    return cnst

def save_hdr_ini2(cfg, path_ini_file):
    config = ConfigParser()
    config["acquisition"] = {}
    config["acquisition"]["out-foldername"] = cfg["histo_folder"]
    config["acquisition"]["with-signal"] = cfg["with-signal"]
    config["acquisition"]["time-acq"] = cfg["time-acq"]
    config["en-wins-1"] = {}
    config["en-wins-2"] = {}
    for i in range(1, 3):
        for j in range(1, 5):
            config["en-wins-" + str(i)]["det-" + str(j)] = " ".join(map(str, cfg["en_range"][j-1]))
    #config["en-wins-1"]["det-1"] = " ".join(map(str, cfg["en_range"][0]))
    #config["en-wins-1"]["det-2"] = " ".join(map(str, cfg["en_range"][1]))
    #config["en-wins-1"]["det-3"] = " ".join(map(str, cfg["en_range"][2]))
    #config["en-wins-1"]["det-4"] =  " ".join(map(str, cfg["en_range"][3]))
    
    with open(path_ini_file, "w") as ini_file_fd:
        config.write(ini_file_fd)
    
def set_ini(path_ini_file, **kwargs):
    with open(path_ini_file, "w") as fp_ini:
        ini = json.load(fp_ini)

        for path_name in ini.keys():
            if path_name in kwargs.keys():
                ini[path_name] = kwargs[path_name]

def parse_ini(path_ini_file):
    ini = {}

    with open(path_ini_file, "r") as fp_ini:
        ini = json.load(fp_ini)

    return ini


def parse_cfg(path_cfg_file):
    cfg = {}
    
    with open(path_cfg_file, 'r') as cfg_file:
        #if the size of cfg.json will be larger than 2048 bytes, please change it
        config = cfg_file.read(2048)
        config_vals = json.loads(config)
        
        cfg["porog"] = config_vals["porog for detectors"]#[90, 90, 90, 90]
        cfg["delay"] = config_vals["time of coincidence [ch]"]#100
        cfg["coinc"] = 1 if config_vals["coinc mode?"] == "True" else 0#1
        cfg["acq_time"] = config_vals["time of exposition [s]"]#10
        cfg["histo_folder"] = config_vals["histo folder"]
        cfg["en_range"] = []
        for i in range(det_num):
            cfg["en_range"].append(config_vals["en_range D" + str(i + 1)])

    return cfg


def save_cfg(cfg, path_cfg_file):
    with open(path_cfg_file, 'w') as cfg_file:
        def dict_to_true_cfg_str(d):
            res = "{\n"
            res += "\t\"porog for detectors\": " + str(d["porog"]) + ",\n"
            res += "\t\"time of coincidence [ch]\": " + str(d["delay"]) + ",\n"
            res += "\t\"coinc mode?\": " + "\"" + str(d["coinc"] == 1) + "\"" + ",\n"
            res += "\t\"time of exposition [s]\": " + str(d["acq_time"]) + ",\n"
            res += "\t\"histo folder\": " + "\"" + str(d["histo_folder"]) + "\"" + ",\n"
            for i in range(0, 3):
                res += "\t\"en_range D" + str(i + 1) + "\": " + str(d["en_range"][i])  + ",\n"
            i = 3
            res += "\t\"en_range D" + str(i + 1) + "\": " + str(d["en_range"][i])  + "\n"
            res += "}"

            return res

        cfg_file.write(dict_to_true_cfg_str(cfg))


def new_save_cfg(cfg, path_cfg_file):
    with open(path_cfg_file, 'w') as cfg_file:
        json.dump(cfg,
                  cfg_file,
                  indent="    ")     


def parse_constants(path_constants_file):
    cnst = {}
    
    with open(path_constants_file, 'r') as const_file:
        #if the size of cfg.json will be larger than 2048 bytes, please change it
        consts = const_file.read(2048)
        consts_vals = json.loads(consts)
        cnst["T_SCALE"] = consts_vals["T_SCALE"]

    return cnst

        
if __name__ == "__main__":
    main()
