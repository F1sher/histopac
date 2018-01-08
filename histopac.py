#!/usr/bin/python3


from json import loads as json_loads
from json import dumps as json_dumps

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk3agg import FigureCanvasGTK3Agg as FigureCanvas
from matplotlib.backends.backend_gtk3 import NavigationToolbar2GTK3 as NavigationToolbar

import numpy as np
from struct import unpack

import gui_name


det_num = 4
t_spk_num = 12

histo_size = 4096 #(chanels)
histo_finisize = 512
histo_fsize = 16384 #(bytes)
histo_en_fnames = ["BUFKA" + str(i) + ".SPK" for i in range(1, 5)]
histo_t_fnames = ["T12", "T21"]

det_clrs = ["r", "g", "b", "m"]

cfg_fname = "./testspk/cfg_.json"


            
def main():
    ui = Create_UI()
    
    en_spk, t_spk = get_histos_from_folder("../sum_spk/res2/")

    ui.set_histos(en_spk, t_spk)
    ui.plot_all_histos()

    cfg = parse_cfg("./testspk/cfg.json")
    save_cfg(cfg, "./testspk/cfg_.json")
    ui.set_cfg(cfg)
    
    ui.show_all()
    Gtk.main()


class Analyze_en():
    def __init__(self, en_spk, x_l, x_r):
        self.en_spk = en_spk
        self.x_l = x_l
        self.x_r = x_r

        self.mean = self.calc_mean()
        self.fwhm = self.calc_fwhm()
        print("mean = {}, calc = {}".format(self.mean, self.fwhm))
        

    def calc_mean(self):
        w_sum = sum([i * y for i, y in zip(range(self.x_l, self.x_r), self.en_spk[self.x_l:self.x_r])])
        
        return w_sum / sum(self.en_spk[self.x_l:self.x_r])

    
    def calc_fwhm(self):
        max_en_spk = max(self.en_spk[self.x_l:self.x_r])
        
        for y in self.en_spk[self.x_l:self.x_r]:
            if y >= max_en_spk / 2:
                i = self.en_spk[self.x_l:self.x_r].index(y)
                break

        k_l = self.en_spk[self.x_l + i] - self.en_spk[self.x_l + i - 1]
        b_l = self.en_spk[self.x_l + i] - k_l * i

        for y in self.en_spk[self.x_l + i:self.x_r]:
            if y <= max_en_spk / 2:
                i = self.en_spk[self.x_l:self.x_r].index(y)
                break

        k_r = self.en_spk[self.x_l + i] - self.en_spk[self.x_l + i - 1]
        b_r = self.en_spk[self.x_l + i] - k_r * i

        self.fwhm_ch_l = (max_en_spk / 2 - b_l) / k_l
        self.fwhm_ch_r = (max_en_spk / 2 - b_r) / k_r
        #print("max_en_spk = {}".format(max_en_spk))
        #print("fwhm_l = {}, fwhm_r = {}".format(self.x_l + fwhm_l, self.x_l + fwhm_r))
        
        return (self.fwhm_ch_r - self.fwhm_ch_l)
        

class Calibr_en():
    def __init__(self):
        self.ch = [-1, -1]
        self.en = [-1.0, -1.0]

        self.set_ch_en_from_file("./testspk/calibr_en.json")

        self.calc_k_b()

        
    def calc_k_b(self):
        try:
            self.k = (self.en[0] - self.en[1]) / (self.ch[0] - self.ch[1])
        except ZeroDivisionError:
            self.k = 1.0
        self.b = self.en[0] - self.k * self.ch[0]

        
    def set_ch_en_from_file(self, fname):
        with open(fname, "r") as f:
            data = f.read(2048)
            vals = json_loads(data)

            for i in range(2):
                self.ch[i] = int(vals["channels"][i])
                self.en[i] = float(vals["energies"][i])


    def set_ch_en_from_input(self, ch, en):
        self.ch = ch
        self.en = en


    def save_file_ch_en(self, fname):
        with open(fname, "w+") as f:
            data = "{\n"
            data += "\t\"channels\": ["
            data +=  str(self.ch[0])
            data += ", "
            data += str(self.ch[1])
            data += "],\n"
            data += "\t\"energies\": ["
            data += str(self.en[0])
            data += ", "
            data += str(self.en[1])
            data += "]\n"
            data += "}"
            f.write(data)
        

class Dialog_calibr_en(Gtk.Dialog):
    def __init__(self, parent, calibr_en):
        Gtk.Dialog.__init__(self, "Energy calibration", parent, 0,
                            (Gtk.STOCK_OK, Gtk.ResponseType.OK,
                             Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL))

        self.set_default_size(200, 100)

        label_ch = Gtk.Label("Channels:")
        label_en = Gtk.Label("Energies:")
        
        self.entry_ch = [Gtk.Entry(), Gtk.Entry()]
        self.entry_en = [Gtk.Entry(), Gtk.Entry()]

        for i in range(2):
            self.entry_ch[i].set_text(str(calibr_en.ch[i]))
            self.entry_en[i].set_text(str(calibr_en.en[i]))
            
        self.grid = Gtk.Grid()

        self.grid.attach(label_ch, 0, 0, 1, 1)
        self.grid.attach(label_en, 1, 0, 1, 1)
        self.grid.attach(self.entry_ch[0], 0, 1, 1, 1)
        self.grid.attach(self.entry_en[0], 1, 1, 1, 1)
        self.grid.attach(self.entry_ch[1], 0, 2, 1, 1)
        self.grid.attach(self.entry_en[1], 1, 2, 1, 1)
        
        box = self.get_content_area()
        box.add(self.grid)

        self.show_all()
        

class Create_UI(Gtk.Window):
    def __init__(self):
        Gtk.Window.__init__(self, title="histopac")
        self.connect("delete-event", Gtk.main_quit)

        box_main = Gtk.Box(orientation = Gtk.Orientation.VERTICAL)
        self.add(box_main)
        
        stack = Gtk.Stack()
        stack_switch = Gtk.StackSwitcher()
        stack_switch.set_stack(stack)
        box_main.pack_start(stack_switch, False, False, 0)
        box_main.pack_start(stack, True, True, 0)
        
        hbox_en = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL)
        hbox_t = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL)

        stack.add_titled(hbox_en, "EN stack", "EN histo")
        stack.add_titled(hbox_t, "T stack", "T histo")
        
        self.fig_en = Figure(figsize = (5, 4), dpi = 100, tight_layout = True)
        
        self.canvas_en = FigureCanvas(self.fig_en)
        self.canvas_en.mpl_connect("motion_notify_event", self.motion_mpl_en)
        self.canvas_en.mpl_connect("button_press_event", self.press_btn_mpl_en)
        
        self.ax_en = self.fig_en.add_subplot(111)
        self.ax_en.autoscale(False, "both", True)
        
        nav_toolbar_en = NavigationToolbar(self.canvas_en, self)
        
        vbox_mpl_en = Gtk.Box(orientation = Gtk.Orientation.VERTICAL)
        vbox_mpl_en.pack_start(self.canvas_en, True, True, 0)
        vbox_mpl_en.pack_start(nav_toolbar_en, False, False, 0)
        ###

        self.fig_t = Figure(figsize=(5, 4), dpi=100)
        self.canvas_t = FigureCanvas(self.fig_t)

        self.calibr_en = Calibr_en()
        
        grid_en = Gtk.Grid()
        grid_en.set_row_spacing(10)
        
        #Check buttons for EN
        vbox_en_spk_chooser = Gtk.Box(orientation = Gtk.Orientation.VERTICAL)
        grid_en.attach(vbox_en_spk_chooser, 0, 0, 1, 1)

        vbox_en_spk_chooser.pack_start(Gtk.Label("EN spk"), False, False, 0)
        self.check_btn_en = []
        for i in range(0, det_num):
            self.check_btn_en.append(Gtk.CheckButton.new_with_label(gui_name.en[i]))
            self.check_btn_en[-1].set_active(True)
            self.check_btn_en[-1].connect("toggled", self.toggle_check_btn_en, gui_name.en[i])
            vbox_en_spk_chooser.pack_start(self.check_btn_en[-1], False, False, 0)

        vbox_calc_en = Gtk.Box(orientation = Gtk.Orientation.VERTICAL)
        textview_calc_en = Gtk.TextView()
        self.textbuf_calc_en = textview_calc_en.get_buffer()
        self.textbuf_calc_en.set_text("Calc results:\n - 1\n - 2")
        vbox_calc_en.pack_start(textview_calc_en, False, False, 0)
        grid_en.attach(vbox_calc_en, 0, 1, 1, 1)
        
        ##add buttons_cntrl
        grid_btns_cntrl = Gtk.Grid()
        btn_calibr_en = Gtk.Button("EN Calibraion")
        btn_calibr_en.connect("clicked", self.click_btn_calibr_en)
        btn_analyze_en = Gtk.Button("Analyze")
        btn_analyze_en.connect("clicked", self.click_btn_analyze_en)
        btn_clr_analyze_en = Gtk.Button("Clear Analyze")
        btn_clr_analyze_en.connect("clicked", self.click_btn_clr_analyze_en)
        btn_set_lwin_en = Gtk.Button("Set Left Win")
        btn_set_lwin_en.connect("clicked", self.click_btn_set_lwin_en)
        btn_set_rwin_en = Gtk.Button("Set Right Win")
        btn_set_rwin_en.connect("clicked", self.click_btn_set_rwin_en)
        btn_show_wins_en = Gtk.Button("Show Wins")
        btn_show_wins_en.connect("clicked", self.click_btn_en_show_wins)
        
        grid_btns_cntrl.attach(btn_calibr_en, 0, 0, 1, 1)
        grid_btns_cntrl.attach(btn_analyze_en, 0, 1, 1, 1)
        grid_btns_cntrl.attach(btn_clr_analyze_en, 1, 1, 1, 1)
        grid_btns_cntrl.attach(btn_set_lwin_en, 0, 2, 1, 1)
        grid_btns_cntrl.attach(btn_set_rwin_en, 1, 2, 1, 1)
        grid_btns_cntrl.attach(btn_show_wins_en, 0, 3, 1, 1)
        
        grid_en.attach(grid_btns_cntrl, 0, 2, 1, 1)
        
        ##add entry_pointer
        vbox_ptr_en = Gtk.Box(orientation = Gtk.Orientation.VERTICAL)
        self.entry_ptr_en = Gtk.Entry()
        self.entry_ptr_en.editable = False
        vbox_ptr_en.pack_start(self.entry_ptr_en, False, False, 0)
        grid_en.attach(vbox_ptr_en, 0, 3, 1, 1)
                       
        hbox_en.pack_start(vbox_mpl_en, True, True, 0)
        hbox_en.pack_start(grid_en, False, False, 5)
        
        hbox_t.pack_start(self.canvas_t, True, True, 0)

        
    def count_act_check_btns(self):
        num_act_btns = 0
        btn_ind = -1
        for i in range(det_num):
            if self.check_btn_en[i].get_active():
                num_act_btns += 1 
                btn_ind = i

        return num_act_btns, btn_ind

        
    def motion_mpl_en(self, event):
        txt = ""
        num_act_btns, btn_ind = self.count_act_check_btns()
                
        if num_act_btns == 1:
            try:
                x = int(float(event.xdata))
                y = self.en_spk[btn_ind][x]
                txt = "{:d} ({:.0f}) {:d}".format(x, self.calibr_en.k * x + self.calibr_en.b, y)
            except TypeError:
                txt = "Out of range"

        self.entry_ptr_en.set_text(txt)


    def press_btn_mpl_en(self, event):
        r_mouse_btn = 3
        num_act_btns, btn_ind = self.count_act_check_btns()

        if num_act_btns == 1:
            if event.button == r_mouse_btn:
                try:
                    self.vlines_en

                    if len(self.vlines_en) == 2:
                        self._clr_vlines_en()
                except AttributeError:
                    self.vlines_en = []
                    self.x_vlines_en = []
                     
                x = int(event.xdata)
                y = self.en_spk[btn_ind][x]

                self.x_vlines_en.append(x)
                
                self.vlines_en.append(self.ax_en.vlines(x = x,
                                                        ymin = 0,
                                                        ymax = y,
                                                        color = "black"))
                self.canvas_en.draw()
        
                                   
    def toggle_check_btn_en(self, btn, name):        
        if name in gui_name.en:
            ind = gui_name.en.index(name)
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
            ch = [int(dialog.entry_ch[0].get_text()), int(dialog.entry_ch[1].get_text())]
            en = [float(dialog.entry_en[0].get_text()), float(dialog.entry_en[1].get_text())]
            #save data to file
            self.calibr_en.set_ch_en_from_input(ch, en)
            self.calibr_en.save_file_ch_en("./testspk/calibr_en.json")
            #calc k, b
            self.calibr_en.calc_k_b()
        elif response == Gtk.ResponseType.CANCEL:
            None

        dialog.destroy()


    def click_btn_analyze_en(self, btn):
        num_act_btns, btn_ind = self.count_act_check_btns()

        try:
            self._clr_analyze_en()
        except AttributeError:
            print("AttributeError: in click_btn_analyze_en()")
            None
            
        if num_act_btns == 1:
            x_l = min(self.x_vlines_en)
            x_r = max(self.x_vlines_en)
            res = Analyze_en(self.en_spk[btn_ind], x_l, x_r)
            self.analyze_curve_en = self.ax_en.fill_between(range(x_l, x_r), self.en_spk[btn_ind][x_l:x_r],
                                                            color="#42f4ee")
            self.analyze_mean_en = self.ax_en.vlines(res.mean, 0, self.en_spk[btn_ind][int(res.mean)])
            self.analyze_fwhm_en = self.ax_en.hlines(max(self.en_spk[btn_ind][x_l:x_r]) / 2,
                                                     x_l + res.fwhm_ch_l,
                                                     x_l + res.fwhm_ch_r)
            self.canvas_en.draw()


    def click_btn_clr_analyze_en(self, btn):
        try:
            self._clr_analyze_en()
        except AttributeError:
            print("AttributeError: in click_btn_clr_analyze_en()")
            None

        self.canvas_en.draw()
            
                
    def click_btn_set_lwin_en(self, btn):
        if len(self.x_vlines_en) != 2:
            None
        else:
            num_act_btns, btn_ind = self.count_act_check_btns()
            
            self.x_vlines_en.sort()
            #cfg lef wins set
            self.cfg["en_range"][btn_ind][0] = self.x_vlines_en[0]
            self.cfg["en_range"][btn_ind][1] = self.x_vlines_en[1]
            #save cfg to file
            save_cfg(self.cfg, cfg_fname)

            self._clr_vlines_en()
            self.canvas_en.draw()

            
    def click_btn_set_rwin_en(self, btn):
        if len(self.x_vlines_en) != 2:
            None
        else:
            num_act_btns, btn_ind = self.count_act_check_btns()
            
            self.x_vlines_en.sort()
            #cfg lef wins set
            self.cfg["en_range"][btn_ind][2] = self.x_vlines_en[0]
            self.cfg["en_range"][btn_ind][3] = self.x_vlines_en[1]
            #save cfg to file
            save_cfg(self.cfg, cfg_fname)
            
            self._clr_vlines_en()
            self.canvas_en.draw()

            
    def click_btn_en_show_wins(self, btn):
        if btn.get_label() == "Show Wins":
            num_act_btns, btn_ind = self.count_act_check_btns()

            if num_act_btns == 1:
                btn.set_label("Hide Wins")
                
                lwin_l, lwin_r = self.cfg["en_range"][btn_ind][0], self.cfg["en_range"][btn_ind][1]
                rwin_l, rwin_r = self.cfg["en_range"][btn_ind][2], self.cfg["en_range"][btn_ind][3]
                
                self.lines_lwin, = self.ax_en.bar(lwin_l,
                                                  self.y_en_max0,
                                                  lwin_r - lwin_l,
                                                  color = "#f4be41",
                                                  alpha = 0.4)
            
                self.lines_rwin, = self.ax_en.bar(rwin_l,
                                                  self.y_en_max0,
                                                  rwin_r - rwin_l,
                                                  color = "#d3f441",
                                                  alpha = 0.4)
                
        else:
            btn.set_label("Show Wins")
                        
            self.lines_lwin.remove()
            self.lines_rwin.remove()

        self.canvas_en.draw()
        
        
    def set_histos(self, en_spk, t_spk):
        self.en_spk = en_spk
        self.t_spk = t_spk

    def set_cfg(self, cfg):
        self.cfg = cfg

        
    def plot_histo_name(self, name):
        #for line in self.ax_en.lines:
        #    self.ax_en.lines.remove(line)
        x = range(histo_size)
        
        if name in gui_name.en:
            ind = gui_name.en.index(name)
            #self.ax_en.plot(x, self.en_spk[ind], det_clrs[ind])
            self.en_lines[ind].set_ydata(self.en_spk[ind])
            
            self.canvas_en.draw()
            
    def clr_histo_name(self, name):
        x = range(histo_size)
        y_zeros = [0] * histo_size
        #self.ax_en.cla()
        
        for n in gui_name.en:
            ind = gui_name.en.index(n)
            btn = self.check_btn_en[ind]

            if btn.get_active() is False:
                self.en_lines[ind].set_ydata(y_zeros)
                #self.ax_en.plot(x, self.en_spk[ind], det_clrs[ind])

        self.canvas_en.draw()
        
    def plot_all_histos(self):
        self.en_lines = []
        x = range(histo_size)
        
        for i in range(det_num):
            self.en_lines.append( self.ax_en.plot(x, self.en_spk[i], det_clrs[i])[0] )
        self.set_lim_vals_en(0)
            
        self.canvas_en.draw()


    def set_lim_vals_en(self, flag):
        if flag == 0:
            self.x_en_min0, self.x_en_max0 = 0, histo_size
            self.y_en_min0, self.y_en_max0 = 0, 1.05 * max([max(spk) for spk in self.en_spk])
            self.ax_en.set_xlim(self.x_en_min0, self.x_en_max0)
            self.ax_en.set_ylim(self.y_en_min0, self.y_en_max0)


    def _clr_analyze_en(self):
        try:
            self.analyze_curve_en.remove()
            self.analyze_mean_en.remove()
            self.analyze_fwhm_en.remove()
        except ValueError:
            print("ValueError: in _clr_analyze_en()")
            
        
    def _clr_vlines_en(self):
        self.vlines_en[0].remove()
        self.vlines_en[1].remove()

        self.vlines_en = []
        self.x_vlines_en = []
            

def get_histos_from_folder(foldername = "./testspk/"):
    foldername = foldername if foldername[-1] == "/" else foldername + "/"
    
    en_spk = []
    i = 0
    for fname in histo_en_fnames:
        with open(foldername + fname, "rb") as fd:
            fd.seek(histo_finisize, 0)
            en_spk.append(fd.read(histo_fsize))
            en_spk[i] = unpack(str(histo_size) + "i", en_spk[i])
            
        i += 1

    t_spk = []
    i = 0

    return en_spk, t_spk
            

def parse_cfg(path_cfg_file):
    cfg = {}
    
    with open(path_cfg_file, 'r') as cfg_file:
        #if the size of cfg.json will be larger than 2048 bytes, please change it
        config = cfg_file.read(2048)
        config_vals = json_loads(config)
        
        cfg["porog"] = config_vals["porog for detectors"]#[90, 90, 90, 90]
        cfg["delay"] = config_vals["time of coincidence [ch]"]#100
        cfg["coinc"] = 1 if config_vals["coinc mode?"] == "True" else 0#1
        cfg["acq_time"] = config_vals["time of exposition [s]"]#10
        cfg["histo_folder"] = config_vals["histo folder"]
        cfg["en_range"] = []
        for i in range(det_num):
            cfg["en_range"].append( config_vals["en_range D" + str(i + 1)] )

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
            res +="}"

            return res

        cfg_file.write(dict_to_true_cfg_str(cfg))
    

if __name__ == "__main__":
    main()

