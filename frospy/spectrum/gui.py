# import numpy as np
from frospy.spectrum.app import check_input
from frospy.spectrum.spectrum_gui import set_defaults_init_main, print_gui
from frospy.spectrum import spectrum_gui
from frospy.spectrum.controllers import process_input
from frospy.core.spectrum.plot import init_plot_spectrum
from obspy.core import AttribDict

import sys
from tkinter import (
                     Frame, Text, END, StringVar, OptionMenu,
                     Button, Entry, Label
                      )
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.gridspec as gridspec
# from matplotlib.figure import Figure
import matplotlib as mpl
mpl.use('TkAgg')


class frospyGUI(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.master.geometry("1080x806")

        self.master.title('frospy')

        # self.UIFrame = Frame(self.master, relief=RAISED, borderwidth=1)
        # self.UIFrame.pack(side=RIGHT, fill=BOTH, expand=YES)
        # self.plotFrame = Frame(self.master)
        # self.plotFrame.pack(side=LEFT)

        self.plotFrame = Frame(self.master, width=700, height=500,
                               padx=3)
        self.UIFrame = Frame(self.master, width=200, height=220,
                             padx=3, pady=3)
        self.RespondFrame = Frame(self.master, width=200, height=480,
                                  padx=3, pady=3)
        self.QuitFrame = Frame(self.master, width=200, height=80,
                               padx=3)

        self.master.grid_rowconfigure(0, weight=1)
        self.master.grid_columnconfigure(0, weight=1)

        self.plotFrame.grid(row=0, column=0, rowspan=3,
                            sticky="nsew")
        self.UIFrame.grid(row=0, column=1, sticky="new")
        self.RespondFrame.grid(row=1, column=1, sticky="n")
        self.QuitFrame.grid(row=2, column=1, sticky="s")

        self.master.bind('<Return>', self.func)

        self.init_plot()
        self.create_widgets()

    def init_plot(self):
        # Create plotFrame
        # gs = gridspec.GridSpec(12, 6)
        # self.toplevel.bind("<Configure>", resize)
        gs = gridspec.GridSpec(3, 1)
        self.fig, self.ax, self.seg_ax = init_plot_spectrum(gs, 'pretty')
        self.fig.set_size_inches((6, 8))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plotFrame)
        self.canvas.get_tk_widget().grid(row=0, column=0)

    def create_widgets(self):
        _dpath = '/Users/Simon/mounts/eejit/splitting/alldatafiles/060994A.ahx'

        # # Create UIFrame
        self.label_datapath = Label(self.UIFrame, text="Data path")
        self.label_datapath.grid(row=0, column=1, sticky='w')
        self.box_datapath = Entry(self.UIFrame)
        self.box_datapath.insert(END, _dpath)
        self.box_datapath.grid(row=0, column=2, columnspan=2)

        self.label_synpath = Label(self.UIFrame, text="Syn path")
        self.label_synpath.grid(row=1, column=1, sticky='w')
        self.box_synpath = Entry(self.UIFrame)
        self.box_synpath.insert(END, _dpath)
        self.box_synpath.grid(row=1, column=2, columnspan=2)

        self.label_tw = Label(self.UIFrame, text='Time Window (hrs)')
        self.label_tw.grid(row=4, column=1, sticky='w')
        self.box_tw1 = Entry(self.UIFrame, width=9)
        self.box_tw1.insert(END, '5')
        self.box_tw1.grid(row=4, column=2, sticky='ew')

        self.box_tw2 = Entry(self.UIFrame, width=9)
        self.box_tw2.grid(row=4, column=3, sticky='ew')
        self.box_tw2.insert(END, '60')
        self.button_get_tw = Button(self.UIFrame, text="update time window",
                                    fg="black", state="disabled",
                                    command=self.get_field_tw)
        self.button_get_tw.grid(row=5, column=2, sticky='w')

        self.label_fw = Label(self.UIFrame, text='Freq Window (mHz)')
        self.label_fw.grid(row=6, column=1, sticky='w')
        self.box_fw1 = Entry(self.UIFrame, width=9, text='2')
        self.box_fw1.insert(END, '2')
        self.box_fw1.grid(row=6, column=2, sticky='ew')

        self.box_fw2 = Entry(self.UIFrame, width=9, text='2.5')
        self.box_fw2.insert(END, '2.5')
        self.box_fw2.grid(row=6, column=3, sticky='ew')
        self.button_get_fw = Button(self.UIFrame, text="update freq window",
                                    fg="black", state="disabled",
                                    command=self.get_field_fw)
        self.button_get_fw.grid(row=7, column=2, sticky='w')

        # Create previous/next button
        self.label_fw = Label(self.UIFrame, text='command')
        self.label_fw.grid(row=9, column=1, sticky='w')
        # self.box_cmd = Entry(self.UIFrame, state='disabled')
        # self.box_cmd.grid(row=9, column=2, columnspan=3, sticky='ew')

        self.button_cmd = Button(self.UIFrame, text='go', state="disabled",
                                 command=self.get_field_cmd)
        self.button_cmd.grid(row=9, column=3, sticky='ew')

        self.button_prev = Button(self.UIFrame, text="Prev", state="disabled",
                                  command=self.get_prev)
        self.button_prev.grid(row=11, column=2, sticky="w")

        self.button_init_spectrum = Button(self.UIFrame, text="Plot",
                                           fg="black",
                                           command=self.init_spectrum)
        self.button_init_spectrum.grid(row=11, column=3,
                                       sticky='w')

        self.button_next = Button(self.UIFrame, text="Next", state="disabled",
                                  command=self.get_next)
        self.button_next.grid(row=11, column=4, sticky="w")

        cmdvar = StringVar(self.UIFrame)
        choices = {'Show Modes', 'More options', 'even mode'}
        cmdvar.set('Show Modes')
        self.popupMenu = OptionMenu(self.UIFrame, cmdvar, *choices)
        self.popupMenu.grid(row=9, column=2)

        # Create RespondFrame
        self.respond = Text(self.RespondFrame, width=55, height=40,
                            state='disabled')
        self.respond.grid(row=1, column=1, columnspan=3, sticky='ew')

        # Create Quit button
        self.quit = Button(self.QuitFrame, text="Quit", fg="black",
                           command=self.quit)
        self.quit.grid(row=0, column=0)

    def get_field_tw(self):
        self.main.tw = [float(self.box_tw1.get()), float(self.box_tw2.get())]
        self.main = spectrum_gui.run(self)

    def get_field_fw(self):
        self.main.fw = [float(self.box_fw1.get()), float(self.box_fw2.get())]
        self.main = spectrum_gui.run(self)

    def get_field_cmd(self):
        # print('Go go go')
        self.cmd = self.box_cmd.get()
        msg = self.cmd
        out = process_input(self.main.spec, self.main.pick, self.main, msg=msg)
        self.main.spec, self.main.pick, self.main = out
        print_gui(self, '\n\n not ready yet', False)
        # self.clear_response()
        # self.main = spectrum_gui.run(self)

    def get_prev(self):
        self.main.i -= 1
        self.main = spectrum_gui.run(self)

    def get_next(self):
        self.clear_response()
        self.main.i += 1
        self.main = spectrum_gui.run(self)

    def quit(self):
        self.master.destroy
        sys.exit()

    def func(self, event):
        # print("You hit return.")
        self.get_field_cmd()

    def clear_response(self):
        self.respond.configure(state="normal")
        self.respond.replace("0.0", END, '')
        self.respond.configure(state="disabled")
        self.master.update()

    def init_spectrum(self):
        self.datapath = self.box_datapath.get()
        self.synpath = self.box_synpath.get()

        if self.datapath == '':
            self.datapath = None
        if self.synpath == '':
            self.synpath = None

        self.datapath, self.iargs = check_input(self.datapath,
                                                {'syn': self.synpath})

        self.main = AttribDict(self.iargs)
        self.main = set_defaults_init_main(self.main)
        self.main.minispec = 'pretty'
        # self.main.verbose = False
        self.main.tw = [float(self.box_tw1.get()), float(self.box_tw2.get())]
        self.main.fw = [float(self.box_fw1.get()), float(self.box_fw2.get())]

        msg = "Input\n\nData: %s" % self.datapath
        msg += "\nSyn: %s" % self.synpath
        msg += "\ntw: %s %s" % (self.main.tw[0], self.main.tw[1])
        msg += "\nfw: %s %s\n\n" % (self.main.fw[0], self.main.fw[1])
        print_gui(self, msg)

        self.main.rfig = self.fig
        self.main.rax = self.ax
        self.main.seg_ax = self.seg_ax

        self.main.pick = None
        self.main.spec = None
        self.main = spectrum_gui.run(self)

        # self.box_cmd.configure(state="normal")
        self.button_cmd.configure(state="normal")
        self.button_next.configure(state="normal")
        self.button_prev.configure(state="normal")
        self.button_get_fw.configure(state="normal")
        self.button_get_tw.configure(state="normal")
        self.button_init_spectrum.configure(state="disabled")
        self.button_init_spectrum.configure(state="disabled")


root = tk.Tk()
app = frospyGUI(master=root)
app.mainloop()
