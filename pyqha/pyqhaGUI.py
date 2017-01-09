#encoding: UTF-8
# Copyright (C) 2016 Mauro Palumbo
# This file is distributed under the terms of the # MIT License. 
# See the file `License' in the root directory of the present distribution.


"""
A tentative GUI for pyqha. Jut started...

It uses :py:mod:`wxPython`. Note that at the time of writing this library is 
available for Python 2.x only. Porting of the library to Python 3.x is ongoing.
""" 

import wx
#matplotlib.use('WXAgg')
import matplotlib.pyplot as plt

import time, sys, os
import numpy as np
from pyqha.read import read_EtotV, read_Etot
#from pyqha import fit_Murn, print_data, calculate_fitted_points, fit_anis, find_min
from pyqha import fitEtot, fitEtotV, fitFvib, fitFvibV, fitCT, compute_alpha_gruneisein
from pyqha.eos import fit_Murn, print_eos_data

from pyqha.plotutils import plot_EV


class RedirectText:
    def __init__(self,aWxTextCtrl):
        self.out=aWxTextCtrl
 
    def write(self,string):
        self.out.WriteText(string)


class MyApp(wx.App):
    def OnInit(self):
        frame = MyFrame("pyQHA", (50, 60), (680, 480))
        frame.Show()
        self.SetTopWindow(frame)
        return True
    
class MyFrame(wx.Frame):
    def __init__(self, title, pos, size):
                
        # Initialize some flags
        self.IsEtotRead = False
        
        wx.Frame.__init__(self, None, -1, title, pos, size)                
        panel = wx.Panel(self, -1)
        panel.SetBackgroundColour("White")
        
        text = wx.TextCtrl(self, wx.ID_ANY, size=size, style = wx.TE_MULTILINE|wx.TE_READONLY|wx.HSCROLL)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(text, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.Fit()
        
        self.Bind(wx.EVT_CLOSE, self.OnQuit)                             
        self.createMenuBar()       
        self.CreateStatusBar()        
        self.SetStatusText("Welcome to pyQHA")
        
        # Redirect stout to the TextCtrl in the main panel
        self.redir=RedirectText(text)
        sys.stdout=self.redir
        sys.stderr=self.redir

    ############################################################################
    #
    # Menu functions
    #
    
    def menuData(self):
        """
        This function returns the menu items, including shortkeys and help strings.
        To modify the menus of the application you need to change them here.
        """
        return (("&File",
            ("Open Etot(V) file...\tctrl+O", "", self.OnOpenEtotV), 
            ("", "", ""),
            ("Quit...\tctrl+Q", "Exit the program", self.OnQuit)),
                ("&Compute",
            ("Fit iso E_tot\tctrl+E", "Fit the total energies as a function of volume", self.OnFitisoEtot),
            ("Fit aniso E_tot\tctrl+F", "Fit the total energies as a function of lattice parameters", self.OnFitanisoEtot)),
            
                ("&Help",
            ("&Help contents\tctrl+H", "Compute the bare potential from the electronic charge", self.OnHelp),
            ("&About\tctrl+A", "Compute the bare plus the Hartree potential from the electronic charge", self.OnAbout)))
        
        
    def createMenuBar(self):
        menuBar = wx.MenuBar()
        for eachMenuData in self.menuData():
            menuLabel = eachMenuData[0]
            menuItems = eachMenuData[1:]
            menuBar.Append(self.createMenu(menuItems), menuLabel)
        self.SetMenuBar(menuBar)
        
        
    def createMenu(self, menuData):
        menu = wx.Menu()
        for eachLabel, eachStatus, eachHandler in menuData:
            if not eachLabel:
                menu.AppendSeparator()
                continue
            menuItem = menu.Append(-1, eachLabel, eachStatus)
            self.Bind(wx.EVT_MENU, eachHandler, menuItem)
        return menu


    ############################################################################
    #
    # Handling methods, rather self explaining
    #
    
    def OnOpenEtotV(self, event):
        wildcard =  "dat file (*.dat)|*.dat|" \
                    "All files (*)|*"
        dialog = wx.FileDialog(None, "Choose a file with total energies", os.getcwd(),"", wildcard)
        if dialog.ShowModal() == wx.ID_OK:
            try:
                fname = dialog.GetPath()
                self.V, self.E = read_EtotV(fname)   
                self.IsEtotVRead = True
                self.SetStatusText("E_tot(V) file "+fname+" read")
            except:
                wx.MessageBox("Something wrong while opening the E_tot file... not loaded.",
                "", wx.OK | wx.ICON_EXCLAMATION, self)  

        dialog.Destroy()
        
    
        
    def OnFitisoEtot(self, event):
        try:
            
            self.a, self.cov, self.chi = fit_Murn(self.V,self.E)
            print_eos_data(self.V,self.E,self.a,self.chi,"Etot")
    
            fig1 = plot_EV(self.V,self.E,self.a)                  	# plot the E(V) data and the fitting line
            fig1.savefig("figure_1.png")

        except:
            wx.MessageBox("Something wrong while fitting total energies...",
                "", wx.OK | wx.ICON_EXCLAMATION, self)  
                
    def OnFitanisoEtot(self, event):
        try:
            wx.MessageBox("Not implemented yet...",
                "", wx.OK | wx.ICON_EXCLAMATION, self)  

        except:
            wx.MessageBox("Something wrong while fitting total energies...",
                "", wx.OK | wx.ICON_EXCLAMATION, self)  
                
    
    def OnQuit(self, event):
        plt.close("all")  # close all matplotlib figures
        self.Destroy()

                
    def OnHelp(self, event):
        wx.MessageBox("A program to perform quasi-harmonic calculations based on the pyqha module",
        "Help Contents", wx.OK | wx.ICON_INFORMATION, self)


    def OnAbout(self, event):
        wx.MessageBox("pyqhaGUI: a program to perform quasi-harmonic calculations",
        "About pyQHA", wx.OK | wx.ICON_INFORMATION, self)      
        

if __name__ == "__main__":
    app = MyApp(False)
    app.MainLoop()
