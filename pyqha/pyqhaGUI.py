#encoding: UTF-8
# Copyright (C) 2016 Mauro Palumbo
# This file is distributed under the terms of the # MIT License. 
# See the file `License' in the root directory of the present distribution.


import wx
#import matplotlib
#matplotlib.use('WXAgg')
import matplotlib.pyplot as plt

import time, sys, os
import numpy as np
from read import read_EtotV, read_Etot
from eos import fit_Murn, print_data, calculate_fitted_points, write_Etotfitted
from fitutils import fit_anis
from minutils import find_min
from fitEtot import fitEtot, fitEtotV
from fitFvib import fitFvib, fitFvibV
from fitC import fitCT
from alphagruneisenp import compute_alpha_gruneisein


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
            ("&Open E_tot file...\tctrl+E", "Load total energies file", self.OnOpenEtot),
            ("&Open DOS file(s)...\tctrl+D", "Load DOS file(s)", self.OnOpenDOS),
            ("", "", ""),
            ("Quit...\tctrl+Q", "Exit the program", self.OnQuit)),
                ("&Compute",
            ("Fit E_tot\tctrl+F", "Fit the total energies", self.OnFitEtot)),

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
    
    def OnOpenEtot(self, event):
        wildcard =  "dat file (*.dat)|*.dat|" \
                    "All files (*)|*"
        dialog = wx.FileDialog(None, "Choose a file with total energies", os.getcwd(),"", wildcard)
        if dialog.ShowModal() == wx.ID_OK:
            try:
                fname = dialog.GetPath()
                self.V, self.E = read_EtotV(fname)   
                self.IsEtotRead = True
                self.SetStatusText("Read E_tot file: "+fname)
            except:
                wx.MessageBox("Something wrong while opening the E_tot file... not loaded.",
                "", wx.OK | wx.ICON_EXCLAMATION, self)  

        dialog.Destroy()
        
        
    def OnOpenDOS(self, event):
        wildcard =  "dos file (*.dos)|*.dos|" \
                    "All files (*)|*"
        dialog = wx.FileDialog(None, "Choose one or more files with phonon DOS", os.getcwd(),"", wildcard, style= wx.FD_MULTIPLE)
        if dialog.ShowModal() == wx.ID_OK:
            try:
                fnames = dialog.GetPaths()
                print (fnames)
            except:
                wx.MessageBox("Something wrong while opening the DOS file(s)... not loaded.",
                "", wx.OK | wx.ICON_EXCLAMATION, self)  

        dialog.Destroy()
    
        
    def OnFitEtot(self, event):
        try:
            a, cov, chi = fit_Murn(self.V,self.E)
            print_data(self.V,self.E,a,chi,"Etot")
            #write_Etotfitted(outputfileEtot,V,E,a,chi,"Etot")
                        
            # Plotting using matplotlib
            Vdense, Edensefitted = calculate_fitted_points(self.V,a)
    
            plt.plot(self.V, self.E, 'o', label='Etot data', markersize=10)
            plt.plot(Vdense, Edensefitted, 'r', label='Fitted EOS')
            plt.legend()
            plt.xlabel('V (a.u.^3)')
            plt.ylabel('E (a.u.) ')
            plt.show()
        except:
            wx.MessageBox("Something wrong while fitting total energies...",
                "", wx.OK | wx.ICON_EXCLAMATION, self)  
    
    def OnQuit(self, event):
        plt.close("all")  # close all matplotlib figures
        self.Destroy()

                
    def OnHelp(self, event):
        wx.MessageBox("A program to perform quasi-harmonic calculations",
        "Help Contents", wx.OK | wx.ICON_INFORMATION, self)


    def OnAbout(self, event):
        wx.MessageBox("pyQHA: a program to perform quasi-harmonic calculations",
        "About pyQHA", wx.OK | wx.ICON_INFORMATION, self)      
        

if __name__ == "__main__":
    app = MyApp(False)
    app.MainLoop()
