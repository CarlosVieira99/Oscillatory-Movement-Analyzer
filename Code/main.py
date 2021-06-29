import sys
import warnings
from PyQt4 import QtGui, uic
from PyQt4.QtGui import QFileDialog, QTextEdit, QPixmap
from PyQt4.Qwt5.Qwt import QwtPlot, QwtPlotCurve
from PyQt4 import QtCore as QtCore

# Libraries
from visual import * #import vPython library
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
import matplotlib.animation as animation
from matplotlibwidget import MatplotlibWidget
from time import sleep
import os
import sympy as sym
import urllib2
import time
import csv
from scipy.signal import find_peaks_cwt
from scipy.optimize import curve_fit

x = []
y = []
    
def _translate(context, text, disambig):
    return QtGui.QApplication.translate(context, text, disambig, _encoding)
def _translate(context, text, disambig):
    return QtGui.QApplication.translate(context, text, disambig)
    self.retranslateUi(MainWindow)

class MyWindow(QtGui.QMainWindow):

    global d
    global t
    global dPlot
    global tPlot
    global h
    global speed_sound
    global spring_k
#    
    def __init__(self):
        super(MyWindow, self).__init__()
        uic.loadUi('App_Designer.ui', self)
        
        #self.calibrationButton.clicked.connect(self.calibration_handle)
        self.simulation_button.clicked.connect(self.simulation)
        self.tabWidget.setCurrentIndex(0)
    
        self.calibrationPlot.hide() 
        self.hookePlot.hide()
        self.dampedPlot.hide()
        self.dataTable.cellChanged.connect(self.table_handle)
        self.dataTable2.cellChanged.connect(self.table2_handle)
        self.actionCSV_File.triggered.connect(self.save_datacsv) 
        self.actionCSV_File_2.triggered.connect(self.open_datacsv)
        self.actionCSV_File_of_Dumping.triggered.connect(self.open_dampedcsv)
        self.dataTable3.cellChanged.connect(self.w0_calc)
        
        self.startDamped_Button.clicked.connect(self.startDamped)
        self.stopDamped_Button.clicked.connect(self.damped_stop)
        self.simDamped_Button.clicked.connect(self.damped_sim)
        
        global d
        d = zeros(10)
        global t
        t = zeros(10)
        global dPlot
        dPlot = []
        global tPlot
        tPlot = []
        global mass
        mass = zeros(10)
        global d2
        d2 = zeros(10)
        global massPlot
        massPlot = []
        global d2Plot
        d2Plot = []
        global curve1
        
        global xs, ys       
        xs = []
        ys = []
                
        curve1 = QwtPlotCurve("Curve 1")
        #curve1.setData(rand, np.zeros(30))
        curve1.attach(self.qwtPlot)
        #QwtPlot.setAxisTitle(axisId, title)[source]
        #Set title
        self.qwtPlot.show()
        
        #app = QtGui.QApplication(sys.argv)
        self.show()
        
        #Take screenshot
        #filename = 'Screenshot.jpg'
        #p = QPixmap.grabWindow(self.Tab_Intro.winId())
        #p.save(filename, 'jpg')
        #QPixmap.grabWindow(self.Tab_Intro.winId()).save(filename, 'jpg')
        
        #self.simulationPlot = MatplotlibWidget(title='Simulation', xlabel='sal', ylabel='salgrosso')
        
        self.simulationPlot.hide()        
        self.simulationPlot.figure.set_facecolor("white")        
        self.simulationPlot.axes.set_title('GRAPH') # Nothing Happens
        self.simulationPlot.axes.plot(x,y) 
        self.simulationPlot.axes.set_xlabel('Sal')
        
        #global ax2
        #ax2 = self.calibrationPlot.axes.twinx()        
        
        self.show()
        
    def save_datacsv(self):
        #Open Dialog Window
        filename = QtGui.QFileDialog.getSaveFileName(self, 'Save File',
                                                     '~/Desktop',
                                                     'Data files (*.csv)')
        ## Write Data        
        with open(filename, mode='wb') as csv_file: #default open as read mode
            fieldnames = ['Distance_Calibration', 'Time_Calibration',
                          'Temperature_Calibration', 'Humidity_Calibration',
                          'Notes_Calibration',
                          'Mass_Hooke', 'Distance_Hooke',
                          'MassSpring_Hooke', 'Temperature_Hooke',
                          'Length_Hooke','Notes_Hooke',
                          'Temperature_Damped', 'MassOfOscillator_Damped',
                          'OscillationPeriod_Damped', 'Damping_Damped',
                          'Frequency0_Damped', 'OscillationFrequency_Damped',
                          'OscillationPeriod_Damped', 'QualityFactor_Damped',
                          'DecayAmplitude_Damped', 'Notes_Damped']
                          
            writer = csv.DictWriter(csv_file, fieldnames = fieldnames)
            writer.writeheader()
            for i in range(0, 10):
                writer.writerow({'Distance_Calibration': self.dataTable.item(i, 0).text(), 
                                 'Time_Calibration': self.dataTable.item(i, 1).text(),
                                 'Temperature_Calibration': self.resultsTable.item(0, 1).text(),
                                 'Humidity_Calibration': self.resultsTable.item(0, 2).text(),
                                 'Notes_Calibration': self.textEdit.toPlainText(),
                                 'Mass_Hooke': self.dataTable2.item(i, 0).text(),
                                 'Distance_Hooke': self.dataTable2.item(i, 1).text(),
                                 'MassSpring_Hooke': self.resultsTable2.item(0, 1).text(),
                                 'Temperature_Hooke': self.resultsTable2.item(0, 2).text(),
                                 'Length_Hooke': self.resultsTable2.item(0, 4).text(),
                                 'Notes_Hooke': self.textEdit_2.toPlainText(),
                                 'Temperature_Damped': self.dataTable3.item(0, 0).text(),
                                 'MassOfOscillator_Damped': self.dataTable3.item(1,0).text(),
                                 'OscillationPeriod_Damped': self.dataTable3.item(2,0).text(),
                                 'Damping_Damped': self.dataTable3.item(3, 0).text(),
                                 'Frequency0_Damped': self.dataTable3.item(4, 0).text(),
                                 'OscillationFrequency_Damped': self.dataTable3.item(5, 0).text(),
                                 'OscillationPeriod_Damped': self.dataTable3.item(6, 0).text(),
                                 'QualityFactor_Damped': self.dataTable3.item(7, 0).text(),
                                 'DecayAmplitude_Damped': self.dataTable3.item(8, 0).text(),
                                 'Notes_Damped': self.textEdit_3.toPlainText()
                                 })
                                 
        
    def open_datacsv(self):
        global speed_sound
        global spring_k
        
        tPlotOpen = np.zeros(10)
        dPlotOpen = np.zeros(10)
        massPlotOpen = np.zeros(10, dtype=float)
        d2PlotOpen = np.zeros(10)
        
        self.dataTable.cellChanged.disconnect(self.table_handle)
        self.dataTable2.cellChanged.disconnect(self.table2_handle)
        
        filename = QFileDialog.getOpenFileName(self,
                                               'Open File', 
                                               '~/Desktop',
                                               'Data files (*.csv)')
                                               
        i = 0
        
        with open(filename, mode='rb') as csv_file:
            csv_reader = csv.DictReader(csv_file)
            for row in csv_reader:
                self.dataTable.item(i, 0).setText(_translate("MyWindow", row['Distance_Calibration'], None))
                self.dataTable.item(i, 1).setText(_translate("MyWindow", row['Time_Calibration'], None))
                self.resultsTable.item(0, 1).setText(row['Temperature_Calibration'])
                self.resultsTable.item(0, 2).setText(_translate("MyWindow", row['Humidity_Calibration'], None))
                self.textEdit.setPlainText(_translate("MyWindow", row['Notes_Calibration'], None))
                
                self.dataTable2.item(i, 0).setText(_translate("MyWindow", row['Mass_Hooke'], None))
                self.dataTable2.item(i, 1).setText(_translate("MyWindow", row['Distance_Hooke'], None))
                self.resultsTable2.item(0, 1).setText(row['MassSpring_Hooke'])
                self.resultsTable2.item(0, 2).setText(row['Temperature_Hooke'])
                self.resultsTable2.item(0, 4).setText(row['Length_Hooke'])
                self.textEdit_2.setPlainText(_translate("MyWindow", row['Notes_Hooke'], None))
                
                self.dataTable3.item(0, 0).setText(_translate("MyWindow", row['Temperature_Damped'], None))
                self.dataTable3.item(1, 0).setText(_translate("MyWindow", row['MassOfOscillator_Damped'], None))
                self.dataTable3.item(2, 0).setText(_translate("MyWindow", row['OscillationPeriod_Damped'], None))
                self.dataTable3.item(3, 0).setText(_translate("MyWindow", row['Damping_Damped'], None))
                self.dataTable3.item(4, 0).setText(_translate("MyWindow", row['Frequency0_Damped'], None))
                self.dataTable3.item(5, 0).setText(_translate("MyWindow", row['OscillationFrequency_Damped'], None))
                self.dataTable3.item(6, 0).setText(_translate("MyWindow", row['OscillationPeriod_Damped'], None))
                self.dataTable3.item(7, 0).setText(_translate("MyWindow", row['QualityFactor_Damped'], None))
                self.dataTable3.item(8, 0).setText(_translate("MyWindow", row['DecayAmplitude_Damped'], None))
                self.textEdit_3.setPlainText(_translate("MyWindow", row['Notes_Damped'], None))                
                
                i=i+1
        
        for i in range(0, 10):
            tPlotOpen[i] = int(self.dataTable.item(i, 0).text())
            dPlotOpen[i] = float(self.dataTable.item(i, 1).text())
            
        for i in range(0, 10):
            massPlotOpen[i] = float(self.dataTable2.item(i, 0).text()) * 9.80665 / 1000
            d2PlotOpen[i] = float(self.dataTable2.item(i, 1).text())
            
        initDistance = float(d2PlotOpen[0]/1000)
            
        for i in range(len(d2PlotOpen)-1, -1, -1):
            d2PlotOpen[i] = d2PlotOpen[0] - d2PlotOpen[i]
        
        self.calibrationPlot.figure.set_facecolor("white")        
        self.calibrationPlot.axes.plot(tPlotOpen, dPlotOpen,'go--')
        self.calibrationPlot.axes.set_title('SONAR Calibration') # Nothing Happens
        self.calibrationPlot.axes.set_xlabel('Distance (mm)')
        self.calibrationPlot.axes.set_ylabel('Time (us)')
        self.calibrationPlot.axes.grid()
        self.calibrationPlot.axes.set_xlim(xmin=0)
        self.calibrationPlot.axes.set_ylim(ymin=0)
        self.calibrationPlot.show()
        
        self.hookePlot.figure.set_facecolor("white")        
        self.hookePlot.axes.plot(d2PlotOpen, massPlotOpen,'go--')
        self.hookePlot.axes.set_title("Hooke's Law") # Nothing Happens
        self.hookePlot.axes.set_xlabel('Stretching (mm)')
        self.hookePlot.axes.set_ylabel('Force (N)')
        self.hookePlot.axes.grid()
        self.hookePlot.axes.set_xlim(xmin=0)
        self.hookePlot.axes.set_ylim(ymin=0)
        self.hookePlot.show()
        
        self.dataTable.cellChanged.connect(self.table_handle)
        self.dataTable2.cellChanged.connect(self.table2_handle)
        
        m, b = np.polyfit(tPlotOpen, dPlotOpen, 1);  
        m = (1/m)*1000*2
       
        speed_sound = round(m, 2) 
        
        #m, b = np.polyfit(massPlotOpen, d2PlotOpen, 1);  
        m, b = np.polyfit(d2PlotOpen, massPlotOpen, 1);  
        spring_k = round(m*1000, 2)
        
        ref_value = 331.3 + (0.6*float(self.resultsTable.item(0, 1).text()[:3]))
        
        self.resultsTable.item(0, 0).setText(_translate("MyWindow", str(speed_sound) + ' +/- 1', None))
        self.resultsTable.item(0, 3).setText(_translate("MyWindow", str(ref_value), None))
        
        self.resultsTable2.item(0, 0).setText(_translate("MyWindow", str(spring_k) + ' +/- 1', None))
        self.resultsTable2.item(0, 3).setText(_translate("MyWindow", str(initDistance) + ' +/- 1', None))
        
    def w0_calc(self):
        #m = float(self.dataTable3.item(1, 0).text())
        #spring_k = float(self.resultsTable2.item(0, 0).text())
        m = 0.220
        spring_k = 4.39
        w0 = round(np.sqrt(spring_k/m), 3)
        self.dataTable3.item(4, 0).setText(_translate("MyWindow", str(w0), None))
        
        ts = round(2*pi/w0, 3)

        self.dataTable3.item(2, 0).setText(_translate("MyWindow", str(ts), None))
        
    def open_dampedcsv(self):
        def Read_Two_Column_File(file_name):
            with open(file_name, 'r') as data:
                x = []
                y = []
                for line in data:
                    p = line.split()
                    x.append(float(p[0]))
                    y.append(float(p[1]))
        
            return x, y

        x, y = Read_Two_Column_File('matrix.dat')
        
        t = [i / 10000 for i in x]
        d = [i * speed_sound / 10000 for i in y]
        
        t = np.array(t)
        d = np.array(d)
        
        peakind = find_peaks_cwt(d, np.arange(0.001/max(d), 20))
        peakind2 = np.array(peakind)
        
        periodo = np.diff(t[peakind2])
        periodo = np.mean(periodo)
        
        print("Period between peaks=", periodo)
        
        def func(x, a, b):
            return a*np.exp(-b*x)
            
        
        popt, pcov = curve_fit(func, t[peakind2], d[peakind2])
        
        print(popt, pcov)
        
        print("A=", popt[0]/100)
        
        print("Gama=", popt[1]*2)
        
        gama = popt[1]*2
        
        self.dataTable3.item(3, 0).setText(_translate("MyWindow", str(round(gama, 3)), None))
        
        qf = 4.467 / gama
        
        self.dataTable3.item(7, 0).setText(_translate("MyWindow", str(round(qf, 3)), None))
        
        e = func(t, popt[0], popt[1])
            
        self.dampedPlot.figure.set_facecolor("white")        
        self.dampedPlot.axes.plot(t, d,'-')
        self.dampedPlot.axes.hold(True)
        self.dampedPlot.axes.plot(t[peakind2], d[peakind2], 'o')
        self.dampedPlot.axes.plot(t, e, 'r-')
        self.dampedPlot.axes.set_title('Damped Oscillator') # Nothing Happens
        self.dampedPlot.axes.set_xlabel('Time (s)')
        self.dampedPlot.axes.set_ylabel('Distance (cm)')
        self.dampedPlot.axes.grid()
        self.dampedPlot.axes.set_xlim(xmin=0)
        self.dampedPlot.axes.set_ylim(ymin=0)
        self.dampedPlot.show()
        self.dampedPlot.axes.hold(False)
        
    def getData(self):
        resp = urllib2.urlopen('http://192.168.4.1')
        thePage = resp.read()
        thePageArray = thePage.split(',')
        time = int(thePageArray[1])
        temperature = float(thePageArray[2])
        temperature = int(round(temperature, 0))
        humidity = float(thePageArray[3])
        humidity = int(round(humidity, 0))
        return time, temperature, humidity

    def table_handle(self):
        global d
        global t
        global dPlot
        global tPlot
        global h
        global speed_sound
        
        tData, tempData, humData = self.getData()
        tSonar = float(tData)/2.0  
        
        cur = self.dataTable.currentRow()
        self.dataTable.cellChanged.disconnect(self.table_handle)
        self.dataTable.item(cur, 1).setText(_translate("MyWindow", str(tSonar), None))
        self.dataTable.cellChanged.connect(self.table_handle)
        
        d[cur] = (float(self.dataTable.item(cur, 0).text())/100)
        t[cur] = float(tSonar)
        
        if(len(dPlot) < len(d) and len(dPlot) < (cur+1)):
            for j in range(len(dPlot), len(d)):
                if d[j] != 0:
                    dPlot.append(d[j])
                if t[j] != 0:
                    tPlot.append(t[j])
        else:
            dPlot[cur] = d[cur]
            tPlot[cur] = t[cur]
        
        m, b = np.polyfit(tPlotOpen, dPlotOpen, 1);  
        m = (1/m)*1000*2
        speed_sound = round(m, 2)
        
        ref_value = 331.3 + (0.6*tempData)
        
        self.resultsTable.item(0, 0).setText(_translate("MyWindow", str(speed_sound) + ' +/- 1', None))
        self.resultsTable.item(0, 1).setText(_translate("MyWindow", str(tempData)  + ' +/- 1', None))
        self.resultsTable.item(0, 2).setText(_translate("MyWindow", str(humData)  + ' +/- 1', None))
        self.resultsTable.item(0, 3).setText(_translate("MyWindow", str(ref_value), None))
        
        h = zeros(size(dPlot))
        
        for k in range(0, len(dPlot)):
            h[k] = (dPlot[k]*m+b)
        
        def animate_cal(i):
            global dPlot
            global tPlot
            global h
            
            dPlot2 = [x*100 for x in dPlot]
            
            self.calibrationPlot.figure.set_facecolor("white")        
            self.calibrationPlot.axes.plot(dPlot2, tPlot,'go--')
            self.calibrationPlot.axes.set_title('SONAR Calibration') # Nothing Happens
            self.calibrationPlot.axes.set_xlabel('Distance (mm)')
            self.calibrationPlot.axes.set_ylabel('Time (us)')
            self.calibrationPlot.axes.set_xlim(xmin=0)
            self.calibrationPlot.axes.set_ylim(ymin=0)
            self.calibrationPlot.axes.grid()
            #ax2.plot(dPlot, h, 'b-', label="Linear Regression")              
            
        ani = animation.FuncAnimation(self.calibrationPlot.figure, animate_cal, interval=100)
        self.calibrationPlot.show()
            
    def table2_handle(self):
        global d2
        global mass
        global d2Plot
        global massPlot
        global h2
        global spring_k
        global speed_sound
        
        tData, tempData, humData = self.getData()    
        
        tSonar = float(tData)/2.0
        dSonar = round((tSonar / 1000000)* speed_sound * 100, 2)
        
        cur = self.dataTable2.currentRow()
        self.dataTable2.cellChanged.disconnect(self.table2_handle)
        self.dataTable2.item(cur, 1).setText(_translate("MyWindow", str(dSonar), None))
        if cur == 0:
            self.resultsTable2.item(0, 3).setText(_translate("MyWindow", str(dSonar) + ' +/- 1', None))
        self.dataTable2.cellChanged.connect(self.table2_handle)
        
        mass[cur] = (float(self.dataTable2.item(cur, 0).text()))
        d2[cur] = dSonar
        
        if(len(massPlot) < len(mass) and len(massPlot) < (cur+1)):
            for j in range(len(massPlot), len(mass)):
                if mass[j] != 0:
                    massPlot.append(mass[j])
                if d2[j] != 0:
                    d2Plot.append(d2[j])
        else:
            massPlot[cur] = mass[cur]
            d2Plot[cur] = d2[cur]
        
        m, b = np.polyfit(d2Plot, massPlot, 1);  
        spring_k = round(m, 2)
        
        self.resultsTable2.item(0, 0).setText(_translate("MyWindow", str(spring_k) + ' +/- 1', None))
        self.resultsTable2.item(0, 2).setText(_translate("MyWindow", str(tempData) + ' +/- 1', None))
        
        def animate_hooke(i):
            global d2Plot
            global massPlot
            
            #d2Plot2 = [x*100 for x in dPlot]
            
            self.hookePlot.figure.set_facecolor("white")        
            self.hookePlot.axes.plot(d2Plot, massPlot,'go--')
            self.hookePlot.axes.set_title("Hooke's Law") # Nothing Happens
            self.hookePlot.axes.set_xlabel('Distance (mm)')
            self.hookePlot.axes.set_ylabel('Mass (g)')
            self.hookePlot.axes.grid()
            
        anie = animation.FuncAnimation(self.hookePlot.figure, animate_hooke, interval=100)
        self.hookePlot.show()
        
    def startDamped(self):
        # set up the QwtPlot (pay attention!)
        self.timer = QtCore.QTimer() #start a timer (to call replot events)
        self.timer.start(5.0) #set the interval (in ms)
        self.qwtPlot.connect(self.timer, QtCore.SIGNAL('timeout()'), self.damped_start)
        global tempo1
        tempo1 = time.time() # Time reference        
        
    def damped_start(self):
        yD, _, _ = self.getData()  
        
        global tempo, tempo1, curve1, xs, ys
        global speed_sound
        
        tempo = time.time() - tempo1
        # Limit to 30 samples displayed
        if(len(xs) == 50):
            xs.pop(0)
        if(len(ys) == 50):
            ys.pop(0)
            
        xs.append(tempo)
        ys.append(yD*speed_sound/100000)
        
        print(yD*speed_sound/100000)
        
        curve1.setData(xs, ys)
        self.qwtPlot.setTitle('Damping Visualization')
        self.qwtPlot.setAxisTitle(0, 'Distance (cm)')
        self.qwtPlot.setAxisTitle(2, 'Time (s)')
        self.qwtPlot.replot()  
    
    def damped_stop(self):
        self.qwtPlot.disconnect(self.timer, QtCore.SIGNAL('timeout()'), self.damped_start)
        
        _, temp, _ = self.getData()
        
        self.dataTable3.item(0, 0).setText(_translate("MyWindow", str(temp) + ' +/- 1', None))
        
    def damped_sim(self):
        xs = []
        ys = []
        
        global speed_sound
        
        # Scene Configuration

        MyScene=display(title='Real Time 3D Spring Simulation') #Create your scene and give it a title.
        MyScene.width=800  #We can set the dimension of your visual box. 800X800 pixels works well on my screen
        MyScene.height= 800
        
        
        #MyScene.autoscale=False #We want to set the range of the scene manually for better control. Turn autoscale off
        #MyScene.range = (12,12,12) #Set range of your scene to be 12 inches by 12 inches by 12 inches. 
        
        #Structure Creation
        
        support_base = box(pos=(0, -1, 0), size=(6, 0.1, 6), material = materials.wood)
        
        #support_base.visible = False
        
        support_tall = box(pos=(  support_base.pos.x, 1.5 + support_base.pos.y -(support_base.size.y/2)+18.5, 
                                        support_base.pos.z - support_base.size.z/2 - support_base.size.y/2), 
                                                size=(6, 40, 0.1), material = materials.wood)
                                                
        support_ceil = box(pos=(support_base.pos.x, support_tall.size.y + support_base.pos.y, support_base.pos.z - support_tall.size.z),
                                                size=(6, 0.1, 6), material = materials.wood)
                                                
        spring = helix(pos=(support_ceil.pos.x, support_ceil.pos.y, support_ceil.pos.z), axis=(0, -1, 0), 
               radius=1, coils=8, thickness = 0.01, stiffness=1, length = 4.28)
        
        eq_pos = spring.pos + spring.axis
                
        endSpring  =    box(pos=eq_pos, color=color.yellow,  size=(0.75, 0.01, 0.5) ,
                            mass=1.0, velocity=vector(0, 0.1, 0),
                            force=(0, 0, 0), make_trail=True, retain = 10, trail_color = color.green, opacity = 0.5)
        
        #ball = sphere(  pos=eq_pos, color=color.red, size=(0.75, 0.1, 0.5))
                       # mass=1.0, velocity=(0, 0.1, 0),
                       # force=(0, 0, 0)                   )
                        
        sonar_base  =    box(pos=(support_base.pos.x,support_base.pos.y+support_base.size.y, support_base.pos.z-0.1), size=(4.5, 0.1, 2), color = color.blue)
        
        sonar_echo  =    cylinder(color=color.green, pos=(sonar_base.pos.x - 1, sonar_base.pos.y, sonar_base.pos.z), radius=0.75, length=1, axis=(0, 1, 0) )
        
        sonar_trigger   =  cylinder(color=color.green, pos=(sonar_base.pos.x + 1, sonar_base.pos.y, sonar_base.pos.z), radius=0.75,length=1, axis=(0, 1, 0))
                        
        labl =          label(pos=(1, 1, 0), text="Simulation 3D") # Stays always in the same place
        
        dist = endSpring.pos.y - sonar_base.pos.y
        
        dist_label =    label(pos=(1, 0.5, 0)) # Stays always in the same place
        
        ## Camera Controls
        
        #print(scene.camera.pos)
        #print(scene.camera.axis)
        
        #scene. = vector(0, 20, 0)
        scene.camera.axis = vector(0, 0, 20)        
        
        #scene.camera.follow(endSpring)
        
        dt = 0.01 # Time step
        
#        while(True):
#            rate(1000)
#            
#            tData, tempData, humData = self.getData()
#            yD = tData * speed_sound / 1000
#            
#            endSpring.pos = (0, yD/1000, 0)
#            spring.axis = endSpring.pos - spring.pos
#            dist = endSpring.pos.y - sonar_base.pos.y
#            dist_label.text = ("Distance: " + "{0:.2f}".format(dist) + "cm")
        
    def simulation(self):
        scene = display(title='Simple Harmonic Oscillator', autoscale=0,
                width=500, height=500, range=7, forward=vector(0,0,-1))
        chao = box(pos=vector(0,-6.1,0),size=vector(20,0.1,10),material=materials.wood)
        parede = box(pos=vector(0,3.85,-5.05),size=vector(20,20,0.1), color=vector(0.7,0.7,0.7))
        #sitio = text(pos=vector(0,3.85,-5), text='Carlos Silva', color=color.blue, align='center', depth=0)
        
        # Mola e cilindro
        bar0 = curve(radius=0.03)
        for ang in arange(pi,-pi/2.,-0.1):
            bar0.append(vector(0, 5.2+0.23*sin(ang), 0.23*cos(ang)))
        bar0.append(pos=vector(0,4.5,0))
        bar0.append(pos=vector(0,4.5,0.3))
        mola1 = helix(pos=vector(0,4.5,0), radius=0.3, thickness=0.05, coils=30,
                      axis = vector(0,-5.8,0))
    
        bar1 = curve(radius=0.03, pos=[vector(0,-1.6,0),vector(0,-1.3,0),vector(0,-1.3,0.3)])
        c1 = cylinder(pos=vector(0,-2,0),radius=0.5,axis=vector(0,0.4,0),
                      color=vector(0.3,0.3,0.3))
    
        # Barras do suporte
        s1 = cylinder(pos=vector(3,5.2,0),radius=0.2, axis=vector(-4,0,0))
        s2 = cylinder(pos=vector(2.5,5.2,0),radius=0.6, axis=vector(-1,0,0),color=vector(0.5,0.5,0.6))
        s3 = cylinder(pos=vector(2,-5,0.4),radius=0.2, axis=vector(0,11,0))
        s4 = cylinder(pos=vector(2,-6.05,0.4),radius=0.8, axis=vector(0,1,0),color=vector(0.9,0.9,0.6))
    
        # Funcaoo que alonga/comprime a mola y1 unidades
        def alongar_mola(y1):
            c1.pos.y = y1 - 2
            bar1.origin = vector(0, y1, 0)
            mola1.axis.y = y1 - 5.8
    
        def uf(rf):          # velocidade no espaco de fase
            a = (-k1*rf.y-b1*rf.z)/m1
            return vector(0,rf.z, a)
    
        rf = vector(0, 0.9, 0)
        m1 = 0.3; k1 = 16; b1 = 0.02
        dt = 0.01
        alongar_mola(rf.y)
    
        while True:
            rate(100)
            uf1 = uf(rf)
            uf2 = uf(rf + dt*uf1 / 2.)
            uf3 = uf(rf + dt*uf2 / 2.)
            uf4 = uf(rf + dt*uf3)
            ufm = (uf1 + 2*uf2 + 2*uf3 + uf4)/6.
            rf += dt*ufm
            alongar_mola(rf.y)
            
if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    window = MyWindow()
    sys.exit(app.exec_())