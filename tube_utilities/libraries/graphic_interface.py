from Tkinter import *
import tkFileDialog
import matplotlib.pyplot as plot
import matplotlib.image as mpimg
from matplotlib import rcParams
from matplotlib import rc
import pickle
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import scipy.interpolate


rootWindowHandler = []
imageFilename = []
curvesCorners = []
curvesData = []
curvesVoltages = []
outputDataVar = []
max_voltage = 0
max_current = 0
yesNoReturnedValue = 0


class YesNoDialog:

    def __init__(self, parent, questionText, button1Text, button2Text):

        top = self.top = Toplevel(parent)

        Label(top, text=questionText).pack()

        b_ok = Button(top, text=button1Text, command=self.buttonYesPressedEvent)
        b_ok.pack(pady=5)

        b_cancel = Button(top, text=button2Text, command=self.buttonNoPressedEvent)
        b_cancel.pack(pady=5)

    def buttonYesPressedEvent(self):
        global yesNoReturnedValue
        yesNoReturnedValue = 1
        self.top.destroy()

    def buttonNoPressedEvent(self):
        global yesNoReturnedValue
        yesNoReturnedValue = 0
        self.top.destroy()


class InputDialog:

    def __init__(self, parent, title):

        top = self.top = Toplevel(parent)

        Label(top, text=title).pack()

        self.e = Entry(top)
        self.e.pack(padx=5)

        b = Button(top, text="OK", command=self.okButtonEvent)
        b.pack(pady=5)

    def okButtonEvent(self):

        receivedValue = self.e.get()
        try:
            parsedValue = float(receivedValue)
            global outputDataVar
            outputDataVar = parsedValue
            #curvesVoltages.extend([parsedValue])
            self.top.destroy()
        except:
            pass


def startGraphicInterface():

    global rootWindowHandler
    rootWindowHandler = Tk()

    # Allocates the buttons needed
    Button(rootWindowHandler, text="Capture curves from image", command=captureCurvesCallback).pack()
    Button(rootWindowHandler, text="Read curves from file", command=readCurvesFromFileCallback).pack()
    Button(rootWindowHandler, text="Computations from curves", command=computationsFromCurves).pack()
    
    # Gives title to the main window
    rootWindowHandler.title("Tube utilities - I. D. Simon, 2015.")

    # Shows the main window
    rootWindowHandler.mainloop()


def displayCurvesCapturing(filename):

    # Extracts the filename from the complete path
    onlyFilename = filename.split("/")[-1]

    global imageFilename
    imageFilename = filename

    # Interactive mode on to allow the rest of the program to continue
    plot.ion()
    updatePlots()

    # Sets title to the window
    fig = plot.gcf()
    fig.canvas.set_window_title("Acquiring curves from file: '" + onlyFilename + "'.")

    # Asks for maximum axes values
    global max_voltage
    global max_current
    d = InputDialog(rootWindowHandler, "Max. value on X-axis (V):")
    rootWindowHandler.wait_window(d.top)
    max_voltage = outputDataVar
    d = InputDialog(rootWindowHandler, "Max. value on Y-axis (mA):")
    rootWindowHandler.wait_window(d.top)
    max_current = outputDataVar

    plot.title("Step 1 - Select the four corners of the plot. [4 left]")

    # Wire the event for the click
    fig.canvas.mpl_connect('button_press_event', onCapturePlotClickedEvent)
    fig.canvas.mpl_connect('close_event', onCapturePlotWindowClosedEvent)

    # Initialize output variables
    global curvesCorners
    global curvesData
    global curvesVoltages
    curvesCorners = []
    curvesData = []
    curvesVoltages = []


def onCapturePlotWindowClosedEvent(event):

    # Cancel if there's nothing saved
    if curvesVoltages == []:
        return

    # Ask for filename and save
    _outputFile = tkFileDialog.asksaveasfilename(title="Save file to ...", defaultextension=".cur")
    if _outputFile != "":
        saveDataToFile(_outputFile, curvesData, curvesVoltages, curvesCorners, max_voltage, max_current)
        print "Data saved to: ", _outputFile


def onCapturePlotClickedEvent(event):

    if event.xdata is None or event.ydata is None:
        print "Position outside boundaries."
        return

    updatedPlotText = "Step 3 - Insert curves with mouse click.\nLeft: new point. Right: delete last point, Double click: finish curve."

    # Add the boundaries of the plot
    if len(curvesCorners) < 4:
        curvesCorners.extend([[event.xdata, event.ydata]])
        if len(curvesCorners) < 4:
            updatedPlotText = "Step 1 - Select the four corners of the plot. [" + str(4 - len(curvesCorners)) + " left]"

    # Add the curves
    else:

        # Left click, keep adding points
        if event.button == 1:
            if len(curvesData) == 0:
                curvesData.extend([[[event.xdata, event.ydata]]])
            else:
                # If it's double click add a new curve
                if event.dblclick == 1:
                    if not len(curvesData) == 0:
                        if not len(curvesData[-1]) == 0:
                            curvesData.extend([[]])

                            # Asks for the corresponding voltage
                            d = InputDialog(rootWindowHandler, "Voltage of inserted curve:")
                            rootWindowHandler.wait_window(d.top)
                            curvesVoltages.extend([outputDataVar]);

                # Otherwise just add the new point
                else:
                    curvesData[-1].extend([[event.xdata, event.ydata]])

        # Right click, delete last point
        if event.button == 3:
            if len(curvesData) > 0:
                # If the current curve is not closed
                if len(curvesData[-1]) > 0:
                    curvesData[-1].remove(curvesData[-1][-1])

                # Otherwise, remove the point from the previous curve
                else:
                    curvesData.remove(curvesData[-1])
                    if len(curvesData) > 0:
                        curvesData[-1].remove(curvesData[-1][-1])
                        curvesVoltages.remove(curvesVoltages[-1])

    updatePlots(updatedPlotText)


def xyCoordinatesToTwoVectors(xyCoordinates):

    xCoordinates = []
    yCoordinates = []

    for newPairOfCoordinates in xyCoordinates:
        xCoordinates.extend([newPairOfCoordinates[0]])
        yCoordinates.extend([newPairOfCoordinates[1]])

    return [xCoordinates, yCoordinates]


def updatePlots(updatedPlotText=""):

    plot.clf()

    _image = mpimg.imread(imageFilename)
    # Flip upside-down image to compensate for the inverted Y-axis
    _image = np.flipud(_image)
    imgplot = plot.imshow(_image)

    # Invert vertical axis to have the "normal" coordinates
    plot.gca().set_ylim(plot.gca().get_ylim()[::-1])

    # Removes the toolbar from the plot window
    rcParams['toolbar'] = 'None'

    # Removes the axis
    plot.axis('off')
    savedAxis = plot.axis()

    # Plot corners if ready
    if len(curvesCorners) == 4:
        [xCoords, yCoords] = xyCoordinatesToTwoVectors(curvesCorners)
        xCoords.extend([xCoords[0]])
        yCoords.extend([yCoords[0]])
        plot.plot(xCoords, yCoords, 'blue')
        xCoords.remove([xCoords[-1]])
        yCoords.remove([yCoords[-1]])

    # Plot all curves that are available
    if len(curvesData) > 0:
        for newCurve in curvesData:
            [xCoords, yCoords] = xyCoordinatesToTwoVectors(newCurve)
            plot.plot(xCoords, yCoords, marker='x', color='red')

    # Restores margins and axis
    plot.subplots_adjust(left=0.02, right=0.98, top=0.9, bottom=0.02)
    plot.axis(savedAxis)

    # Sets the title
    plot.title(updatedPlotText)

    # Redraws the plot
    plot.draw()


def askForInputFile(fileFilter="*.*"):
    _filename = tkFileDialog.askopenfilename(filetypes = [("Open file ...", fileFilter)])
    return _filename


def loadDataFromFile(inputFilename):
    '''
        Loads a saved file with captured curves. Returns
        data as two vectors of: ([[curves]], [voltages])
    '''

    with open(inputFilename, 'rb') as input:
        _loadedData1 = pickle.load(input)
        _loadedData2 = pickle.load(input)

    return _loadedData1, _loadedData2


def readCurvesFromFileCallback():

    # Ask for input file
    _filename = askForInputFile(fileFilter="*.cur")
    if _filename == "":
        return

    # Load data from file
    [_readCurves, _readVoltages] = loadDataFromFile(_filename)
    for i in range(0, len(_readVoltages)):
        _voltages = []
        _currents = []
        for j in range(0, len(_readCurves[i])):
            _voltages.extend([_readCurves[i][j][0]])
            _currents.extend([_readCurves[i][j][1]])

    # Plot read data
    onlyFilename = _filename.split("/")[-1]
    plotCurves(_readCurves, _readVoltages, onlyFilename, interpolate=True)

    # Ask if wants to save plot to PDF file
    d = YesNoDialog(rootWindowHandler, 'Save plot to PDF file?', 'Yes', 'No')
    rootWindowHandler.wait_window(d.top)
    if not yesNoReturnedValue:
        return
    _fileTypes = '*.*'
    _filename = askForOutputFilename(_fileTypes)
    if _filename == "":
        return
    _filename = fixExtensionOfFilename(_filename, 'pdf')
    pp = PdfPages(_filename)
    plot.savefig(pp, format='pdf')
    pp.close()
    print 'Saved curves to PDF file: "%s"' % _filename


def frange(start, end, increment):
    
    _stepsCount = int((end - start)/float(increment))
    _steps = [0] * _stepsCount

    i = 1
    _steps[0] = start
    while i < _stepsCount:
        _steps[i] = _steps[i-1] + increment
        i += 1
    
    return _steps


def plotCurves(curves, associatedVoltages, filename="", interpolate=False, interpolateDelta=0.5):

    # Removes the toolbar from the plot window
    rcParams['toolbar'] = 'None'
    fig = plot.gcf()
    plot.ion()
    plot.rc('text', usetex=True)
    plot.rc('font', family='serif')
    
    _interpolatedCurves = []

    for i in range(0, len(associatedVoltages)):
        _voltages = []
        _currents = []
        for j in range(0, len(curves[i])):
            _voltages.extend([curves[i][j][0]])
            _currents.extend([curves[i][j][1]])

        if interpolate:
            # Interpolation, obtaining (by default) values every 0.5 volts of Vak
            [_voltages, _currents] = interpolateCurve(_voltages, _currents, interpolateDelta, extendToZero=True)
            plot.plot(_voltages, _currents, '-', markersize=2, color='blue')
        else:
            plot.plot(_voltages, _currents, 'o-', markersize=2, color='blue')

        plot.text(_voltages[-1], _currents[-1], r'$V_g=' + str(associatedVoltages[i]) + '$', horizontalalignment='left', verticalalignment='bottom', rotation=45, fontsize=13)

        _interpolatedCurves.extend([[_voltages, _currents]])

    plot.grid()
    plot.xlabel(r'$V_{AK}$ [$V$]', fontsize=14)
    plot.ylabel(r'$I_{AK}$ [$mA$]', fontsize=14)
    plot.title(r'Set of curves $Vak/Ia=f(V_{grid})$')

    limits = plot.gca().get_ylim()
    plot.ylim([limits[0], limits[1]*1.2])
    limits = plot.gca().get_xlim()
    plot.xlim([limits[0], limits[1]*1.1])

    # Connect the close event of the plot window
    fig.canvas.mpl_connect('close_event', onFilePlotWindowClosedEvent)

    if filename != "":
        fig.canvas.set_window_title("Curves read from file: '%s'" % filename)
    plot.draw()
    plot.show()

    return _interpolatedCurves

def fixExtensionOfFilename(filename, extension):

    if len(filename) < len(extension)+1:
        _correctedFilename = filename + '.' + extension
    else:
        if filename[-(len(extension)+1):] != '.' + extension:
            _correctedFilename = filename + '.' + extension
        else:
            _correctedFilename = filename

    return _correctedFilename



def onFilePlotWindowClosedEvent(event):
    print "[Debug] Read curves from file plot window closed."


def askForOutputFilename(fileExtension="*.*"):
    _outputFilename = tkFileDialog.asksaveasfilename(title="Save file to ...", defaultextension=fileExtension)
    return _outputFilename


def askForInputFilename(_fileTypes=None):
    if _fileTypes == None:
        _fileTypes = [('All files', '*.*')]

    _openDialogHandler = tkFileDialog.Open(rootWindowHandler, filetypes=_fileTypes)
    
    _filename = _openDialogHandler.show()
    return _filename


def captureCurvesCallback():

    _fileTypes = [('Png image files', '*.png'), ('All files', '*')]
    _filename = askForInputFilename(_fileTypes)
    if _filename == "":
        return

    displayCurvesCapturing(_filename)


def saveDataToFile(outputFilename, curvesData, curvesVoltages, curvesCorners, maxVoltage, maxCurrent):
    '''
        Saves the captured data to a text file. Data is
        scaled from the captured pixels to volts / amps.
        '''
    
    # Find the origin and remove it from the list of corners
    _positionOrigin = curvesCorners[0]
    for i in range(0, len(curvesCorners)):
        if np.linalg.norm(curvesCorners[i]) < np.linalg.norm(_positionOrigin):
            _positionOrigin = curvesCorners[i]
    curvesCorners.remove(_positionOrigin)

    # Find the maximum and remove it
    _maxCorner = curvesCorners[0]
    for i in range(0, len(curvesCorners)):
        if np.linalg.norm(curvesCorners[i]) > np.linalg.norm(_maxCorner):
            _maxCorner = curvesCorners[i]
    curvesCorners.remove(_maxCorner)
    
    # Find limit of X axis
    _xLim = curvesCorners[0]
    for i in range(0, len(curvesCorners)):
        if curvesCorners[i][0] > _xLim[0]:
            _xLim = curvesCorners
    curvesCorners.remove(_xLim)
    _xLim = _xLim[0]

    # Find limit of Y axis
    _yLim = curvesCorners[0]
    for i in range(0, len(curvesCorners)):
        if curvesCorners[i][1] > _yLim[1]:
            _yLim = curvesCorners
    curvesCorners.remove(_yLim)
    _yLim = _yLim[1]

    # Save all data to a file
    _allCurves = []
    _allVoltages = []
    for _curveIndex in range(0, len(curvesVoltages)):
        _newCurve = []
        for i in range(len(curvesData[_curveIndex])):
            # Add new V/I point
            _v = (curvesData[_curveIndex][i][0] - _positionOrigin[0]) / (_xLim - _positionOrigin[0]) * maxVoltage
            _i = (curvesData[_curveIndex][i][1] - _positionOrigin[1]) / (_yLim - _positionOrigin[1]) * maxCurrent
            _newCurve.extend([[_v, _i]])
        _allCurves.extend([_newCurve])
        _allVoltages.extend([curvesVoltages[_curveIndex]])


    # Save all data to a file
    if len(outputFilename) < 4:
        _filenameToSave = outputFilename + '.cur'
    else:
        if outputFilename[-4:] != '.cur':
            _filenameToSave = outputFilename + '.cur'
        else:
            _filenameToSave = outputFilename
    with open(_filenameToSave, 'wb') as output:
        pickle.dump(_allCurves, output, pickle.HIGHEST_PROTOCOL)
        pickle.dump(_allVoltages, output, pickle.HIGHEST_PROTOCOL)


def generateLoadLine(Vmax, Imax, points=1000):

    # Obtains the set of X values (Voltages)
    _dataRange = frange(0, Vmax, float(Vmax)/float(points))

    _x_coordinates = [0]*len(_dataRange)
    _y_coordinates = [0]*len(_dataRange)

    # Slope of the load line
    m = Imax
    k = -(float(Imax)/float(Vmax))

    # Generate associated currents
    for i in range(0, len(_dataRange)):
        _x_coordinates[i] = _dataRange[i]
        _y_coordinates[i] = (k * _dataRange[i]) + m

    return [_x_coordinates, _y_coordinates]


def interpolateCurve(xValues, yValues, interpolateDelta, extendToZero=False):

    if extendToZero:
        # Need to add a first point that makes the curve converge to zero current
        # 1 - Find slope of first two points
        slope = (yValues[0] - yValues[1]) / (xValues[0] - xValues[1]);
        offset = yValues[0] - (slope * xValues[0]);
        newPoint_x = -offset / slope;
        newPoint_y = 0

        # 2 - Create two new lists with the new point at the beginning
        new_xValues = [0]*len(xValues)
        new_yValues = [0]*len(yValues)
        new_xValues[1:] = xValues
        new_yValues[1:] = yValues
        new_xValues[0] = newPoint_x
        new_yValues[0] = newPoint_y
        xValues = new_xValues
        yValues = new_yValues

    # Interpolation, obtaining (by default) values every 0.5 volts of Vak
    _interpolatedFunction = scipy.interpolate.UnivariateSpline(xValues, yValues, s=0)
    _voltages = frange(min(xValues), max(xValues), interpolateDelta)
    _currents = list(_interpolatedFunction(_voltages))

    return _voltages, _currents


def generateLoadLineFunction(Vmax, Imax, pointVoltage):
    # Slope of the load line
    m = Imax
    k = -(float(Imax)/float(Vmax))
    return k * float(pointVoltage) + m







def computationsFromCurves():
    _fileTypes = [('Captured curves files', '*.cur'), ('All files', '*')]
    _filename = askForInputFilename(_fileTypes)
    if _filename == "":
        return

    # Read data from file
    [_readCurves, _gridVoltages] = loadDataFromFile(_filename)

    # Working conditions of the circuit
    Vb = 400
    IBias = 1.1 / 1000
    VBias = -2.5
    Ra=150000
    
    # Computations of load line
    Imax= float(Vb) / float(Ra) * 1000 # in mA

    _intersectionPoints = []
    for i in range(len(_gridVoltages)):

        _currentCurveV, _currentCurveI = xyCoordinatesToTwoVectors(_readCurves[i])

        # Interpolate to increase the resolution and extend the data range to the cross with zero current
        [_currentCurveV, _currentCurveI] = interpolateCurve(_currentCurveV, _currentCurveI, interpolateDelta=1, extendToZero=True)

        # Obtain the interpolator function
        _interpolatorFunction = scipy.interpolate.UnivariateSpline(_currentCurveV, _currentCurveI, s=0)

        # Obtains the set of X values (Voltages)
        _dataRange = frange(min(_currentCurveV), max(_currentCurveV), float(max(_currentCurveV))/float(1000))

        # Find the closest "intersection" point
        _difference = 100000
        for _x in _dataRange:
            _newDifference = np.sqrt((float(_interpolatorFunction(_x)) - generateLoadLineFunction(Vb, Imax, _x))**2)
            if _newDifference < _difference:
                _difference = _newDifference
                _intersectionPoint = [_x, float(_interpolatorFunction(_x))]
        _intersectionPoints.extend([_intersectionPoint])


    # Debug:
    plotCurves(_readCurves, _gridVoltages, filename=_filename, interpolate=True, interpolateDelta=1)
    for i in range(len(_intersectionPoints)):
        plot.plot(_intersectionPoints[i][0], _intersectionPoints[i][1], 'o', color='green')
    _loadLine = generateLoadLine(Vb, Imax, points=3000)
    plot.plot(_loadLine[0], _loadLine[1], color='red')


    # Find the Gm values from the intersection points
    _valuesGm = []
    for i in range(len(_intersectionPoints)):
        if i != 0 and i!= len(_intersectionPoints)-1:
            _currentVgrid = _gridVoltages[i]
            _currentV = _intersectionPoints[i][0]

            # Previous curve - same V
            [_prevCurveV, _prevCurveI] = xyCoordinatesToTwoVectors(_readCurves[i-1])
            [_prevCurveV, _prevCurveI] = interpolateCurve(_prevCurveV, _prevCurveI, interpolateDelta=0.1, extendToZero=True)
            difference_prev = 1000000
            for j in range(len(_prevCurveV)):
                if np.sqrt((_prevCurveV[j] - _currentV)**2) < difference_prev:
                    difference_prev = np.sqrt((_prevCurveV[j] - _currentV)**2)
                    _cutI_prev = _prevCurveI[j]

            # Next curve - same V
            [_nextCurveV, _nextCurveI] = xyCoordinatesToTwoVectors(_readCurves[i+1])
            [_nextCurveV, _nextCurveI] = interpolateCurve(_nextCurveV, _nextCurveI, interpolateDelta=0.1, extendToZero=True)
            difference_next = 1000000
            for j in range(len(_nextCurveV)):
                if np.sqrt((_nextCurveV[j] - _currentV)**2) < difference_next:
                    difference_next = np.sqrt((_nextCurveV[j] - _currentV)**2)
                    _cutI_next = _nextCurveI[j]

            # Compute the Gm associated
            if (difference_prev < 10) and (difference_next < 10):
                _newGm = (_cutI_prev - _cutI_next) / (_gridVoltages[i-1] - _gridVoltages[i+1]);
                _valuesGm.extend([[_gridVoltages[i], _newGm]])

    # Find the Ri values from the slope of every curve on the intersection point
    _valuesRi = []
    for i in range(len(_intersectionPoints)):
        [_voltages, _currents] = xyCoordinatesToTwoVectors(_readCurves[i])
        [_voltages, _currents] = interpolateCurve(_voltages, _currents, interpolateDelta=0.2, extendToZero=True)
        for j in range(len(_voltages)):
            if np.abs(_voltages[j] - _intersectionPoints[i][0]) < 0.5:
                if j > 5 and j < len(_voltages)-5:
                    _newRi = (_voltages[j+5] - _voltages[j-5]) / (_currents[j+5] - _currents[j-5]) * 1000;
                    _valuesRi.extend([[_gridVoltages[i], _newRi]])
                break


    print _valuesRi
    print _valuesGm



    return

















