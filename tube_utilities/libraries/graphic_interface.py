from Tkinter import *
import tkFileDialog
import matplotlib.pyplot as plot
import matplotlib.image as mpimg
from matplotlib import rcParams

import numpy as np


rootWindowHandler = []
imageFilename = []
curvesCorners = []
curvesData = []
curvesVoltages = []


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
            curvesVoltages.extend([parsedValue])
            self.top.destroy()
        except:
            pass


def startGraphicInterface():

    global rootWindowHandler
    rootWindowHandler = Tk()

    # Allocates the buttons needed
    Button(rootWindowHandler, text="Import curves from image", command=captureCurvesCallback).pack()

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

    plot.title("Step 1 - Select the four corners of the plot. [4 left]")

    # Wire the event for the click
    fig.canvas.mpl_connect('button_press_event', onPlotClickedEvent)
    fig.canvas.mpl_connect('close_event', onPlotWindowClosedEvent)


def onPlotWindowClosedEvent(event):

    # TODO: ask user for the output filename
    # TODO: save this to a file
    print "Collected data:"
    print curvesData
    print curvesVoltages


def onPlotClickedEvent(event):

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
    imgplot = plot.imshow(_image)

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



def captureCurvesCallback():

    _fileTypes = [('Png image files', '*.png'), ('All files', '*')]
    _openDialogHandler = tkFileDialog.Open(rootWindowHandler, filetypes=_fileTypes)

    _filename = _openDialogHandler.show()
    if _filename == "":
        return

    displayCurvesCapturing(_filename)
