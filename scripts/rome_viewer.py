'''
@auther Yongbei.Ma
@date 2016.04.03
@version 1.1
@usage : make sure you have install numpy,scipy and PIL library
         python rome_viewer.py
@description

framework :               DataHandler           ************************
                              /\             \  **                    **
                             /||\ ------------\ **                    **
                              ||  ------------/ **    ImagesViewer    **
                             \||/            /  **                    **
                              \/                **                    **
**************************************************                    **********************
*****                         TopMenu           ************************            ********
********************************************************************************************
**                                                                                        **
**                                                                                        **
**                                  ClassesViewer                                         **
**                                                                                        **
**                                                                                        **
********************************************************************************************
*****                         StatusViewer                                          ********
********************************************************************************************

'''
import sys,os
import math
from PIL import Image, ImageTk
import numpy
from scipy import ndimage

try:
    # Python2
    import Tkinter as tk
    import tkFileDialog as filedialog
    import tkMessageBox
except ImportError:
    # Python3
    import tkinter as tk
    import tkinter.filedialog as filedialog
    from tkinter import messagebox as tkMessageBox

################### handle the mrcs and star data ############
### the data flow is : raw_data => filter_data => page_data => image_data ###
class DataHandler :
    # constructor
    def __init__(self):
        #---- flag whether the data is loaded or not ----#
        self.initFile()
        #---- mrcs head  ----#
        self.initMrcshead()
        #---- star data    ----#
        self.initMetadata()
        #---- ClassesViewer page information   ----#
        self.initClassPagedata()
        #---- ImageViewer page information ----#
        self.initImagePagedata()
        #---- filtered data,may remove the no-empty class and sorted by image number   ----#
        self.initFilteredClassdata()
        # ---- tk photo image class  ----#
        self.initClassdata()
        # ---- tk photo image image  ----#
        self.initImagedata()
        # ---- tk photo image image  ----#
        self.initFilteredImagedata()
        #---- state data   ----#
        self.initStatedata()

    ###################  initialization function #####################
    def initFile(self):
        self.fileHasStar            = False
        self.fileHasMrcs            = False
        self.fileFnStar             = "None"
        self.fileFnMrcs             = "None"
        self.fileFnPath             = "None"

    def initMrcshead(self):
        self.mrcshead_nx            = 0
        self.mrcshead_ny            = 0
        self.mrcshead_N             = 0
        self.mrcshead_mode          = 0

    def initMetadata(self):
        self.metaHead               = []
        self.metaData               = [] # star metadata line
        self.metaImageIndex         = []
        self.metaImageName          = []
        self.metaProb               = []
        self.metaClass              = []
        self.metaShiftx             = []
        self.metaShifty             = []
        self.metaRotation           = []

    def initClassPagedata(self):
        self.classPageNumber        = 0
        self.classPageIndex         = 0
        self.classPageClassesNumber = 0
        self.classPageClasses       = []
        self.classPageIndex_start   = 0
        self.classPageIndex_end     = 0

    def initImagePagedata(self):
        self.imagePageNumber        = 0
        self.imagePageIndex         = 0
        self.imagePageImagesNumber  = 0
        self.imagePageImages        = 0
        self.imagePageIndex_start   = 0
        self.imagePageIndex_end     = 0

    def initFilteredClassdata(self):
        self.filteredClasses        = [] # sorted class or no-empty class

    def initFilteredImagedata(self):
        self.filteredImages         = []  # filtered image by current selected class,the metadata index

    def initClassdata(self):
        self.classes                = []
        self.classesNumber          = 0
        self.classSize              = 0

    def initImagedata(self):
        self.images                 = []
        self.imagesNumber           = []
        self.imageSize              = []

    def initStatedata(self):
        self.stateClassSelectOrNot  = []  # flag the class selected or not
        self.stateClassImagesNumber = []  # the number of images each class has
        self.stateImageSelectOrNot  = []  # the image selected state,metadata index

    ############# load the mrcs or star data,if they have same name,load both ################
    def loadFile(self,fn):
        #---- set the file infomation ----#
        self.initFile()
        self.fileFnPath = os.path.dirname(fn)
        #---- open star file and its mrcs file to particle picking -----#
        if ".star" in fn :
            self.fileFnStar = fn
            self.fileFnMrcs = fn[0:-4] + "mrcs"
            # load star and mrcs file
            if os.path.exists(self.fileFnStar) and os.path.exists(self.fileFnMrcs) :
                # TODO : need check deeply whether read successly
                self.loadStarFile()
                self.fileHasStar = True
                self.loadMrcsFile()
                self.fileHasMrcs = True
            else :
                print("(DataHandler.loadFile)cannot find the ",self.fileFnStar," and ",self.fileFnMrcs," file.")
        #---- only open mrcs file to classes view ----#
        elif ".mrcs" in fn :
            # TODO : clear the star file
            self.fileFnMrcs = fn
            # load mrcs file
            if os.path.exists(self.fileFnMrcs) :
                # TODO : need check deeply whether read successly
                self.loadMrcsFile()
                self.fileHasMrcs = True
            else :
                print("(DataHandler.loadFile)cannot find the ",self.fileFnMrcs," file.")
        else :
            print("(DataHandler.loadFile)file name extension is incorrect.")
        #return flag
        self.resetStatedata()
        print("(DataHandler.loadFile)reset state data.")
        return (self.fileHasStar,self.fileHasMrcs)

    ############# load the star data to starData   ##################
    def loadStarFile(self):
        # NOTE : clear the list or dict,maybe not good to free memory
        self.initMetadata()
        #---- each metadata(class,prob,shiftx,shifty,rotation,image) column index ----#
        classIndex      = 0
        probIndex       = 0
        imageIndex      = 0
        shiftxIndex     = 0
        shiftyIndex     = 0
        rotationIndex   = 0
        #---- open and read data ----#
        starfile = open(self.fileFnStar,'r')
        lineIndex = 0
        for line in starfile:
            if '#' in line or lineIndex < 7:# set the starHead
                self.metaHead.append(line)
                if "_rlnclassnumber" in line.lower():
                    classIndex      = int(line[line.find("#")+1:])
                elif "_rlnmaxvalueprobdistribution" in line.lower():
                    probIndex       = int(line[line.find("#")+1:])
                elif "_rlnimagename" in line.lower():
                    imageIndex      = int(line[line.find("#") + 1:])
                elif "_rlnoriginx" in line.lower():
                    shiftxIndex     = int(line[line.find("#")+1:])
                elif "_rlnoriginy" in line.lower():
                    shiftyIndex     = int(line[line.find("#") + 1:])
                elif "_rlnanglepsi" in line.lower():
                    rotationIndex   = int(line[line.find("#") + 1:])
            elif len(line) > 3:# escape empty line,set the starData
                self.metaData.append(line)
                col = line.split()
                self.metaClass      .append(int(col[classIndex-1]))
                self.metaProb       .append(float(col[probIndex-1]))
                self.metaShiftx     .append(float(col[shiftxIndex - 1]))
                self.metaShifty     .append(float(col[shiftyIndex - 1]))
                self.metaRotation   .append(float(col[rotationIndex - 1]))
                self.metaImageIndex .append(int(col[imageIndex - 1].split('@')[0]))
                self.metaImageName  .append(col[imageIndex - 1].split('@')[1])
            lineIndex = lineIndex + 1
        print("load ",self.fileFnStar," file")
        starfile.close()

    ############# load the mrcs data file to mrcsData   ##################
    def loadMrcsFile(self):
        # NOTE : clear the mrcsData list,maybe not good to free memory
        self.initMrcshead()
        # open and read data
        mrcsfile = open(self.fileFnMrcs, 'rb')
        mrcsHead = numpy.fromfile(mrcsfile, dtype=numpy.int32, count=256)  # read the mrcs head data into numpy
        mrcsfile.close()
        self.mrcshead_nx    = mrcsHead[0]
        self.mrcshead_ny    = mrcsHead[1]
        self.mrcshead_N     = mrcsHead[2]
        self.mrcshead_mode  = mrcsHead[3]
        if self.mrcshead_nx != self.mrcshead_ny or self.mrcshead_mode != 2:
            print("wrong mrcs data file.")
        print("load mrcs data successfully")
        print("image dimesion x : ", self.mrcshead_nx)
        print("image dimesion y : ", self.mrcshead_ny)
        print("    image number : ", self.mrcshead_N)
        print(" image data mode : ", self.mrcshead_mode)
        return (self.mrcshead_nx,self.mrcshead_ny,self.mrcshead_N,self.mrcshead_mode)

    ############# get metadata for image   ##################
    def getMetadata(self,iimage):
        return (self.metaImageIndex[iimage],self.metaImageName[iimage],self.metaProb[iimage])

    ############   reset class filtered data    ############
    def resetFilteredClassdata(self,filterType):
        self.initFilteredClassdata()
        if filterType == "None":
            self.filteredClasses = [iclass+1 for iclass in range(self.mrcshead_N)]
        elif filterType == "NonEmpty":
            if self.fileHasStar : # only in particle select mode
                for iclass in range(self.mrcshead_N):
                    iclassplus1 = iclass + 1
                    if (iclassplus1 in self.metaClass) and (iclassplus1 not in self.filteredClasses):
                        self.filteredClasses.append(iclassplus1)
            else:
                print("(DataHandler.resetFilteredClassdata)can not no-empty class in gerneral view mode.")
                self.filteredClasses = [iclass + 1 for iclass in range(self.mrcshead_N)]
        elif filterType == "ByImageDist":
            if self.fileHasStar: # only in particle select mode
                imageDist = sorted(zip(range(self.mrcshead_N),self.stateClassImagesNumber),
                                   key=lambda x: x[1],reverse=True)
                for (iclass,N) in imageDist:
                    if N != 0:
                        self.filteredClasses.append(iclass+1)
            else:
                print("(DataHandler.resetFilteredClassdata)can not sort class by image number distribution in gerneral view mode.")
                self.filteredClasses = [iclass + 1 for iclass in range(self.mrcshead_N)]

    ###########  filter the image by selected class  ##############
    def resetFilteredImagedata(self,filterType):
        if filterType == "ByProb":
            probDist = sorted(zip(self.filteredImages,[self.metaProb[iimage] for iimage in self.filteredImages]),
                                key=lambda x: x[1], reverse=True)
            self.initFilteredImagedata()
            for (iimage, Prob) in probDist:
                self.filteredImages.append(iimage)
        else:
            iclassplus1 = filterType
            self.initFilteredImagedata()
            for iimage in range(len(self.metaData)):
                if self.metaClass[iimage] == iclassplus1:
                    self.filteredImages.append(iimage)

    ############        page infomation        ############
    def resetClassPagedata(self,classPageIndex,classPageClassNumber):
        self.initClassPagedata()
        self.classPageClassesNumber = classPageClassNumber
        self.classPageNumber = int(int(len(self.filteredClasses)-1) / int(self.classPageClassesNumber)) + 1
        if classPageIndex < self.classPageNumber:
            self.classPageIndex = classPageIndex
        else:
            print("(DataHandler.resetClassPagedata)wrong resetClassPagedata,classPageIndex out of range.",classPageIndex,self.classPageNumber)
        self.classPageIndex_start = self.classPageIndex*self.classPageClassesNumber
        self.classPageIndex_end   = min((self.classPageIndex+1)*self.classPageClassesNumber,len(self.filteredClasses))
        self.classPageClasses = self.filteredClasses[self.classPageIndex_start:self.classPageIndex_end]

    def getClassPageIndo(self):
        return (self.classPageIndex,self.classPageNumber,self.classPageClassesNumber)

    def getClassPageClasses(self):
        return self.classPageClasses

    def resetImagePagedata(self,imagePageIndex,imagePageImageNumber):
        self.initImagePagedata()
        self.imagePageImagesNumber = imagePageImageNumber
        self.imagePageNumber = int(int(len(self.filteredImages)-1) / int(self.imagePageImagesNumber)) + 1
        if imagePageIndex < self.imagePageNumber:
            self.imagePageIndex = imagePageIndex
        else:
            print("(DataHandler.resetImagePagedata)wrong resetImagePagedata,imagePageIndex out of range.", imagePageIndex,self.imagePageNumber)
        self.imagePageIndex_start = self.imagePageIndex * self.imagePageImagesNumber
        self.imagePageIndex_end = min((self.imagePageIndex + 1) * self.imagePageImagesNumber,len(self.filteredImages))
        print ("dataHander.resetImagePagedata,imagePageIndex_start = "+str(self.imagePageIndex_start)+",imagePageIndex_end = "+str(self.imagePageIndex_end)+",ImagePageNumber = "+str(self.imagePageNumber))
        self.imagePageImages = self.filteredImages[self.imagePageIndex_start:self.imagePageIndex_end]

    def getImagePageInfo(self):
        return (self.imagePageIndex,self.imagePageNumber,self.imagePageImagesNumber)

    def getImagePageImages(self):
        return self.imagePageImages

    ############ reset the state of class(selected ro not)  #######
    def resetStatedata(self):
        self.initStatedata()
        self.stateClassSelectOrNot = [False]*self.mrcshead_N
        self.stateClassImagesNumber = [0]*self.mrcshead_N
        self.stateImageSelectOrNot = [False]*len(self.metaClass)
        for iclassplus1 in self.metaClass:
            self.stateClassImagesNumber[iclassplus1-1] += 1

    def getStateClassSelectOrNot(self,iclassplus1):
        return self.stateClassSelectOrNot[iclassplus1-1]

    def setStateClassSelectOrNot(self,iclassplus1,state):
        self.stateClassSelectOrNot[iclassplus1 - 1] = state
        for iimage in range(len(self.stateImageSelectOrNot)):
            if(self.metaClass[iimage] == iclassplus1):
                self.stateImageSelectOrNot[iimage] = state

    def getStateImageSelectOrNot(self,iimage):
        return self.stateImageSelectOrNot[iimage]

    def setStateImageSelectOrNot(self,iimage,state):
        self.stateImageSelectOrNot[iimage]  =state

    def getStateClassImageNumber(self,iclassplus1):
        allimagesNumber = sum(self.stateClassImagesNumber)
        classImageNumber = self.stateClassImagesNumber[iclassplus1 - 1]
        if allimagesNumber != 0:
            return (classImageNumber,float(classImageNumber) / float(allimagesNumber))
        else:
            return (0,0.0)

    ############ reset and get Tkinter Photo Image  ############
    def resetClassImages(self):
        self.initClassdata()
        self.classSize = self.mrcshead_nx
        self.classesNumber = len(self.classPageClasses)
        classSize2 = self.classSize*self.classSize
        #---- open file and read data ----#
        mrcsfile = open(self.fileFnMrcs, 'rb')
        for iclassplus1 in self.classPageClasses:
            offset = 256*4 + (numpy.int64(iclassplus1)-1)*classSize2*4
            mrcsfile.seek(offset,os.SEEK_SET)
            # read the data into numpy
            classData = numpy.fromfile(mrcsfile, dtype=numpy.float32, count=classSize2).reshape(self.classSize,self.classSize)
            if classData.max() == 0 and classData.min() == 0:
                normalize = 0.0
            else:
                normalize = 255.0 / (classData.max() - classData.min())
            classDataInt = (normalize * (classData - classData.min())).astype(numpy.uint8)
            self.classes.append(Image.fromarray(classDataInt))
        mrcsfile.close()

    def getClassImages(self):
        return (self.classes,self.classSize,self.classesNumber)

    def resetImageImages(self):
        self.initImagedata()
        self.imageSize = self.mrcshead_nx
        self.imagesNumber = len(self.imagePageImages)
        imageSize2 = self.imageSize*self.imageSize
        #---- open mrcs one by one to read image data ----#
        for iimage in self.imagePageImages:
            mrcsfile = open(self.fileFnPath+"/"+self.metaImageName[iimage],'rb')
            offset = 256*4 + (numpy.int64(self.metaImageIndex[iimage])-1)*imageSize2*4
            mrcsfile.seek(offset,os.SEEK_SET)
            # read the data into numpy
            imageData = numpy.fromfile(mrcsfile, dtype=numpy.float32, count=imageSize2).reshape(self.imageSize,self.imageSize)
            # NOTE :  shift and rotate the image
            imageData = ndimage.shift(imageData,(self.metaShifty[iimage],self.metaShiftx[iimage]),order=2,mode='wrap')
            #imageData = ndimage.shift(imageData, (0, 0), order=1,mode='wrap')
            imageData = ndimage.rotate(imageData, -self.metaRotation[iimage], reshape=False, order=3, mode='wrap')
            # donot consider about zero case
            normalize = 255.0 / (imageData.max() - imageData.min())
            imageDataInt = (normalize * (imageData - imageData.min())).astype(numpy.uint8)
            self.images.append(Image.fromarray(imageDataInt))
            mrcsfile.close()

    def getImageImages(self):
        return (self.images,self.imageSize,self.imagesNumber)

    ###########   save selected images to star file  ##########
    def saveStarFile(self,fn_star):
        starfile = open(fn_star,'w')
        starfile.write("".join(self.metaHead))
        selectMetaData = filter(lambda e:e[0] == True,zip(self.stateImageSelectOrNot,self.metaData))
        starfile.write("".join([e[1] for e in selectMetaData]))
        print("(DataHandler.saveStarFile)save star "+fn_star)
        starfile.close()

    ############  save mrcs data file as tiff file
    def saveTiffFile(self,fn_tiff):
        for classplus1_index in range(len(self.filteredClasses)):
            classplus1 = self.filteredClasses[classplus1_index]
            if self.getStateClassSelectOrNot(classplus1):
                save_fn = fn_tiff+"_class_"+str(classplus1)+".tiff"
                self.classes[classplus1_index].save(save_fn)
                print("save tiff file name : "+save_fn)

############# status viewer class,each time do some event it will show some info in bottom status bar #############
class StatusViewer:
    def __init__(self,root,defalutText):
        self.root = root
        self.statusBar = tk.Label(root,text = defalutText,bd = 1,relief = tk.SUNKEN,anchor = tk.W)
        self.statusBar.pack(side = tk.BOTTOM,fill = tk.X)

    #update status
    def showMessage(self,status):
        self.statusBar.configure(text = status,fg = "black")
    def showError(self,status):
        self.statusBar.configure(text = status,fg = "red")

    # NOTE : using call back to implement the refresh bar
    # increase from 0 to maximum
    def setupProgress(self,increasedFunc,maximum):
        # for grogress bar
        self.increasedFunc = increasedFunc
        self.maximum = maximum
        self.increase = 0

    def refreshProgress(self):
        progress_percentage = int(math.floor(float(self.increase+1)/float(self.maximum)*100.))
        progress_bar_length = int(progress_percentage*0.8)
        status = "( Load : "+str(progress_percentage)+"%)"+"||"*progress_bar_length
        self.statusBar.configure(text = status,fg = "blue")

        self.increasedFunc(self.increase)
        self.increase += 1

        if self.increase < self.maximum :
            self.root.after(1, self.refreshProgress)

############# top Menu #############
class TopMenu:
    # constructor
    def __init__(self, root, classesViewer, imagesViewer, statusViewer, dataHandler):
        #---- initilaize some class varibale ----#
        self.classesViewer = classesViewer
        self.imagesViewer = imagesViewer
        self.statusViewer = statusViewer
        self.dataHandler = dataHandler
        #---- create top menu ----#
        topMenu = tk.Menu(root)
        root.config(menu = topMenu)
        #----  create file menu  ----#
        fileMenu = tk.Menu(topMenu,tearoff=0)
        topMenu.add_cascade(label = "File",menu = fileMenu)
        fileMenu.add_command(label="Open STAR File(particle selection)",command = self.fileOpenStarFile)
        fileMenu.add_command(label="Open MRCS File(general view)",command = self.fileOpenMrcsFile)
        fileMenu.add_separator()
        fileMenu.add_command(label="Save Selected Classes into STAR File",command = self.fileSaveStarFile)
        fileMenu.add_command(label="Save Selected Classes to TIFF File",command = self.fileSaveTiffFile)

        #---- create edit menu  ----#
        editMenu = tk.Menu(topMenu,tearoff=0)
        topMenu.add_cascade(label = "Edit",menu = editMenu)
        editMenu.add_command(label = "Select All",command = self.editSelectAll)
        editMenu.add_command(label = "Clear Selection",command = self.editUnselectAll)
        editMenu.add_command(label = "Reverse Selection",command = self.editReverseSelect)

        #----  crate view menu  ----#
        viewMenu = tk.Menu(topMenu,tearoff=0)
        topMenu.add_cascade(label = "View",menu = viewMenu)
        #----  for flitering the data  ----#
        self.filteredMode = tk.StringVar()
        self.filteredMode.set("None")
        viewMenu.add_radiobutton(label = "Show All Classes",value = "None",variable = self.filteredMode,
                                command = self.viewShowAll)
        viewMenu.add_radiobutton(label = "Show Non-Empty Classes",value = "NonEmpty",variable = self.filteredMode,
                                command = self.viewShowNonEmptyClass)
        viewMenu.add_radiobutton(label = "Sort Classes by Images Population",
                                value = "ByImageDist",variable = self.filteredMode,
                                command = self.viewSortByImagesDistribution)
        viewMenu.add_separator()
        #----  page control  ----#
        viewMenu.add_command(label = "Next Page",command = self.viewNextPage)
        viewMenu.add_command(label = "Previous Page",command = self.viewPreviousPage)
        viewMenu.add_command(label = "First Page", command=self.viewFirstPage)
        viewMenu.add_command(label = "Last Page", command=self.viewLastPage)
        classNumberMenu = tk.Menu(topMenu,tearoff=0)
        self.classPageClassesNumber = tk.IntVar()
        self.classPageClassesNumber.set(100)
        viewMenu.add_cascade(label = "Each Page Classes Number",menu = classNumberMenu)
        classNumberMenu.add_radiobutton(label = "100",value = 100,variable = self.classPageClassesNumber,
                                            command = self.viewClassNumberPage)
        classNumberMenu.add_radiobutton(label = "300",value = 300,variable = self.classPageClassesNumber,
                                            command = self.viewClassNumberPage)
        classNumberMenu.add_radiobutton(label = "500",value = 500,variable = self.classPageClassesNumber,
                                            command = self.viewClassNumberPage)
        classNumberMenu.add_radiobutton(label = "1000",value = 1000, variable=self.classPageClassesNumber,
                                        command=self.viewClassNumberPage)
        viewMenu.add_separator()
        #----  scalable and scrolled view  ----#
        self.viewingMode = tk.IntVar()
        self.viewingMode.set(0)
        #----  scalable view  ----#
        viewMenu.add_radiobutton(label = "Scalable Viewing Mode",value = 0,variable = self.viewingMode,
                                command = self.viewScalableView)
        # scrolled view
        # NOTE : value from 1~5 take as the class scale
        scrolledViewMenu = tk.Menu(topMenu,tearoff=0)
        viewMenu.add_cascade(label = "Scrolled Viewing Mode",menu = scrolledViewMenu)
        scrolledViewMenu.add_radiobutton(label = "1 / 1 Scale",value = 1,variable = self.viewingMode,
                                            command = self.viewScrolledView)
        scrolledViewMenu.add_radiobutton(label = "1 / 2 Scale",value = 2,variable = self.viewingMode,
                                            command = self.viewScrolledView)
        scrolledViewMenu.add_radiobutton(label = "1 / 3 Scale",value = 3,variable = self.viewingMode,
                                            command = self.viewScrolledView)
        scrolledViewMenu.add_radiobutton(label = "1 / 4 Scale",value = 4,variable = self.viewingMode,
                                            command = self.viewScrolledView)
        scrolledViewMenu.add_radiobutton(label = "1 / 5 Scale",value = 5,variable = self.viewingMode,
                                            command = self.viewScrolledView)
        scrolledViewMenu.add_radiobutton(label = "1 / 6 Scale",value = 6,variable = self.viewingMode,
                                            command = self.viewScrolledView)

        #----  crate windows menu  ----#
        winMenu = tk.Menu(topMenu,tearoff=0)
        topMenu.add_cascade(label = "Windows",menu = winMenu)
        winMenu.add_command(label = "ImageViewer",command = self.winImageViewer)

        #----  crate help menu  ----#
        helptMenu = tk.Menu(topMenu,tearoff=0)
        topMenu.add_cascade(label = "Help",menu = helptMenu)
        helptMenu.add_command(label = "Help",command = self.helpHelp)
        helptMenu.add_command(label = "About",command = self.helpAbout)

    ###################       file bound event         #####################

    #############  tell datahander and classes viewer to reset  #############
    def tellDataHanderReset(self):
        filteredMode = self.filteredMode.get()
        self.dataHandler.resetFilteredClassdata(filteredMode)
        # NOTE : the page class is base on the filtered class data
        classPageClassesNumber = self.classPageClassesNumber.get()
        self.dataHandler.resetClassPagedata(0, classPageClassesNumber)
        # set up and running the generate tk photo images for class
        self.dataHandler.resetClassImages()

    def tellClassViewerReset(self):
        # reset the class viewer
        viewingMode = self.viewingMode.get()
        if viewingMode == 0:
            self.classesViewer.resetClassesViewer("scalable")
        else:
            self.classesViewer.resetClassesViewer("scrolled",viewingMode)

    #############  open star file  #############
    def fileOpenStarFile(self):
        filetypes = (("star files","*.star"),("all files","*.*"))
        fn = filedialog.askopenfilename(initialdir = "./",title = "Choose your file",filetypes = filetypes)
        # if fn != "":
        (loadStar,loadMrcs) = self.dataHandler.loadFile(fn)
        if (loadStar,loadMrcs) == (True,True) :
            # TODO : these message donnot work
            self.statusViewer.showMessage("Load Data.....waiting.........")
            # reset the state,class select and class image number distribution
            self.tellDataHanderReset()
            # reset class data and class viewer
            self.tellClassViewerReset()
            self.statusViewer.showMessage("Load "+fn[0:-5]+".star and .mrcs data successly.")
        else :
            self.statusViewer.showError("Cannot open file(make sure your star and mrcs have same name and directory).")

    #############  open mrcs file  #############
    def fileOpenMrcsFile(self):
        filetypes = (("mrcs files","*.mrcs"),("all files","*.*"))
        fn = filedialog.askopenfilename(initialdir = "./",title = "choose your file",filetypes = filetypes)
        (loadStar,loadMrcs) = self.dataHandler.loadFile(fn)
        # if fn != "":
        if (loadStar,loadMrcs) == (False,True) :
            # TODO : these message donnot work
            self.statusViewer.showMessage("Load Data.....waiting.........")
            # reset the state,class select and class image number distribution
            self.tellDataHanderReset()
            # reset class data and class viewer
            self.tellClassViewerReset()
            self.statusViewer.showMessage("Load "+fn[0:-5]+".mrcs data successly.")
        else :
            self.statusViewer.showError("Open mrcs file failed.")

    # save the selected class to star file
    def fileSaveStarFile(self):
        filetypes = (("star files","*.star"),("all files","*.*"))
        fn = filedialog.asksaveasfilename(initialdir = "./",title = "Save your file",filetypes = filetypes)
        # print(fn)
        if fn != "":
            self.dataHandler.saveStarFile(fn)
        self.statusViewer.showMessage("Save star file.")

    def fileSaveTiffFile(self):
        filetypes = (("tiff files", "*.tiff"), ("all files", "*.*"))
        fn = filedialog.asksaveasfilename(initialdir="./", title="Save your file", filetypes=filetypes)
        # print(fn)
        if fn != "":
            self.dataHandler.saveTiffFile(fn)
        self.statusViewer.showMessage("Save tiff file.")

    ###################  edit bound event  #####################
    def editSelectAll(self):
        classPageClasses = dataHandler.getClassPageClasses()
        for iclassplus1 in classPageClasses:
            dataHandler.setStateClassSelectOrNot(iclassplus1,True)
        self.classesViewer.updateLabelsBG()
        self.statusViewer.showMessage("Select all classes (images)")

    def editUnselectAll(self):
        classPageClasses = dataHandler.getClassPageClasses()
        for iclassplus1 in classPageClasses:
            dataHandler.setStateClassSelectOrNot(iclassplus1,False)
        self.classesViewer.updateLabelsBG()
        self.statusViewer.showMessage("Unselect all classes (images)")

    def editReverseSelect(self):
        classPageClasses = dataHandler.getClassPageClasses()
        for iclassplus1 in classPageClasses:
            selectOrNot = dataHandler.getStateClassSelectOrNot(iclassplus1)
            dataHandler.setStateClassSelectOrNot(iclassplus1,not selectOrNot)
        self.classesViewer.updateLabelsBG()
        self.statusViewer.showMessage("Reverse select all classes (images)")

    ###################  view bound event  #####################
    def viewShowAll(self):
        self.tellDataHanderReset()
        self.tellClassViewerReset()
        self.statusViewer.showMessage("Show all classes (images).")

    # show all not empty classes
    def viewShowNonEmptyClass(self):
        self.tellDataHanderReset()
        self.tellClassViewerReset()
        self.statusViewer.showMessage("Remove all empty classes.")

    # sort the classes by images distribution.
    def viewSortByImagesDistribution(self):
        self.tellDataHanderReset()
        self.tellClassViewerReset()
        self.statusViewer.showMessage("Sort the classes by images populations.")

    # update the page #
    def viewSetPage(self,newClassPageIndex,classPageNumber):
        if newClassPageIndex < classPageNumber and newClassPageIndex >= 0:
            self.dataHandler.resetClassPagedata(newClassPageIndex, self.classPageClassesNumber.get())
            self.dataHandler.resetClassImages()
            self.tellClassViewerReset()
            self.statusViewer.showMessage("Show next page : "+str(newClassPageIndex+1)+" / "+str(classPageNumber))
        else:
            self.statusViewer.showError("At end of the pages : "+str(newClassPageIndex)+" / "+str(classPageNumber))

    def viewNextPage(self):
        (classPageIndex,classPageNumber,classPageClassesNumber) = self.dataHandler.getClassPageIndo()
        classPageIndex += 1
        self.viewSetPage(classPageIndex,classPageNumber)

    def viewPreviousPage(self):
        (classPageIndex, classPageNumber, classPageClassesNumber) = self.dataHandler.getClassPageIndo()
        classPageIndex -= 1
        self.viewSetPage(classPageIndex, classPageNumber)

    def viewFirstPage(self):
        (classPageIndex, classPageNumber, classPageClassesNumber) = self.dataHandler.getClassPageIndo()
        classPageIndex = 0
        self.viewSetPage(classPageIndex, classPageNumber)

    def viewLastPage(self):
        (classPageIndex, classPageNumber, classPageClassesNumber) = self.dataHandler.getClassPageIndo()
        classPageIndex = classPageNumber - 1
        self.viewSetPage(classPageIndex, classPageNumber)

    def viewClassNumberPage(self):
        classPageClassesNumber = self.classPageClassesNumber.get()
        self.dataHandler.resetClassPagedata(0, classPageClassesNumber)
        self.dataHandler.resetClassImages()
        self.tellClassViewerReset()
        self.statusViewer.showMessage("Show "+str(classPageClassesNumber)+" classes each page.")

    def viewScalableView(self):
        self.tellClassViewerReset()
        self.statusViewer.showMessage("Change to scalable mode(your can adjust the class size by changing the windows size).")

    def viewScrolledView(self):
        self.tellClassViewerReset()
        viewingMode = self.viewingMode.get()
        self.statusViewer.showMessage("Scrolled mode,class downsize to 1/"+str(viewingMode)+".")

    def winImageViewer(self):
        self.imagesViewer.setupIamgeViewer()

    def helpHelp(self):
        print("TODO,jump to web")

    def helpAbout(self):
        tkMessageBox.showinfo("About", "ROME 1.0 ")

#*** the adjustable class thumb viewer ***#
class ClassesViewer:
    ###################       constructor                        #######################################
    def __init__(self, root, imagesViewer, statusViewer, dataHandler):
        self.root = root
        self.statusViewer = statusViewer
        self.imagesViewer = imagesViewer
        self.dataHandler = dataHandler
        # two mode scalable or scrolled
        self.classesThumbLabel = []
        # create the classes show Frame
        self.FrameHeight = self.root.winfo_screenheight() / 2
        self.FrameWidth = self.root.winfo_screenwidth() / 2
        self.borderwidth = 2
        print("(classesViewer.__init__)classesViewerFrame : Width = ", self.FrameWidth, ",Height = ", self.FrameHeight)
        # default is scalable mode
        self.viewingMode = "scalable"
        self.currentSelectClass = 0
        self.setupScalableLayout()

    #############  tell datahander and imageViewer to reset  #############
    def tellDataHanderAndImageViewerReset(self,selectedClass):
        # when the imageviewer is turn on
        if self.imagesViewer.getImageViewerState() != False :
            self.dataHandler.resetFilteredImagedata(selectedClass)
            # reset the labels
            self.imagesViewer.resetLabels()

    ###################     viewingMode = "scalable"(defalut),"scrolled&&scale"   #####################
    def resetClassesViewer(self,newViewingMode, newScale=None):
        # 1)clear all the labels
        self.clearLables()
        # 2)reset layout if needed
        if newViewingMode == "scalable" and self.viewingMode == "scalable":
            print("(ClassesViewer.resetClassesViewer)set scalable view mode.")
        elif newViewingMode == "scalable" and self.viewingMode == "scrolled":
            self.viewingMode = newViewingMode
            self.clearScrolledLayout()
            self.setupScalableLayout()
        elif newViewingMode == "scrolled" and self.viewingMode == "scalable":
            self.viewingMode = newViewingMode
            self.scale = newScale
            self.clearScalableLayout()
            self.setupScrolledLayout()
        # NOTE : this means in scrolled mode but change the scale
        elif newViewingMode == "scrolled" and newScale != self.scale:
            self.scale = newScale
            # do not need to clear and reset the main layout
            # NOTE : this will not urefresh the windows,when I scroll bottom down
            # add clear and reset back
            self.clearScrolledLayout()
            self.setupScrolledLayout()
        elif newViewingMode == "scrolled" and self.viewingMode == "scrolled":
            print("(ClassesViewer.resetClassesViewer)set scrolled view mode.")
        else:
            print("(ClassesViewer.resetClassesViewer)not implement this viewing mode yet.")
            print("(ClassesViewer.resetClassesViewer) : " + newViewingMode)
        # 3)setup classes thunmb,donot need to clear,old one will be replace
        # TODO : keep track of filtered mode,if it is not need to change
        self.setupClassesThumb()
        # 4)set up labels
        self.setupLabels()

    ####################################################################################################
    #############                      set up different layout                     #####################
    ####################################################################################################
    #############  set up the scalable layout,prepare Main required tk weights for these Layout  ######
    def setupScalableLayout(self):
        # add the scalable Frame
        self.classesFrame = tk.Frame(self.root, bd=self.borderwidth, bg="white",
                                    height=self.FrameHeight, width=self.FrameWidth)
        # fill and expand the Frame in whole root windows
        self.classesFrame.pack(fill=tk.BOTH, expand=True)
        # set up pop out menu,TODO : how to remove this in clearLayout?
        self.popMenu = tk.Menu(self.root,tearoff=0)
        self.popMenu.add_command(label="Unselect", command=self.eventUnselectClass)
        # add event
        self.classesFrame.bind("<Configure>", self.eventAdjustScalableLayout)

    # NOTE : call this after clearLabels()
    def clearScalableLayout(self):
        self.classesFrame.pack_forget()
        self.classesFrame.destroy()

    #############  set up the scrolled layout,prepare Main required tk weights for these Layout  #######
    def setupScrolledLayout(self):
        # setup canvas and put classes frame on it
        # self.scrollCanvas = tk.Canvas(self.root, borderwidth = 0, background = "red")
        # self.classesFrame = tk.Frame(self.scrollCanvas, background = "blue")
        self.scrollCanvas = tk.Canvas(self.root, borderwidth=0, background="white")
        self.classesFrame = tk.Frame(self.scrollCanvas, background="white")
        # set up pop out menu,TODO : how to remove this in clearLayout?
        self.popMenu = tk.Menu(self.root,tearoff=0)
        self.popMenu.add_command(label="Unselect", command=self.eventUnselectClass)
        # add scrollbar and bind this with canvas
        self.scrollBar = tk.Scrollbar(self.root, orient="vertical", command=self.scrollCanvas.yview)
        self.scrollCanvas.configure(yscrollcommand=self.scrollBar.set)
        # pack the canvas and scrollbar
        self.scrollBar.pack(side="right", fill="y")
        self.scrollCanvas.pack(side="left", fill="both", expand=True)
        self.scrollCanvas.create_window((0, 0), window=self.classesFrame, anchor="nw")
        # add event
        self.classesFrame.bind("<Configure>", self.eventAdjustScrolledLayout)
        self.scrollCanvas.bind_all("<MouseWheel>", self.eventMouseWheelScrolledLayout)

    # NOTE : call this after clearLabels()
    def clearScrolledLayout(self):
        self.classesFrame.destroy()
        self.scrollBar.pack_forget()
        self.scrollBar.destroy()
        self.scrollCanvas.pack_forget()
        self.scrollCanvas.destroy()

    ####################################################################################################
    #############                      label configuration                         #####################
    ####################################################################################################
    #############             set up the classes thumbs data                               #############
    def setupClassesThumb(self):
        # 1) get the original data from different data filtered mode
        (self.classes, self.classSize,self.classesNumber) = self.dataHandler.getClassImages()
        self.classesThumb = [None] * self.classesNumber
        # 2)get the appropriate class thumb size and grid dimension info
        if self.viewingMode == "scrolled":
            self.getScrolledLayoutInfo()
        else:
            self.getScalableLayoutInfo()
        # 3)set the classes thumbs by scale
        self.updateClassesThumb()

    #############  update classes thumb      #############
    def updateClassesThumb(self):
        # need to get the copy of class data
        if self.upORdown == "down":
            newSize = int(self.classSize / self.scale)
            for iclass in range(self.classesNumber):
                # NOTE : in Tkinter use subsample,in PIL use resize
                self.classesThumb[iclass] = ImageTk.PhotoImage(
                                            self.classes[iclass].resize((newSize,newSize),resample=Image.NEAREST))
        elif self.upORdown == "up":
            newSize = int(self.classSize * self.scale)
            for iclass in range(self.classesNumber):
                self.classesThumb[iclass] = ImageTk.PhotoImage(
                                            self.classes[iclass].resize((newSize, newSize), resample=Image.NEAREST))

    #############  set up all the class thumbs labels  #############
    def setupLabels(self):
        # create class thumb label
        classPageClasses = self.dataHandler.getClassPageClasses()
        for iclass in range(self.classesNumber):
            selectOrNot = self.dataHandler.getStateClassSelectOrNot(classPageClasses[iclass])
            if selectOrNot:
                self.classesThumbLabel.append(tk.Label(self.classesFrame, image=self.classesThumb[iclass],
                                                     text=classPageClasses[iclass], bd=self.borderwidth, bg="red"))
            else:
                self.classesThumbLabel.append(tk.Label(self.classesFrame, image=self.classesThumb[iclass],
                                                      text=classPageClasses[iclass],bd=self.borderwidth, bg="white"))
        # update label layout
        self.updateLabelsGrid()
        # add event
        for label in self.classesThumbLabel:
            label.bind("<Button-1>", self.eventSelectClassThumb)
            # NOTE : the right clicked Menu is 2 on mac and 3 on unix and windows
            label.bind("<Button-2>", self.eventPopMenu)
            label.bind("<Button-3>", self.eventPopMenu)

    #############  clear all labels    #############
    def clearLables(self):
        # TODO : is this efficient to remove all the label?
        # print("remove all label.")
        for label in self.classesThumbLabel:
            label.grid_forget()  # just make invisiable
            label.destroy()
        self.classesThumbLabel = []

    #############  update labels' classess  #############
    def updateLabelsClass(self):
        for iclass in range(0, self.classesNumber):
            self.classesThumbLabel[iclass].configure(image=self.classesThumb[iclass])

    #############  put or update the label grid layout on Frame  #############
    def updateLabelsGrid(self):
        iclass = 0
        for label in self.classesThumbLabel:
            Icol = iclass % int(self.gridColNum)
            Irow = int((iclass - Icol) / int(self.gridColNum))
            # NOTE : sticky=tk.W+tk.E not work here,padx = self.gridRowPad,pady = self.gridColPad
            label.grid(row=Irow, column=Icol)
            iclass = iclass + 1

    #############  update label background color           #############
    def updateLabelBG(self, label):
        if label.cget('bg') == "white":
            label.configure(bg="red")
        else:
            label.configure(bg="white")

    #############     update all labels background color  #############
    ### NOTE : this is used for TopMenu call
    def updateLabelsBG(self):
        classPageClasses = self.dataHandler.getClassPageClasses()
        for iclass in range(self.classesNumber):
            selectOrNot = self.dataHandler.getStateClassSelectOrNot(classPageClasses[iclass])
            if selectOrNot:
               self.classesThumbLabel[iclass].configure(bg="red")
            else:
                self.classesThumbLabel[iclass].configure(bg="white")

    ####################################################################################################
    #############                      all events                                  #####################
    ####################################################################################################
    #############       adjust the Frame Layout                                            #############
    def eventAdjustScalableLayout(self, event):
        print("((ClassesViewer.eventAdjustscalableLayout)AdjustLayout,Classes viewer Frame width = ", event.width, ",height = ", event.height)
        # NOTE : if we just open main windows and not open any file,return this
        # these mean we donnot get layout info and set label yet
        if (not hasattr(self, 'gridRowNum')):
            self.FrameHeight = event.height
            self.FrameWidth = event.width
            return
        # NOTE : each 21 step update,donot update frequence,high frequence will cause problem
        if event.width % 5 == 0 or event.height % 5 == 0:
            # record the old layout info
            old_gridRowNum = self.gridRowNum
            old_gridColNum = self.gridColNum
            old_scale = self.scale
            old_upORdown = self.upORdown
            # get new layout info
            self.FrameHeight = event.height
            self.FrameWidth = event.width
            self.getScalableLayoutInfo()
            # update classes thumb
            if old_scale != self.scale or old_upORdown != self.upORdown:
                print("(ClassesViewer.eventAdjustScalableLayout)AdjustLayout,update class thumb size")
                for label in self.classesThumbLabel:
                    label.grid_forget()
                # update label class
                self.updateClassesThumb()
                # update label classes
                self.updateLabelsClass()
                # update label layout
                self.updateLabelsGrid()
            elif old_gridRowNum != self.gridRowNum or old_gridColNum != self.gridColNum:
                print("(ClassesViewer.eventAdjustScalableLayout)AdjustLayout,update classes grid.")
                for label in self.classesThumbLabel:
                    label.grid_forget()
                self.updateLabelsGrid()

    def eventAdjustScrolledLayout(self, event):
        self.scrollCanvas.configure(scrollregion=self.scrollCanvas.bbox("all"))

    def eventMouseWheelScrolledLayout(self, event):
        self.scrollCanvas.yview_scroll(-1 * (event.delta / 120), "units")

    # select class event
    def eventSelectClassThumb(self, event):
        # get the class thumb state first
        iclassplus1 = event.widget.cget("text")
        selectOrNot = self.dataHandler.getStateClassSelectOrNot(iclassplus1)
        #--------  if image viewer is not open   -----------#
        if self.imagesViewer.getImageViewerState() == False :
            self.updateLabelBG(event.widget)
            # if the class already has been selected, unselect the class,else select the class
            self.dataHandler.setStateClassSelectOrNot(iclassplus1, not selectOrNot)
        else : # if the image viewer is open
            if not selectOrNot :
                self.updateLabelBG(event.widget)
                self.dataHandler.setStateClassSelectOrNot(iclassplus1, True)
            self.tellDataHanderAndImageViewerReset(iclassplus1)
        self.currentSelectClass = iclassplus1
        (imagesDistribution, probability) = self.dataHandler.getStateClassImageNumber(iclassplus1)
        statusMessage = "the " + str(iclassplus1 ) + "th class."
        statusMessage += "class Info : images = " + str(imagesDistribution) + ",probability = " + str(probability)
        self.statusViewer.showMessage(statusMessage)

    def eventPopMenu(self,event):
        self.popMenu.post(event.x_root,event.y_root)
        self.currentSelectClass = event.widget.cget("text")

    def eventUnselectClass(self):
        iclassplus1 = self.currentSelectClass
        selectOrNot = self.dataHandler.getStateClassSelectOrNot(iclassplus1)
        # when you have select this class
        if selectOrNot:
            self.dataHandler.setStateClassSelectOrNot(iclassplus1, False)
            self.updateLabelsBG()
        (imagesDistribution, probability) = self.dataHandler.getStateClassImageNumber(iclassplus1)
        statusMessage = "the " + str(iclassplus1) + "th class."
        statusMessage += "class Info : images = " + str(imagesDistribution) + ",probability = " + str(probability)
        self.statusViewer.showMessage(statusMessage)

    ###################################################################################################
    ############# get the layout info ,how to put the class thumb in class viewer  ####################
    ###################################################################################################
    #############  get the info for scalable frame grid when user adjust the Frame size  #############
    def getScalableLayoutInfo(self):
        # enumerate the appropriate square size to fill win rectangle
        # NOTE : the zoom and subsample provided by tk.PhotoImage only can do integral scale
        self.gridRowNum = math.floor(self.FrameHeight / self.classSize)
        self.gridColNum = math.floor(self.FrameWidth / self.classSize)
        gridNum = self.gridRowNum * self.gridColNum
        self.scale = 1
        self.upORdown = "none"
        # decide to downsize or upsize of the class
        if gridNum < self.classesNumber:  # cannot show all classes,downsize the class
            self.upORdown = "down"
            # searching for appropriate thumb size
            while gridNum < self.classesNumber:
                self.scale = self.scale + 1
                # NOTE : tk photo classes subsample() only scale integer size
                (self.thumbClassSize, self.gridRowNum, self.gridColNum) = self.getLayoutInfo(self.upORdown, self.scale)
                gridNum = self.gridRowNum * self.gridColNum
            # show status viewer message
            statusMessage = "LayoutInfo : downsize class,scale = 1/" + str(self.scale)
            statusMessage += ",gridRowNum = " + str(self.gridRowNum)
            statusMessage += ",gridColNum = " + str(self.gridColNum)
            statusMessage += ",Thumb classSize = " + str(self.thumbClassSize)
            self.statusViewer.showMessage(statusMessage)
        elif gridNum > self.classesNumber:  # show more classes,upsize the class to fit
            self.upORdown = "up"
            while gridNum > self.classesNumber:
                self.scale = self.scale + 1
                # NOTE : tk photo classes zoom() only scale integer size
                (self.thumbClassSize, self.gridRowNum, self.gridColNum) = self.getLayoutInfo(self.upORdown, self.scale)
                gridNum = self.gridRowNum * self.gridColNum
            if gridNum == self.classesNumber:
                self.scale = self.scale - 0
            else:  # gridNum < self.classesNumber
                self.scale = self.scale - 1
                (self.thumbClassSize, self.gridRowNum, self.gridColNum) = self.getLayoutInfo(self.upORdown, self.scale)
            # show status viewer message
            statusMessage = "LayoutInfo : downsize class,scale = " + str(self.scale)
            statusMessage += ",gridRowNum = " + str(self.gridRowNum)
            statusMessage += ",gridColNum = " + str(self.gridColNum)
            statusMessage += ",Thumb classSize = " + str(self.thumbClassSize)
            self.statusViewer.showMessage(statusMessage)
            # Aslo need to update grid padding,beacase sticky= not work in grid() method
            # self.gridRowPad = int((self.FrameHeight-self.borderwidth*2 - self.gridRowNum*self.thumbClassSize)/self.gridRowNum/2)
            # self.gridColPad = int((self.FrameWidth-self.borderwidth*2 - self.gridColNum*self.thumbClassSize)/self.gridColNum/2)

    #############  NOTE : in scrolled viewing mode,we fixed the class thumbs size ,it depends on user  #############
    def getScrolledLayoutInfo(self):
        self.upORdown = "down"
        # NOTE : the photoImage subsample() using ceil size,add the label borderwidth
        (self.thumbClassSize, self.gridRowNum, self.gridColNum) = self.getLayoutInfo(self.upORdown, self.scale)

    #############  get the layout info : grid row and column number,and thumb class size  #############
    def getLayoutInfo(self, mode, scale):
        # NOTE : the photoImage zoom() and subsample() using ceil size,add the label borderwidth
        thumbClassSize  = -1.0
        if mode == "down":
            thumbClassSize = math.ceil(float(self.classSize) / float(scale)) + self.borderwidth * 2
        elif mode == "up":
            thumbClassSize = math.ceil(float(self.classSize) * float(scale)) + self.borderwidth * 2
        else:
            print("(ClassesViewer.getLayoutInfo)wrong mode for get layout info.")
        # NOTE : compute the grid dimension,minus Frame borderwidth
        gridRowNum = math.floor(float(self.FrameHeight - self.borderwidth * 2) / float(thumbClassSize))
        gridColNum = math.floor(float(self.FrameWidth - self.borderwidth * 2) / float(thumbClassSize))
        return (thumbClassSize, gridRowNum, gridColNum)

#*** Images Viewer to show images for each class***#
class ImagesViewer:
    def __init__(self, root, dataHandler):
        self.root = root
        self.dataHandler = dataHandler
        self.initLayoutData()

    def tellDataHandlerReset(self):
        self.layoutPageIndex = 0
        self.dataHandler.resetImagePagedata(self.layoutPageIndex, self.layoutGridNum)
        # set up and running the generate tk photo images for images
        self.dataHandler.resetImageImages()
        (imagePageIndex, self.layoutPageNumber, imagePageClassesNumber) = self.dataHandler.getImagePageInfo()

    #############           create and show th top level windows           #############
    def setupIamgeViewer(self):
        #-----  create top win   -----------#
        self.imageViewer = tk.Toplevel(bg='white')
        self.imageViewer.title('ROME2D Images Viewer')
        self.imagesThumbLabel = []
        self.FrameHeight = self.imageViewer.winfo_screenheight()/2
        self.FrameWidth = self.imageViewer.winfo_screenwidth()/2
        self.borderwidth = 2
        #------ set up viewer   -----------#
        self.setupLayout()

    def initLayoutData(self):
        self.layoutGridRowNum = 0
        self.layoutGridColNum = 0
        self.layoutGridNum    = 0
        self.layoutThumbScale = 1 # downsize the image thumb
        self.layoutPageIndex  = 0
        self.layoutPageNumber = 0

    #############  get the state of image viewer(windows show on or not)    #############
    def getImageViewerState(self):
        try:
            state = self.imageViewer.state()
        except:
            state = False
        return state

    #############        set up ,clear ,update, adjust the layout           #############
    def setupLayout(self):
        # add the scalable Frame
        self.imagesFrame = tk.Frame(self.imageViewer, bd=self.borderwidth, bg="white",
                                    height=self.FrameHeight, width=self.FrameWidth)
        # fill and expand the Frame in whole root windows
        self.imagesFrame.pack(fill=tk.BOTH, expand=True)
        # add bottom status viewer
        self.statusViewer = StatusViewer(self.imageViewer,"Page : 0 / 0 | ")
        # add menu
        self.popMenu = tk.Menu(self.root,tearoff=0)
        self.popMenu.add_command(label="Next Page",command = self.eventPopMenuNextPage)
        self.popMenu.add_command(label="Previous Page",command = self.eventPopMenuPreviousPage)
        self.popMenu.add_command(label="First Page",command = self.eventPopMenuFirstPage)
        self.popMenu.add_command(label="Last Page",command = self.eventPopMenuLastPage)
        self.popMenu.add_separator()
        self.popMenu.add_command(label="Refresh Page", command=self.eventRefreshPage)
        self.popMenu.add_separator()
        #-----    add image scale  ---------#
        imageScaleMenu = tk.Menu(self.popMenu,tearoff=0)
        self.imageScaleMode = tk.IntVar()
        self.imageScaleMode.set(1)
        # NOTE : value from 1~5 take as the image scale
        self.popMenu.add_cascade(label = "Image Scale",menu = imageScaleMenu)
        imageScaleMenu.add_radiobutton(label = "1 / 1 Scale",value = 1,variable = self.imageScaleMode,
                                            command = self.eventPopMenuImageScale)
        imageScaleMenu.add_radiobutton(label = "1 / 2 Scale",value = 2,variable = self.imageScaleMode,
                                            command = self.eventPopMenuImageScale)
        imageScaleMenu.add_radiobutton(label = "1 / 3 Scale",value = 3,variable = self.imageScaleMode,
                                            command = self.eventPopMenuImageScale)
        self.popMenu.add_separator()
        self.popMenu.add_command(label="Select All", command=self.eventPopMenuSelectAll)
        self.popMenu.add_command(label="Unselect All", command=self.eventPopMenuUnselectAll)
        self.popMenu.add_separator()
        self.popMenu.add_command(label="Sort by Probality", command=self.eventPopMenuSortProb)

        # add event
        self.imagesFrame.bind("<Configure>", self.eventAdjustLayout)
        # NOTE : the right clicked Menu is 2 on mac and 3 on unix and windows
        self.imagesFrame.bind("<Button-2>", self.eventPopMenu)
        self.imagesFrame.bind("<Button-3>", self.eventPopMenu)

    def clearLayout(self):
        self.imagesFrame.pack_forget()
        self.imagesFrame.destroy()

    def updateLayout(self):
        print ("TODO")

    def resetLayoutInfo(self):
        (temp1,imageSize,temp2) = self.dataHandler.getClassImages()
        print("(imageViewer.resetLayoutInfo),image size get from class size = "+str(imageSize))
        thumbImageSize = math.ceil(float(imageSize) / float(self.layoutThumbScale)) + self.borderwidth * 2
        self.layoutGridRowNum = int(math.floor(float(self.FrameHeight - self.borderwidth * 2) / float(thumbImageSize)))
        self.layoutGridColNum = int(math.floor(float(self.FrameWidth  - self.borderwidth * 2) / float(thumbImageSize)))
        self.layoutGridNum = self.layoutGridRowNum*self.layoutGridColNum
        print("(imageViewer.resetLayoutInfo),Row:"+str(self.layoutGridRowNum)+",Col:"+str(self.layoutGridColNum)+",Num:"+str(self.layoutGridNum))

    def eventAdjustLayout(self,event):
        # NOTE : each 21 step update,donot update frequence,high frequence will cause problem
        #if event.width % 3 == 0 or event.height % 3 == 0:
        self.FrameHeight = event.height
        self.FrameWidth = event.width
        # if there are label in frame,update the layout

    def eventRefreshPage(self):
        self.resetLabels()

    #############           set up , reset, clear ,update the labels             #############
    def setupLabels(self):
        #--------   convert the image to to TK Photo image  ----------#
        (self.images, self.imageSize, self.imagesNumber) = self.dataHandler.getImageImages()
        newSize = int(self.imageSize / self.layoutThumbScale)
        self.imagesThumb = [ImageTk.PhotoImage(image.resize((newSize,newSize),resample=Image.NEAREST)) for image in self.images]
        #--------    clear labels      -------#
        self.clearLabels()
        # ---------   set up all label    ------------#
        self.imagesThumbLabel = []
        imagePageImage = self.dataHandler.getImagePageImages()
        for iimage in range(self.imagesNumber):
            selectOrNot = self.dataHandler.getStateImageSelectOrNot(imagePageImage[iimage])
            self.imagesThumbLabel.append(tk.Label(self.imagesFrame, image=self.imagesThumb[iimage],
                                                text=imagePageImage[iimage], bd=self.borderwidth, bg=("white","red")[selectOrNot]))
        #--------            add event         ------------#
        for label in self.imagesThumbLabel:
            label.bind("<Button-1>", self.eventSelectImageThumb)
            # NOTE : the right clicked Menu is 2 on mac and 3 on unix and windows
            label.bind("<Button-2>", self.eventPopMenu)
            label.bind("<Button-3>", self.eventPopMenu)
        #--------    put all label on layout  -------------#
        for iimage in range(self.imagesNumber):
            Icol = iimage % int(self.layoutGridColNum)
            Irow = int((iimage - Icol) / int(self.layoutGridColNum))
            # NOTE : sticky=tk.W+tk.E not work here,padx = self.gridRowPad,pady = self.gridColPad
            self.imagesThumbLabel[iimage].grid(row=Irow, column=Icol)
        # show status
        self.eventStatusViewerShow("")

    def resetLabels(self):
        #------- get layout info   ------------#
        self.resetLayoutInfo()
        if self.layoutGridNum == 0: # return when your layout is bad
            return
        #------- prepare data first     -----------#
        self.tellDataHandlerReset()
        self.setupLabels()

    def clearLabels(self):
        #--------    remove all label  on layout  -------------#
        ## TODO : is this efficient to remove all the label?
        for label in self.imagesThumbLabel:
            label.grid_forget()  # just make invisiable
            label.destroy()
        self.imagesThumbLabel = []

    def updateLabels(self):
        (imagePageIndex, imagePageNumber, imagePageClassesNumber) = self.dataHandler.getImagePageInfo()
        # if current page is last page(self.layoutPageIndex == imagePageNumber-1) or
        # go to the last page (imagePageIndex == imagePageNumber-1)
        if self.layoutPageIndex == imagePageNumber-1 or imagePageIndex == imagePageNumber-1:
            #--------  clear label first     -------#
            self.clearLabels()
            #--------  setup all labels      -------#
            self.setupLabels()
        else : #just configure the label,update image,text,and color
            # --------   convert the image to to TK Photo image  ----------#
            (self.images, self.imageSize, self.imagesNumber) = self.dataHandler.getImageImages()
            newSize = int(self.imageSize / self.layoutThumbScale)
            imagesThumb = [ImageTk.PhotoImage(image.resize((newSize, newSize), resample=Image.NEAREST)) for image in self.images]
            imagePageImage = self.dataHandler.getImagePageImages()
            for iimage in range(self.imagesNumber):
                selectOrNot = self.dataHandler.getStateImageSelectOrNot(imagePageImage[iimage])
                # update the background
                self.imagesThumbLabel[iimage].configure(bg = ("white","red")[selectOrNot])
                self.imagesThumbLabel[iimage].configure(text = imagePageImage[iimage])
                self.imagesThumbLabel[iimage].configure(image = imagesThumb[iimage])
            # NOTE : keep the image reference
            self.imagesThumb = imagesThumb
        self.layoutPageIndex = imagePageIndex
        self.layoutPageNumber = imagePageNumber
        # show status
        self.eventStatusViewerShow("")


    #############  update label background color           #############
    def updateLabelBG(self, label):
        if label.cget('bg') == "white":
            label.configure(bg="red")
        else:
            label.configure(bg="white")

    ####################################################################
    ##############           define all events         #################
    ####################################################################
    #############    select the image thumb                #############
    def eventSelectImageThumb(self,event):
        # get the image thumb state first
        iimage = event.widget.cget("text")
        selectOrNot = self.dataHandler.getStateImageSelectOrNot(iimage)
        self.updateLabelBG(event.widget)
        self.dataHandler.setStateImageSelectOrNot(iimage,not selectOrNot)
        (imageMetaIndex,imageMetaName,imageMetaProb) = self.dataHandler.getMetadata(iimage)
        message = "Image : "+str(imageMetaIndex)+"@"+str(imageMetaName)
        message += " | maxProbality : "+str(imageMetaProb)
        self.eventStatusViewerShow(message)

    def eventPopMenu(self,event):
        self.popMenu.post(event.x_root,event.y_root)

    ###########         set the page          ###############
    def eventPopMenuSetPage(self,newImagePageIndex,imagePageNumber):
        if newImagePageIndex < imagePageNumber and newImagePageIndex >= 0:
            self.dataHandler.resetImagePagedata(newImagePageIndex, self.layoutGridNum)
            # set up and running the generate tk photo images for images
            self.dataHandler.resetImageImages()
            self.updateLabels()
        else:
            print("(imageViewer.eventPopMenuSetPage,out of range)"+str(newImagePageIndex)+","+str(imagePageNumber))

    def eventPopMenuNextPage(self):
        (imagePageIndex, imagePageNumber, imagePageClassesNumber) = self.dataHandler.getImagePageInfo()
        imagePageIndex += 1
        self.eventPopMenuSetPage(imagePageIndex,imagePageNumber)

    def eventPopMenuPreviousPage(self):
        (imagePageIndex, imagePageNumber, imagePageClassesNumber) = self.dataHandler.getImagePageInfo()
        imagePageIndex -= 1
        self.eventPopMenuSetPage(imagePageIndex, imagePageNumber)

    def eventPopMenuFirstPage(self):
        (imagePageIndex, imagePageNumber, imagePageClassesNumber) = self.dataHandler.getImagePageInfo()
        imagePageIndex = 0
        self.eventPopMenuSetPage(imagePageIndex, imagePageNumber)

    def eventPopMenuLastPage(self):
        (imagePageIndex, imagePageNumber, imagePageClassesNumber) = self.dataHandler.getImagePageInfo()
        imagePageIndex = imagePageNumber - 1
        self.eventPopMenuSetPage(imagePageIndex, imagePageNumber)

    def eventPopMenuImageScale(self):
        imageScaleMode = self.imageScaleMode.get()
        self.layoutThumbScale = imageScaleMode
        self.resetLabels()

    def eventPopMenuSelectAll(self):
        imagePageImage = self.dataHandler.getImagePageImages()
        for iimage in range(len(imagePageImage)):
            self.dataHandler.setStateImageSelectOrNot(imagePageImage[iimage], True)
            selectOrNot = self.dataHandler.getStateImageSelectOrNot(imagePageImage[iimage])
            self.imagesThumbLabel[iimage].configure(bg=("white","red")[selectOrNot])

    def eventPopMenuUnselectAll(self):
        imagePageImage = self.dataHandler.getImagePageImages()
        for iimage in range(len(imagePageImage)):
            self.dataHandler.setStateImageSelectOrNot(imagePageImage[iimage], False)
            selectOrNot = self.dataHandler.getStateImageSelectOrNot(imagePageImage[iimage])
            self.imagesThumbLabel[iimage].configure(bg=("white","red")[selectOrNot])

    def eventPopMenuSortProb(self):
        self.dataHandler.resetFilteredImagedata("ByProb")
        self.resetLabels()

    def eventStatusViewerShow(self,message):
        self.statusViewer.showMessage(
            "Page : " + str(self.layoutPageIndex + 1) + " / " + str(self.layoutPageNumber) + " | "+message)

# handle the Mrcs,STAR,TIFF datafile
dataHandler = DataHandler()

# create the root window
root = tk.Tk()
root.title("ROME2D Viewer")
# add bottom status bar
statusViewer = StatusViewer(root,"Preparing,you should open file first.")
# add assist viewer
imagesViewer = ImagesViewer(root,dataHandler)
# add images viewer
classesViewer = ClassesViewer(root,imagesViewer,statusViewer,dataHandler)
# add top Menu
topMenu = TopMenu(root,classesViewer,imagesViewer,statusViewer,dataHandler)

root.mainloop()

