from decimal import Decimal

import sys
#----------------------------------------------------------------------
# supported categorizations
#----------------------------------------------------------------------

categorizationToNumCategories = {
    "1r91eta":  1,
    "1r92eta":  1 * 2,
    "2r92eta":  2 * 2,
    "2r92eta2pTh" : 2 * 2 * 2,

    # Chris' analysis
    "VBF": 1,
    }

# first index is the categorization, the second
# is the category number
categoryToName = {
    "2r92eta2pTh": [
        "high pT(gg), barrel, high R9",
        "high pT(gg), barrel, low R9",
        "high pT(gg), endcap, high R9",
        "high pT(gg), endcap, low R9",
        "low pT(gg), barrel, high R9",
        "low pT(gg), barrel, low R9",
        "low pT(gg), endcap, high R9",
        "low pT(gg), endcap, low R9",
        ]
    }

#----------------------------------------------------------------------

def getGraphYvalues(graph):
    
    size = graph.GetN()
    yvalues = graph.GetY()
    
    return [ yvalues[i] for i in range(size)]

#----------------------------------------------------------------------

def makeGraphFromArrays(xvalues, yvalues):
    import array, ROOT
    
    assert(len(xvalues) == len(yvalues))

    return ROOT.TGraph(len(xvalues),
                     array.array('d', xvalues),
                     array.array('d', yvalues)
                     )

#----------------------------------------------------------------------

def getGraphMaximumY(graph):
    return max(getGraphYvalues(graph))

def getGraphsMaximumY(graphs):
    return max([ getGraphMaximumY(graph) for graph in graphs])

#----------------------------------------------------------------------

def optError(msg):
    print >> sys.stderr,msg + " Run with --help to get the list of options."
    sys.exit(1)

#----------------------------------------------------------------------

def readCSV(csvFile, includeExpected):
    """ @param csvFile must be an open file object""" 

    import csv
    
    reader = csv.reader(csvFile)

    # read the header
    header = reader.next()

    #--------------------
    # 2011-09-21 added some cross checks on format
    # mass,xsect_lim_obs,xsect_lim_expected

    massPos                = header.index('mass')
    xsect_lim_obs_pos      = header.index('xsect_lim_obs')
    xsect_lim_expected_pos = header.index('xsect_lim_expected')

    # the following are optional
    xsect_lim_expected_minus_2_sigma_pos = listFind(header,'xsect_lim_expected_minus_2_sigma')
    xsect_lim_expected_minus_1_sigma_pos = listFind(header,'xsect_lim_expected_minus_1_sigma')

    xsect_lim_expected_plus_1_sigma_pos = listFind(header,'xsect_lim_expected_plus_1_sigma')
    xsect_lim_expected_plus_2_sigma_pos = listFind(header,'xsect_lim_expected_plus_2_sigma')

    #--------------------
    
    masses = []
    observedValues = []
    if includeExpected:
        expectedValues = []
        xsect_lim_expected_minus_2_sigma_values = []
        xsect_lim_expected_minus_1_sigma_values = []
        xsect_lim_expected_plus_1_sigma_values  = []
        xsect_lim_expected_plus_2_sigma_values  = []

    else:
        expectedValues = None
        xsect_lim_expected_minus_2_sigma_values = None
        xsect_lim_expected_minus_1_sigma_values = None
        xsect_lim_expected_plus_1_sigma_values =  None
        xsect_lim_expected_plus_2_sigma_values =  None
        
    for row in reader:
        
        if not row:
            # skip empty lines
            continue
        masses.append(float(row[massPos]))

        # observed cross section limit
        observedValues.append(float(row[xsect_lim_obs_pos]))

        if includeExpected:
            expectedValues.append(float(row[xsect_lim_expected_pos]))

            for pos, values in (
                (xsect_lim_expected_minus_2_sigma_pos, xsect_lim_expected_minus_2_sigma_values),
                (xsect_lim_expected_minus_1_sigma_pos, xsect_lim_expected_minus_1_sigma_values),
                (xsect_lim_expected_plus_1_sigma_pos,  xsect_lim_expected_plus_1_sigma_values ),
                (xsect_lim_expected_plus_2_sigma_pos,  xsect_lim_expected_plus_2_sigma_values ),
                ):
                if pos != -1:
                    values.append(float(row[pos]))


    # sort them by increasing xvalue order
    # (necessary for plotting and for the filled
    indices = range(len(masses))
    indices.sort(lambda i,j: cmp(masses[i],masses[j]))

    masses = [ masses[i] for i in indices ]
    observedValues = [ observedValues[i] for i in indices ]
    
    if expectedValues:
        expectedValues = [ expectedValues[i] for i in indices ]

    return masses, observedValues, expectedValues, \
           xsect_lim_expected_minus_2_sigma_values, \
           xsect_lim_expected_minus_1_sigma_values, \
           xsect_lim_expected_plus_1_sigma_values, \
           xsect_lim_expected_plus_2_sigma_values 

#----------------------------------------------------------------------

# some functions used by the cross section limit plotting scripts

#----------------------------------------------------------------------

def getXsectTimesBR(mass, fermiophobic):

    # Marco's cross sections
    # import CrossSectionData 
    # observedValues[index] *= CrossSectionData.getInterpolatedXsect(mass)

    import LHCxsect7TeV

    if fermiophobic:
        import LHCFermioPhobicBR

        # this does not work yet with mass steps smaller than one GeV.
        xsectTimesBR = LHCxsect7TeV.getFermiophobicXsect(mass) * LHCFermioPhobicBR.getBR(mass,"gammaGamma")

    else:
        import LHCsmHiggsBR
        xsectTimesBR = LHCxsect7TeV.getXsect(mass) * LHCsmHiggsBR.getBR(mass,"gammaGamma")
    

    return xsectTimesBR

#----------------------------------------------------------------------

def listFind(list, item):

    if item in list:
        return list.index(item)
    else:
        return -1

#----------------------------------------------------------------------

def writeCSV(csvFile, masses, observedValues, expectedValues):
    """ @param csvFile must be an open file object""" 

    assert(len(masses) == len(observedValues))

    headers = ['mass', 'xsect_lim_obs']

    if expectedValues != None:
        assert(len(masses) == len(expectedValues))
        headers.append('xsect_lim_expected')
    else:
        expectedValues = [ None ] * len(masses)

    # print the header
    print >> csvFile,",".join(headers)

    for mass, obs, exp in zip(masses, observedValues, expectedValues):

        values = [ mass, obs]
        if exp != None:
            values.append(exp)

        print >> csvFile,",".join([str(x) for x in values ])


#----------------------------------------------------------------------

def __applyOperatorOnArrayWithXsectAndBR(masses, origValues, fermiophobic, opFunc):
    """ helper function for multiplyArrayByXsectAndBR and divideArrayByXsectAndBR
    """

    if origValues == None:
        return None
        
    assert(len(masses) == len(origValues))
        
    retval = []
         
    for mass, value in zip(masses, origValues):
        # assume the inputs were normalized to the LHC working group
        # cross sections

        # apply the operator function (multiply/divide) by the absolute cross section
        retval.append(opFunc(value,getXsectTimesBR(mass, fermiophobic)))

    return retval

#----------------------------------------------------------------------

def multiplyArrayByXsectAndBR(masses,xsectRatios, fermiophobic):
    """ this is used to multiply a set of cross section ratios
        by cross section times branching ration into gamma gamma 
        
        For convenience, xsectRatios may be None in which case
        this function returns None (this is useful e.g. when 
        the expected values should not be included and thus
        are set to None instead of a list in the calling function)
        
        @return the list with the normalized values or None
          if xsectRatios was None.
        """

    import operator
    return __applyOperatorOnArrayWithXsectAndBR(masses, xsectRatios, fermiophobic, operator.mul)

#----------------------------------------------------------------------

def divideArrayByXsectAndBR(masses, values, fermiophobic):
    """ similar to multiplyArrayByXsectAndBR but divides
        'values' element-wise by cross section * branching ratio.

        Useful e.g. to go from absolute limits on the cross section
        to relative limits (w.r.t. the predicted cross section)
        """

    import operator
    return __applyOperatorOnArrayWithXsectAndBR(masses, values, fermiophobic, operator.div)

#----------------------------------------------------------------------





def defaultParameters():
    return dict(
        name = "norm",
        inputoutputfile = "may26_schif4_CIC4fix_ST_ST_json13_smear8_shift8_EBEE_300pb_rered",
        # categorization = "1r91eta",
        # bins = "250MeV",
        # fit = "mpexp",
        whichcat = 0,         # type of categorization
        data = 1,
        soverb = 0.1,         # S/B cut
        scalesyst = 1.12,
        correlatederror = 0.02,
        smearsig = 0,
        nexperiments = 100000,
        # shiftmass = -2,
        shiftmass = 0,
        lumi = 1,
        scalesig = 1
        )

#----------------------------------------------------------------------

def getOutputFileNames(name,
          inputoutputfile,
          # categorization,
          # bins,
          # fit,
          whichcat,         # type of categorization
          data,
          soverb,         # S/B cut
          scalesyst,
          correlatederror,
          smearsig,
          nexperiments,
          shiftmass,        # not used anymore
          lumi,
          scalesig,
          signalIntMass,
                       ):
    """ returns the names of the output files we think clfast_marco produces
    for the given parameters """

    # currently defined in ClRun_marco.C (left are parameter names in ClRun_marco.C, right
    # are those used in runClRun.C)
    smearbg     = scalesyst
    thisntries  = nexperiments
    shiftsignal = shiftmass
    scalsig     = scalesig

    # old convention (with signal mass shift but not absolute signal mass)
    # fname = "%s_%s_%s_smear_%s_%s.root_%d_%f_%f_%d_%d_%f" % (inputoutputfile, categorization, bins, fit, name,
    # 
    #                                                          # after .root
    #                                                          whichcat,soverb,smearbg,thisntries,shiftsignal,scalsig)

    # until 2011-07-04
    ## fname = "%s_%s_%s_smear_%s_%s.root_%d_%f_%f_%d_%.1f_%.5f" % (inputoutputfile, categorization, bins, fit, name,
    ##                                                          
    ##                                                          # after .root
    ##                                                          whichcat,soverb,smearbg,thisntries,signalIntMass / 10.0,scalsig)

    fname = "%s_%d_%f_%f_%d_%.1f_%.5f" % (inputoutputfile, 
                                          whichcat,soverb,smearbg,thisntries,signalIntMass / 10.0,scalsig)


    fname += "_output.dat"

    return (fname, fname + ".gif")

#----------------------------------------------------------------------

def getPrepareHistosErrorOutputFileName(params):
    """ name of the file produced by GlobePrepareHistos containing the uncertainties of the background fit
        per category """

    fname = "%s_%s_%s_smear_%s_%s.root_errors.dat" % (params['inputoutputfile'], params['categorization'], params['bins'], params['fit'], params['name'])
    return fname

#----------------------------------------------------------------------

def parseOutputFileName(fname):
    """ returns the list of parameters given an output file name """

    origFname = fname
    
    if fname.endswith(".gif"):
        fname = fname[:-4]

    if fname.endswith(".dat"):
        fname = fname[:-4]

    parts = fname.split('_')
    if parts[-1] != "output":
        raise Exception("file name " + origFname + " does not have _output_ at the expected position")
    parts.pop(-1)

    # example file name:
    #   may26_schif4_CIC4fix_ST_ST_json13_smear8_shift8_EBEE_300pb_rered_1r91eta_250MeV_smear_mpexp_norm.root_0_0.000000_1.120000_100000_115_1.000000_output.dat.gif

    # example file name (2011-08-10):
    #   jul30_all_2r92eta2pTh_250MeV_smear_mpexp_120.root_80_0.011000_0.000000_200000_120.0_1.10000_output.dat
    
    retval = {}

    retval['scalesig']        = Decimal(parts.pop(-1))


    # changed convention 2011-06-09
    # retval['shiftmass']       = int(parts.pop(-1))
    retval['signalIntMass']   = int(float(parts.pop(-1)) * 10 + 0.5)
    
    retval['nexperiments']    = int(parts.pop(-1))
    retval['scalesyst']       = Decimal(parts.pop(-1))
    retval['soverb']          = Decimal(parts.pop(-1))
    retval['whichcat']        = int(parts.pop(-1))

    # before .root
    retval['inputoutputfile'] = "_".join(parts)

    value = parts.pop(-1)
    assert(value.lower().endswith('.root'))
    retval['name']            = value[:-5]

    retval['fit']             = parts.pop(-1)

    assert(parts.pop(-1) == 'smear')

    retval['bins']            = parts.pop(-1)
    retval['categorization']  = parts.pop(-1)

    # assume the rest is the name of the original file
    # retval['inputoutputfile'] = "_".join(parts)

    return retval
    
#----------------------------------------------------------------------
def testParseOutputFileName():
    """ a test for the parseOutputFileName(..) function """

    # needs to be updated for the new format (after 2011-06-08)

    fname = "may26_schif4_CIC4fix_ST_ST_json13_smear8_shift8_EBEE_300pb_rered_1r91eta_250MeV_smear_mpexp_norm.root_0_0.000000_1.120000_100000_115_1.000000_output.dat.gif"

    params = parseOutputFileName(fname)

    # pprint(params)

    #--------------------
    # these parameters are not reflected
    # in the output file name
    defparam = defaultParameters()
    for name in ['data', 'correlatederror', 'smearsig', 'lumi' ]:
        params[name] = defparam[name]

    #--------------------
    
    # print len(params)
    outputFnames = getOutputFileNames(**params)

    assert(fname in outputFnames)

    
#----------------------------------------------------------------------
def makeBinarySearchOrder(values):
    
    if len(values) < 3:
        return values
    
    values = values[:]
    
    retval = []
    
    
    retval.append(values.pop(0))
    retval.append(values.pop(-1))
    # retval.append(values.pop(len(values) / 2))
    
    queue = []
    queue.append(values)
    
    # not the most efficient implementation but should work
    while queue:
        thisValues = queue.pop(0)
        
        if len(thisValues) < 2:
            retval.extend(thisValues)
            continue

        mid = len(thisValues) / 2
        
        left = thisValues[:mid]
        retval.append(thisValues[mid])
        right = thisValues[mid+1:]
        
        queue.append(left)
        queue.append(right)

        
    return retval


#----------------------------------------------------------------------
def binsizeFromString(binningString):

    # convert it to GeV
    binningUnit = binningString.lower()[-3:]
    binningValue = float(binningString[:-3])
    if binningUnit == "mev":
        return binningValue * 1e-3
    elif binningUnit == "gev":
        return binningValue
    else:
        raise Exception("unknown binning unit " + binningUnit)

    
def readResultFile(fname):

    fin = open(fname)

    retval = {}

    for line in fin.readlines():
        line = line.split('\n')[0].strip()

        if line == '':
            continue

        parts = line.split('=',1)

        parts = [ x.strip() for x in parts ]

        # assume these are all floats
        retval[parts[0]] = float(parts[1])


    fin.close()
    return retval


#----------------------------------------------------------------------

def runIt(name,
          inputoutputfile,
          # categorization,
          # bins,
          # fit,
          whichcat,         # type of categorization
          data,
          soverb,         # S/B cut
          scalesyst,
          correlatederror,
          smearsig,
          nexperiments,
          shiftmass,
          lumi,
          scalesig,
          signalIntMass
          ):

    # TString name, TString inputoutputfile, TString categorization, TString bins, TString fit, int whichcat, int data, float soverb, float scalesyst, float correlatederror, float smearsig, int nexperiments, int shiftmass, float lumi, float scalesig) 


    arguments = [
        ('name', str),
        ('inputoutputfile', str),
        # ('categorization', str),
        # ('bins', str),
        # ('fit', str),
        ('whichcat', int),
        ('data', int),
        ('soverb', float),
        ('scalesyst', float),
        ('correlatederror', float),
        ('smearsig', float),
        ('nexperiments', int),
        ('shiftmass', int),
        ('lumi', float),
        ('scalesig', float),
        ('signalIntMass', int),
        ]


    # put string together for the command
    argParts = []

    for argName, argType in arguments:
        value = locals()[argName]
        if not isinstance(value,argType):
            if argType == float and isinstance(value,int):
                pass
            else:
                raise Exception("parameter " + argName + " is not of type " + str(argType))

        if argType == str:
            argParts.append('"%s"' % value)
        else:
            argParts.append(str(value))

    cmd = "root -l -q -x 'runClRun.C(%s)'" % ",".join(argParts)

    # print "cmd=",cmd

    import tempfile
    fout = tempfile.NamedTemporaryFile(delete = False)

    cmd += " > " + fout.name + " 2>&1"
    res = os.system(cmd)

    if res != 0:
        import sys
        print >> sys.stderr,"error running cl calculation, log file is " + fout.name
        assert(res == 0)

    return fout.name

#----------------------------------------------------------------------
import threading

class RunnerThread(threading.Thread):
    """ a wrapper around a function which should be run in a
    thread making sure that not more than a given number of threads run
    at the same time."""

    #--------------------

    def __init__(self, threadLimiter, func):
        threading.Thread.__init__(self)
        self.threadLimiter = threadLimiter

        self.func = func

    #--------------------
    def run(self):
        self.threadLimiter.acquire()
        try:
            self.func()
        finally:
            self.threadLimiter.release()
            
#----------------------------------------------------------------------
class ThreadLimiterScheduler:

    def __init__(self, maxNumThreads):
        self.threadLimiter = threading.BoundedSemaphore(maxNumThreads)
        self.threads = []

    def run(self, func):
        """ queues func in a thread and calls
        func() when there is a slot available to run """

        self.threads.append(RunnerThread(self.threadLimiter, func))
        self.threads[-1].start()

    def runAll(self, funcs):
        for func in funcs:
            self.run(func)

    def waitForAll(self):
        # waits for all threads to complete
        for thread in self.threads:
            thread.join()

#----------------------------------------------------------------------
import os, shutil

class ClRunner:

    #----------------------------------------
    def __init__(self, params, outputDir):
        # make a copy
        self.params = dict(params)

        self.outputDir = outputDir

    #----------------------------------------
    def __call__(self):
        # run the external program
        logFname = runIt(**self.params)

        # get the output files
        self.outputfiles = getOutputFileNames(**self.params)

        datFname = None

        for fname in self.outputfiles:
            if not os.path.exists(fname):
                raise Exception("expected output file " + fname + " does not exist")

            if fname.endswith(".dat"):
                datFname = fname

        # for fname in self.outputfiles:
        #     shutil.move(fname, self.outputDir + "/" + fname)

        # move log file
        basedir = os.path.dirname(datFname)
        shutil.move(logFname, basedir + "/" + os.path.basename(datFname) + ".log")


#----------------------------------------------------------------------

def readResultFile(fname):

    fin = open(fname)

    retval = {}

    for line in fin.readlines():
        line = line.split('\n')[0].strip()

        if line == '':
            continue

        # do not read the parameter section
        # (it's formatted differently)
        import re
        if re.match("-+$",line):
            break

        parts = line.split('=',1)

        parts = [ x.strip() for x in parts ]

        # assume these are all floats
        retval[parts[0]] = float(parts[1])


    fin.close()
    return retval

#----------------------------------------------------------------------


def findParamSetsFromFileNames(fnames):
    paramSets = []

    filesToParamSets = {}

    varyingParamNames = set([])

    for fname in fnames:

        # do not look at incomplete files
        if os.path.getsize(fname) == 0:
            continue

        thisParams = parseOutputFileName(fname)

        # see how this parameter set differs from the previous parameter sets
        # (we assume that the set of parameter NAMES always stays the same)

        for prevParams in paramSets:

            for name, value in prevParams.items():

                if thisParams[name] != value:
                    varyingParamNames.add(name)


        paramSets.append(thisParams)
        filesToParamSets[fname] = thisParams

    return filesToParamSets, varyingParamNames

#----------------------------------------------------------------------

def getNumberOfCategories(params):
    
    retval = params['whichcat']
    
    if retval == 0:
        # this is data, one category
        return 1
    else:
        if retval >= 20:
            # data
            return retval / 10
        else:
            # MC
            return retval

#----------------------------------------------------------------------
#def scaleGraphY(graph, scaleFactor):
#    # multiplies all y values of the graph by the given scaleFactor 
#    
#  
#  for (int i = 0; i < graph->GetN(); ++i)
#  {
#    double x,y;
#    graph->GetPoint(i, x, y);
#    graph->SetPoint(i, x, y * scaleFactor);
#  }
#
#}

#----------------------------------------------------------------------
def getWhichCat(numCategories, bgFitFromData):

    # 1..16: background fit from MC
    if not bgFitFromData:
        return numCategories

    # 10..160: background fit from data but instead of 10 use 0
    if numCategories == 1:
        return 0 # instead of 10 (which is used for MC, 10 categories), special value

    return 10 * numCategories
    
#----------------------------------------------------------------------

def countProcessors():
    """ determines the number of 'processors' (could include hyperthreads)
    according to /proc/cpuinfo"""

    retval = 0

    for line in open("/proc/cpuinfo").readlines():
        line = line.strip()

        if line == "":
            continue

        left, right = line.split(':',1)
        left = left.strip()

        if left == 'processor':
            retval += 1

    return retval
    
#----------------------------------------------------------------------

def addCommonCommandLineOptions(parser):
    """ for runCL.py and runXSectLimits.py: adds common command
    line options """

    parser.add_option(
        "--categorization",
        help="Categorization. Supported values are " + ", ".join(categorizationToNumCategories.keys()),
        default=None,
        )

    parser.add_option(
        "--num-threads",
        type=int,
        help="maximum number of threads (mass points) to run in parallel",
        default= countProcessors(),
        )

    parser.add_option(
        "--bgfit",
        type=str,
        help="'mc' or 'data' depending on whether to take the background from the fit to data or background MC",
        default='data',
        )

    parser.add_option(
        "--masses",
        type=str,
        help="specification of mass range: can be either a single value, or start,stop (assumes step of 1 GeV) or of the form start,stop,step (all masses in GeV). If not given, all masses available in the input file are processed",
        default=None
        )

    parser.add_option(
        "--output",
        type=str,
        dest="outputDir",
        help="manually specify the output directory. Must not yet exist.",
        default=None
        )

    parser.add_option(
        "--smearsig",
        type=float,
        dest="smearsig",
        help="signal systematic uncertainty (normalization). e.g. specify 0.1 for 10%",
        default=0.0,
        )

    parser.add_option(
        "--correlatederror",
        type=float,
        dest="correlatederror",
        help="background uncertainty (normalization) correlated across all channels. E.g. specify 0.01 for 1%",
        default=0.0,
        )

    parser.add_option(
        "--num-trials",
        type=int,
        dest="numTrials",
        help="number of pseudo-experiments to perform",
        default=100000,
        )

    parser.add_option(
        "--scalesyst",
        type=float,
        dest="scalesyst",
        help="by how much to scale the (background ?!) systematic uncertainties",
        default=0,
        )

#----------------------------------------------------------------------
def makeMassRange(startMass, stopMass, massStep = 0.1):
    retval = []

    mass = startMass
    while mass <= stopMass:

        retval.append(mass)
        mass += massStep


    return retval
    
#----------------------------------------------------------------------

def parseMassRangeString(massRangeString):
    """ takes a comma separated list of values (maximum three)
        and returns the corresponding range"""
    
    tmp = massRangeString.split(',')

    massStep = 1

    if len(tmp) == 1:
        startMass = float(tmp[0])
        stopMass = startMass
    elif len(tmp) == 2:
        startMass = float(tmp[0])
        stopMass = float(tmp[1])
    elif len(tmp) == 3:
        startMass = float(tmp[0])
        stopMass = float(tmp[1])
        massStep = float(tmp[2])
    else:
        optError("Illegal mass range specification '%s'." % options.masses)


    return makeMassRange(startMass, stopMass, massStep)

#----------------------------------------------------------------------


def getMassesFromInputFileOrCommandLineArguments(options, inputRootFile):
    """ determines the masses to run on either from the given command line options
    or from the masses available in the input file.

    @return a list of the masses to be processed
    """
    if options.masses == None:
        # not specified on the command line, check which masses
        # are available in the input file

        masses = findMassHypothesesFromGlobePrepareHistosOutput(inputRootFile)
        print >> sys.stderr,"found the following masses in the input file:",masses

        return masses

    else:
        return parseMassRangeString(options.masses)


#----------------------------------------------------------------------
def processCommonCommandLineOptions(options, ARGV, params):
    """ sets params """
    
    if len(ARGV) != 1:
        print >> sys.stderr,"no input file specified. Run with --help to get the list of options."
        sys.exit(1)

    options.input_root_file = ARGV[0]

    if options.categorization != None and not options.categorization in categorizationToNumCategories.keys():
        optError("invalid categorization '%s'." % options.categorization)

    if not os.path.exists(options.input_root_file):
        print >> sys.stderr,"input file " + options.input_root_file + " does not exist"
        sys.exit(1)

    if not options.bgfit.lower() in ('data', 'mc'):
        optError("invalid background fit type '" + options.bgfit + "'.")

    #--------------------
    # parse mass range
    #--------------------

    options.masses = getMassesFromInputFileOrCommandLineArguments(options, options.input_root_file)

    #----------------------------------------
    # set parameters
    #----------------------------------------

    params['name'] = "norm"
    #--------------------
    # params['categorization'] = options.categorization


    if options.categorization != None:
        # specified on the command line, convert it to the number
        # of categories, combined with which background fit (to MC or to data)
        # we should use
        params['whichcat'] = getWhichCat(categorizationToNumCategories[options.categorization], options.bgfit == 'data')
    else:
        # determine the number of categories from the given input file
        categories = findCategoryNumbersFromGlobePrepareHistosOutput(options.input_root_file)

        # ensure they are consecutive
        if not categories:
            print >> sys.stderr,"no categories found in input file %s" % options.input_root_file
            sys.exit(1)

        if categories[-1] != len(categories) - 1:
            print >> sys.stderr,"categories are not consecutive in input file %s" % options.input_root_file
            sys.exit(1)

        params['whichcat'] = getWhichCat(len(categories), options.bgfit == 'data')

    #--------------------

    if options.bgfit.lower() == "data":
        params['data'] = 1                 # for picking the right error from the fit (1 = data, 0 = MC)
    elif options.bgfit.lower() == "mc":
        params['data'] = 0
    else:
        raise Exception("internal error")

    # params['scalesyst'] = 0 # should be 1.12 if not excluding the fitted range 
    # params['scalesyst'] = 1.12 # if not excluding the fitted range 1.12
    params['scalesyst'] = options.scalesyst

    #--------------------
    # original 
    #--------------------
    # params['correlatederror'] = 0.02 # what to set this to ?
    # params['smearsig'] = 0.1 

    #--------------------
    # new
    #--------------------
    params['correlatederror'] = options.correlatederror
    params['smearsig'] = options.smearsig 

    #--------------------

    params['nexperiments'] = options.numTrials

    # already applied in prepare histos ?!
    params['lumi'] = 1

    # 1 for SM cross section
    # Marco's suggestion 2011-06-09
    if hasattr(options,'scalesig'):
        params['scalesig'] = options.scalesig
    else:
        params['scalesig'] = 10

    params['soverb'] = params['scalesig'] * 0.1
    # params['soverb'] = 0


    #----------------------------------------------------------------------
    # DEBUGGING ONLY
    # params['smearsig'] = 0
    # params['scalesyst'] = 0   # background uncertainty
    #----------------------------------------------------------------------



#----------------------------------------------------------------------

def getTreeV1(tree):

    numSelectedRows = tree.GetSelectedRows()

    values = tree.GetV1()

    return [ values[i] for i in range(numSelectedRows) ]

#----------------------------------------------------------------------
def getTreeV2(tree):

    numSelectedRows = tree.GetSelectedRows()

    values = tree.GetV2()

    return [ values[i] for i in range(numSelectedRows) ]

#----------------------------------------------------------------------

def findMassHypothesesFromGlobePrepareHistosOutput(fin):
    """ given an output file of GlobePrepareHistos (an input to the confidence
    level calculation), returns the mass hypotheses found.

    fin can either be a TFile or a string (in which case the TFile is opened
    and closed again)
    """

    import ROOT
    import re

    

    if isinstance(fin,str):
        fin = ROOT.TFile.Open(fin)
        mustCloseFin = True
    else:
        mustCloseFin = False

    retval = []

    for key in fin.GetListOfKeys():

        obj = fin.Get(key.GetName())

        if not isinstance(obj, ROOT.TDirectoryFile):
            continue

        mo = re.match("mass(\d+\.\d)$", key.GetName())
        if not mo:
            continue

        retval.append(float(mo.group(1)))


    if mustCloseFin:
        fin.Close()
        
    return sorted(retval)

#----------------------------------------------------------------------

def findCategoryNumbersFromGlobePrepareHistosOutput(fin):
    """ given an output file of GlobePrepareHistos (an input to the confidence
    level calculation), returns the categories found.

    This is just obtained from the data histograms (it is assumed
    that the categories are the same for the background and signal
    histograms...)

    fin can either be a TFile or a string (in which case the TFile is opened
    and closed again)
    """

    import ROOT
    import re

    if isinstance(fin,str):
        fin = ROOT.TFile.Open(fin)
        mustCloseFin = True
    else:
        mustCloseFin = False

    retval = set()

    for key in fin.GetListOfKeys():

        obj = fin.Get(key.GetName())

        if not isinstance(obj, ROOT.TH1):
            continue

        mo = re.match("dat_cat(\d+)$", key.GetName())
        if not mo:
            continue

        retval.add(int(mo.group(1)))


    if mustCloseFin:
        fin.Close()
        
    return list(retval)

#----------------------------------------------------------------------

def makeOutputDir(template):

    import time, errno

    while True:

        # see also http://stackoverflow.com/questions/273192/
        outputDir = "" + time.strftime(template, time.localtime(time.time()))
        try:
            # just try to create the directory
            # if already existed, we get an exception
            os.mkdir(outputDir)
            return outputDir
        
        except OSError,ex:
            if ex.errno != errno.EEXIST:
                raise ex

            # wait a bit before trying again
            time.sleep(1)


#----------------------------------------------------------------------

def cumulativeNormalFunction(x):
    """ returns the integral of - infinity to x of the normal distribution """

    # see http://stackoverflow.com/a/457475/288875 (and one of the comments)
    from scipy.special import erf
    import math

    # see http://en.wikipedia.org/wiki/Error_function#The_name_.22error_function.22
    return 0.5*(1+erf(x/math.sqrt(2)))

#----------------------------------------------------------------------

def makeBandGraphFromArrays(xvalues, lowerYvalues, upperYvalues):
    xvaluesForFilledGraph = xvalues[:]
    yvaluesForFilledGraph = lowerYvalues[:]

    xvaluesForFilledGraph.extend(xvalues[::-1])
    yvaluesForFilledGraph.extend(upperYvalues[::-1])

    # just to be sure: duplicate first point
    xvaluesForFilledGraph.append(xvaluesForFilledGraph[0])
    yvaluesForFilledGraph.append(yvaluesForFilledGraph[0])
    
    grFilled = makeGraphFromArrays(xvaluesForFilledGraph, yvaluesForFilledGraph)

    return grFilled



#----------------------------------------------------------------------

def getQuantile(sortedValues, quantile):
    """
        returns the given quantile (e.g. 0.5 for the median).
        
        takes care of underflow/overflow effects. sortedValues must
        be sorted already.
        """

    if not sortedValues:
        # no values, no quantile
        return None

    pos_in_array = int(round(quantile * len(sortedValues)))
    if pos_in_array < 0:
        return sortedValues[0]
    elif pos_in_array >= len(sortedValues):
        return sortedValues[-1]
    else:
        return sortedValues[pos_in_array]

#----------------------------------------------------------------------
