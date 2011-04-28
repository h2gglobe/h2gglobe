# Python Configuration Handler for Util
# Original Author - Nicholas Wardle

PYDEBUG = 0

# System python imports
import sys,os
import array

# Local python imports
from datBlocks import *
from makeFilelist import *

class configProducer:

  def __init__(self,Ut,conf_filename,Type):

    print "h2gglobe: step %d, with Config %s" \
	 %(Type,conf_filename)

    self.ut_   = Ut;
    self.type_ = Type;
    
    self.conf_filename = str(conf_filename)
    self.lines_ = []

    self.conf_    = confBlock()
    self.plotvar_ = datBlock()
  
    if self.type_ == 1:
      self.init_reduce()
      self.init_cuts()

    elif self.type_ == 2:
      self.init_loop()
      self.init_cuts()
      self.init_histos()
      self.init_counters()

    elif self.type_ == 3:
      self.init_loop()

    else: 
      sys.exit("No Such Type As: %d"%self.type_)
      
    
  def read_file(self,conf_filename,lines=None):
    if lines == None:
      self.lines_ = [ ]
      lines = self.lines_
    if not os.path.isfile(conf_filename):
      sys.exit("No Configuration file named - %s"%conf_filename)
    else:
      cf = open(conf_filename,"r")
      comment_status = False
      for l in cf.readlines():
        line = l.strip(" ").lstrip(" ")
        if len(line) < 2: continue
        if line[0:2] == '->':	
          if comment_status: comment_status = False
          else:		  comment_status = True
          
        elif not comment_status:
          lines.append(line)
        else:
          self.conf_.comments+=line

  def init_cuts(self):
    self.read_dat_cuts('cuts.dat')
    for dum in self.plotvar_.vardef:
       self.ut_.AddCut(
         dum['cutname'	  ]
         ,dum['ncat' 	  ]
         ,dum['dir'  	  ]
         ,dum['fin'  	  ]
         ,dum['cutValuel' ]
         ,dum['cutValueh' ]
         )
       
  def init_counters(self):
    self.read_dat_counters('counters.dat')
    self.ut_.InitCounters()
    for dum in self.plotvar_.vardef:
      self.ut_.AddCounter(
			  dum['ncat']
        ,dum['countername' ]
        ,dum['denomname1'  ]
        ,dum['denomname2'  ]
        ,dum['denomname3'  ]
        )
      
  def init_histos(self):
    self.read_dat_plotvariables('plotvariables.dat')
    self.ut_.InitHistos()
    for dum in self.plotvar_.vardef:
       self.ut_.BookHisto(
         dum['htyp'	]
         ,dum['plot'	]
         ,dum['default'	]
         ,dum['ncat'	]
         ,dum['xbins'	]
         ,dum['ybins'	]
         ,dum['xmin'	]
         ,dum['xmax'	]
         ,dum['ymin'	]
         ,dum['ymax'	]
         ,dum['name'	]
			 )

  def init_reduce(self):
    self.read_config_reduce(self.conf_filename)
    if PYDEBUG: self.conf_.print_conf()
    self.add_files()
    self.ut_.SetTypeRun(self.type_,self.conf_.outfile)
    # self.ut_.ReadInput(self.type_)
    ## self.read_struct_from_file(self.ut_.loops.vtxAlgoParams,)

  def init_loop(self):
    self.read_config_loop(self.conf_filename)
    if PYDEBUG: self.conf_.print_conf()
    self.add_files()
    self.ut_.SetTypeRun(self.type_,self.conf_.histfile)
    for dum in self.conf_.confs:
      self.ut_.DefineSamples(dum['Nam'  ] 
			    ,dum['typ'  ] 
			    ,dum['ind'  ] 
			    ,dum['draw' ] 
			    ,dum['red'  ] 
			    ,dum['tot'  ] 
			    ,dum['intL' ] 
			    ,dum['lum'  ] 
			    ,dum['xsec' ] 
			    ,dum['kfac' ] 
			    ,dum['scal' ]
			    )
 
  def add_files(self):
    for t_f in self.conf_.files:
      self.ut_.AddFile(t_f[0],t_f[1])

# File Parsers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  def read_dat_plotvariables(self,f): 
     "Parsing of the plotting variables dat"
     self.plotvar_.clear()
     self.read_file(f)
     map_c   = {}
     default = 0
     for line in self.lines_:       
      if len(line) < 2: continue
      if (len(line.split()) < 1): continue
       # Decide whether this is a define line or a file line:
      if "default" in line:   
        split_line = line.split()
	# We have the definiion line
        for sp in split_line:
          val = sp.split("=")
	  if val[0] == "default":
           default = int(val[1])
          else: sys.exit("Unrecognised Argument:\n ' %s ' in line:\n ' %s '"
			%(val[0],line))

      elif "name" in line:
        # We have one of the file def lines
        split_line = line.split()
        for sp in split_line:
	 val = sp.split("=")
         if val[0] == "htyp":
           map_c["htyp"] = int(val[1])
         elif val[0] == "plot":
           map_c["plot"] = int(val[1])
         elif val[0] == "ncat":
           map_c["ncat"] = int(val[1])
         elif val[0] == "xbins":
           map_c["xbins"] = int(val[1])
         elif val[0] == "ybins":
           map_c["ybins"] = int(val[1])
         elif val[0] == "xmin":
           map_c["xmin"] = float(val[1])
         elif val[0] == "xmax":
           map_c["xmax"] = float(val[1])
         elif val[0] == "ymin":
           map_c["ymin"] = float(val[1])
         elif val[0] == "ymax":
           map_c["ymax"] = float(val[1])
         elif val[0] == "name":
           map_c["name"] = str(val[1])
         

         else: sys.exit("Unrecognised Argument:\n ' %s ' in line:\n ' %s '"
			%(val[0],line))

        self.plotvar_.vardef.append(map_c.copy())
      else: sys.exit("Config Line Unrecognised:\n ' %s '"%line)
     for cc in self.plotvar_.vardef:
	cc["default"] = default

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def read_dat_counters(self,f): 
     print "Parsing of the counters dat"
     self.plotvar_.clear()
     self.read_file(f)
     map_c   = {}
     for line in self.lines_:       
      if len(line) < 2: continue
      if (len(line.split()) < 1): continue

      if "countername" in line:
        # We have one of the file def lines
        split_line = line.split()
        for sp in split_line:
	 val = sp.split("=")
         if val[0] == "ncat":
           map_c["ncat"] = int(val[1])
         elif val[0] == "countername":
           map_c["countername"] = str(val[1])
         elif val[0] == "denomname1":
           map_c["denomname1"] = str(val[1])
         elif val[0] == "denomname2":
           map_c["denomname2"] = str(val[1])
         elif val[0] == "denomname3":
           map_c["denomname3"] = str(val[1])

         else: sys.exit("Unrecognised Argument:\n ' %s ' in line:\n ' %s '"
			%(val[0],line))

        self.plotvar_.vardef.append(map_c.copy())
      else: sys.exit("Config Line Unrecognised:\n ' %s '"%line)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  def read_dat_cuts(self,f): 
     "Parsing of the cuts dat"
     self.plotvar_.clear()
     self.read_file(f)
     map_c   = {}
     for line in self.lines_:       
      if len(line) < 2: continue
      if (len(line.split()) < 1): continue

      if "cutname" in line:
        # We have one of the file def lines
        split_line = line.split()

        val_arr    = []
        minval_arr = []
        maxval_arr = []

        for sp in split_line:
	 val = sp.split("=")
         if val[0] == "cutname":
           map_c["cutname"] = str(val[1])
         elif val[0] == "ncat":
           map_c["ncat"] = int(val[1])
         elif val[0] == "dir":
           map_c["dir"] = int(val[1])
         elif val[0] == "fin":
           map_c["fin"] = int(val[1])
         elif val[0] == "minval":
           minval_arr.append(float(val[1]))
	 elif val[0] == "maxval":
           maxval_arr.append(float(val[1]))
         elif val[0] == "val":
	   val_arr.append(float(val[1]))

         else: sys.exit("Unrecognised Argument:\n ' %s ' in line:\n ' %s '"
			%(val[0],line))
        if map_c["dir"] == 2 : 
	   map_c["cutValuel"] = array.array('f',minval_arr)
	   map_c["cutValueh"] = array.array('f',maxval_arr)
	else: 
	   map_c["cutValuel"] = array.array('f',val_arr)
	   map_c["cutValueh"] = array.array('f',[float(-9999.)\
					         for v in val_arr])
	
        self.plotvar_.vardef.append(map_c.copy())
      else: sys.exit("Config Line Unrecognised:\n ' %s '"%line)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def read_config_reduce(self,f):
     "Parsing of the reduction configuration"        
     self.read_file(f)
     comment_status = False 

     for line in self.lines_:
       # Decide whether this is a define line or a file line:
       if "output" in line:   
         split_line = line.split()
         # We have the definition line
         for sp in split_line:
           val = sp.split("=")
           if val[0] == "output":
             self.conf_.outfile = str(val[1])
           else: sys.exit("Unrecognised Argument:\n ' %s ' in line:\n ' %s '"
                          %(val[0],line))
           
       elif "Fil=" in line or "Dir=" in line:
         values = { "CaDir" : "", "Dir" : "", "typ" : -1, "Fil" : ""  }; 
         # We have one of the file def lines
         split_line = line.split()
         for sp in split_line:
           name,val = sp.split("=")
           if name in values:
             values[name] = type(values[name])(val)
           else:
             sys.exit( "Unrecognised Argument:\n ' %s ' in line:\n ' %s '"% (name,line) )
         directory = values["Dir"]
         cas_directory = values["CaDir"]
         fi_name   = values["Fil"]
         fi_type   = values["typ"]
         
         if fi_name != '':
	   if not os.path.isfile(fi_name): 
	     sys.exit("No Input File Named: %s"%fi_name)
             tuple_n = fi_name, fi_type
             self.conf_.files.append(tuple_n)
             
         if cas_directory != '':
           ca_files = makeCaFiles(cas_directory)
           for file_s in ca_files:
             self.conf_.files.append((str(directory+'/'+file_s),fi_type))
             
         if os.path.isdir(directory): 
           di_files = os.listdir(directory)
 	   di_files = filter(lambda x: ".root" in x, di_files)
           for file_s in di_files:
             self.conf_.files.append((str(directory+'/'+file_s),fi_type))

       # read a generic member of the LoopAll class
       else:
         split_line = line.split()
         struct_name = split_line.pop(0)
         try:
           struct = self.ut_.loops.__getattribute__(struct_name)            
           if os.path.isfile(split_line[0]):
             self.read_struct_from_file(split_line[0],struct)
           else:
             self.read_struct(split_line,struct)
         except:
           sys.exit("Unkown data member %s at line:\n\n%s\n"% (struct_name, line) )
           
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def read_config_loop(self,f):
     "Parsing of the looper configuration"        
     self.read_file(f)
     comment_status = False 
     intL	    = 0
     map_c	    = {}

     for line in self.lines_:
      if len(line) < 2: continue
      if line[0:2] == '->':	
       if comment_status: comment_status = False
       else:		  comment_status = True

      else:
       if not comment_status:
        # Decide whether this is a define line or a file line:
        if "histfile" in line:  
        # We have the definition line
        
         split_line = line.split()
         for sp in split_line:
          val = sp.split("=")
	  if val[0] == "output":
	   self.conf_.outfile = str(val[1])
	  elif val[0] == 'histfile':
	   self.conf_.histfile = str(val[1])
	  elif val[0] == 'intL':
           intL = float(val[1])
          else: sys.exit("Unrecognised Argument:\n ' %s ' in line:\n ' %s '"
			%(val[0],line))

        
        elif "Fil=" in line or "Dir=" in line:
         directory = ''
         cas_directory = ''
         fi_name   = ''
         fi_type   = ''
        # We have one of the file def lines
         split_line = line.split()
         for sp in split_line:
	  val = sp.split("=")
          if val[0] == "Fil":
           fi_name = str(val[1])
          elif val[0]== "Dir":
           directory=str(val[1])
          elif val[0]== "CaDir":
           cas_directory=str(val[1])
	  elif val[0] == "typ":
	   fi_type = int(val[1])
           map_c["typ"] = int(val[1])
          elif val[0] == "Nam":
           map_c["Nam"] = str(val[1])
	  elif val[0] == "draw":
	   map_c["draw"] = int(val[1])
	  elif val[0] == "ind":
	   map_c["ind"] = int(val[1])
	  elif val[0] == "tot":
	   map_c["tot"] = int(val[1])
	  elif val[0] == "red":
	   map_c["red"] = int(val[1])
	  elif val[0] == "lum":
	   map_c["lum"] = float(val[1])
	  elif val[0] == "xsec":
	   map_c["xsec"] = float(val[1])
	  elif val[0] == "kfac":
	   map_c["kfac"] = float(val[1])
	  elif val[0] == "scal":
	   map_c["scal"] = float(val[1])
          else: sys.exit("Unrecognised Argument:\n ' %s ' in line:\n ' %s '"
			%(val[0],line))

	 if fi_name != '':
	   if not os.path.isfile(fi_name): 
	     sys.exit("No Input File Named: %s"%fi_name)
	   tuple_n = fi_name, fi_type
	   self.conf_.files.append(tuple_n)
	   self.conf_.confs.append(map_c.copy())

         if cas_directory != '':
           ca_files = makeCaFiles(cas_directory)
           for file_s in ca_files:
             self.conf_.files.append((str(directory+'/'+file_s),fi_type))
             self.conf_.confs.append(map_c.copy())
	
         if os.path.isdir(directory): 
           di_files = os.listdir(directory)
 	   di_files = filter(lambda x: ".root" in x, di_files)
           for file_s in di_files:
             self.conf_.files.append((str(directory+'/'+file_s),fi_type))
             self.conf_.confs.append(map_c.copy())
        

	else: sys.exit("Config Line Unrecognised:\n ' %s '"%line)
       else: self.conf_.comments+=line
       for cc in  self.conf_.confs:
         cc["intL"] = intL 

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def read_struct(self,lines,struct):
    print "Reading structure %s" % ( str(struct) )
    for line in lines:
      split_line = line.split()
      for sp in split_line:
        name,val = [ s.lstrip(" ").rstrip(" ") for s in sp.split("=") ]
        try:
          var = struct.__getattribute__(name)
          var = type(var)(val)
        except AttributeError:
          sys.exit( "Error parsing file %s. Unkown attribute: %s\nline: %s" % ( f, name, line ) )
        except ValueError, e:
          sys.exit( "Error parsing file %s. Syntax error in line: %s\n%s" % ( f, line, e ) )
    

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def read_struct_from_file(self,f,struct):
    lines = []
    print "Reading structure from file %s %s " % ( str(struct), f )
    self.read_file(f,lines)
    
    self.read_struct(lines,struct)
    
    
#--------------------------------------------------------//
#EOF 
