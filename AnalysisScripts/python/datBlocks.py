
class confBlock:
  def __init__(self):
   self.outfile  = ""
   self.histfile = ""
   self.files = []
   self.confs = []
   self.comments = ''
  def nfiles(self):
    return len(self.files)
  def print_conf(self):
    if len(self.comments) > 0:
      print "**************************************"
      print self.comments
      print "**************************************"

class datBlock:
  def __init__(self,):
   self.vardef  = []
  def clear(self):
   self.vardef  = []
