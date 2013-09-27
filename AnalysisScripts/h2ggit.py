#!/bin/env python

import os, sys, json, commands
from optparse import OptionParser, make_option
from copy import copy

def run(cmd,msg,expect=0):
    sta,out = commands.getstatusoutput(cmd)
    ### print cmd, "[%d]" % sta
    ### print out
    check(sta,msg,out,expect)
    return out

def check(sta,msg,extra,expect=0):
    if sta != expect:
        exit(sta,msg,extra)
        
def exit(sta,msg,extra=None):
    print
    print msg
    if extra:
        print
        print extra
    print
    os.sys.exit(sta)

def trim(docstring):
    if not docstring:
        return ''
    # Convert tabs to spaces (following the normal Python rules)
    # and split into a list of lines:
    lines = docstring.expandtabs().splitlines()
    # Determine minimum indentation (first line doesn't count):
    indent = sys.maxint
    for line in lines[1:]:
        stripped = line.lstrip()
        if stripped:
            indent = min(indent, len(line) - len(stripped))
            # Remove indentation (first line is special):
    trimmed = [lines[0].strip()]
    if indent < sys.maxint:
        for line in lines[1:]:
            trimmed.append(line[indent:].rstrip())
            # Strip off trailing and leading blank lines:
    while trimmed and not trimmed[-1]:
        trimmed.pop()
    while trimmed and not trimmed[0]:
        trimmed.pop(0)
    # Return a single string:
    return '\n'.join(trimmed)
                                                                                                            

def wrapcommand(method):
    def inner_func(*args, **kwargs):
        for a in args:
            if a == "-h" or a == "--help":
                helpstr = """Help for %s:

%s
"""        % ( method.__name__, trim(method.__doc__) )
                sys.exit( helpstr )
        method(*args)
        
    return inner_func

  
class GitHelper:

    _commands = ["ls-issues","tagme","tagit"]
    
    @staticmethod
    def commands():
        commands = GitHelper._commands[0]
        for a in GitHelper._commands[1:]:
            commands += ", %s" %a
        return commands
    
    def __init__(self,gh,settings,options):
        self._gh = gh
        self._options = copy(settings)
        for nam,val in options.__dict__.iteritems():
            self._options[nam] = val
            
        self._upstream = self._gh.repository(self._options['organization'],self._options['repo'])
        self._fork = self._gh.repository(self._options['user'],self._options['repo'])
        self._upstreamtags = set()
        self._forktags = set()
        self._localtags = set()


    def command(self,name,args):
        if name not in self._commands:
            raise NameError('No command %s defined' % name)
        getattr(self,name.replace("-","_"))(*args)
    
    def gettags(self):
        out = run("git ls-remote %s" % self._upstream.ssh_url,"Unable to list tags from %s " % self._upstream.ssh_url)
        for l in out.split("\n"):
            try:
                sha, ref = l.split("\t")
            except:
                continue
            if "tag" in ref:
                self._upstreamtags.add(ref.split("/",3)[-1])

        out = run("git ls-remote %s" % self._fork.ssh_url,"Unable to list tags from %s " % self._fork.ssh_url)
        for l in out.split("\n"):
            try:
                sha, ref = l.split("\t")
            except:
                continue
            if "tag" in ref:
                self._forktags.add(ref.split("/",3)[-1])

        out = run("git tag -l","Unable to list local tags.")
        for l in out.split("\n"):
            self._localtags.add(l)

    @wrapcommand
    def tagme(self,tagname,title):
        """
        Prepare a tagme request.
           usage: tagme <tagname> <title>
        """
        self.gettags()
        print ""
        print "Preparing tagme request %s" %tagname
        print "----------------------------------------------------"
        print title
        print "----------------------------------------------------"
        print
        
        print
        print " * Getting status of local repository"
        run("git status --porcelain | grep -v h2ggit.py | egrep -v '^\\?'",
            "There are uncommitted changes in your workarea. Commit them first.",256)

        localbr = run("git branch | grep '*'","").rsplit(" ",1)[-1]
        print "Local branch is %s" % localbr

        print
        print " * Checking if the tag already exists"        
        if tagname in self._upstreamtags:
            exit(-1,"Tag %s already in upstream repository." % tagname)
        if tagname in self._forktags:
            print "Warning: tag %s already in forked repository." % tagname            
        if not tagname in self._localtags:
            print "Tag not found. Creating it."
            run("git tag -a %s -m '%s'" % (tagname, title), "Unable to tag local repository" )
        else:
            print
            print " * Tag already exists. Checking differences."
            print run("git diff %s  | filterdiff h2ggit.py | grep [+-]" % tagname,
                      "Tag %s exists but does not correspond to your workarea." % tagname,256)

        self.gettags()
        if not tagname in self._localtags:
            exit(-1,"Failed to tag local repository.")

        print
        print " * Checking if tagme request already exists."
        head = "%s:tagme_%s" % ( self._options['user'], tagname )
        tagmes = self.gettagmes()
        for t in tagmes:
            if t.head == head or t.tagme == tagname:
                self.ls_issues()
                exit(1,"Error: tagme request for %s alerady exists.\nNothing done." % tagname)
                
        print
        print " * Pushing to personal fork"
        run("git push -u %s %s:tagme_%s" % (self._fork.ssh_url, localbr, tagname),"")
        run("git push %s %s" % (self._fork.ssh_url, tagname),"")

        print
        print " * Opening an issue for the tagme request."
        issue = self._upstream.create_issue(title=title,
                                            labels=["tagme"],
                                            body="head: %s\nbase: master\ntagme: %s"
                                            % (head, tagname)
                                            )
        self.ls_issues()

    def gettagmes(self):
        issues = self._upstream.iter_issues(state='open')
        tagmes = []
        
        for i in issues:
            for l in  i.labels:
                if str(l) == "tagme":
                    for l in i.body.split("\n"):
                        if "head:" in l:
                            i.head = l.rsplit(" ",1)[-1]
                        if "tagme:" in l:
                            i.tagme = l.rsplit(" ",1)[-1]
                        if "base:" in l:
                            i.base = l.rsplit(" ",1)[-1]
                    tagmes.append(i)
                    break

        return tagmes
    
    @wrapcommand
    def ls_issues(self):
        """List all open issues.
        """
        pulls  = self._upstream.iter_pulls(state='open')
        issues = self._upstream.iter_issues(state='open')

        pullsdic = {}
        for p in pulls:
            pullsdic[p.number]=p 

        print "Open issues"
        print "-----------"
        print
        for i in issues:
            i.tagme = False
            print ("#%d" % i.number).rjust(5), str(i.title).ljust(10),
            print " [",
            for l in  i.labels:
                print str(l),
                if str(l) == "tagme":
                    i.tagme = True
            print "]"
            if i.number in pullsdic:
                print "      %s" % str(pullsdic[i.number].html_url)
            else:
                print "      %s" % str(i.html_url)
            if i.tagme:
                for l in i.body.split("\n"):
                    print "      %s" % l
            print

    @wrapcommand
    def tagit(self,tagname):
        """
        Accept a tagme request.

        Usage: tagit <tagname>
        """
        self.gettags()
        issues = self._upstream.iter_issues(state='open')

        issue = None
        tagmes = self.gettagmes()
        for t in tagmes:
            if t.tagme == tagname:
                issue = t
                break            

        if not issue:
            exit("No tagme issue found %s" % tagname,1)

        user,branch = issue.head.split(":")
        remote = self._gh.repository(user,self._options['repo'])

        print
        print "Adding remotes"
        run("(git remote | grep fork_%s) || (git remote add fork_%s %s)" % (user,user,remote.ssh_url), "" )
        run("(git remote | grep upstream) || (git remote add upstream %s)" % (self._upstream.ssh_url), "" )

        if tagname in self._upstreamtags:
            print
            print "Tag is already in upstream repository."
        else:
            print
            print "Checking out tag to be pushed"
            run("git fetch fork_%s" % user, "" )
            run("git checkout %s" % (tagname), "" )
        
            print
            print "Pushing the tag"
            out = run("git push upstream %s" % (tagname), "" )
            print out
            issue.create_comment(out)
        
        print
        print "Checking out remote branch"
        run("git checkout -b tagme_%s fork_%s/tagme_%s" % (tagname,user,tagname), "" )

        print
        print "Comparing remote branch with upstream/master"
        run("git fetch upstream", "")

        commits = run("git log --no-merges --oneline tagme_%s ^upstream/master" % tagname, "")
        if commits != "":
            print
            print "There are commits to be merged\n%s" % commits

            if not issue.pull_request:
                print "Making pull request"
                self._upstream.create_pull_from_issue(issue.number,base=issue.base,head=issue.head)
            
        print

def install(repo,tag):
    cwd = os.getcwd()
    scratch = os.path.expanduser("~/.scratch")
    if not os.path.exists(scratch):
        os.mkdir(scratch)
    os.chdir(scratch)
    local = os.path.basename(repo).replace(".git","")
    sta,out = commands.getstatusoutput("rm -rf %s" % (local) )
    print out
    sta,out = commands.getstatusoutput("git clone -b %s %s" % (tag, repo) )
    if sta != 0:
        print "Failed to clone %s %s\n%s" % ( repo, tag, out )
        sys.exit(sta)
    print out
    os.chdir(local)
    sta,out = commands.getstatusoutput("python ./setup.py install --user")
    if sta != 0:
        print "Failed to install %s\n%s" % ( local, out )
        sys.exit(sta)
    print out
    os.chdir(cwd)

def setuplib():
    print
    print 
    print "This is the first time you run this script."
    print "I will now install the dependencies."
    deps = [("https://github.com/kennethreitz/requests.git","v1.2.3"),
            ("git://github.com/musella/github3.py.git","0.7.0_patch_issues_url"),
            ]
    for repo,tag in deps:
        print repo, tag
        install(repo,tag)
    print "done"
    print

def setupssett(settings):
    from github3 import login, authorize
    from getpass import getpass
    user = ''
    password = ''
    while not user:
        print "GitHub user name (will be asked only the first time): ",
        user = raw_input()

    while not password:
        password = getpass('GitHub password for {0}: '.format(user))

    note = 'h2ggit.py'
    note_url = 'http://cern.ch/cms'
    scopes = ['user', 'repo']

    auth = authorize(user, password, scopes, note, note_url)

    fout = open(settings,"w+")

    settings = { "token" : auth.token, "id" : auth.id, "user" : user }
    fout.write( json.dumps(settings)  )
    fout.close()
    
    return login(token=auth.token), settings

def main(options,args):
    settings=os.path.expanduser("~/.h2ggit.json")
    try:
        from github3 import login
    except Exception, e:
        setuplib()
        sys.exit("All dependencies installed. Please run again the script")
        
    if os.path.exists(settings):
        pin = open(settings)
        settings = json.loads( pin.read() )
        pin.close()
        
    if 'token' in settings:
        gh = login( token=settings['token'] )
    else:
        gh,settings = setupssett(settings)
        
    command = args[0]
    helper = GitHelper(gh,settings,options)

    helper.command(command,args[1:])
    
    
if __name__ == "__main__":
    
    cmds = GitHelper.commands()
    parser = OptionParser(
        usage="""%%prog [options] <command> [command args]

The main purpose of this script is to handle `tagme` requests, i.e. requests to
propagate tags created in a personal fork to the upstream repository.

The workflow is the following:
* Request creation:
  * Issued by the tag author with the `tagme` command.
  * A tag is created locally and committed to a dedicated branch on the presonal
  fork.
  * A `tagme` issue is opened on github with the information needed to propagate
  the tag an possible other changes.
  
* Request acceptance:
  * Issued by a librarian with the `tagit` command.
  * The tagme issue is retrieved from github and import the tag is imported
  'as is' from the personal fork to the upstream repository.
  * If the tag contains extra commits not present in the upstrem repository, the
  issue is then turned into a pull request.
  
Commands:
%s

See %%prog <command> --help for the syntax of specific commands
        """ % cmds,
        option_list=[
            make_option("-o", "--organization",
                        action="store", type="string", dest="organization",
                        default="h2gglobe",
                        help="organization name",
                        ),
            make_option("-r", "--repository",
                        action="store", type="string", dest="repo",
                        default="h2gglobe",
                        help="repository name",
                        ),
            ])
    parser.disable_interspersed_args()
    
    (options, args) = parser.parse_args()

    main(options,args)
