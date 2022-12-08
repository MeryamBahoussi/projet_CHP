#! /usr/bin/env python
# $Id: depmod.py,v 1.3 2004/12/15 09:26:29 gayp Exp $
# Pierre Gay
# generation de dependances four les fichiers sources Fortran (77/90)

import string
import sys
import os.path

def get_svn_info():
    remsh = "svn info"
    for ligne in os.popen(remsh).readlines():
        pos = string.find(ligne, "Rev:")
        if (pos >= 0):
            return ligne[pos+4:]
    return None

def get_svn_version():
    i=0
    try:
        for ligne in open("./.svn/entries").readlines():
            i=i+1
            if (i == 11):
                return ligne
    except:
        return "0"

def get_svn_depot():
    i=0
    for ligne in open("./.svn/entries").readlines():
        i=i+1
        if (i == 6):
            return ligne

def get_svn_look():
    # "svnlook youngest " + get_svn_depot()
    remsh = "svnlook youngest /stock/geniaero/SVN_SOURCES/kiss3d"
    for ligne in os.popen(remsh).readlines():
        return ligne

def search_file(filename,paths=['.']):
    res = 0
    for p in paths:
        newpath=os.path.join(p,filename)
        if os.path.exists(newpath):
            return newpath
    return None

def get_modules(filename):
    "returns a tuple containing the list of defined modules and the list of used modules of a fortran source file"

    f=open(filename,'r')
    lines = f.readlines()
    f.close()

    defined=[]
    used=[]
    included=[]

    for lsave in lines:
        l=string.expandtabs(string.lower(lsave)[:-1],1)
        words=string.split(string.lstrip(l))
        if len(words) > 0:
            if words[0] == 'use':
                used.append(string.split(words[1],',')[0])
            if words[0] == 'module':
                if len (words) == 2 or words[1] != "procedure":
                    defined.append(words[1])
            if words[0] == 'include':
                newstring = string.replace(words[1],'\'','')
                newstring = string.replace(newstring,'\"','')
                included.append(newstring)
        l=string.expandtabs(lsave[:-1],1)
        words=string.split(string.lstrip(l))
        if len(words) > 0:
            if words[0] == '#include':
                newstring = string.replace(words[1],'\"','')
                included.append(newstring)

    return defined,used,included

def set_dependencies(filename,modfile,moddep,incdep={}, defd = {}):
    "adds dependencies for 'filename' in dictionnaries"
    
    defined,used,included=get_modules(filename)

#    sys.stderr.write("file=%s includes=%s\n"%(filename,`included`))

    defd[filename] = defined

    for module in defined:
        modfile[module]=filename

    if len(used) > 0:
        moddep[filename]=used

    if len(included) > 0:
        incdep[filename]=included

def dependency_list(filenames,excludes,upperinc=0,incpath=['.'], prefix = ""):
    "builds dependency string for 'filename'"
    
    res = []
    modfile={}
    moddep={}
    incdep={}
    defd = {}

    for exc in excludes:
        filenames.remove(exc)

    for filename in filenames:
        set_dependencies(filename,modfile,moddep,incdep, defd)

    for filename in filenames:
        objname = os.path.join(prefix, filename)
	comment = "# %s modules definis: %s" % (objname, `defd[filename]`)
        depstring = objname + ' :'
        deptoadd=0

        if moddep.has_key(filename):
            deps={}
            for mod in moddep[filename]:
                try:
                    deps[modfile[mod]]=1
                except:
                    sys.stderr.write("warning: module "+mod+" non defini...\n")
                    #return []

            if len(deps.keys()) > 0:
                for dep in deps.keys():
                    depstring = depstring + ' ' + os.path.join(prefix, dep)
                deptoadd = 1

        if incdep.has_key(filename):
            incs={}
            for inc in incdep[filename]:
                incs[inc]=1

            if len(incs.keys()) > 0:
                for inc in incs.keys():
                    if upperinc: inc = string.upper(inc)
                    cmpfile = search_file(inc,incpath)
                    if cmpfile:
                        depstring = depstring + ' ' + cmpfile
                deptoadd = 1

        res.append (comment)
        if deptoadd:
            res.append(depstring)

        res.append("")
               

    return res

if __name__ == '__main__':
    args = sys.argv[1:]
    excludes=[]
    incpath=[]
    prefix = ""
    files = []
    for arg in args:
        if arg[0:5] == '-exc:':
            excludes.append(arg[5:])
            continue
        if arg[0:2] == '-I':
            incpath.append(arg[2:])
            continue
        if arg[0:3] == '-o:':
            prefix = arg[3:]
            continue
        files.append(arg)

    deplist = dependency_list(files,excludes,0,incpath, prefix)
    for dep in deplist:
        res = string.replace(dep,'.f90','.o')
        res = string.replace(res,'.F90','.o')
        res = string.replace(res,'.F','.o')
        print res
    print "SVN=" + get_svn_version()
