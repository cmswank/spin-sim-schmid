#!/usr/bin/env python

import sys
import os
import re
import getopt
import time

def strip_comments(lines):
    match = re.compile('(.*)((\/\/|\/\*)(.*))').match
    lines = iter(lines)
    line = lines.next()
    while line:
        line = line.strip()
        m = match(line)
        if m is None:
            yield line
            line = lines.next()
        elif m.group(3) in ('//', None):
            yield m.group(1).strip()
            try:
                line = lines.next()
            except StopIteration:
                line = None
        else:
            pass
            # got /* find closing */
            ret = m.group(1)
            temp = m.group(4)
            while True:
                try:
                    i = temp.index('*/')
                except ValueError:
                    if ret:
                        yield ret.strip()
                        ret = None
                    temp += lines.next()
                else:
                    line = temp[i+2:]
                    break

def find_classes(filename, verbose=False):
    class_search = re.compile('(struct|class)\s+(\w+)\s*(\S)?').search
    f = open(filename)
    classes = []
    for i,line in enumerate(strip_comments(f)):
        m = class_search(line)
        if m is not None:
            if verbose:
                print "possible candidate: %s" % (m.group(0),)
            # make sure definition of class
            # (must continue with either '{' or ':'
            while m.group(3) is None:
                line += f.next().strip()
                m = class_search(line)
            if m.group(3) not in ('{', ':'):
                args = (m.group(3), m.group(2), i, filename)
                print "W: invalid char %s for possible class %s at line %d in %s" % args
                continue
            if verbose:
                print "FOUND: %s" % (m.group(2),)
            classes.append(m.group(2))
    return classes

def format_classes(filename, classes):
    ret = '';
    if len(classes) == 0:
        ret += '// no classes in header %s\n' % (filename,)
    else:
        ret += '// Parsed from header %s\n' % (filename,)
        for c in classes:
            ret += '#pragma link C++ class %s;\n' % (c,)
    ret += '\n'
    return ret

def format_header(header_files):
    ret = """/***************************
 * Autogenerated by %s on %s
 * input header files:
""" % (__file__, time.ctime())
    for f in header_files:
        ret += ' * - %s\n' % (f,)
    ret += "*/\n"
    ret += """
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

"""
    return ret

def format_LinkDef(header_files, verbose=False):
    header = format_header(header_files)
    file_classes = [find_classes(f, verbose=verbose) for f in header_files]
    body = ''.join([format_classes(f,c) for f,c in zip(header_files, file_classes)])
    footer = '\n#endif\n'
    return header+body+footer

def main(opts, args):
    if '-n' in opts:
        LinkDef_output = '%s_LinkDef.h' % (opts['-n'],)
    else:
        LinkDef_output = 'LinkDef.h'
    linkDefText = format_LinkDef(args, verbose=('-v' in opts))
    with open(LinkDef_output, 'w') as fout:
        fout.write(linkDefText)

if __name__ == '__main__':
    opts, args = getopt.getopt(sys.argv[1:], 'vn:')
    main(dict(opts), args)