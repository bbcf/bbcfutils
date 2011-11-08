#!/usr/bin/env python

'''
parsegtf.py
by Julien Delafontaine, <julien.delafontaine@epfl.ch>
010.08.2011

Parser for .gtf (Gene Transfer Format) files
GFF structure: <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
In GTF format, [attributes] have the following form: <desc> "value";
'''

import sys, os, pickle, time, getopt

usage = """Usage: parse_gtf.py input.gtf [-o output] [-n nlines]
>>> parse_gtf.py Ciona_intestinalis.JGI2.63.gtf -o parsegtftest -n 20

input.gtf  A file of genomic features in GTF format
Options:
-h     Print this message and exit
-o     Name of the output file (default: species name)
-n     First n lines of the input file to consider (for testing purposes)
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def timer(f):
    """ A decorator that make every decorated function return its execution time """
    def wrapper(*args, **kwargs):
        t1 = time.time()
        result = f(*args, **kwargs)
        t2 = time.time()
        print "Execution time of function", f.__name__, ":", str(t2-t1), "s."
        return result
    return wrapper

#@timer
def parse_gtf(input, output=None, nlines=None):
    """
    Parse a .gtf file and returns a list of dictionaries, each one of them
    corresponding to a line (sequence element) of the file.
    *input*: a gtf file;
    *ouput*: name of the output pickle file;
    *nlines*: select only the *nlines* first lines of the input file for testing purposes.
    """
    with open(input,"rb") as gtf:
        contents = []; l=0;
        for line in gtf:
            l+=1;
            if l>(nlines or float("infinity")):break;
            if len(line) >1:
                sequence_name, source, feature, start, end, score, strand, frame = line.split()[:8]   #common to GFF format
                content = {"sequence_name": sequence_name, "source": source, "feature": feature,
                           "start": start, "end": end, "score": score, "strand": strand, "frame": frame}
                gtf_attributes = line.split('\t')[-1].split(';')[:-1]   #specific to GTF format; last one is '\n'
                for att in gtf_attributes:
                    content[att.split()[0]] = att.split()[1].strip("\"")
                contents.append(content)
    output_filename = "parsed_"+os.path.basename(input)+".pickle"
    with open(output_filename,"wb") as f:
        pickle.dump(contents,f)
    return output_filename

#@timer
def load_pickle(pickle_file):
    """
    Open a pickle file and returns its content as a Python object
    (so that you don't even have to know what a 'pickle' is)
    """
    with open(pickle_file,"rb") as f:
        contents = pickle.load(f)
    return contents

def search_gtf(data, request):
    """
    From a list of dictionaries *data* (as returned by load_pickle),
    retrieves the elements corresponding to a *request* in the form of another dictionary.
    E.g.:
    >>> data = [{1:'a', 2:'b'}, {1:'c', 2:'b'}, {1:'a',2:'b',3:'c'}]
    >>> search_gtf(data, {1:'a',2:'b'})
    [{1:'a', 2:'b'}, {1:'a',2:'b',3:'c'}]
    >>> search_gtf(data, {3:'c'})
    [{1:'a',2:'b',3:'c'}]
    """
    result = []
    for c in contents:
        if all([request[k] == c.get(k) for k in request.keys()]):
            result.append(c)
    return result


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    nlines = None
    output_filename = None

    try:
        try:
            opts, args = getopt.getopt(argv, "h", ["help"])
        except getopt.error, msg:
            raise Usage(msg)
        if len(args) < 1:
            raise Usage("parsegtf.py takes at least 1 argument.")

        input_filename = str(args[0])

        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                print usage
                sys.exit(0)
            if o in ("-n",):
                nlines = a
            if o in ("-o",):
                output_filename = a
            else:
                raise Usage("Unhandled option: " + o)

        ### Program body ###
        contents = parse_gtf(input_filename, output_filename, nlines)

        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2


if __name__ == '__main__':
    sys.exit(main())
