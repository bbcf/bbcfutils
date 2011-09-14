#!/bin/env python
"""
Creates an `MA-plot` to compare transcription levels of a set of genes
(or other features) in two different conditions, from a CSV file.
One can enter several datasets (CSV files), each of which will be plotted
in a different color, and be annotated if requested.

The class AnnoteFinder is used to create interactive - clickable - plots.
"""
import sys, os, pickle, json, pysam, urllib, math, time
import numpy
import csv
import getopt
from bein import *
from bein.util import unique_filename_in
from scipy import stats
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
from bbcflib.common import timer, results_to_json
import matplotlib
matplotlib.use('Agg', warn=False) #force backend (trying to avoid problems with -X ssh sessions)
import matplotlib.pyplot as plt

usage = """maplot.py [-h] [-l limspath] [-u via] [-m mode] [-d deg] [-b bins] [-a aid]
                     [-q quantiles] [-a annotate] data_1 .. data_n

Creates an `MA-plot` to compare transcription levels of a set of genes
(or other features) in two different conditions.

**Input**: CSV files with rows of the form
(feature_name, expression_in_condition1, expression_in_condition2).
**Output**: (str,str) - the name of the .png file produced, and the name of a json
containing enough information to reconstruct the plot using Javascript.

Options:
-h, --help   print this message and exit.
-l, --lims   name of or path to the bein's MiniLIMS to receive output files.
-u  --via    protocol, may be 'local' or 'lsf'.
-m, --mode   'normal' - static .pgn output -, 'interactive' - clic to display gene names,
             or 'json' - .json output for Javascript web interface.
-d, --deg    degree of the interpolant percentile splines.
-b, --bins   number of divisions of the x axis to calculate percentiles.
-a, --aid    identifier for the Genrep assembly (e.g. 'hg19') - used to add more information
             about features into the json output
-q, --quantiles   if other than True, quantile splines wont't be drawn. Thus may
                  improve speed and lisibility in some cases.
-n, --annotate    indication of which datasets to annotate (if 'normal' mode). Must be a
                  binary string of the same lenght as the number of datasets, 1 indicating
                  to annotate the corresponding set, 0 not to annotate. E.g. For a list of
                  datasets d1,d2,d3, if you want to annotate only d3, type --annotate 001.
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def MAplot(dataset, annotate=None, mode="normal", deg=4, bins=30, assembly_id=None, quantiles=True):
    """
    Creates an "MA-plot" to compare transcription levels of a set of genes
    in two different conditions. It returns the name of the .png file produced,
    and the name of a json containing enough information to reconstruct the plot using Javascript.

    * *data*:  list or string, containing names of CSV files with rows of
    the form (feature_name, mean_expression, fold_change)
    * *mode*:  string, display mode:
    - if `interactive`, click on a point to display its name
    - if `normal`, name of genes over 99%/under 1% quantile are displayed
    - if `json`, a .json file is produced that allows to reproduce th graph
    in a web interface using Javascript.
    * *deg*:  int, the degree of the interpolating polynomial splines
    * *bins*:  int, the number of divisions of the x axis for quantiles estimation
    * *assembly_id*:  string of integer. If an assembly ID is given,
    the json output will provide links to information on genes.
    """

    # Extract data from CSV
    if isinstance(dataset,str): dataset = [dataset]
    names=[]; means=[]; ratios=[]; pvals=[]; points=[]; delimiter=None; groups={}
    for data in dataset:
        with open(data,'r') as f:
            header = f.readline()
            for d in ['\t',',',' ',':','-']:
                if len(header.split(d)) in [3,4]:
                    delimiter = d; break;
            if not delimiter:
                print """\n Each line of the CSV file must be of the form \n
                         Feature_name    Expression_cond1    Expression_cond2 \n
                         Accepted delimiters: (space) , : - \\t    \n"""
            csvreader = csv.reader(f, delimiter=delimiter, quoting=csv.QUOTE_NONE)
            n=[]; m=[]; r=[]; p=[]
            for row in csvreader:
                #numpy.seterr(all='raise') # testing
                c1 = float(row[1]); c2 = float(row[2])
                if c1>=1 and c2>=1 and c1<=10e8 and c2<=10e8:
                    """ Counts of 1 may become slightly inferior after normalization processes,
                    thus features with count 1 may be excluded - they shouldn't matter anyway.
                    Another threshold of 10^8 was added in case operations on counts suffered
                    of numerical instability or bad conversions of NA values.
                    You may want to modify these values. """
                    n.append(row[0])
                    m.append(numpy.log10(numpy.sqrt(c1*c2)))
                    r.append(numpy.log2(c1/c2))
                if len(row)==4:
                    p.append(float(row[3]))
                else: p.append(None)
        groups[data] = zip(n, m, r, p)
        names.extend(n); means.extend(m); ratios.extend(r); pvals.extend(p)
    points = zip(names, means, ratios, pvals)
    xmin = min(means); xmax = max(means); ymin = min(ratios); ymax = max(ratios)
    xlen = abs(xmax-xmin); ylen = abs(ymax-ymin)

    # Figure initialization
    fig = plt.figure(figsize=[14,9])
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.98)

    # Points
    colors = iter(["black","red","green","cyan","magenta","yellow"])
    datacolors = {}
    for data in dataset:
        pts = zip(*groups[data])
        datacolors[data] = colors.next()
        ax.plot(pts[1], pts[2], ".", color=datacolors[data])

    # Lines (best fit of percentiles)
    if quantiles:
        # Create bins
        """ If counts < (2,3), their mean < log10(sqrt(2*3)). You may want to
        modify this value to obtaine nicer splines for your data. """
        meaningless = math.log10(math.sqrt(6))
        intervals = numpy.linspace(meaningless, xmax, bins+1) #xmin?

        points_in = []
        for i in range(len(intervals)-2,-1,-1): #from l-1 to 0, decreasing
            p_in_i = [p for p in points if p[1]>=intervals[i] and p[1]<intervals[i+1]]
            if len(p_in_i) < 10:
                intervals = numpy.delete(intervals,i); bins-=1
            else:
                points_in.append(p_in_i)
        points_in.reverse()
        x = intervals[:-1]+(intervals[1:]-intervals[:-1])/2. #the middle point of each bin

        # Compute percentiles in each bin
        spline_annotes=[]; spline_coords={};
        percentiles = [1,5,25,50,75,90,99]
        for k in percentiles:
            h=[]
            for b in range(bins):
                h.append( stats.scoreatpercentile([p[2] for p in points_in[b]], k) )

            #xs = numpy.concatenate((x,[4])) # add a factice point (10,0) - corresponding to zero
            #hs = numpy.concatenate((h,[0]))  # features with expression level of 10^10.
            coeffs = numpy.polyfit(x, h, deg)
            x_spline = numpy.array(numpy.linspace(x[0], 0.85*x[-1], 10*bins))
            y_spline = numpy.polyval(coeffs, x_spline)

            ax.plot(x_spline, y_spline, "-", color="blue")
            #ax.plot(x, h, "o", color="blue") # testing
            spline_annotes.append((k,x_spline[0],y_spline[0])) #quantile percentages
            spline_coords[k] = zip(x_spline,y_spline)

        # Annotation of splines (percentage)
        for sa in spline_annotes:
            ax.annotate(str(sa[0])+"%", xy=(sa[1],sa[2]), xytext=(-33,-5), textcoords='offset points',
                    bbox=dict(facecolor="white",edgecolor=None,boxstyle="square,pad=.4"))

    # Annotation of points
    annotes={}
    if not annotate:
        if mode == "normal": annotate = [0 for data in dataset]
        elif mode == "json": annotate = [1 for data in dataset]
    if mode == "interactive":
        af = AnnoteFinder( means, ratios, names )
        plt.connect('button_press_event', af)
        plt.draw()
        plt.show()
    elif mode == "normal" or mode == "json":
        for data,annote in zip(dataset,annotate):
            if annote==1:
                annotes[data] = []
                for p in groups[data]:
                    ax.annotate(p[0], xy=(p[1],p[2]))
                    annotes[data].append(p[0])

    # Output figure
    ax.set_xlabel("Log10 of sqrt(x1*x2)")
    ax.set_ylabel("Log2 of x1/x2")
    plt.xlim([xmin-0.1*xlen,xmax+0.1*xlen])
    plt.ylim([ymin-0.1*ylen,ymax+0.1*ylen])
    figname = unique_filename_in()+".png"
    fig.savefig(figname)


    # Output for Javascript
    if mode == "json":
        def rgb_to_hex(rgb):
            return '#%02x%02x%02x' % rgb
        jsdata = []
        # - data points
        for data in dataset:
            name = os.path.basename(data)
            gdata = zip(*groups[data])
            pvals = gdata[3]
            datapts = zip(*gdata[1:3])
            jsdata.append({"label": "Data points from file "+name,
                           "data": datapts,
                           "labels": annotes.get(data),
                           "points": {"symbol":"circle", "show":True},
                           "color": datacolors[data],
                           "pvals": pvals})
        # - mean
        jsdata.append({"label": "Mean", "data": spline_coords[50],
                       "lines": {"show":True}, "color": rgb_to_hex((255,0,255)) })

        # - splines
        transparencies = iter([0.12, 0.18, 0.3, 0.3, 0.18, 0.12])
        percentiles.pop(percentiles.index(50))
        for i in percentiles:
            jsdata.append({"id": str(i)+" % Quantile",
                           "data": spline_coords[i],
                           "lines": {"show":True, "lineWidth":0, "fill": transparencies.next()},
                           "color": rgb_to_hex((255,0,255)) })
        splinelabels = {"id": "Spline labels",
                        "data": [spline_coords[k][0] for k in percentiles],
                        "points": {"show":False},
                        "lines": {"show":False},
                        "labels": ["1%","5%","25%","50%","75%","95%","99%"]}
        jsdata = "var data = " + json.dumps(jsdata) + ";\n" \
                 + "var splinelabels = " + json.dumps(splinelabels) + ";\n"

        # - url for more info on features
        if assembly_id:
            nr_assemblies = urllib.urlopen("http://bbcftools.vital-it.ch/genrep/nr_assemblies.json").read()
            nr_assemblies = json.loads(nr_assemblies)
            for a in nr_assemblies:
                if a['nr_assembly']['name'] == assembly_id or a['nr_assembly']['id'] == assembly_id:
                        assembly_id = a['nr_assembly']['id'];
                        md5 = a['nr_assembly']['md5']; break
            url_template = urllib.urlopen("http://bbcftools.vital-it.ch/genrep/nr_assemblies/" +
                            "get_links/" + str(assembly_id) + ".json?gene_name=%3CName%3E&md5=" + md5)
            jsdata = jsdata + "var url_template = " + url_template.read() + ";"
        else: jsdata = jsdata + "var url_template = null;"
        print >> sys.stdout, jsdata

    return figname#!/bin/env python
"""
Creates an `MA-plot` to compare transcription levels of a set of genes
(or other features) in two different conditions, from a CSV file.
One can enter several datasets (CSV files), each of which will be plotted
in a different color, and be annotated if requested.

The class AnnoteFinder is used to create interactive - clickable - plots.
"""
import sys, os, pickle, json, pysam, urllib, math, time
import numpy
import csv
import getopt
from bein import *
from bein.util import unique_filename_in
from scipy import stats
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
from bbcflib.common import timer, results_to_json
import matplotlib
matplotlib.use('Agg', warn=False) #force backend (trying to avoid problems with -X ssh sessions)
import matplotlib.pyplot as plt

usage = """maplot.py [-h] [-l limspath] [-u via] [-m mode] [-d deg] [-b bins] [-a aid]
                     [-q quantiles] [-a annotate] data_1 .. data_n

Creates an `MA-plot` to compare transcription levels of a set of genes
(or other features) in two different conditions.

**Input**: CSV files with rows of the form
(feature_name, expression_in_condition1, expression_in_condition2).
**Output**: (str,str) - the name of the .png file produced, and the name of a json
containing enough information to reconstruct the plot using Javascript.

Options:
-h, --help   print this message and exit.
-l, --lims   name of or path to the bein's MiniLIMS to receive output files.
-u  --via    protocol, may be 'local' or 'lsf'.
-m, --mode   'normal' - static .pgn output -, 'interactive' - clic to display gene names,
             or 'json' - .json output to stdout for Javascript web interface.
-d, --deg    degree of the interpolant percentile splines.
-b, --bins   number of divisions of the x axis to calculate percentiles.
-a, --aid    identifier for the Genrep assembly (e.g. 'hg19') - used to add more information
             about features into the json output
-q, --quantiles   if other than True, quantile splines wont't be drawn. Thus may
                  improve speed and lisibility in some cases.
-n, --annotate    indication of which datasets to annotate (if 'normal' mode). Must be a
                  binary string of the same lenght as the number of datasets, 1 indicating
                  to annotate the corresponding set, 0 not to annotate. E.g. For a list of
                  datasets d1,d2,d3, if you want to annotate only d3, type --annotate 001.
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

@timer
def MAplot(dataset, annotate=None, mode="normal", deg=4, bins=30, assembly_id=None, quantiles=True):
    """
    Creates an "MA-plot" to compare transcription levels of a set of genes
    in two different conditions. It returns the name of the .png file produced,
    and the name of a json containing enough information to reconstruct the plot using Javascript.

    * *data*:  list or string, containing names of CSV files with rows of
    the form (feature_name, mean_expression, fold_change)
    * *mode*:  string, display mode:
    - if `interactive`, click on a point to display its name
    - if `normal`, name of genes over 99%/under 1% quantile are displayed
    - if `json`, a .json file is produced that allows to reproduce th graph
    in a web interface using Javascript.
    * *deg*:  int, the degree of the interpolating polynomial splines
    * *bins*:  int, the number of divisions of the x axis for quantiles estimation
    * *assembly_id*:  string of integer. If an assembly ID is given,
    the json output will provide links to information on genes.
    """

    # Extract data from CSV
    if isinstance(dataset,str): dataset = [dataset]
    names=[]; means=[]; ratios=[]; pvals=[]; points=[]; delimiter=None; groups={}
    for data in dataset:
        with open(data,'r') as f:
            header = f.readline()
            for d in ['\t',',',' ',':','-']:
                if len(header.split(d)) in [3,4]:
                    delimiter = d; break;
            if not delimiter:
                print """\n Each line of the CSV file must be of the form \n
                         Feature_name    Expression_cond1    Expression_cond2 \n
                         Accepted delimiters: (space) , : - \\t    \n"""
            csvreader = csv.reader(f, delimiter=delimiter, quoting=csv.QUOTE_NONE)
            n=[]; m=[]; r=[]; p=[]
            for row in csvreader:
                #numpy.seterr(all='raise') # testing
                c1 = float(row[1]); c2 = float(row[2])
                if c1>=1 and c2>=1 and c1<=10e8 and c2<=10e8:
                    """ Counts of 1 may become slightly inferior after normalization processes,
                    thus features with count 1 may be excluded - they shouldn't matter anyway.
                    Another threshold of 10^8 was added in case operations on counts suffered
                    of numerical instability or bad conversions of NA values.
                    You may want to modify these values. """
                    n.append(row[0])
                    m.append(numpy.log10(numpy.sqrt(c1*c2)))
                    r.append(numpy.log2(c1/c2))
                if len(row)==4:
                    p.append(float(row[3]))
                else: p.append(None)
        groups[data] = zip(n, m, r, p)
        names.extend(n); means.extend(m); ratios.extend(r); pvals.extend(p)
    points = zip(names, means, ratios, pvals)
    xmin = min(means); xmax = max(means); ymin = min(ratios); ymax = max(ratios)
    xlen = abs(xmax-xmin); ylen = abs(ymax-ymin)

    # Figure initialization
    fig = plt.figure(figsize=[14,9])
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.98)

    # Points
    colors = iter(["black","red","green","cyan","magenta","yellow"])
    datacolors = {}
    for data in dataset:
        pts = zip(*groups[data])
        datacolors[data] = colors.next()
        ax.plot(pts[1], pts[2], ".", color=datacolors[data])

    # Lines (best fit of percentiles)
    if quantiles:
        # Create bins
        """ If counts < (2,3), their mean < log10(sqrt(2*3)). You may want to
        modify this value to obtaine nicer splines for your data. """
        meaningless = math.log10(math.sqrt(6))
        intervals = numpy.linspace(meaningless, xmax, bins+1) #xmin?

        points_in = []
        for i in range(len(intervals)-2,-1,-1): #from l-1 to 0, decreasing
            p_in_i = [p for p in points if p[1]>=intervals[i] and p[1]<intervals[i+1]]
            if len(p_in_i) < 10:
                intervals = numpy.delete(intervals,i); bins-=1
            else:
                points_in.append(p_in_i)
        points_in.reverse()
        x = intervals[:-1]+(intervals[1:]-intervals[:-1])/2. #the middle point of each bin

        # Compute percentiles in each bin
        spline_annotes=[]; spline_coords={};
        percentiles = [1,5,25,50,75,90,99]
        for k in percentiles:
            h=[]
            for b in range(bins):
                h.append( stats.scoreatpercentile([p[2] for p in points_in[b]], k) )

            #xs = numpy.concatenate((x,[4])) # add a factice point (10,0) - corresponding to zero
            #hs = numpy.concatenate((h,[0]))  # features with expression level of 10^10.
            coeffs = numpy.polyfit(x, h, deg)
            x_spline = numpy.array(numpy.linspace(x[0], 0.85*x[-1], 10*bins))
            y_spline = numpy.polyval(coeffs, x_spline)

            ax.plot(x_spline, y_spline, "-", color="blue")
            #ax.plot(x, h, "o", color="blue") # testing
            spline_annotes.append((k,x_spline[0],y_spline[0])) #quantile percentages
            spline_coords[k] = zip(x_spline,y_spline)

        # Annotation of splines (percentage)
        for sa in spline_annotes:
            ax.annotate(str(sa[0])+"%", xy=(sa[1],sa[2]), xytext=(-33,-5), textcoords='offset points',
                    bbox=dict(facecolor="white",edgecolor=None,boxstyle="square,pad=.4"))

    # Annotation of points
    annotes={}
    if not annotate:
        if mode == "normal": annotate = [0 for data in dataset]
        elif mode == "json": annotate = [1 for data in dataset]
    if mode == "interactive":
        af = AnnoteFinder( means, ratios, names )
        plt.connect('button_press_event', af)
        plt.draw()
        plt.show()
    elif mode == "normal" or mode == "json":
        for data,annote in zip(dataset,annotate):
            if annote==1:
                annotes[data] = []
                for p in groups[data]:
                    ax.annotate(p[0], xy=(p[1],p[2]))
                    annotes[data].append(p[0])

    # Output figure
    ax.set_xlabel("Log10 of sqrt(x1*x2)")
    ax.set_ylabel("Log2 of x1/x2")
    plt.xlim([xmin-0.1*xlen,xmax+0.1*xlen])
    plt.ylim([ymin-0.1*ylen,ymax+0.1*ylen])
    figname = unique_filename_in()+".png"
    fig.savefig(figname)


    # Output for Javascript
    if mode == "json":
        def rgb_to_hex(rgb):
            return '#%02x%02x%02x' % rgb
        jsdata = []
        # - data points
        for data in dataset:
            name = os.path.basename(data)
            gdata = zip(*groups[data])
            pvals = gdata[3]
            datapts = zip(*gdata[1:3])
            jsdata.append({"label": "Data points from file "+name,
                           "data": datapts,
                           "labels": annotes.get(data),
                           "points": {"symbol":"circle", "show":True},
                           "color": datacolors[data],
                           "pvals": pvals})
        # - mean
        jsdata.append({"label": "Mean", "data": spline_coords[50],
                       "lines": {"show":True}, "color": rgb_to_hex((255,0,255)) })

        # - splines
        transparencies = iter([0.12, 0.18, 0.3, 0.3, 0.18, 0.12])
        percentiles.pop(percentiles.index(50))
        for i in percentiles:
            jsdata.append({"id": str(i)+" % Quantile",
                           "data": spline_coords[i],
                           "lines": {"show":True, "lineWidth":0, "fill": transparencies.next()},
                           "color": rgb_to_hex((255,0,255)) })
        splinelabels = {"id": "Spline labels",
                        "data": [spline_coords[k][0] for k in percentiles],
                        "points": {"show":False},
                        "lines": {"show":False},
                        "labels": ["1%","5%","25%","50%","75%","95%","99%"]}
        jsdata = "var data = " + json.dumps(jsdata) + ";\n" \
                 + "var splinelabels = " + json.dumps(splinelabels) + ";\n"

        # - url for more info on features
        if assembly_id:
            nr_assemblies = urllib.urlopen("http://bbcftools.vital-it.ch/genrep/nr_assemblies.json").read()
            nr_assemblies = json.loads(nr_assemblies)
            for a in nr_assemblies:
                if a['nr_assembly']['name'] == assembly_id or a['nr_assembly']['id'] == assembly_id:
                        assembly_id = a['nr_assembly']['id'];
                        md5 = a['nr_assembly']['md5']; break
            url_template = urllib.urlopen("http://bbcftools.vital-it.ch/genrep/nr_assemblies/" +
                            "get_links/" + str(assembly_id) + ".json?gene_name=%3CName%3E&md5=" + md5)
            jsdata = jsdata + "var url_template = " + url_template.read() + ";"
        else: jsdata = jsdata + "var url_template = null;"
        print >> sys.stdout, jsdata

    return figname


#----------------------------------------------------------#

class AnnoteFinder:
  """
  callback for matplotlib to display an annotation when points are clicked on.  The
  point which is closest to the click and within xtol and ytol is identified.

  Use this function like this:

  plot(xdata, ydata)
  af = AnnoteFinder(xdata, ydata, annotes)
  connect('button_press_event', af)
  """

  def __init__(self, xdata, ydata, annotes, axis=None, xtol=None, ytol=None):
    self.data = zip(xdata, ydata, annotes)
    self.xrange = max(xdata) - min(xdata)
    self.yrange = max(ydata) - min(ydata)
    if xtol is None:
      xtol = len(xdata)*(self.xrange/float(len(xdata)))/10
    if ytol is None:
      ytol = len(ydata)*(self.yrange/float(len(ydata)))/10
    self.xtol = xtol
    self.ytol = ytol
    if axis is None:
      self.axis = plt.gca()
    else:
      self.axis= axis
    self.drawnAnnotations = {}
    self.links = []

  def distance(self, x1, x2, y1, y2):
    """distance between two points"""
    return math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

  def __call__(self, event):
    if event.inaxes:
      clickX = event.xdata
      clickY = event.ydata
      if self.axis is None or self.axis==event.inaxes:
        annotes = []
        for x,y,a in self.data:
          if  clickX-self.xtol < x < clickX+self.xtol and  clickY-self.ytol < y < clickY+self.ytol :
            annotes.append((self.distance(x/self.xrange,clickX/self.xrange,y/self.yrange,clickY/self.yrange),x,y, a) )
        if annotes:
          annotes.sort()
          distance, x, y, annote = annotes[0]
          self.drawAnnote(event.inaxes, x, y, annote)
          for l in self.links:
            l.drawSpecificAnnote(annote)

  def drawAnnote(self, axis, x, y, annote):
    """
    Draw the annotation on the plot
    """
    if (x,y) in self.drawnAnnotations:
      markers = self.drawnAnnotations[(x,y)]
      for m in markers:
        m.set_visible(not m.get_visible())
      self.axis.figure.canvas.draw()
    else:
      t = axis.text(x,y, annote)
      m = axis.scatter([x],[y], marker='d', c='r', zorder=100)
      self.drawnAnnotations[(x,y)] =(t,m)
      self.axis.figure.canvas.draw()

  def drawSpecificAnnote(self, annote):
    annotesToDraw = [(x,y,a) for x,y,a in self.data if a==annote]
    for x,y,a in annotesToDraw:
      self.drawAnnote(self.axis, x, y, a)


#---------------------------- MAIN ------------------------------#


def main(argv=None):
    limspath = None
    via = "lsf"
    mode = "normal"
    deg = 4
    bins = 30
    assembly_id = None
    annotate = None
    quantiles = 'True'

    if argv is None:
        argv = sys.argv
    try:
        try:
            opts,args = getopt.getopt(sys.argv[1:],"hl:u:m:d:b:a:n:q:",
                         ["help","lims=","via=","mode=","deg=","bins=","aid=","annotate=","quantiles="])
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                print usage
                sys.exit(0)
            elif o in ("-l", "--lims"):
                limspath = a
            elif o in ("-u", "--via"):
                if a=="local": via = "local"
                elif a=="lsf": via = "lsf"
                else: raise Usage("Via (-u) can only be \"local\" or \"lsf\", got %s." % (a,))
            elif o in ("-m", "--mode"):
                if a=="interactive": mode = "interactive"
                if a=="json": mode = "json"
                else: mode = "normal"
            elif o in ("-d", "--deg"):
                deg = int(a)
            elif o in ("-b", "--bins"):
                bins = int(a)
            elif o in ("-a", "--aid"):
                assembly_id = a
            elif o in ("-n", "--annotate"):
                annotate = a
            elif o in ("-q", "--quantiles"):
                quantiles = a
            else: raise Usage("Unhandled option: " + o)

        if len(args) < 1:
            raise Usage("maplot.py needs at least one argument (CSV file).")

        if annotate: annotate = [eval(a) for a in annotate] # 0100101 -> [0,1,0,0,1,0,1]
        if quantiles != 'True': quantiles = False
        args = [os.path.abspath(a) for a in args]

        # Program body #
        if limspath:
            M = MiniLIMS(limspath)
            with execution(M) as ex:
                figname = MAplot(args, mode=mode, deg=deg, bins=bins, assembly_id=assembly_id,
                                         annotate=annotate, quantiles=quantiles)
                if figname: ex.add(figname, description="png:MAplot of data")
            results_to_json(M, ex.id)
        else:
            figname = MAplot(args, mode=mode, deg=deg, bins=bins, assembly_id=assembly_id,
                                     annotate=annotate, quantiles=quantiles)
            print "png:", figname
        # End of program body #

        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())



#----------------------------------------------------------#

class AnnoteFinder:
  """
  callback for matplotlib to display an annotation when points are clicked on.  The
  point which is closest to the click and within xtol and ytol is identified.

  Use this function like this:

  plot(xdata, ydata)
  af = AnnoteFinder(xdata, ydata, annotes)
  connect('button_press_event', af)
  """

  def __init__(self, xdata, ydata, annotes, axis=None, xtol=None, ytol=None):
    self.data = zip(xdata, ydata, annotes)
    self.xrange = max(xdata) - min(xdata)
    self.yrange = max(ydata) - min(ydata)
    if xtol is None:
      xtol = len(xdata)*(self.xrange/float(len(xdata)))/10
    if ytol is None:
      ytol = len(ydata)*(self.yrange/float(len(ydata)))/10
    self.xtol = xtol
    self.ytol = ytol
    if axis is None:
      self.axis = plt.gca()
    else:
      self.axis= axis
    self.drawnAnnotations = {}
    self.links = []

  def distance(self, x1, x2, y1, y2):
    """distance between two points"""
    return math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

  def __call__(self, event):
    if event.inaxes:
      clickX = event.xdata
      clickY = event.ydata
      if self.axis is None or self.axis==event.inaxes:
        annotes = []
        for x,y,a in self.data:
          if  clickX-self.xtol < x < clickX+self.xtol and  clickY-self.ytol < y < clickY+self.ytol :
            annotes.append((self.distance(x/self.xrange,clickX/self.xrange,y/self.yrange,clickY/self.yrange),x,y, a) )
        if annotes:
          annotes.sort()
          distance, x, y, annote = annotes[0]
          self.drawAnnote(event.inaxes, x, y, annote)
          for l in self.links:
            l.drawSpecificAnnote(annote)

  def drawAnnote(self, axis, x, y, annote):
    """
    Draw the annotation on the plot
    """
    if (x,y) in self.drawnAnnotations:
      markers = self.drawnAnnotations[(x,y)]
      for m in markers:
        m.set_visible(not m.get_visible())
      self.axis.figure.canvas.draw()
    else:
      t = axis.text(x,y, annote)
      m = axis.scatter([x],[y], marker='d', c='r', zorder=100)
      self.drawnAnnotations[(x,y)] =(t,m)
      self.axis.figure.canvas.draw()

  def drawSpecificAnnote(self, annote):
    annotesToDraw = [(x,y,a) for x,y,a in self.data if a==annote]
    for x,y,a in annotesToDraw:
      self.drawAnnote(self.axis, x, y, a)


#---------------------------- MAIN ------------------------------#


def main(argv=None):
    limspath = None
    via = "lsf"
    mode = "normal"
    deg = 4
    bins = 30
    assembly_id = None
    annotate = None
    quantiles = 'True'

    if argv is None:
        argv = sys.argv
    try:
        try:
            opts,args = getopt.getopt(sys.argv[1:],"hl:u:m:d:b:a:n:q:",
                         ["help","lims=","via=","mode=","deg=","bins=","aid=","annotate=","quantiles="])
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                print usage
                sys.exit(0)
            elif o in ("-l", "--lims"):
                limspath = a
            elif o in ("-u", "--via"):
                if a=="local": via = "local"
                elif a=="lsf": via = "lsf"
                else: raise Usage("Via (-u) can only be \"local\" or \"lsf\", got %s." % (a,))
            elif o in ("-m", "--mode"):
                if a=="interactive": mode = "interactive"
                if a=="json": mode = "json"
                else: mode = "normal"
            elif o in ("-d", "--deg"):
                deg = int(a)
            elif o in ("-b", "--bins"):
                bins = int(a)
            elif o in ("-a", "--aid"):
                assembly_id = a
            elif o in ("-n", "--annotate"):
                annotate = a
            elif o in ("-q", "--quantiles"):
                quantiles = a
            else: raise Usage("Unhandled option: " + o)

        if len(args) < 1:
            raise Usage("maplot.py needs at least one argument (CSV file).")

        if annotate: annotate = [eval(a) for a in annotate] # 0100101 -> [0,1,0,0,1,0,1]
        if quantiles != 'True': quantiles = False
        args = [os.path.abspath(a) for a in args]

        # Program body #
        if limspath:
            M = MiniLIMS(limspath)
            with execution(M) as ex:
                figname = MAplot(args, mode=mode, deg=deg, bins=bins, assembly_id=assembly_id,
                                         annotate=annotate, quantiles=quantiles)
                if figname: ex.add(figname, description="png:MAplot of data")
            results_to_json(M, ex.id)
        else:
            figname = MAplot(args, mode=mode, deg=deg, bins=bins, assembly_id=assembly_id,
                                     annotate=annotate, quantiles=quantiles)
            print "png:", figname
        # End of program body #

        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())
