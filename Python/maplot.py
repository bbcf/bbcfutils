"""
Creates an MA-plot from a CSV file.

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
from bbcflib.common import timer
#import matplotlib
#matplotlib.use('Agg') #trying to avoid problems with -X ssh sessions (force backend)
import matplotlib.pyplot as plt

usage = """maplot.py [-h] [-l limspath] [-u via] [-m mode] [-d deg] [-b bins] [-a aid] CSV_data

Creates an "MA-plot" to compare transcription levels of a set of genes
(or other features) in two different conditions.

**Input**: str - a CSV file with rows of the form (feature_name, mean_expression, fold_change).
**Output**: str - the name of the .png file produced, and the name of a json
containing enough information to reconstruct the plot using Javascript.

Options:
-h, --help   print this message and exit.
-l, --lims   str - name of or path to the bein's MiniLIMS to receive output files.
-u  --via    str - protocol, may be 'local' or 'lsf'.
-m, --mode   str - 'normal' - static .pgn output - or 'interactive' - clic to display gene names.
-d, --deg    int - degree of the interpolant percentile splines.
-b, --bins   int - number of divisions of the x axis to calculate percentiles.
-a, --aid  int or str - identifier for the Genrep assembly (e.g. 'hg19', or 76)
Note: the assembly_id is used to add more information on features into the json output.
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg


def MAplot(data, mode="normal", deg=4, bins=30, assembly_id=None):
    """
    Creates an "MA-plot" to compare transcription levels of a set of genes
    in two different conditions. It returns the name of the .png file produced,
    and the name of a json containing enough information to reconstruct the plot using Javascript.

    * *data*:  string, name of a CSV file with rows (feature_name, mean_expression, fold_change)
    * *mode*:  string, display mode:
    - if `interactive`, click on a point to display its name
    - if `normal`, name of genes over 99%/under 1% quantile are displayed
    * *deg*:  int, the degree of the interpolating polynomial splines
    * *bins*:  int, the number of divisions of the x axis for quantiles estimation
    * *assembly_id*:  string of integer. If an assembly ID is given,
    the json output will provide links to information on genes.
    """

    # Extract data from CSV
    names=[]; means=[]; ratios=[]; points=[]
    with open(data,'r') as f:
        header = f.readline()
        for d in ['\t',',',' ',':','-']:
            if len(header.split(d)) == 3:
                delimiter = d; break;
            else:
                print """Each line of the CSV file must be of the form \n
                              Feature_name    Mean    fold_change \n
                           Accepted delimiters: (space) , : - \t    """
        csvreader = csv.reader(f, delimiter=delimiter, quoting=csv.QUOTE_NONE)
        for row in csvreader:
            if float(row[1])!=0:
                names.append(row[0])
                means.append(numpy.log10(float(row[1])))
                ratios.append(numpy.log2(float(row[2])))
    points = zip(names, means, ratios)
    points = [p for p in points if (p[1]!=0 and p[2]!=0)]
    xmin = min(means); xmax = max(means); ymin = min(ratios); ymax = max(ratios)

    # Create bins
    N = len(points); dN = N/bins #points per bin
    rmeans = numpy.sort(means)
    intervals = []
    for i in range(bins):
        intervals.append(rmeans[i*dN])
    intervals.append(xmax)
    intervals = numpy.array(intervals)

    points_in = {}; perc = {}
    for b in range(bins):
        points_in[b] = [p for p in points if p[1]>=intervals[b] and p[1]<intervals[b+1]]
        perc[b] = [p[2] for p in points_in[b]]

    # Figure
    fig = plt.figure(figsize=[14,9])
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.98)
    figname = None

    # Points
    points = zip(*points)
    ax.plot(points[1], points[2], ".", color="black")
    points = zip(*points)

    # Lines (best fit of percentiles)
    annotes=[]; spline_annotes=[]; spline_coords={}
    percentiles = [1,5,25,50,75,95,99]
    for k in percentiles:
        h = numpy.ones(bins)
        for b in range(bins):
            if points_in[b] != []:
                h[b] = stats.scoreatpercentile(perc[b], k)
                if k==1:
                    for p in points_in[b]:
                        if p[2]<h[b]:
                            annotes.append(p[0])
                if k==99:
                    for p in points_in[b]:
                        if p[2]>h[b]:
                            annotes.append(p[0])
            else: h[b] = h[b-1]
        x = intervals[:-1]+(intervals[1:]-intervals[:-1])/2.
        spline = UnivariateSpline(x, h, k=deg)
        xs = numpy.array(numpy.linspace(xmin, xmax, 10*bins)) #to increase spline smoothness
        ys = numpy.array(spline(xs))
        l = len(xs)
        xi = numpy.arange(l)[numpy.ceil(l/6):numpy.floor(8*l/9)]
        x_spline = xs[xi]
        y_spline = ys[xi]
        ax.plot(x_spline, y_spline, "-", color="blue") #ax.plot(x, h, "o", color="blue")
        spline_annotes.append((k,x_spline[0],y_spline[0])) #quantile percentages
        spline_coords[k] = zip(x_spline,y_spline)

    # Decoration
    ax.set_xlabel("Log10 of sqrt(x1*x2)")
    ax.set_ylabel("Log2 of x1/x2")
    for sa in spline_annotes:
        ax.annotate(str(sa[0])+"%", xy=(sa[1],sa[2]), xytext=(-33,-5), textcoords='offset points',
                    bbox=dict(facecolor="white",edgecolor=None,boxstyle="square,pad=.4"))
    if mode == "interactive":
        af = AnnoteFinder( means, ratios, names )
        plt.connect('button_press_event', af)
        plt.draw()
        plt.show()
    else:
        for p in annotes:
            ax.annotate(p[0], xy=(p[1],p[2]) )
    figname = unique_filename_in()+".png"
    fig.savefig(figname)

    # Output for Javascript
    def rgb_to_hex(rgb):
        return '#%02x%02x%02x' % rgb
    jsdata = [{"label": "Data points",
               "data": points,
               "labels": annotes,
               "points": {"symbol":"circle", "show":True},
               "color": "black"},
              ## {"label": "Genes with p-value < " + str(alpha),
              ##  "data": redpoints,
              ##  "labels": annotes_red,
              ##  "pvals": pvals_red,
              ##  "points": {"symbol":"circle", "show":True},
              ##  "color": "red"},
              ## {"label": "Genes with p-value > " + str(alpha),
              ##  "data": blackpoints,
              ##  "labels": annotes_black,
              ##  "pvals": pvals_black,
              ##  "points": {"symbol":"circle", "show":True},
              ##  "color": "black"},
              {"label": "Mean", "data": spline_coords[50],
               "lines": {"show":True}, "color": rgb_to_hex((255,0,255)) },
              {"id": "1% Quantile", "data": spline_coords[1], "lines": {"show":True, "lineWidth":0, "fill": 0.12},
               "color": rgb_to_hex((255,0,255))},
              {"id": "5% Quantile", "data": spline_coords[5], "lines": {"show":True, "lineWidth":0, "fill": 0.18},
               "color": rgb_to_hex((255,0,255))},
              {"id": "25% Quantile", "data": spline_coords[25], "lines": {"show":True, "lineWidth":0, "fill": 0.3},
               "color": rgb_to_hex((255,0,255))},
              {"id": "75% Quantile", "data": spline_coords[75], "lines": {"show":True, "lineWidth":0, "fill": 0.3},
               "color": rgb_to_hex((255,0,255))},
              {"id": "95% Quantile", "data": spline_coords[95], "lines": {"show":True, "lineWidth":0, "fill": 0.18},
               "color": rgb_to_hex((255,0,255))},
              {"id": "99% Quantile", "data": spline_coords[99], "lines": {"show":True, "lineWidth":0, "fill": 0.12},
               "color": rgb_to_hex((255,0,255))}
             ]
    splinelabels = {"id": "Spline labels",
                    "data": [spline_coords[k][0] for k in percentiles],
                    "points": {"show":False}, "lines": {"show":False},
                    "labels": ["1%","5%","25%","50%","75%","95%","99%"]}
    jsdata = "var data = " + json.dumps(jsdata) + ";\n" \
             + "var splinelabels = " + json.dumps(splinelabels) + ";\n"
    if assembly_id:
        nr_assemblies = urllib.urlopen("http://bbcftools.vital-it.ch/genrep/nr_assemblies.json").read()
        nr_assemblies = json.loads(nr_assemblies)
        if isinstance(assembly_id,str):
            for a in nr_assemblies:
                if a['nr_assembly']['name'] == assembly_id:
                    assembly_id = a['nr_assembly']['id']; break
        url_template = urllib.urlopen("http://bbcftools.vital-it.ch/genrep/nr_assemblies/" \
                                      + str(assembly_id) + "/get_links.json?gene_name=%3CName%3E")
        jsdata = jsdata + "var url_template = " + url_template.read() + ";"
    jsname = unique_filename_in()+".js"
    #jsname = "data.js"
    with open(jsname,"w") as js:
        js.write(jsdata)
    
    return figname, jsname

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
      self.axis = pylab.gca()
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


def results_to_json(lims, exid):
    """Create a JSON string describing the results of execution *exid*.

    The execution is sought in *lims*, and all its output files and
    their descriptions are written to the string.
    """
    produced_file_ids = lims.search_files(source=('execution',exid))
    d = dict([(lims.fetch_file(i)['description'], lims.path_to_file(i))
              for i in produced_file_ids])
    j = json.dumps(d)
    return j


#---------------------------- MAIN ------------------------------#


def main(argv=None):
    limspath = None 
    via = "lsf"
    mode = "normal"
    deg = 4
    bins = 30
    assembly_id = None

    if argv is None:
        argv = sys.argv
    try:
        try:
            opts,args = getopt.getopt(sys.argv[1:],"hu:m:d:b:a:",
                         ["help","via","mode","deg","bins","aid"])
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
                else: mode = "normal"
            elif o in ("-d", "--deg"):
                if isinstance(a,int): deg = a
                else: raise Usage("The polynom degree must be an integer got %s." % (a,))
            elif o in ("-b", "--bins"):
                if isinstance(a,int): bins = a
                else: raise Usage("The number of bins must be an integer got %s." % (a,))
            elif o in ("-a", "--aid"):
                assembly_id = a
            else: raise Usage("Unhandled option: " + o)

        if len(args) < 1:
            raise Usage("maplot.py needs at least one argument (CSV file).")

        data = str(args[0])

        # Program body #
        if limspath:
            M = MiniLIMS(limspath)
            with execution(M) as ex:
                figname, jsname = MAplot(data, mode, deg, bins, assembly_id)
                ex.add(figname, description="png:MAplot of data from file: "+data)
                ex.add(jsname, description="json:json output for file: "+data)
            results_to_json(M, ex.id)
        else:
            figname, jsname = MAplot(data, mode, deg, bins, assembly_id)
            print "png:", figname, "; json:", jsname
        # End of program body #

        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())
