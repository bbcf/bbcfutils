#!/usr/bin/env python

"""
Creates an `MA-plot` to compare transcription levels of a set of genes
(or other features) in two different conditions, from a CSV file.
One can enter several datasets (CSV files) in the same format, each of which
will be plotted in a different color, and be annotated if requested.

The class AnnoteFinder is used to create interactive - clickable - plots.
"""

import sys, os, json, urllib, math, csv, optparse
import numpy
from scipy import stats
import matplotlib.pyplot as plt

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def rstring(len=20):
    import string
    import random
    return "".join([random.choice(string.letters+string.digits) for x in range(len)])

def MAplot(dataset, cols=[2,3],annotate=None, mode="normal", data_format="counts", limits=[None,None,None,None],
           slimits=[None,None], deg=3, bins=30, assembly_id=None, quantiles=True, title="MA-plot"):
    """
    Creates an "MA-plot" to compare transcription levels of a set of genes
    in two different conditions. It returns the name of the .png file produced,
    and the name of a json containing enough information to reconstruct the plot using Javascript.

    :param cols: list containing the numbers of the two columns containing the numeric data to compare.
    :param dataset: list or string, containing names of up to six CSV files with rows
    of the form (feature_name, sample1, sample2, ...).
    :param annotate: in 'normal' mode, choose which for which datasets you want the
    points to be labeled. In a list, enter 1 to annotate, 0 not to annotate, in the
    same order as datasets were entered. E.g. [0,0,1] to annotate only the third of 3
    datasets.
    :param mode: string, display mode:
    * If `normal`, name of genes over 99%/under 1% quantile are displayed.
    * If `interactive`, click on a point to display its name.
    * If `json`, a .json file is produced that allows to reproduce th graph
    in a web interface using Javascript.
    :param data_format: `counts` or `rpkm`.
    :param deg: int, the degree of the interpolating polynomial splines
    :param bins: int, the number of divisions of the x axis for quantiles estimation
    :param assembly_id: string of integer. If an assembly ID is given,
    the json output will provide links to information on genes.
    """
    # Constants:
    if data_format == "counts":
        lower = 1
        spline_xmin = math.log10(math.sqrt(10)) #counts of 2,3 -> log(sqrt(2*3))
        spline_xmax = None
        slimits[0] = slimits[0] or 1
    elif data_format == "rpkm":
        lower = 0
        spline_xmin = math.log10(0.1)
        spline_xmax = None
        slimits[0] = slimits[0] or -1
    min_pts_per_bin = 20
    extreme_ratios = False # print extremes ratios

    # Extract data from CSV
    if isinstance(dataset,str): dataset = [dataset]
    names=[]; means=[]; ratios=[]; pvals=[]; points=[]; groups={}; counts={}
    for data in dataset:
        with open(data,'r') as f:
            dialect = csv.Sniffer().sniff(f.read(1024)); f.seek(0)
            f.readline()
            csvreader = csv.reader(f, dialect=dialect, quoting=csv.QUOTE_NONE)
            n=[]; m=[]; r=[]; p=[]
            for row in csvreader:
                #numpy.seterr(all='raise') # testing
                c1 = float(row[cols[0]-1]); c2 = float(row[cols[1]-1])
                if (c1*c2 > lower):
                    counts[row[0]] = (c1,c2)
                    n.append(row[0])
                    m.append(numpy.log10(numpy.sqrt(c1*c2)))
                    r.append(numpy.log2(c1/c2))
                p.append(None) # future p-values
        groups[data] = zip(n, m, r, p)
        names.extend(n); means.extend(m); ratios.extend(r); pvals.extend(p)
    points = zip(names, means, ratios, pvals)
    xmin = min(means); xmax = max(means); ymin = min(ratios); ymax = max(ratios)

    # Figure initialization
    fig = plt.figure(figsize=[14,9])
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.95)

    # Points
    colors = iter(["black","red","green","cyan","magenta","yellow"])
    datacolors = {}
    for data in dataset:
        pts = zip(*groups[data])
        datacolors[data] = colors.next()
        ax.plot(pts[1], pts[2], ".", color=datacolors[data])

    # Lines (best fit of percentiles)
    if quantiles:
        # Create bins; consider only a subset [spline_xmin, spline_xmax] of meaningful points.
        spline_xmin = spline_xmin or xmin
        spline_xmax = spline_xmax or xmax
        intervals = numpy.linspace(spline_xmin, spline_xmax, bins+1)

        points_in = []
        for i in range(len(intervals)-2,-1,-1): #from l-1 to 0, decreasing
            p_in_i = [p for p in points if p[1]>=intervals[i] and p[1]<intervals[i+1]]
            if len(p_in_i) < min_pts_per_bin:
                intervals = numpy.delete(intervals,i); bins-=1
            else:
                points_in.append(p_in_i)
        points_in.reverse()
        x = intervals[:-1]+(intervals[1:]-intervals[:-1])/2. #the middle point of each bin

        # Compute percentiles in each bin
        spline_annotes=[]; spline_coords={}; extremes=[]
        percentiles = [1,5,25,50,75,95,99]
        for k in percentiles:
            h=[]
            for b in range(bins):
                score = stats.scoreatpercentile([p[2] for p in points_in[b]], k)
                h.append(score)
                if k==1: extremes.extend([p for p in points_in[b] if p[2] < score])
                if k==99: extremes.extend([p for p in points_in[b] if p[2] > score])

            coeffs = numpy.polyfit(x, h, deg)
            smin = slimits[0] or x[0]
            smax =  slimits[1] or 0.85*x[-1]
            x_spline = numpy.array(numpy.linspace(smin, smax, 10*bins))
            y_spline = numpy.polyval(coeffs, x_spline)

            ax.plot(x_spline, y_spline, "-", color="blue")
            #ax.plot(x, h, "o", color="blue") # testing
            spline_annotes.append((k,x_spline[0],y_spline[0])) #quantile percentages
            spline_coords[k] = zip(x_spline,y_spline)

        if extreme_ratios:
            with open("extremes_ratios"+rstring(5),"w") as f:
                c = csv.writer(f,delimiter="\t")
                c.writerow(["Name","countsC1","countsC2","log10Mean","log2Fold"])
                for p in extremes:
                    c.writerow([p[0], str(counts[p[0]][0]), str(counts[p[0]][1]), str(p[1]), str(p[2])])

        # Annotation of splines (percentage)
        for sa in spline_annotes:
            ax.annotate(str(sa[0])+"%", xy=(sa[1],sa[2]), xytext=(-33,-5), textcoords='offset points',
                    bbox=dict(facecolor="white",edgecolor=None,boxstyle="square,pad=.4"))

    # Decoration
    ax.set_xlabel("Log10 of sqrt(x1*x2)")
    ax.set_ylabel("Log2 of x1/x2")
    ax.set_title(title)
    if limits[0] is not None: xmin = limits[0]
    if limits[1] is not None: xmax = limits[1]
    if limits[2] is not None: ymin = limits[2]
    if limits[3] is not None: ymax = limits[3]
    xlen = abs(xmax-xmin); ylen = abs(ymax-ymin)
    plt.xlim([xmin-0.1*xlen,xmax+0.1*xlen])
    plt.ylim([ymin-0.1*ylen,ymax+0.1*ylen])

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
    if mode == "normal":
        figname = rstring()+".png"
        fig.savefig(figname)
        print "png:", figname

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
        splinelabels = {"id": "Spline labels",
                        "data": [spline_coords[k][0] for k in percentiles],
                        "points": {"show":False},
                        "lines": {"show":False},
                        "labels": ["1%","5%","25%","50%","75%","95%","99%"]}
        transparencies = iter([0.12, 0.2, 0.35, 0.35, 0.2, 0.12])
        percentiles.pop(percentiles.index(50))
        for i in percentiles:
            jsdata.append({"id": str(i)+" % Quantile",
                           "data": spline_coords[i],
                           "lines": {"show":True, "lineWidth":0, "fill": transparencies.next()},
                           "color": rgb_to_hex((255,0,255)) })
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

    return 0


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

usage = """

maplot.py [-h] [-c --cols] [-m --mode] [-f --format] [-d --deg] [-b --bins] [-a --assembly]
          [-q --noquantiles] [-n --annotate] [--xmin --xmax --ymin --ymax] [--smin --smax]
          data_1 .. data_n

**Input**: CSV files with rows of the form
            (feature_name, expression_in_condition1, expression_in_condition2).
**Output**: (str,str) - the name of the .png file produced,
            and the name of a json containing enough information to reconstruct
            the plot using Javascript. """

help = iter([
"""The numbers of the two columns containing the numeric data to compare, separated by
commas. E.g. --cols 3,5.""",
"""Display mode: 'normal' for static .pgn output,
'interactive' - clic to display gene names, or
'json' - json output to stdout for Javascript web interface.""",
"""Data type: 'counts' for raw count data (default), 'rpkm' for normalized data.""",
"""Degree of the interpolant percentile splines.""",
"""Number of divisions of the x axis to calculate percentiles.""",
"""Identifier for the Genrep assembly used to add more
information about features into the json output. E.g. -a hg19.""",
"""Don't draw quantile splines. This may improve speed and lisibility in some cases.""",
"""Indication of which datasets to annotate (if 'normal' mode).
Must be a binary string of the same lenght as the number of datasets, \
1 indicating to annotate the corresponding set, 0 not to annotate. \
E.g. For a list of datasets d1,d2,d3, if you want to annotate only d3, \
type --annotate 001.""",
"""Minimum x value to be displayed on the output graph.""",
"""Maximum x value to be displayed on the output graph.""",
"""Minimum y value to be displayed on the output graph.""",
"""Maximum y value to be displayed on the output graph.""",
"""Left bound to draw splines.""",
"""Right bound to draw splines.""",
"""Title of the graph"""
])

def main():
    try:
        parser = optparse.OptionParser(usage=usage, description="Creates an `MA-plot` to \
                                       compare transcription levels of a set of genes \
                                       (or other features) in two different conditions.")

        parser.add_option("-c", "--cols", default='2,3', help = help.next())
        parser.add_option("-m", "--mode", default='normal', help = help.next())
        parser.add_option("-f", "--format", default='counts', help = help.next())
        parser.add_option("-d", "--deg", default=4, type="int", help = help.next())
        parser.add_option("-b", "--bins", default=30, type="int", help = help.next())
        parser.add_option("-a", "--assembly", default=None, help = help.next())
        parser.add_option("-q", "--noquantiles", default=True, action= "store_false", help = help.next())
        parser.add_option("-n", "--annotate", default=None, help = help.next())
        parser.add_option("--xmin", default=None, type="float", help = help.next())
        parser.add_option("--xmax", default=None, type="float", help = help.next())
        parser.add_option("--ymin", default=None, type="float", help = help.next())
        parser.add_option("--ymax", default=None, type="float", help = help.next())
        parser.add_option("--smin", default=None, type="float", help = help.next())
        parser.add_option("--smax", default=None, type="float", help = help.next())
        parser.add_option("-t","--title", default="MA-plot", help = help.next())

        (opt, args) = parser.parse_args()
        args = [os.path.abspath(a) for a in args]

        annotate = None
        if opt.mode not in ["normal","interactive","json"]:
            parser.error("--mode must be one of 'normal','interactive', or 'json'.")
        if opt.format not in ["counts","rpkm"]:
            parser.error("--format must be one of 'counts' or 'rpkm'.")
        if opt.annotate: annotate = [eval(a) for a in opt.annotate] # 0100101 -> [0,1,0,0,1,0,1]
        limits = [None or opt.xmin, None or opt.xmax, None or opt.ymin, None or opt.ymax]
        slimits = [None or opt.smin, None or opt.smax]
        cols = [int(c) for c in opt.cols.split(",")]

        # Program body #
        MAplot(args, cols=cols, mode=opt.mode, data_format=opt.format, limits=limits, slimits=slimits,
                         deg=opt.deg, bins=opt.bins, assembly_id=opt.assembly,
                         annotate=annotate, quantiles=opt.noquantiles, title=opt.title)
        # End of program body #

        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())
