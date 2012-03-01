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
    import string, random
    return "".join([random.choice(string.letters+string.digits) for x in range(len)])

def guess_file_format(f, sep=None):
    """Guess format of file object *f*"""
    assert isinstance(f,file), "f must be a file object"
    header = 'None'
    f.readline()
    dialect = csv.Sniffer().sniff(f.readline()); f.seek(0)
    if sep: dialect.delimiter = sep
    if csv.Sniffer().has_header(f.readline()):
        f.seek(0); header=f.readline()
    return dialect, header

def name_or_index(cols, dialect, header):
    """Given an array *cols*, detect if elements are indices of *header* or elements of *header*."""
    if all([c in header.split(dialect.delimiter) for c in cols]):
        cols = [header.split(dialect.delimiter).index(c)+1 for c in cols]
    else:
        try: cols = [int(c) for c in cols]
        except ValueError:
            print "\nError: --cols must contain column names or indices (got %s)." % cols
            print "Detected header: %s" % header
    return cols


def MAplot(dataset, cols=[2,3], labels=[1], annotate=None, mode="normal", data_format="counts",
           sep=None, limits=[None,None,None,None], slimits=[None,None], deg=3, bins=30,
           assembly_id=None, quantiles=True, title="MA-plot", extremes=False):
    """
    Creates an "MA-plot" to compare transcription levels of a set of genes
    in two different conditions. It returns the name of the .png file produced,
    and the name of a json containing enough information to reconstruct the plot using Javascript.

    :param dataset: (list or string) names of up to six CSV files with rows
    of the form (feature_name, sample1, sample2, ...).
    :param cols: (list) indices of the two columns with the numeric data to compare.
    :param labels: (list) indices of the columns used as labels. The first column
    specified must contain unique elements.
    :param annotate: (list) in 'normal' mode, choose which for which datasets you want the
    points to be labeled. Enter 1 to annotate, 0 not to annotate, in the
    same order as datasets were entered. E.g. [0,0,1] to annotate only the third of 3
    datasets.
    :param mode: (str) display mode:
    * If `normal`, name of genes over 99%/under 1% quantile are displayed.
    * If `interactive`, click on a point to display its name.
    * If `json`, a .json file is produced that allows to reproduce th graph.
    in a web interface using Javascript.
    :param data_format: (str) `counts` or `rpkm`.
    :param sep: (str) character delimiting the columns.
    :param limits: (list[4]) bounds of the region displayed on the output graph: [minx,maxx,miny,maxy].
    :param slimits: (list[2]) left and right bounds of the section of the splines to be displayed.
    :param deg: (int) the degree of the interpolating polynomial splines.
    :param bins: (int) the number of divisions of the x axis for quantiles estimation.
    :param assembly_id: (str or int) if an assembly ID is given,
    the json output will provide links to information on genes.
    :param quantiles: (bool) if False, no quantile splines are drawn.
    :param title: (str) title to be written on top of the graph.
    :param extremes: (int) create an output file containing features for which ratios were outside the specified
    percentile (two-sided). For the moment, must be 1 or 5. The file is named *extreme_ratios_xxxxx* .
    """
    # Constants:
    if data_format == "counts":
        lower = 1
        #bounds of what is taken into account for the computation of the splines
        spline_xmin = math.log10(math.sqrt(10)) #counts of 2,3 -> log(sqrt(2*3))
        spline_xmax = None
        #bounds of the section of the splines that is displayed
        slimits[0] = slimits[0] or 1
    elif data_format == "rpkm":
        lower = 0
        spline_xmin = math.log10(0.1)
        spline_xmax = None
        slimits[0] = slimits[0] or -1
    min_pts_per_bin = 20
    output_filename = rstring()

    # Extract data from CSV
    if isinstance(dataset,str): dataset = [dataset]
    names=[]; means=[]; ratios=[]; pvals=[]; points=[]; groups={}; counts={}
    for data in dataset:
        with open(data,'r') as f:
            # Guess file format
            dialect,header = guess_file_format(f,sep)
            try: csvreader = csv.reader(f, dialect=dialect, quoting=csv.QUOTE_NONE)
            except TypeError: csvreader = csv.reader(f, dialect='excel-tab', quoting=csv.QUOTE_NONE)
            cols = name_or_index(cols, dialect, header)
            labels = name_or_index(labels, dialect, header)
            # Read the file
            n=[]; m=[]; r=[]; p=[]
            for row in csvreader:
                try: c1 = float(row[cols[0]-1]); c2 = float(row[cols[1]-1])
                except ValueError: continue # Skip line if contains NA, nan, etc.
                if (c1*c2 > lower):
                    counts[row[labels[0]-1]] = (c1,c2)
                    n.append(' | '.join([row[l-1] for l in labels]))
                    m.append(numpy.log10(numpy.sqrt(c1*c2)))
                    r.append(numpy.log2(c1/c2))
                p.append(None) # future p-values, not used yet
        groups[data] = zip(n, m, r, p)
        names.extend(n); means.extend(m); ratios.extend(r); pvals.extend(p)
    points = zip(names, means, ratios, pvals)
    try: xmin = min(means); xmax = max(means); ymin = min(ratios); ymax = max(ratios)
    except ValueError:
        print "\nError: non-numeric columns selected. Try to specify more suitable columns with [-c].\n"; raise

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
        spline_annotes=[]; spline_coords={}; extreme_ratios=[]
        percentiles = [1,5,25,50,75,95,99]
        for k in percentiles:
            h=[]
            for b in range(bins):
                score = stats.scoreatpercentile([p[2] for p in points_in[b]], k)
                h.append(score)
                if k == extremes: extreme_ratios.extend([p for p in points_in[b] if p[2] < score])
                if k == 100-extremes: extreme_ratios.extend([p for p in points_in[b] if p[2] > score])

            coeffs = numpy.polyfit(x, h, deg)
            smin = slimits[0] or x[0]
            smax =  slimits[1] or 0.85*x[-1]
            x_spline = numpy.array(numpy.linspace(smin, smax, 10*bins))
            y_spline = numpy.polyval(coeffs, x_spline)

            ax.plot(x_spline, y_spline, "-", color="blue")
            #ax.plot(x, h, "o", color="blue") # testing
            spline_annotes.append((k,x_spline[0],y_spline[0])) #quantile percentages
            spline_coords[k] = zip(x_spline,y_spline)

        if extremes:
            extremes_filename = "extreme_ratios_"+output_filename
            with open(extremes_filename,"w") as f:
                c = csv.writer(f,delimiter="\t")
                c.writerow(["Name","countsC1","countsC2","log10Mean","log2Fold"])
                for p in extreme_ratios:
                    c.writerow([p[0], str(counts[p[0]][0]), str(counts[p[0]][1]), str(p[1]), str(p[2])])
            print "Ratios r<"+str(extremes)+"% "+str(100-extremes)+"%<r :", extremes_filename

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

    # Annotation of points, draw
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
        figname = "maplot_"+output_filename +".png"
        fig.savefig(figname)
        print "Figure:", figname

    # Output for Javascript
    if mode == "json":
        def rgb_to_hex(rgb):
            return '#%02x%02x%02x' % rgb
        jsdata = []
        # - data points
        for data in dataset:
            gdata = zip(*groups[data])
            pvals = gdata[3]
            datapts = zip(*gdata[1:3])
            jsdata.append({"label": "Data points",
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
        transparencies = iter([0.15, 0.2, 0.35, 0.35, 0.2, 0.15])
        percentiles.pop(percentiles.index(50))
        for i in percentiles:
            jsdata.append({"id": str(i)+"%",
                           "data": spline_coords[i],
                           "lines": {"show":True, "lineWidth":0, "fill": transparencies.next()},
                           "color": rgb_to_hex((255,0,255))})

        jsdata = "var data = " + json.dumps(jsdata) + ";\n" \
                 + "var splinelabels = " + json.dumps(splinelabels) + ";\n"

        # - url for more info on features
        if assembly_id:
            assemblies = urllib.urlopen("http://bbcftools.vital-it.ch/genrep/assemblies.json").read()
            assemblies = json.loads(assemblies)
            try: assembly_id = int(assembly_id)
            except: pass
            for a in assemblies:
                if a['assembly']['name'] == assembly_id or a['assembly']['id'] == assembly_id:
                    nr_assembly_id =  a['assembly']['nr_assembly_id']
                    md5 = a['assembly']['md5']; break
            url_template = urllib.urlopen("http://bbcftools.vital-it.ch/genrep/nr_assemblies/" +
                            "get_links/" + str(nr_assembly_id) + ".json?gene_name=%3CName%3E&md5=" + md5)
            jsdata = jsdata + "var url_template = " + url_template.read() + ";"
        else: jsdata = jsdata + "var url_template = null;"
        print >> sys.stdout, jsdata

    return 0


#----------------------------------------------------------#

class AnnoteFinder:
  """
  Callback for matplotlib to display an annotation when points are clicked on.  The
  point which is closest to the click and within xtol and ytol is identified.

  Typical usage::

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
    """Distance between two points"""
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

maplot.py [-h] [-c --cols] [-l --labels] [-m --mode] [-f --format] [-s --sep]
          [-d --deg] [-b --bins] [-a --assembly] [-q --noquantiles]
          [--annotate] [--xmin --xmax --ymin --ymax] [--smin --smax]
          [-t --title] [-e --extremes]
          data_1 .. data_n

**Input**: CSV files containing at least two numeric columns representing the two
           samples to compare, and one with labels for each pair of points. Files must all
           be in the same format. By default numeric columns are 2 and 3,
           whilst column 1 contains the labels. Else, precise which columns to use with the
           `--cols` argument.

**Output**: depending on the chosen `--mode`, prints to stdout the name of the .png file produced,
            prints to stdout a json containing enough information to reconstruct the plot using
            Javascript, or produces an interactive matplotlib figure. """

help = iter([
"""The numbers or names of the two columns containing the numeric data to compare, separated by
commas. E.g. --cols 3,5.""",
"""The numbers or names of the columns containing the labels of the points, separated by
commas. The first element must contain unique labels; others will be concatenated.
E.g. --labels 1,6,7 may produce `Id | name | desc` labels.""",
"""Display mode: 'normal' for static .pgn output,
'interactive' - clic to display gene names, or
'json' - json output to stdout for Javascript web interface.""",
"""Data type: 'counts' for raw count data (default), 'rpkm' for normalized data.""",
"""The character delimiting the columns of the file. If not specified, the program tries
to detect it automatically. Use 'C^V' or '\t' for a <tab> delimiter.""",
"""Degree of the interpolant percentile splines.""",
"""Number of divisions of the x axis to calculate percentiles.""",
"""Identifier for the Genrep assembly (e.g. 'hg19' or 7) used to add more
information about features into the json output.""",
"""Don't draw quantile splines. This may improve speed and lisibility in some cases.""",
"""(In 'normal' mode) Indication of which datasets to annotate.
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
"""Title of the graph""",
"""Create an output file containing features for which ratios were outside the specified
 percentile (two-sided). For the moment, must be 1 or 5. The file is named *extreme_ratios_xxxxx* ."""
])

description = """Creates an `MA-plot` to compare transcription levels of a set of
genes (or other genomic features) in two different conditions."""

def main():
    try:
        parser = optparse.OptionParser(usage=usage, description=description)

        parser.add_option("-c", "--cols", default='2,3', help = help.next())
        parser.add_option("-l", "--labels", default='1', help = help.next())
        parser.add_option("-m", "--mode", default='normal', help = help.next())
        parser.add_option("-f", "--format", default='counts', help = help.next())
        parser.add_option("-s", "--sep", default=None, help = help.next())
        parser.add_option("-d", "--deg", default=4, type="int", help = help.next())
        parser.add_option("-b", "--bins", default=30, type="int", help = help.next())
        parser.add_option("-a", "--assembly", default=None, help = help.next())
        parser.add_option("-q", "--noquantiles", default=True, action="store_false", help = help.next())
        parser.add_option("--annotate", default=None, help = help.next())
        parser.add_option("--xmin", default=None, type="float", help = help.next())
        parser.add_option("--xmax", default=None, type="float", help = help.next())
        parser.add_option("--ymin", default=None, type="float", help = help.next())
        parser.add_option("--ymax", default=None, type="float", help = help.next())
        parser.add_option("--smin", default=None, type="float", help = help.next())
        parser.add_option("--smax", default=None, type="float", help = help.next())
        parser.add_option("-t","--title", default="MA-plot", help = help.next())
        parser.add_option("-e","--extremes", default=False, type="int", help=help.next())

        (opt, args) = parser.parse_args()
        args = [os.path.abspath(a) for a in args]

        annotate = None
        if opt.annotate:
            annotate = [eval(a) for a in opt.annotate] # 0100101 -> [0,1,0,0,1,0,1]
            assert len(args) == len(annotate), "There must be one digit per dataset in --annotate."
        limits = [None or opt.xmin, None or opt.xmax, None or opt.ymin, None or opt.ymax]
        slimits = [None or opt.smin, None or opt.smax]
        cols = opt.cols.split(",")
        labels = opt.labels.split(",")
        if len(cols) != 2:
            parser.error("--cols must be *two* integers or strings separated by commas (got %s)." % opt.cols)
        if opt.mode not in ["normal","interactive","json"]:
            parser.error("--mode must be one of 'normal','interactive', or 'json' (got %s)." % opt.mode)
        if opt.format not in ["counts","rpkm"]:
            parser.error("--format must be one of 'counts' or 'rpkm' (got %s)." % opt.format)

        MAplot(args, cols=cols, labels=labels, mode=opt.mode, data_format=opt.format, sep=opt.sep,
                limits=limits, slimits=slimits,deg=opt.deg, bins=opt.bins, assembly_id=opt.assembly,
                annotate=annotate, quantiles=opt.noquantiles, title=opt.title, extremes=opt.extremes)

        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())


