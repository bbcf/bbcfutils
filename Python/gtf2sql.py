from bbcflib.btrack import track, ensembl_to_ucsc
from bbcflib.bFlatMajor.common import split_field, map_chromosomes
from bbcflib import genrep
import os, re

def _split_attributes(stream,field='attributes'):
    fi = stream.fields.index(field)
    for x in stream:
        for y in x[fi].split(';'):
            header = re.search(r'\s*(\S+) \S*',y)
            if header:
                yield header.groups()[0].strip('"')


def _frame(x):
    if x == '.' or x is None: x='-1'
    return int(x)

def _exon_number(x):
    if x is None: x='1'
    return int(x)

basepath = "/db/genrep/nr_assemblies/gtf/"
std_outfields = ['gene_id','gene_name','transcript_id','transcript_name','exon_id','exon_number']
params = {'info': {'datatype':'relational'},
          'outtypes': {'exon_number': 'integer'}}
gtf_read_fields = ['chr','source','name','start','end','strand','frame','attributes']
sql_fields = ['chr','biotype','type','start','end','strand','frame']+std_outfields

for _a,info in genrep.GenRep().assemblies_available():
    assembly = genrep.Assembly(_a)
    gtf_path = os.path.join(basepath,assembly.md5+".gtf.gz")
    if not(assembly.bbcf_valid and os.path.exists(gtf_path)): continue
    print info, gtf_path
    gtf = track(gtf_path, intypes={'frame': _frame,
                                   'exon_number': _exon_number},
                chrmeta=assembly.chrmeta)
    all_fields = set(_split_attributes(gtf.read(fields=['attributes'])))
    new_fields = [f for f in all_fields if not(f in std_outfields)]
    xsplit = split_field(ensembl_to_ucsc(
            gtf.read(fields=gtf_read_fields)),
            outfields=std_outfields+new_fields,
            infield='attributes', header_split=' ')
    xsplit.fields[:7] = sql_fields[:7]
    outname = assembly.md5+".sql"
    if os.path.exists(outname): os.unlink(outname)
    outf = track(outname,fields=sql_fields+new_fields,
                 chrmeta=assembly.chrmeta,**params)
    outf.write(map_chromosomes(xsplit,assembly.chromosomes))
