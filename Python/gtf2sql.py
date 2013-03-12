from bbcflib.btrack import track, ensembl_to_ucsc, FeatureStream
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


def _exon_number(x):
    if x is None: x='1'
    return int(x)

def _fix_exon_id(stream,count):
    for x in stream:
        type = x[2]
        exon_id = x[12]
        if type == "exon" and not exon_id:
            exon_id = "gr_exid_%i"%count
            count += 1
            x = x[:12]+(exon_id,)+x[13:]
        yield x

_intypes = {'exon_number': _exon_number}
basepath = "/data/epfl/bbcf/genrep/nr_assemblies/%s"
std_outfields = ['gene_id','gene_name','transcript_id','transcript_name','exon_id','exon_number']
params = {'info': {'datatype':'relational'},
          'outtypes': {'exon_number': 'integer'}}
gtf_read_fields = ['chr','source','name','start','end','strand','frame','attributes']
sql_fields = ['chr','biotype','type','start','end','strand','frame']+std_outfields

for _a,info in genrep.GenRep().assemblies_available():
    for n in range(100):
        assembly = genrep.Assembly(_a)
        gtf_path = os.path.join(basepath%"gtf","%s_%i.gtf.gz"%(assembly.md5,n))
        if not(assembly.bbcf_valid and os.path.exists(gtf_path)): break
        sql_path = os.path.join(basepath%"annot_tracks","%s_%i.sql"%(assembly.md5,n))
        if os.path.exists(sql_path): continue
        print info, gtf_path, sql_path
        gtf = track(gtf_path, intypes=_intypes, chrmeta=assembly.chrmeta)
        all_fields = set(_split_attributes(gtf.read(fields=['attributes'])))
        new_fields = [f for f in all_fields if not(f in std_outfields)]
        xsplit = split_field(ensembl_to_ucsc(gtf.read(fields=gtf_read_fields)),
                             outfields=std_outfields+new_fields, infield='attributes', 
                             header_split=' ', strip_input=True)
        xsplit.fields[:7] = sql_fields[:7]
        exon_count = 1
        outf = track(sql_path,fields=sql_fields+new_fields,chrmeta=assembly.chrmeta,**params)
        outf.write(FeatureStream(_fix_exon_id(map_chromosomes(xsplit,assembly.chromosomes),exon_count),
                                 xsplit.fields))
        outf.close()
