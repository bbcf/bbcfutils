#!/usr/bin/env julia

using MultivariateStats
using DataFrames
using Gadfly


if length(ARGS) > 0; filename = ARGS[1]
else filename = "genes_expression.txt"; end
if length(ARGS) > 1; outprefix = ARGS[2]
else outprefix = "pca_rnaseq"; end
@assert isfile(filename) "File $filename not found"


# PCA itself
df = readtable(filename, header=true, separator='\t')
header = map(string,names(df))
rpkm_cols = filter(i->ismatch(r"^rpkm.*",header[i]), 1:length(header))
rpkm_names = [h[6:end] for h in header[rpkm_cols]]
rpkm = convert(Array,DataArray(df[rpkm_cols]))
trpkm = log(rpkm.+0.5)
p = fit(PCA, trpkm, maxoutdim=10)
proj = trpkm' * p.proj


# Biplot
pl = plot(x=proj[:,1],y=proj[:,2], Geom.point,
          Guide.xlabel("PC1"), Guide.ylabel("PC2"),
          Guide.title("PCA biplot of samples ~ gene expression (log RPKM)"),
          label=rpkm_names, Geom.label(hide_overlaps=false) )
draw(PDF("$(outprefix).pdf", 6inch, 6inch), pl)


# Barplot of variance proportions
vars = p.prinvars / p.tvar
cumvar = string(round(sum(vars[1:2])*100.,2))
barplot = plot(x=1:length(vars), y=vars, Geom.bar(position=:dodge),
               Theme(bar_spacing=0.5cm), Scale.x_discrete,
               Guide.title("PC1+PC2: $cumvar %"),
               Guide.xlabel("PC#"), Guide.ylabel("Proportion of variance") )
draw(PDF("$(outprefix)_sdev.pdf", 6inch, 6inch), barplot)


