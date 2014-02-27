#!/usr/bin/env julia

using DimensionalityReduction
using DataFrames


if length(ARGS) > 0; filename = ARGS[1]
else filename = "genes_expression.tab"; end
if length(ARGS) > 1; outprefix = ARGS[2]
else outprefix = "pca_rnaseq"; end
@assert isfile(filename) "File $filename not found"


# PCA itself
df = readtable(filename, header=true, separator='\t')
header = map(string,names(df))
rpkm_cols = filter(i->ismatch(r"^rpkm.*",header[i]), 1:length(header))
rpkm_names = [h[6:end] for h in header[rpkm_cols]]
rpkm = convert(Array,DataArray(df[rpkm_cols]))
trpkm = transpose(log(rpkm.+0.5))
p = pca(trpkm, scale=false)


# Biplot
using Gadfly
pl = plot(x=p.scores[:,1],y=p.scores[:,2], Geom.point,
          Guide.xlabel("PC1"), Guide.ylabel("PC2"),
          Guide.title("PCA biplot of samples ~ gene expression (log RPKM)"),
          label=rpkm_names, Geom.label(hide_overlaps=false) )
draw(PNG("$(outprefix).png", 6inch, 6inch), pl)
draw(D3("$(outprefix).js", 6inch, 6inch), pl)


# Barplot of variance proportions
cumvar = string(round(sum(p.cumulative_variance[1:2])*100.,2))
barplot = plot(x=1:10, y=p.proportion_of_variance[1:10], Geom.bar(position=:dodge),
               Theme(bar_spacing=0.5cm), Scale.x_discrete,
               label=["PC1+PC2: $cumvar %"], Geom.label,
               Guide.xlabel("PC#"), Guide.ylabel("Proportion of variance") )
draw(PNG("$(outprefix)_sdev.png", 6inch, 6inch), barplot)
draw(D3("$(outprefix)_sdev.js", 6inch, 6inch), barplot)


