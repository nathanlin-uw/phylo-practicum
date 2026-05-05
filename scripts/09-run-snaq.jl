## Use multithreaded Julia session (julia -t 2)

using Distributed
addprocs(4)

@everywhere using PhyloNetworks
@everywhere using SNaQ

## read table of CF
d_sp = readtableCF("../results/09-tableCF_species.csv"); # "DataCF" object for use in snaq!
#read in the species tree from ASTRAL as a starting point
T_sp = readnewick("../results/07-species-tree-astral4.tre")

net = snaq!(T_sp, d_sp, runs=100, Nfail=200, filename= "../results/snaq/09-snaq-h1",seed=8485);

## Plot now``
using PhyloPlots
# net = readnewick("(Ae_sharonensis,Ae_longissima,(Ae_bicornis,(Ae_searsii,((Ae_tauschii,(((Ae_uniaristata,Ae_comosa)1:0.4918206502664954,(Ae_caudata,Ae_umbellulata)1:0.13449338165653227)1:0.00911821493436927,((T_boeoticum,T_urartu)1:1.6460105085783057,(H_vulgare,((Ae_speltoides,Ae_mutica)1:0.07124470208266999)#H26:0.14159810198824307::0.7563186990617421)1:0.1869824746969994)1:0.46325640725730144)0.99651:0.06048455134529327):0.16788568790821257,#H26:0.0::0.2436813009382579):0.45932787904385297)1:0.9296082436533977)1:0.5926597507029276)1;")
plot(net, showedgenumber=true)

## Root on outgroup
rootonedge!(net, 16)
rotate!(net,22)
rotate!(net,23)
rotate!(net,-6)
rotate!(net,12)
rotate!(net,11)
plot(net, showgamma=true)


## Recolor to mimic Figure 5
using DataFrames

tipnodes = [n.number for n in net.node if n.leaf]
tipnames = [n.name for n in net.node if n.leaf]

tipcolors = Dict(
    "T_urartu" => "darkolivegreen",
    "T_boeoticum" => "darkolivegreen",
    "Ae_comosa" => "chocolate",
    "Ae_uniaristata" => "chocolate",
    "Ae_caudata" => "khaki",
    "Ae_umbellulata" => "gold",
    "Ae_tauschii" => "red",
    "Ae_longissima" => "mediumorchid",
    "Ae_sharonensis" => "mediumorchid",
    "Ae_bicornis" => "mediumorchid",
    "Ae_searsii" => "mediumorchid",
    "Ae_mutica" => "dodgerblue",
    "Ae_speltoides" => "navy"
)

colors = [get(tipcolors, name, "black") for name in tipnames]

nodelabel = DataFrame(
    number = tipnodes,
    label = tipnames,
    nodelabelcolor = colors
)

plot(net, nodelabel = nodelabel)