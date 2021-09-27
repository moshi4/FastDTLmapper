# %%

from Bio import Phylo, SeqIO
from Bio.Phylo.BaseTree import Tree

fasta_file = "./example/fasta/MAlvi.faa"
genbank_file = "./example/genbank/MAlvi.gbk"
nwk_file = "./example/Mycoplasma.nwk"
txt_file = "./example/FastDTLmapper_cmd_example.txt"
unrooted_nwk_file = "./example/test_unrooted.nwk"


# SeqIO.read(txt_file, "fasta")
# SeqIO.read(txt_file, "genbank")
# SeqIO.read(fasta_file, "fasta")
# tree = Phylo.read(txt_file, "newick")
# tree = Phylo.read(nwk_file, "newick")
# len(list(tree.find_clades()))

tree: Tree = Phylo.read(unrooted_nwk_file, "newick")
tree.root.is_bifurcating()


# %%
