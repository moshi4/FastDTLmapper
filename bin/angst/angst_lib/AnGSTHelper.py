#!/usr/bin/python
# AnGST

# python libraries
import math
import sys
import time
import pdb
import PyVM
import random
import os

current_path = os.path.dirname(os.path.realpath(__file__))
tree_lib_dir = os.path.join(current_path, '..', 'tree_lib')
tree_lib_dir = os.path.normpath(tree_lib_dir)
sys.path.append(tree_lib_dir)
import multitree as multitree
import node as node
import reconcile as reconcile
import output as output
import AnGSTInput as input

def RunAnGST(input_dict):

    # unpack input variables
    input_info = input_dict['input_info']
    mem_str = input_dict['mem_str']
    write_out = input_dict['write_out']

    # build the reconciliation model #
    model = reconcile.dlt_model(input_info)

    # build the species tree #
    print "* build species tree"
    species_tree_file = open(input_info.species_tree_filename,"r")
    species_tree_string = species_tree_file.read().strip()
    species_tree = multitree.multitree()
    species_tree.Build(species_tree_string)
    species_tree_file.close()
    model.species_tree = species_tree
    mem_str.append(PyVM.MemoryUpdate("build species tree",'return'))

    # make sure no delimiters in the species tree
    for node in species_tree.node_dict:
        if node.count('.') > 0 or node.count('_') > 0:
            print "no delimiters in species tree please."
            sys.exit(1)

    # get all the pairwise distances on the species tree
    species_tree.root.FindDists()
    species_tree.root.GetHeights()

    # get the time constraints
    if input_info.ultra_tree_bool is True:
        ultra_tree_file = open(input_info.species_tree_filename,"r")
        ultra_tree_string = ultra_tree_file.read().strip()
        ultra_tree = multitree.multitree()
        ultra_tree.Build(ultra_tree_string)
        ultra_tree_file.close()
        mem_str.append(PyVM.MemoryUpdate("build ultra tree",'return'))

        ultra_tree.model = model
        ultra_tree.root.FindDists()
        ultra_tree.root.GetHeights()

        print "foo!"

        # link the species and ultrametric tree nodess
        tree_link = species_tree.TreeLink(ultra_tree)
    else:
        ultra_tree = None
        tree_link = None

    # build the gene trees
    print "* build gene trees"
    bootstraps = open(input_info.boots_file,'r').readlines()
    boot_list = []
    for boot in bootstraps:
        boot_list.append(multitree.multitree())
        boot_list[-1].Build(boot)
    mem_str.append(PyVM.MemoryUpdate("build bootstraps",'return'))

    # is there a single rooted bootstrap?
    if len(boot_list) == 1 and boot_list[0].root is not None:
        gene_tree = boot_list[0]
        tree = gene_tree
    else:
        # you've got bootstraps

        # check that all of the trees have the same leaves
        leaf_set = boot_list[0].MakeLeafSet()
        for ind in range(1,len(boot_list)):
            boot_list[ind].TestLeafSet(leaf_set)

        # if doing HR, make sure nodes map one-to-one between species tree
        # and gene tree
        if model.special == 'hr':
            gene_leaves = [leaf.species for leaf in boot_list[0].leaf_dict.values()]
            species_leaves = species_tree.root.leaves.keys()
            for leaf in gene_leaves:
                if leaf in species_leaves:
                    species_leaves.remove(leaf)
                else:
                    print "leaf " + leaf + " in gene tree, but not species tree."
                    print "un-matched nodes; can't do HR"
                    sys.exit(1)
            if len(species_leaves) > 0:
                print "un-matched nodes; can't do HR"
                sys.exit(1)

        # are the bootstraps rooted or unrooted?
        rooted_c = 0
        for boot_tree in boot_list:
            if boot_tree.root is not None:
                rooted_c += 1
            #endif
        #endfor

        if rooted_c == len(boot_list):
            ''' all boots rooted '''

            # reconcile each of the rooted bootstraps
            overall_best_score = sys.maxint
            overall_best_tree = None
            random.shuffle(boot_list)
            for boot_tree in boot_list:
                node_link = reconcile.node_link(species_tree,model,ultra_tree,tree_link)
                node_link.hr_scaling = input_info.hr_scaling
                node_link.learn_events = 1
                rec_res = reconcile.LocalReconcile(boot_tree.root,None,node_link)
                boot_tree.root.lca_lookups = rec_res
                best_score = min(rec_res['scores'].values())
                if best_score < overall_best_score:
                    overall_best_score = best_score
                    overall_best_tree = boot_tree

            # select one that gives best reconciliation
            gene_tree = overall_best_tree

        elif rooted_c == 0:
            ''' all boots unrooted '''

            # link the nodes
            print "* linking nodes"
            node_link = reconcile.node_link(species_tree,model,ultra_tree,tree_link)
            node_link.hr_scaling = input_info.hr_scaling
            for tree in boot_list:
                node_link_dict = {}
                tree.a_node.GetNodeLinkDict(node_link_dict)
                node_link.MergeNodeLinkDicts(node_link_dict)
                del node_link_dict
            mem_str.append(PyVM.MemoryUpdate("link nodes",'return'))

            # optimize ALL leaves
            print "* initial reconciliation"
            for l in node_link.link_dict.keys():
                if model.special == 'hr':
                    node_link.learn_events = 1
                    model.dup_penalty = math.exp(500)
                elif model.event_guide is not None:
                    node_link.learn_events = 1
                else:
                    node_link.learn_events = 0
                node_link.GlobalReconcile(l)
            mem_str.append(PyVM.MemoryUpdate("reconcile",'return'))

            # evaluate all roots
            print "* evaluate roots"
            for leaf_str in node_link.link_dict:
                node_pair = node_link.link_dict[leaf_str][0]
                # only care about one root if you're doing HR
                if model.outgroup is not None:
                    if node_pair[0].name == model.outgroup:
                        node_link.TryRoot(node_pair)
                    if node_pair[1].name == model.outgroup:
                        node_link.TryRoot(node_pair)
                else:
                    node_link.TryRoot(node_pair)

            mem_str.append(PyVM.MemoryUpdate("find root",'return'))
            best_score = node_link.best_res['score']
            best_newick = "(" + node_link.best_res['newick'] + ":0.01337);"

            gene_tree = multitree.multitree()
            gene_tree.Build(best_newick)

        else:
            print "boot root mismatch?!"
            sys.exit(1)
        #endif

    if write_out:
        # perform a final reconciliation on best gene tree to get all data
        print "* final reconciliation"

        species_tree.rec_dict = {}
        species_tree.los_dict = {}

        mem_str.append(PyVM.MemoryUpdate("rebuild gene tree",'return'))
        node_link = reconcile.node_link(species_tree,model,ultra_tree,tree_link)
        node_link.hr_scaling = input_info.hr_scaling
        node_link.learn_events = 2      # 2 == full learning
        rec_res = reconcile.LocalReconcile(gene_tree.root,None,node_link)
        gene_tree.root.lca_lookups = rec_res
        mem_str.append(PyVM.MemoryUpdate("last reconciliation",'return'))
        best_score = min(rec_res['scores'].values())

        # establish output directories
        results_dirs = output.MakeDirectories(input_info.output_dir)

        # collect the data that's going to be printed out
        print "* write output"
        counts = {}
        gene_c = None
        gene_c = len(str(gene_tree.root).split('-'))
        counts["gene"] = gene_c
        counts["spec"] = (len(species_tree.node_dict) + 1)/2
        counts["boots"] = len(bootstraps)

        output_args = {}
        output_args["node_link"] = node_link
        output_args["results"] = results_dirs
        output_args["res_dict"] = rec_res
        output_args["start_time"] = input_dict['start_time']
        output_args["gene_tree"] = gene_tree
        output_args["model"] = model
        output_args["counts"] = counts
        output_args["species_file"] = input_info.species_tree_filename
        output_args["boot_file"] = input_info.boots_file
        output_args["mem_str"] = ''.join(mem_str)
        output_args["event_guide"] = input_info.event_guide_fn
        output.SaveResults(output_args)

    print best_score
    return best_score




