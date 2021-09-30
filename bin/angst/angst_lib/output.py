#!/usr/bin/python
# AnGST

import os
import pdb
import sys
import time
import math
import copy
import shutil

from glob import glob
import output as output

def MakeDirectories(output_dir):
  '''make a new directory for output data.  overwrite if prior directory
  exists '''

  # does write directory exist?
  angst_results_dir = output_dir + "/"
  results_dirs = {}
  results_dirs["angst"] = angst_results_dir
  if os.path.isdir(angst_results_dir):
    # now wipe write directory
    angst_results_dir = angst_results_dir.rstrip("/")
    shutil.rmtree(angst_results_dir)

  os.mkdir(angst_results_dir)
  return results_dirs

def SaveResults(output_args):
  ''' write the results of AnGSTs analysis to file '''

  model = output_args["model"]
  results_dirs = output_args["results"]
  start_time = output_args["start_time"]
  counts = output_args["counts"]
  mem_str = output_args["mem_str"]
  gene_tree = output_args["gene_tree"]
  gene_num = counts["gene"]
  spec_num = counts["spec"]
  node_link = output_args["node_link"]
  #scaling_list = output_args["scaling"]
  scaling_list = []

  hgt_penalty = model.hgt_penalty
  los_penalty = model.los_penalty
  dup_penalty = model.dup_penalty
  spc_penalty = model.spc_penalty
  angst_results_dir = results_dirs["angst"]
  species_tree = model.species_tree
  luca = species_tree.root
  gene_root = gene_tree.root

  # find best rooting
  best_dict = output_args["res_dict"]
  best_node = None
  best_newick = None
  best_events = None

  if model.luca == "Root" or model.special == "hr":
    best_score = best_dict["scores"][luca]
    best_node = luca
  else:
    best_score = min(best_dict["scores"].values())
    for node in best_dict["scores"]:
      if best_dict["scores"][node] == best_score:
        best_node = node
  best_newick = best_dict["newick"][best_node]
  best_events = best_dict["events"][best_node]

  # dump best rootings to output files
  output_newick = open(angst_results_dir + "AnGST.newick",'w')
  output_events = open(angst_results_dir + "AnGST.events",'w')
  output_score = open(angst_results_dir + "AnGST.score",'w')
  output_leaf = open(angst_results_dir + "AnGST.leaf",'w')
  output_counts = open(angst_results_dir + "AnGST.counts",'w')

  # give root branch a "meaningless" length so that rooted gene_tree
  # formatted correctly
  output_newick.write("(" + best_newick + ":0.01337);")
  output_score.write(str(best_score) + "\n")

  # write out a leaf-centric list of events
  LeafBasedEvents(gene_root,best_node)
  PrintLeafEvents(gene_root,output_leaf)

  # print out all events
  output_events.write("[brn]: " + str(best_node) + "\n")
  event_counts = PrintAllEvents(gene_root,best_node,output_events)

  # do counts
  CountGenes(gene_root,species_tree,output_counts)

  # look for paradoxes
  #output_paradox = open(angst_results_dir + "AnGST.paradox",'w')
  #best_paradox = best_dict['paradox'][best_node]
  #paradox_count = ParadoxProcess(best_paradox,output_paradox)
  #output_paradox.close()


  # get branch distances
  #output_branches = open(angst_results_dir + "AnGST.branches",'w')
  #path_dict = {}
  #GetBranchDists(gene_root,best_node,None,path_dict,output_branches)
  #output_branches.close()

  # get distances
  #output_dists = open(angst_results_dir + "AnGST.dists",'w')
  #output_dists.close()

  output_newick.close()
  output_counts.close()
  output_score.close()
  output_events.close()
  output_leaf.close()

  # save information on transfer distances
  #output_hgt_dists = open(angst_results_dir + "AnGST.hgts",'w')
  #PrintHGTDists(gene_root,best_node,species_tree,output_hgt_dists)
  #output_hgt_dists.close()

  # how many possible hgt were there?
  possible_hgt = len(species_tree.possible_hgt)

  # save run statistics and notes
  output_stats = open(angst_results_dir + "AnGST.stats",'w')
  output_stats.write("\n# angst version info #\n")
  output_stats.write("angst v\n")

  output_stats.write("\n# input data #\n")
  output_stats.write("species file:\t")
  output_stats.write(output_args["species_file"] + "\n")
  output_stats.write("species count:\t" + str(spec_num) + "\n")
  output_stats.write("bootstrap file:\t")
  output_stats.write(output_args["boot_file"] + "\n")
  output_stats.write("gene count:\t" + str(gene_num) + "\n")
  output_stats.write("bootstraps:\t" + str(counts["boots"]) + "\n")
  output_stats.write("possible hgt:\t" + str(possible_hgt) + "\n")

  output_stats.write("\n# model settings #\n")
  output_stats.write("event guide:\t")
  output_stats.write(str(output_args["event_guide"]) + "\n")
  output_stats.write("hgt penalty: " + str(hgt_penalty) + "\n")
  output_stats.write("dup penalty: " + str(dup_penalty) + "\n")
  output_stats.write("los penalty: " + str(los_penalty) + "\n")
  output_stats.write("spc penalty: " + str(spc_penalty) + "\n")
  #output_stats.write("cog scaling: " + str(node_link.cog_scaling) + "\n")
  #output_stats.write("hr scaling: " + str(node_link.hr_scaling) + "\n")

  output_stats.write("\n# inference counts #\n")
  output_stats.write("hgt inferred: " + str(event_counts['hgt']) + "\n")
  output_stats.write("dup inferred: " + str(event_counts['dup']) + "\n")
  output_stats.write("los inferred: " + str(event_counts['los']) + "\n")
  output_stats.write("spc inferred: " + str(event_counts['spc']) + "\n")
  #output_stats.write("time paradoxes: " + str(paradox_count) + "\n")

  output_stats.write("\n# running time #\n")
  output_stats.write("run-time: " + str(time.time()-start_time) + " sec\n")

  output_stats.write("\n# memory usage #\n")
  output_stats.write(mem_str)

  output_stats.close()

  return best_score

def PrintEventList(events,out_list):
  '''print out events.  recurse if the list actually has another list
  in it'''
  if type(events) == type(""):
    if len(events) > 0:
      if events.count('[') < 2:
        out_list.append(events)
      else:
        event_parts = events.split('[')
        for event_ind in xrange(1,len(event_parts)):
          out_list.append('[' + event_parts[event_ind])
  else:
    for event in events:
      PrintEventList(event,out_list)
  return out_list


def LeafBasedEvents(this_node,lca):
  ''' print out events in the history of a particular leaf '''
  event_list = this_node.lca_lookups['events'][lca]
  this_node.event_list = PrintEventList(event_list,[])

  if this_node.raw_list is not None:
    print "there should be no raw list here already"
    sys.exit(1)
  this_node.raw_list = copy.copy(this_node.event_list)

  # get the events in your kids
  kids = this_node.kid_nodes.keys()
  for kid in kids:
    kid_lca = this_node.lca_map[lca][kid]
    LeafBasedEvents(kid,kid_lca)
    kid.rec_lca = kid_lca
    # remove kids events from all ancestors
    for event in kid.event_list:
      ParentEventRemove(this_node,event)


def PrintDistances(this_node,lca,dist_file):
  ''' print out distance ratios for all nodes on gene tree '''
  kids = this_node.kid_nodes.keys()
  if len(kids) < 1:
    return

  # break out the distances
  dist_diff = this_node.dist_ratios[lca][2]
  event_str = this_node.dist_ratios[lca][3]
  dist_file.write(str(dist_diff) + "\t" + event_str + "\n")
  for kid in kids:
    PrintDistances(kid,this_node.lca_map[lca][kid],dist_file)
  return

def PrintLeafEvents(this_node,output_file):

  if this_node.leaf_event_list is None:
    this_node.leaf_event_list = []

  # pass to kids the events at this node
  kids = this_node.kid_nodes.keys()
  for kid in kids:
    kid.leaf_event_list = []
    kid.leaf_event_list.extend(this_node.leaf_event_list)
    kid_lca = kid.rec_lca

    # handle events happening at the leaves as well
    all_events = []
    all_events.extend(this_node.event_list)
    if len(kid.branch_list) < 2:
      for kid_event in kid.event_list:
        if 'cur' not in kid_event:
          all_events.append(kid_event)

    for event in all_events:
      if event.count('hgt') > 0:
        spec_spec = event.split(' --> ')[1]
      else:
        if len(event.split(': ')) < 2:
          pdb.set_trace()
        spec_spec = event.split(': ')[1]
      spec_spec = "-" + spec_spec + "-"
      # is this species in there?
      this_spec = "-" + kid_lca.species + "-"
      if spec_spec.count(this_spec) > 0:
        kid.leaf_event_list.append(event)
    PrintLeafEvents(kid,output_file)

  # now, if you're at a leaf, print events
  if len(this_node.branch_list) < 2:
    output_file.write("leaf " + this_node.name + ":\n")
    for event in this_node.leaf_event_list:
      output_file.write("\t" + event + "\n")
    output_file.write("\t" + this_node.event_list[0] + "\n")
    output_file.write("\n")


def ParentEventRemove(this_node,event):
  this_node.event_list.remove(event)
  if this_node is this_node.tree.root:
    return
  p_ends = this_node.parent_branch.ends
  parent = filter(lambda i: i is not this_node,p_ends)[0]
  ParentEventRemove(parent,event)


def PrintHGTDists(this_node,lca,tree,output_file):
  events = this_node.raw_list
  for event in events:
    if 'hgt' in event:
      this_str = event.lstrip('[hgt]: ')
      parts = this_str.split(' --> ')
      donor = tree.node_dict[parts[0]]
      accep = tree.node_dict[parts[1]]
      n1n2_dist = tree.dist_dict[donor][accep]

      output_str = this_str
      output_str += "\t" + str(n1n2_dist)
      output_str += "\t" + str(donor.DistTo(tree.root)[1])
      output_str += "\t" + str(accep.DistTo(tree.root)[1])
      output_str += "\t" + str(tree.height_dict[donor])
      output_str += "\t" + str(tree.height_dict[accep])
      output_file.write(output_str + "\n")


def PrintAllEvents(this_node,lca,output_file):
  event_counts = {'hgt':0, 'dup':0, 'los':0, 'spc':0}
  events = this_node.raw_list
  for event in events:
    output_file.write(event + "\n")
    if event.count('hgt') > 0:
      event_counts['hgt'] += 1
    elif event.count('dup') > 0:
      event_counts['dup'] += 1
    elif event.count('los') > 0:
      event_counts['los'] += 1
    elif event.count('spc') > 0:
      event_counts['spc'] += 1

  return event_counts


def GetLeafList(this_node,leaf_list=[]):
  kids = this_node.kid_nodes.keys()
  if len(kids) < 1:
    leaf_list.append(this_node)
  else:
    for kid in kids:
      GetLeafList(kid,leaf_list)
  return leaf_list


def CountGenes(this_node,species_tree,output_file):
  all_species = species_tree.node_dict.keys()
  species_count_dict = {}
  for species in all_species:
    species_count_dict[species] = 0

  event_list = this_node.raw_list
  for event in event_list:
    if 'cur' in event or 'spc' in event:
      host = event.split(': ')[1]
      species_count_dict[host] += 1

  for species in species_count_dict:
    this_count = species_count_dict[species]
    output_file.write(species + ": " + str(this_count) + "\n")

  return


def ParadoxProcess(paradox_list,paradox_file):
  paradox_count = 0
  for paradox in paradox_list:
    paradox_count += 1
    paradox_file.write(paradox[0] + "\n")
    for ind in range(1,len(paradox[1])):
      paradox_file.write("\t" + paradox[1][ind] + "\n")
  return paradox_count


def GetBranchDists(this_node,this_lca,parent_lca,path_dict,output_f):
  ''' iterate through all nodes on the tree '''
  for kid in this_node.kid_nodes:
    kid_lca = this_node.lca_map[this_lca][kid]
    if kid_lca is not this_lca:
      GetPaths(this_node,kid,this_lca,kid_lca,path_dict,output_f)
    GetBranchDists(kid,kid_lca,this_lca,path_dict,output_f)


def GetPaths(root_node,this_node,root_lca,this_lca,path_dict,output_f):
  last_stop = True
  if this_lca.AreRelated(root_lca):
    for kid in this_node.kid_nodes:
      kid_lca = this_node.lca_map[this_lca][kid]
      if kid_lca is this_lca:
        last_stop = False
        path_dict[kid] = 1
        GetPaths(root_node,kid,root_lca,kid_lca,path_dict,output_f)
  if last_stop == True:
    gene_dist = root_node.DistTo(this_node)[1]
    # hgt
    label_str = "foo"
    if root_lca is this_lca:
      spec_dist = this_lca.parent_branch.length
      label_str = "dup"
    # speciation
    elif root_lca.AreRelated(this_lca):
      spec_dist = root_lca.DistTo(this_lca)[1]
      label_str = "ver"
    # hgt
    else:
      spec_dist = this_lca.parent_branch.length / 2
      label_str = "hgt"

    output_str = str(root_lca) + "\t"
    output_str += str(this_lca) + "\t"
    output_str += str(spec_dist) + "\t"
    output_str += str(gene_dist) + "\t"
    output_str += label_str + "\n"
    output_f.write(output_str)
