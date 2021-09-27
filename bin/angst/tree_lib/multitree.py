# object class for construction of rooted and unrooted binary trees.
#
# lawrence david - ldavid@mit.edu - 2006.

import re
import pdb
import sys
import copy

import node
import branch

class multitree:

    # recursive module for building unrooted trees from newick-notated
    # strings.
    def __init__(self):
	self.build_queue = []
        self.a_node = None
	self.a_branch = None
	self.root = None
        self.model = None
	self.root_branch = None
	self.node_dict = {}
        self.leaf_dict = {}
        self.rec_dict = {}
        self.los_dict = {}
        self.unrooted = True
        self.dist_dict = {}
        self.height_dict = {}
        self.time_dict = {}
        self.possible_hgt = {}
        self.branch_list = []

    def __repr__(self):
        ''' print the tree (intended only for rooted trees '''
	if self.root == None:
	    return "not rooted"
        elif len(self.root.leaves) < 1:
            root = self.root
            newick_s = str(root) + ":" + str(root.parent_branch.length)
            return newick_s
	else:
	    newickString = ""
	    return self.root.treePrint(newickString)


    def Build(self,newick):
        ''' initiate tree object construction from newick string '''

        # remove any bootstraps from newick
        p = re.compile('\)[^:]+:')
        newick = p.sub('):',newick)
        self.RecursiveBuild(newick)


    def RecursiveBuild(self,newick):
        ''' recursively populate tree '''

	# handle the case when you reach a trifurcating node (so
	# unrooted), by looking for when you've got a stretch of 2
	# colons unseparated by a parenthesis
        colon_count = newick.count(':')
	colon_match2 = re.findall(':[^\)]*:[^\)]*:',newick)


        if colon_count == 0:
            # tree is empty.  abdicate.
            print "no nodes in tree"
            return

	# if you reach a root or tree only has 1 node
        elif colon_count == 1:
            root_dist = float(newick.split(':')[1].split(')')[0])
            root_branch = branch.branch(root_dist)

            # handle case of 1 node tree
            if len (self.build_queue) < 1:
                # 1 node tree
                open_paren_i = newick.find('(')
                close_paren_i = newick.find(')')
                node_s = newick[open_paren_i+1:close_paren_i]
                #node_s = newick[1:-2]
                root_node = self.BuildNode(node_s)
            else:
                # multinode tree
                root_node = self.build_queue[0]
                for kid_branches in root_node.branch_list:
                    root_node.child_branches.append(kid_branches)

            root_node.parent_branch = root_branch
	    self.node_dict[root_node.species] = root_node

            # do rooted stuff
            self.unrooted = False
	    self.root = root_node
	    self.root.imposeHierarchy()
            self.labelSubtrees()
            self.Get_Subnodes()

        # handle case of a tree with only 2 nodes:
        elif colon_count == 2:

	    # find an innermost set of parentheses:
            regex = re.search('\([^\(\)]*\)',newick).span()
	    clade = newick[regex[0]+1:regex[1]-1]

	    # split up the clade
	    children = re.split(',',clade)
	    son = children[0]
	    daughter = children[1]

            # build node1
            node1 = self.BuildNode(son)
            node2 = self.BuildNode(daughter)

	    # unite the two nodes
	    center_node = node1.unite(node2)

            # create a fake node:
            dummy_node = self.BuildNode("dummy_node:0.1")
            dummy_node.branch_list = []
            dummy_node.addBranch(0.01)
            dummy_node.myBranch.addNode(center_node)
            self.a_branch = dummy_node.myBranch
            self.a_node = dummy_node
            center_node.UnrootedLeaving()

	# handle the trifurcating node of an unrooted tree
	elif colon_count == 3 and len(colon_match2) > 0:

            # note that this updated block can handle trees with leaf
            # count >= 3

            # pull out the first node
            regex = re.search('\([^,]*,',newick).span()
            first_leaf = newick[regex[0]+1:regex[1]-1]
            newick = re.sub(first_leaf + ',','',newick)
            first_node = self.BuildNode(first_leaf)

            # first_node.branch_list = []
            # don't need to worry about clearing, since we'll be
            # adding branches anyway later down in this block.
            self.build_queue.append(first_node)

	    # join the last 2 nodes:
            regex = re.search('\(.*,[^\)]*\);',newick).span()
            clade = newick[regex[0]+1:regex[1]-2]

            # split up the clade
            children = re.split(',',clade)
            son = children[0]
            daughter = children[1]

            # build two nodes and join
            node1 = self.BuildNode(son)
            node2 = self.BuildNode(daughter)
            center_node = node1.unite(node2)

            # finally, join this new node with the remaining node:
            last_node = self.build_queue[0]
            regex = re.search(':[^,]*,',newick).span()
            last_node.myBranch.addNode(center_node)

            self.a_branch = last_node.myBranch
            self.a_node = last_node

            # do some unrooted stuff, making sure you don't begin on a
            # leaf.
            if len(self.a_branch.ends[0].branch_list) == 3:
                self.a_branch.ends[0].UnrootedLeaving()
            else:
                self.a_branch.ends[1].UnrootedLeaving()

	else:

	    # find an innermost set of parentheses:
	    regex = re.search('\([^\(\)]*\):',newick).span()
	    clade = newick[regex[0]+1:regex[1]-2]
	    # split up the clade
	    children = re.split(',',clade)
	    son = children[0]
	    daughter = children[1]

            # build nodes and unite
            node1 = self.BuildNode(son)
            node2 = self.BuildNode(daughter)
	    center_node = node1.unite(node2)

	    # replace the matched string with the new node:
	    new_newick = newick[0:regex[0]]
	    new_newick += center_node.name
	    new_newick += newick[regex[1]-1:]

	    # add the new node to the queue of extant nodes
	    self.build_queue.append(center_node)

	    # recurse
	    self.RecursiveBuild(new_newick)


    def BuildNode(this_tree,node_string):
        ''' create a new node object from newick string '''
        splitSon = re.split(':',node_string)

        # see if anyone in the queue has the same name
        node1match = False
        for i in range(len(this_tree.build_queue)):
            if this_tree.build_queue[i].name == splitSon[0]:
                node1 = this_tree.build_queue[i]
                this_tree.build_queue.remove(node1)
                node1match = True
                break
        if node1match is not True:
            node1 = node.node(splitSon[0],this_tree)
        node1.addBranch(splitSon[1])
        return node1


    def MakeLeafSet(this_tree):
        ''' make a set of leaves in a tree '''
        leaf_set = set()
        for leaf in this_tree.leaf_dict:
            leaf_set.add(leaf)
        return leaf_set


    def TestLeafSet(this_tree,test_leaf_set):
        ''' test if leaves in tree match reference set '''
        this_leaf_set = set()
        for leaf in this_tree.leaf_dict:
            this_leaf_set.add(leaf)
        # assuming that this is as fast as going through and looking
        # through each item individually
        if test_leaf_set != this_leaf_set:
            print "the leaves in the bootstraps don't all match."
            print "aborting prematurely."
            sys.exit(1)


    def labelSubtrees(self):
        ''' label the subtrees '''
	if self.root == None:
	    print "you're trying to label an unrooted tree"
	else:
	    self.root.subtreeLabel()


    def Get_Subnodes(self):
        ''' help setup datastructures for finding LCAs '''
	self.root.Find_Subnodes()


    def GetSibDists(self,species_tree):
        '''get distances between all siblings on the tree'''
        parent_dict = {}
        dist_list = []
        for leaf in self.leaf_dict.values():

            # get sibling (and make sure to not double count)
            other_nodes = leaf.branch_list[0].ends
            parent = filter(lambda i: i is not leaf,other_nodes)[0]
            if parent in parent_dict:
                continue
            else:
                parent_dict[parent] = 1
            p_other_nodes = parent.other_nodes
            other_leaves = filter(lambda i: len(i.branch_list) < 2,p_other_nodes)

            # handle caterpiller case
            if len(other_leaves) == 1:
                continue
            # handle case of tree of size 3
            if len(other_leaves) == 3:
                return dist_list

            # finally get sibling
            sibling = filter(lambda i: i is not leaf,other_leaves)[0]

            # skip over paralogs
            if sibling.species == leaf.species:
                continue

            # get distances on gene and species trees
            gene_dist = sibling.DistTo(leaf)[1]
            spec_sibling = species_tree.node_dict[sibling.species]
            spec_leaf = species_tree.node_dict[leaf.species]
            spec_dist = spec_leaf.DistTo(spec_sibling)[1]
            dist_list.append((spec_dist,gene_dist))
        return dist_list


    def TreeLink(self,other_tree):
        ''' code to link the nodes of two trees whose newick ordering
        of nodes may be different, but whose overall topology is
        identical. '''

        match_dict = {}
        for self_node_str in self.node_dict:
            self_node = self.node_dict[self_node_str]
            self_leaves = self_node.leaves
            for other_node_str in other_tree.node_dict:
                other_node = other_tree.node_dict[other_node_str]
                if other_node.leaves == self_leaves:
                    match_dict[self_node] = other_node
                    break
        if len(match_dict) != len(self.node_dict):
            print "TreeLink failed"
            sys.exit(1)
        return match_dict


    def UnrootPrint(self):
        ''' take in a rooted tree and return an unrooted one '''

        if self.root is None:
            print "sorry, need to start with a rooted tree"
            sys.exit(1)

        unroot_newick = "("
        root = self.root
        # get kids
        kid_list = []
        for i in root.child_branches:
            for j in i.ends:
                if j is not root:
                    kid_list.append(j)
        # how many kids are there?
        num_kids_left = len(kid_list[0].leaves)
        num_kids_rght = len(kid_list[1].leaves)

        if num_kids_left > 1:
            left_kids = []
            for i in kid_list[0].child_branches:
                for j in i.ends:
                    if j is not kid_list[0]:
                        left_kids.append(j)
            unroot_newick += left_kids[0].treePrint('') + ","
            unroot_newick += left_kids[1].treePrint('') + ","
            unroot_newick += kid_list[1].treePrint('') + ");"
        elif num_kids_rght > 1:
            rght_kids = []
            for i in kid_list[1].child_branches:
                for j in i.ends:
                    if j is not kid_list[1]:
                        rght_kids.append(j)
            unroot_newick += kid_list[0].treePrint('') + ","
            unroot_newick += rght_kids[0].treePrint('') + ","
            unroot_newick += rght_kids[1].treePrint('') + ");"
        else:
            print "tree too small to unroot!"
            sys.exit(1)

        return unroot_newick


    def SPR(self,prune_node,graft_node):

        # find nodes
        prune_parent = prune_node.GetParent()
        graft_parent = graft_node.GetParent()

        # don't bother if they are already siblings
        if prune_parent is graft_parent:
            return self

        # get the sibling
        prune_sibling = prune_node.GetSibling()

        # adjust the pruned subtree
        prune_parent_branch = prune_parent.parent_branch
        prune_parent_branch.length += prune_sibling.parent_branch.length
        prune_sibling.parent_branch = prune_parent_branch
        prune_gparent = prune_parent.GetParent()
        # prune_node.parent_branch.length /= 2

        if prune_gparent is None:
            self.root = prune_sibling
        else:
            del prune_gparent.kid_nodes[prune_parent]
            prune_gparent.kid_nodes[prune_sibling] = 1

        # now, handle the grafting
        new_node = node.node('new',self)
        graft_parent_branch = graft_node.parent_branch
        new_branch_length = graft_parent_branch.length/2
        new_node.parent_branch = branch.branch(new_branch_length)
        graft_node.parent_branch.length = new_branch_length

        del graft_parent.kid_nodes[graft_node]
        graft_parent.kid_nodes[new_node] = 1
        new_node.kid_nodes[graft_node] = 1
        new_node.kid_nodes[prune_node] = 1
        new_tree = multitree()

        tree_str = "(" + str(self) + ");"
        new_tree.Build(tree_str)

        return new_tree



