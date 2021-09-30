import sys
import pdb
import math
import multitree
import branch

class node:

    def __init__(self,raw_name,arbre):
        # remove any leading angled brackets
        raw_name = raw_name.lstrip(">")

        # handle mixed underscore and dot environment
        underscore_index = sys.maxint
        dot_index = sys.maxint
        if raw_name.count('_') > 0:
            underscore_index = raw_name.index('_')
        if raw_name.count('.') > 0:
            dot_index = raw_name.index('.')
        # case of species node
        self.name = raw_name
        if underscore_index == dot_index:
            self.species = raw_name
        # case of gene node
        elif dot_index < underscore_index:
            self.species = raw_name[:dot_index]
        else:
            self.species = raw_name[:underscore_index]

	self.tree = arbre               # tree node belongs to
        self.myBranch = None            # used to construct unrooted trees
	self.branch_list = []           # used to construct unrooted trees
	self.child_branches = []        # list of branches to kids
        self.kid_nodes = {}
	self.parent_branch = None         # should have length 1
	self.subnodes = dict()          # all nodes below this
	self.leaves = None              # leaves below this

        # attach to each node a serial number, to tell nodes apart
        self.serial = hash(self.name)
	self.visited = False            # useful for rooting trees
        self.link_dict_visited = False
        self.unrooted_leaving_visited = False

        # things necessary for UnrootedLeaving
        self.leaf_dict = {}
        self.branch_dict = {}
        self.other_nodes = []

        # things for reconciliation
        self.lca_lookups = None
        self.leaf_event_list = None
        self.rec_lca = None
        self.lca_map = None
        self.event_list = None
        self.raw_list = None

        # distance things
        self.leaf_dists = []
        self.dist_ratios = {}


    def __repr__(self):
        return self.name


    def __eq__(self,other):
	if self is None or other is None:
	    return False
	elif self.serial == other.serial:
	    return True
	else:
	    return False


    def __ne__(self,other):

	if self is None or other is None:
	    return True
	elif self.serial == other.serial:
	    return False
	else:
	    return True


    def __gt__(self,other):
	if self.name > other.name:
	    return True
	else:
	    return False

    def __lt__(self,other):
	if self.name < other.name:
	    return True
	else:
	    return False

    def __hash__(self):
        try:
            return self.serial
        except AttributeError:
            return 0

    def isLeaf(self):

        if len(self.branch_list) < 2:
            return True
        else:
            return False

    def addBranch(self,length):

	self.myBranch = branch.branch(length)
	self.myBranch.addNode(self)

    # recursively print out nodes
    def treePrint(self,newickString):

	count = 0
        # implement sorting so that you can compare trees in a quick
        # and dirty way
        kid_l = self.kid_nodes.keys()
        kid_l.sort()
        for kid_node in kid_l:
            if count == 0:
                newickString += "("
                newickString = kid_node.treePrint(newickString)
            else:
                newickString += ","
                newickString = kid_node.treePrint(newickString)
                newickString += ")"
            count += 1

	if self.parent_branch is not None and len(self.child_branches) == 0:
	    addString = self.name + ":" + str(self.parent_branch.length)
	elif self.parent_branch is not None:
	    addString = ":" + str(self.parent_branch.length)
	else:
	    addString = ""
	return newickString + addString

    # join two nodes in a central node
    def unite(self,node2):
        ''' create a new node '''

        # sort names to ease comparisons between trees
        new_node_name = self.name + "-" + node2.name
	center_node = node(new_node_name,self.tree)
	self.myBranch.addNode(center_node)
	node2.myBranch.addNode(center_node)
	return center_node

    def imposeHierarchy(self):
	'''find all the child nodes that have yet to be visited and
	assign their children'''

	for kid_branch in self.child_branches:
	    for node in kid_branch.ends:
		if not node.visited and node is not self:
                    self.kid_nodes[node] = 1
		    node.visited = True
		    node.parent_branch = kid_branch
                    node.tree.branch_list.append(node.parent_branch)
                    node.tree.node_dict[node.species] = node
		    for child_branch in node.branch_list:
			if child_branch is not kid_branch:
			    node.child_branches.append(child_branch)
			    node.imposeHierarchy()

    # recursively label the roots of subtrees w/ the leaves contained
    # below
    def subtreeLabel(self):

	leaf_vec = []

	if self.leaves is not None:
	    return self.leaves

	# recurse
        for kid_node in self.kid_nodes.keys():
            child_leaves = kid_node.subtreeLabel()
            leaf_vec.extend(child_leaves)

	# once you hit a leaf
	if len(self.child_branches) == 0:
	    leaf_vec.append(self.species)

	# non-duplicates:
	self.leaves = dict(map(lambda i: (i,1),leaf_vec))
	return self.leaves.keys()

    # figure out which species tree node a gene tree node maps to
    def subtreeMap(self,species_node):

	# do the children of the current species node possess the
	# relevant genes?  if so, follow that child.  if not, return
	# the current node

	foundSet = 0
        for kid_node in species_node.kid_nodes.keys():
            s_leaves = kid_node.leaves
            g_leaves = self.leaves

            # count how many elements of the gene set of
            # leaves are not in the species set of leaves
            n = 0
            for i in g_leaves.keys():
                if not i in s_leaves:
                    n += 1
                    break

            # if you can see somewhere to descend, follow it
            if n == 0:
                return self.subtreeMap(kid_node)

	# once you've run out of places to descend, just return
	# wherever you've ended up:
	return species_node

    # recursively store at each node a hash of all nodes below
    def Find_Subnodes(self):

        for kid_nodes in self.kid_nodes.keys():
            new_dict = kid_nodes.Find_Subnodes()
            for k in new_dict.keys():
                if not k in self.subnodes:
                    self.subnodes[k] = 1
	self.subnodes[self.species] = 1
	return self.subnodes


    def FindLCA(self,node1,node2):
        ''' find the last common ancestor of two nodes.  will descend
        from the subroot (self) looking to see which children possess
        both nodes.  if neither children possess both nodes, return
        current node '''
        for kid_nodes in self.kid_nodes.keys():
            if node1.species in kid_nodes.subnodes:
                if node2.species in kid_nodes.subnodes:
                    return kid_nodes.FindLCA(node1,node2)
        return self

    # determine if one node is descended from the other ...
    def AreRelated(self,other_node):
        '''is other_node a child of self?'''
        node1 = self
        node2 = other_node
        if len(node1.leaves) > len(node2.leaves):
            node2 = self
            node1 = other_node
	for e in node2.leaves.keys():
	    if e in node1.leaves:
                return True
        return False

    # method to list all the leaves of this node, keyed by the parent
    def UnrootedLeaving(this_node):

        other_nodes = this_node.GetOtherNodes()
        for parent_node in other_nodes:
            child_nodes = list(set(other_nodes).difference(set([parent_node])))
            this_node.GetChildLeaves(parent_node,child_nodes)

        # when done, move on to neighboring nodes:
        for node in other_nodes:

            if node.isLeaf():
                if node.name in node.tree.leaf_dict:
                    print "duplicate leaves in a bootstrap."
                    print "aborting prematurely."
                    sys.exit(1)
                node.tree.leaf_dict[node.name] = node

            if len(node.leaf_dict) is not len(node.branch_list):
                node.UnrootedLeaving()

    def GetLeaves(this_node,parent_node):

        # if the leaf dict has already been defined:
        if parent_node in this_node.leaf_dict:
            return this_node.leaf_dict[parent_node]

        # if is a leaf
        if len(this_node.branch_list) < 2:

            this_node.leaf_dict[parent_node] = [this_node.name]
            return this_node.leaf_dict[parent_node]

        other_nodes = this_node.GetOtherNodes()
        child_nodes = list(set(other_nodes).difference(set([parent_node])))
        this_node.GetChildLeaves(parent_node,child_nodes)

        return this_node.leaf_dict[parent_node]

    def GetChildLeaves(this_node,parent_node,other_nodes):

        merge_list = []
        for kid_node in other_nodes:
            if parent_node in this_node.leaf_dict:
                if len(this_node.leaf_dict[parent_node]) > 1:
                    continue
            kid_leaves = kid_node.GetLeaves(this_node)

            if type(kid_leaves[0]) is not type(''):
                merge_list.extend(kid_leaves[0])
            else:
                merge_list.extend(kid_leaves)

            # this is a total hack
            this_node.leaf_dict[parent_node] = []

        merge_list.sort()
        this_node.leaf_dict[parent_node].extend(merge_list)

    def GetNodeLinkDict(this_node,node_link_dict):

        this_node.link_dict_visited = True
        leaf_dict = this_node.leaf_dict

        for key in leaf_dict:

            subleaves = leaf_dict[key]
            if type(subleaves[0]) == type([]):
                subleaves = subleaves[0]
            sub_str = repr(subleaves)

            if sub_str in node_link_dict:
                node_link_dict[sub_str].append((this_node,key))
            else:
                node_link_dict[sub_str] = [(this_node,key)]

        # after that's been done, decide where to recurse
        for relative in leaf_dict:

            relative_keys = relative.leaf_dict.keys()
            subleaves = relative.leaf_dict[relative_keys[0]]
            subleaves.sort()

            if not relative.link_dict_visited:
                relative.GetNodeLinkDict(node_link_dict)

    # return a list of all bordering nodes
    def GetOtherNodes(this_node):

        # fill out the leaf dictionary
        other_nodes = []
        for branch in this_node.branch_list:
            for node in branch.ends:
                if node is not this_node:
                    other_nodes.append(node)
                    node.branch_dict[this_node] = branch
                    this_node.branch_dict[node] = branch
        this_node.other_nodes = other_nodes

        return other_nodes

    # get the distance between two nodes:
    # call as follows:
    # found, dist = node1.DistTo(node2)
    def DistTo(this_node,that_node,this_branch=None,dist=0):
        # what to do at the end
        if this_node is that_node:
            return True, dist
        # if not ...
        was_found = False
        found_dist = dist
        for i in this_node.branch_list:
            if i is not this_branch:
                for j in i.ends:
                    if j is not this_node:
                        found, new_dist = j.DistTo(that_node,i,dist+i.length)
                        if found is True:
                            was_found = found
                            found_dist = new_dist

        return was_found, found_dist


    def FindDists(self):
        ''' recursively enumerate all pairwise distances on the
        species tree '''

        dist_dict = self.tree.dist_dict
        kids = self.kid_nodes.keys()

        if self not in dist_dict:
            dist_dict[self] = {}

        # give own distance
        dist_dict[self][self] = 0.0

        # now, to each of those add the kids
        for kid in kids:

            # get a list of everyone 'self' is affiliated with
            parent_partners = dist_dict[self].keys()

            if kid not in dist_dict:
                dist_dict[kid] = {}
            dist_2_parent = kid.parent_branch.length
            for partner in parent_partners:
                this_dist = dist_dict[self][partner] + dist_2_parent
                dist_dict[kid][partner] = this_dist
                dist_dict[partner][kid] = this_dist
            kid.FindDists()


    def ValidTransfer(self,other_node):
        ''' return boolean regarding whether or not the parent
        branches of the queried nodes overlap temporally '''

        # no auto-transfer
        if self is other_node:
            return 0

        donor_node = self
        accep_node = other_node
        root = self.tree.root
        dist_dict = self.tree.dist_dict

        # note that there can be messiness with >= (as opposed to just
        # >), when you consider overlap between parent and child.
        # since we don't want transfer between parent and child,
        # branch lengths will be drawn in slightly and then compared
        min_dist = math.exp(-10)
        d11 = dist_dict[donor_node][root] - min_dist
        d12 = d11 - donor_node.parent_branch.length + min_dist
        d21 = dist_dict[accep_node][root] - min_dist
        d22 = d21 - accep_node.parent_branch.length + min_dist

        if d12 <= d21:
            if d22 <= d11:
                return 1 # direct transfer
            else:
                return 2 # phantom transfer
        else:
            return 0 # no transfer


    def IsAncestor(self,other_node):
        ''' is self an ancestor of other_node (1), the same as
        other_node (2), or neither (0). '''
        if self is other_node:
            return 2
        lca_node = self.tree.root.FindLCA(self,other_node)
        if lca_node is self:
            return 1
        else:
            return 0


    def GetAncestors(self):
        ''' get a list of ancestors from the current node (included in
        list ) to the root of the tree '''

        ancestor_list = [self]
        cur_node = self
        root = self.tree.root
        while cur_node is not root:
            # get the parent
            parent = filter(lambda i: i is not cur_node,cur_node.parent_branch.ends)[0]
            ancestor_list.append(parent)
            cur_node = parent
        return ancestor_list


    def GetHeights(self):
        ''' get MEDIAN distance to leaves at all internal nodes on the
        tree '''

        # get distances from leaves
        kids = self.kid_nodes.keys()
        my_len = self.parent_branch.length
        median = 0.0
        if len(kids) < 1:
            self.leaf_dists = [my_len]
        else:
            # get list of child distances
            child_dists = []
            for kid in kids:
                kid.GetHeights()
                child_dists.extend(kid.leaf_dists)
            # find the median
            child_dists.sort()
            mid_ind = float(len(child_dists))/2
            int_mid_ind = int(mid_ind)
            if mid_ind == int_mid_ind:
                median = sum(child_dists[int_mid_ind-1:int_mid_ind+1])/2
            else:
                median = child_dists[int_mid_ind]
            self.leaf_dists = [val + my_len for val in child_dists]

        # dist_range = [median, median + my_len]
        # dist_range = median + (my_len/2)
        dist_range = median

        self.tree.height_dict[self] = dist_range

        return


    def GetParent(self):
        ''' return the parent node '''
        nodes = self.parent_branch.ends
        parent_list = [node for node in nodes if node is not self]
        if len(parent_list) < 1:
            return None
        else:
            return parent_list[0]


    def AreSiblings(self,other_node):
        ''' evaluate whether or not arguments are sibling nodes '''
        parent1 = self.GetParent()
        parent2 = other_node.GetParent()

        if parent1 is parent2:
            return True
        # if they're not siblings:
        return False


    def GetSibling(self):
        ''' get the node sibling '''
        node_sibling = None
        node_parent = self.GetParent()

        if node_parent is None:
            pdb.set_trace()

        for kid_node in node_parent.kid_nodes.keys():
            if kid_node is not self:
                node_sibling = kid_node
        return node_sibling
