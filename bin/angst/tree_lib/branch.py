class branch:

    def __init__(self,length):
        self.ends = []                  # nodes connected to
        self.length = float(length)	# length (can change during rooting)
        self.visited = False            # used for traversing the tree

    def __repr__(self):

	if len(self.ends) == 2:
	    print_string = "(" + self.ends[0].name + ","
	    print_string += self.ends[1].name + "):" + str(self.length)
	else:
	    print_string = ":" + str(self.length)
	return print_string

    def addNode(self,node):
	self.ends.append(node)
	node.branch_list.append(self)

    # recursion for finding all of the branches in an unrooted tree
    def findBranches(self,all_branches):

	all_branches.append(self)
	self.visited = True

	for node in self.ends:
	    for brch in node.branch_list:
		if not brch.visited:
		    all_branches = brch.findBranches(all_branches)

	return all_branches

