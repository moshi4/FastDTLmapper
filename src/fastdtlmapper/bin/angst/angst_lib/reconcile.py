# stick all the reconciliation-related code in here

import sys
import pdb
import math

import node
import branch

class dlt_model:

    min_val = math.exp(-100)

    def __init__(self,input_info):

        # take input penalties
        penalty_scale = math.pow(10,input_info.penalty_weight)
        if input_info.penalties_filename is not None:
            penalties_file = open(input_info.penalties_filename,'r')
            for line in penalties_file:
                parts = line.strip().split(': ')
                if parts[0].count('hgt') > 0:
                    self.hgt_penalty = float(parts[1])*penalty_scale
                if parts[0].count('spc') > 0:
                    self.spc_penalty = float(parts[1])*penalty_scale
                if parts[0].count('dup') > 0:
                    self.dup_penalty = float(parts[1])*penalty_scale
                if parts[0].count('los') > 0:
                    self.los_penalty = float(parts[1])*penalty_scale
            penalties_file.close()
        self.special = input_info.special
        self.outgroup = input_info.outgroup
        self.luca = input_info.luca
        self.species_tree = None
        self.event_guide = None

        # read in the event guide if present
        if input_info.event_guide_fn is not None:
            self.event_guide = {}
            event_guide_f = open(input_info.event_guide_fn,'r')
            for line in event_guide_f:
                parts = line.strip().split(' | ')
                event = parts[0]
                scale = float(parts[1])
                self.event_guide[event] = scale
            event_guide_f.close()


class node_link:
    ''' a data structure for connecting all of the multitrees  '''

    def __init__(self,species_tree,model,ultra_tree,tree_link):
        self.max_float = math.exp(500)
        self.link_dict = {}
        self.global_lcas = {}
        self.putative_root = {}
        self.dl_eval = {}
        self.tr_eval = {}
        self.best_score = self.max_float
        self.best_res = None
        self.learn_events = None
        self.cog_scaling = -1
        self.hr_scaling = -1

        # this is sort of a hack storing these in here, but it's still
        # better i feel than making these guys global
        self.species_tree = species_tree
        self.model = model
        self.ultra_tree = ultra_tree
        self.tree_link = tree_link

    def __repr__(self):
        for k in self.link_dict:
            print k,
            print "\t",
            print self.link_dict[k]
        return ""


    def Reset(self):
        ''' resets node link for iterative usage over different
        scaling factors '''
        self.global_lcas = {}
        self.putative_root = {}
        self.best_score = self.max_float
        self.dl_eval = {}
        self.tr_eval = {}
        self.best_res = None


    def MergeNodeLinkDicts(self,other_dict):
        '''merge lists of linked nodes'''

        for leaves in other_dict:
            if leaves in self.link_dict:
                for tup in other_dict[leaves]:
                    self.link_dict[leaves].append(tup)
            else:
                self.link_dict[leaves] = other_dict[leaves]

    def GlobalReconcile(node_link,leaves_str):
        '''reconcile a list of linked nodes'''

        node_list = node_link.link_dict[leaves_str]
        global_lcas = node_link.global_lcas

        # don't bother reconciling if done already
        if leaves_str in global_lcas:
            return
        global_lcas[leaves_str] = {}
        merge_global = global_lcas[leaves_str]
        merge_global["scores"] = {}
        merge_global["ranges"] = {}
        merge_global["newick"] = {}
        merge_global["paradox"] = {}
        if node_link.learn_events > 0:
            merge_global["events"] = {}

        # step through each of the nodes and reconcile
        total_branch_len = 0.0
        for node_pair in node_list:
            this_node = node_pair[0]
            parent_node = node_pair[1]
            merge_local = LocalReconcile(this_node,parent_node,node_link)
            MergeLCALookups(merge_global,merge_local,node_link)

            branch_len = this_node.branch_dict[parent_node].length
            total_branch_len += branch_len

        # average over all branch lengths.  arguably, this will
        # prevent overfitting & speed up angst
        avg_len = total_branch_len/len(node_list)
        for node_pair in node_list:
            this_node = node_pair[0]
            parent_node = node_pair[1]
            parent_branch = this_node.branch_dict[parent_node]
            parent_branch.length = avg_len

        # if you don't do this, loss caching will start to return
        # inaccurate results.  (i think this applies only if you're
        # looking at branch lengths)
        #species_tree = node_link.species_tree
        #species_tree.rec_dict = {}
        #species_tree.los_dict = {}

    def TryRoot(node_link,node_tuple):
        ''' try rooting an unrooted tree between the two nodes in the
        given tuple '''
        # pull out the two nodes
        node1 = node_tuple[0]
        node2 = node_tuple[1]
        if node_link.model.special == 'hr':
            TryDL(node2,node1,node1,node2,node_link,None)
        else:
            TryDL(node2,node1,node1,node2,node_link,None)
            TryTransfer(node2,node1,node1,node2,node_link,None)
            TryTransfer(node1,node2,node2,node1,node_link,None)


def LocalReconcile(this_node,parent_node,node_link):
    '''perform a reconciliation at a single node/parent combination
    note that there is no recursion'''

    model = node_link.model
    ultra_tree = node_link.ultra_tree
    tree_link = node_link.tree_link

    # handle differently depending on whether or not what's being
    # reconciled is rooted
    lca_lookups = {}
    lca_lookups["scores"] = lca_scores = {}
    lca_lookups["ranges"] = lca_ranges = {}
    lca_lookups["events"] = lca_events = {}
    lca_lookups["newick"] = lca_newick = {}
    lca_lookups["paradox"] = lca_paradox = {}

    if this_node.tree.unrooted:
        other_nodes = this_node.other_nodes
        child_nodes = filter(lambda i:i is not parent_node,other_nodes)
    else:
        child_nodes = this_node.kid_nodes.keys()
        this_node.lca_map = {}

    # leaves get treated special
    if len(child_nodes) < 1:

        # evaluate all possible ancestors
        species_node = node_link.species_tree.node_dict[this_node.species]
        ancestor_list = species_node.GetAncestors()
        event_guide = model.event_guide

        # run into problems with HR if leaves allowed to have built-in
        # losses
        if model.special == 'hr':
            ancestor_list = [species_node]

        for ancestor in ancestor_list:

            l_str, l_counts = CountLosses(ancestor,species_node.leaves,node_link)
            score  = l_counts['spc']*model.spc_penalty
            score += l_counts['los']*model.los_penalty

            if node_link.learn_events > 0:
                lca_events[ancestor] = ["[cur]: " + this_node.species]
                lca_events[ancestor].extend(l_str)
                # adjust scores based on event guide
                if event_guide is not None:
                    score -= EventGuideAdjust(l_str,event_guide,model)

            # take care of the height
            if node_link.tree_link is not None:
                lower_height = ultra_tree.height_dict[tree_link[ancestor]]
                parent = ancestor.GetParent()
                if parent is None:
                    upper_height = lower_height + 1.0
                else:
                    upper_height = ultra_tree.height_dict[tree_link[parent]]
                lca_ranges[ancestor] = [lower_height, upper_height]
            else:
                lca_ranges[ancestor] = None

            lca_scores[ancestor] = score
            lca_newick[ancestor] = ""
            lca_paradox[ancestor] = []
            if not this_node.tree.unrooted:
                this_node.lca_map[ancestor] = {}
                this_node.lca_map[ancestor][this_node] = ancestor
        return lca_lookups

    kid1 = child_nodes[0]
    kid2 = child_nodes[1]

    if this_node.tree.unrooted:
        for kid in child_nodes:
            leaves_str = repr(kid.leaf_dict[this_node])
            if not leaves_str in node_link.global_lcas:
                node_link.GlobalReconcile(leaves_str)
    else:
        for kid in child_nodes:
            new_table = LocalReconcile(kid,this_node,node_link)
            kid.lca_lookups = new_table

    # evaluate transfer scenarios

    TryTransfer(this_node,this_node,kid1,kid2,node_link,lca_lookups)
    TryTransfer(this_node,this_node,kid2,kid1,node_link,lca_lookups)

    # evaluate non-transfer scenarios
    TryDL(this_node,this_node,kid1,kid2,node_link,lca_lookups)
    return lca_lookups


def TryDL(parent1,parent2,node1,node2,node_link,lca_tables):
    ''' evalute vanilla (just duplication and loss) reconciliation
    scenarios '''

    # load in necessary variables
    model = node_link.model
    species_tree = node_link.species_tree
    event_guide = model.event_guide
    ultra_tree = node_link.ultra_tree
    tree_link = node_link.tree_link

    # if we've computed this scenario before, there's no point
    # repeating ourselves ...
    if node1.tree.unrooted:
        r1 = repr(node1.leaf_dict[parent1])
        r2 = repr(node2.leaf_dict[parent2])
        # avoid duplicating work (don't have to worry about branch
        # lengths, since we've averaged over them already)
        if r1 > r2:
            ID_str = ''.join([r1,r2])
        else:
            ID_str = ''.join([r2,r1])
        if ID_str in node_link.dl_eval:
            return
        node_link.dl_eval[ID_str] = 1
        nl1 = node_link.global_lcas[r1]
        nl2 = node_link.global_lcas[r2]
    else:
        nl1 = node1.lca_lookups
        nl2 = node2.lca_lookups

    node1_lca_scores = nl1["scores"]
    node2_lca_scores = nl2["scores"]
    node1_lca_ranges = nl1["ranges"]
    node2_lca_ranges = nl2["ranges"]
    node1_lca_paradox = nl1["paradox"]
    node2_lca_paradox = nl2["paradox"]

    if node_link.learn_events > 0:
        node1_lca_events = nl1["events"]
        node2_lca_events = nl2["events"]

    for lca1 in node1_lca_scores:
        score1 = node1_lca_scores[lca1]
        for lca2 in node2_lca_scores:
            score2 = node2_lca_scores[lca2]
            base_score = score1 + score2
            lcas_lca = species_tree.root.FindLCA(lca1,lca2)

            # reconcile the two lcas
            eQ, counts = MiniReconcile(species_tree,lca1,lca2,node_link)
            num_dups = counts['dup']
            num_loss = counts['los']
            num_vert = counts['spc']
            if num_dups < 1 and num_loss < 1:
                num_vert += 1

            ancestor_list = lcas_lca.GetAncestors()
            prev_dist_penalty = node_link.max_float
            for ancestor in ancestor_list:

                # is there no chance?
                sscore  = base_score
                if lca_tables is None:
                    if sscore > node_link.best_score:
                        break
                else:
                    cur_lca_scores = lca_tables['scores']
                    if ancestor in cur_lca_scores:
                        if sscore > cur_lca_scores[ancestor]:
                            continue

                # results datastructure
                res = {}
                res['node1'] = node1
                res['node2'] = node2
                res['lca1'] = lca1
                res['lca2'] = lca2
                res['parent1'] = parent1
                res['parent2'] = parent2
                res['node_link'] = node_link
                res['ancestor'] = ancestor

                # pass variables to update score
                dist_penalty = 0
                frac = 0.5
                dist_ratios = None
                sscore += dist_penalty

                if lca_tables is None:
                    if sscore >= node_link.best_score:
                        continue
                else:
                    cur_lca_scores = lca_tables['scores']

                # count scores
                l_str, l_counts = CountLosses(ancestor,lcas_lca.leaves,node_link)
                sscore += num_dups*model.dup_penalty
                sscore += (num_loss+l_counts['los'])*model.los_penalty
                sscore += (num_vert+l_counts['spc'])*model.spc_penalty

                if node_link.learn_events > 0:
                    event0 = []
                    event0.extend(l_str)
                    event0.extend(eQ)
                    if num_dups < 1 and num_loss < 1:
                        event0.extend(["[spc]: " + str(lcas_lca.species)])
                    # solicit input from event guide
                    if event_guide is not None:
                        sscore -= EventGuideAdjust(event0,event_guide,model)

                    event0.extend(node1_lca_events[lca1])
                    event0.extend(node2_lca_events[lca2])
                else:
                    event0 = None

                # handle homologous recombination in a special way.
                # require that you must have an HGT to correspond with
                # every loss
                if node_link.model.special == 'hr':
                    if event_guide is not None:
                        print "can't do both HR and event guiding!"
                        sys.exit(1)

                    # are there any losses proposed?
                    los_list = []
                    for event in event0:
                        if 'los' in event:
                            los_parts = event.split(': ')
                            los_node = los_parts[1].replace('*','')
                            los_list.append(los_node)
                    # get the list of hgt
                    hgt_list = []
                    for event in event0:
                        if 'hgt' in event:
                            hgt_parts = event.split(' --> ')
                            hgt_list.append(hgt_parts[1])

                    # note that you have to look through potentially
                    # suboptimal loss scenarios.  therefore, will have
                    # to "pop open" losses ...
                    raw_los_list = []
                    for los in los_list:
                        raw_los_list.extend(los.split('-'))
                    for hgt in hgt_list:
                        raw_hgt = hgt.split('-')
                        for hgt_part in raw_hgt:
                            if hgt_part in raw_los_list:
                                raw_los_list.remove(hgt_part)
                            else:
                                sscore = node_link.max_float
                                # could speedup here by adding a
                                # double-break statement.

                    if len(raw_los_list) > 0:
                        sscore = node_link.max_float

                    # make sure losses only percolated down, not up
                    for ind in range(len(los_list)):
                        los_list[ind] = "-" + los_list[ind] + "-"
                    los_list.sort(lambda a,b:len(a)-len(b))
                    new_los_list = []
                    for los in hgt_list:
                        matched = False
                        for ind in range(len(los_list)):
                            this_los = "-" + los + "-"
                            if this_los in los_list[ind]:
                                matched = True
                                new_los_list.append(los)
                                this_str = los_list[ind]
                                remaining_str = this_str.replace(this_los,'-')
                                los_list.remove(this_str)
                                if len(remaining_str) > 1:
                                    los_list.append(remaining_str)
                                break
                        los_list.sort(lambda a,b:len(a)-len(b))
                        if not matched:
                            sscore = node_link.max_float
                            break

                    if sscore < node_link.max_float:
                        # since python is insane, i can't directly
                        # remove things from this list.  instead, i'll
                        # just add them to a new list
                        orig_los_count = 0
                        keep_list = []
                        orig_los_list = []
                        for event in event0:
                            if event.count('los') < 1:
                                keep_list.append(event)
                            else:
                                orig_los_count += 1
                                orig_los_list.append(event)

                        # add hgt-inspired ones in
                        new_los_count = 0
                        for los in new_los_list:
                            this_los = '[los]: ' + los
                            keep_list.append(this_los)
                            new_los_count += 1
                        event0 = keep_list
                        # remove loss penalties from score
                        sscore -= orig_los_count*model.los_penalty
                        # add in new loss penalties
                        sscore += new_los_count*model.los_penalty

                paradox_count = []
                paradox1 = node1_lca_paradox[lca1]
                paradox2 = node2_lca_paradox[lca2]
                paradox_count.extend(paradox1)
                paradox_count.extend(paradox2)
                res['paradox'] = paradox_count

                # make a note on the dist ratios that these dists came
                # out of a speciation
                if dist_ratios is not None:
                    if num_dups > 0:
                        dist_ratios.append('yes dup')
                    else:
                        dist_ratios.append('yes spc')

                if node_link.tree_link is not None:
                    lower_height = ultra_tree.height_dict[tree_link[ancestor]]
                    parent = ancestor.GetParent()
                    if parent is None:
                        upper_height = lower_height + 1.0
                    else:
                        upper_height = ultra_tree.height_dict[tree_link[parent]]
                    #endif

                    if num_dups > 0:
                        # need to deal w/ situation where dups erase
                        # information about ranges, which leads to
                        # problems when you're trying to prevent time
                        # paradoxes.  see notes from 05.25.10.
                        if lca1 is lca2 is ancestor:    # see 05.28.10
                            min_age_1 = node1_lca_ranges[lca1][0]
                            min_age_2 = node2_lca_ranges[lca2][0]
                            min_age = max([min_age_1,min_age_2])
                            max_age_1 = node1_lca_ranges[lca1][1]
                            max_age_2 = node2_lca_ranges[lca2][1]
                            max_age = min([max_age_1,max_age_2])
                            ancestor_range = [min_age, max_age]
                        elif lca1 is ancestor:
                            ancestor_range = node1_lca_ranges[lca1]
                        elif lca2 is ancestor:
                            ancestor_range = node2_lca_ranges[lca2]
                        else:
                            ancestor_range = [lower_height, upper_height]
                        #endif
                    else:
                        ancestor_range = [lower_height, upper_height]
                    #endif
                else:
                    ancestor_range = None

                res['score'] = sscore
                res['range'] = ancestor_range
                res['events'] = event0
                res['newick'] = None
                res['dist_ratios'] = dist_ratios

                # what to do when trying to root
                if lca_tables is None:
                    # if you're rooting at the LUCA, make sure to give
                    # infinite score to things not born at the LUCA
                    if model.luca == "Root" or model.special == "hr":
                        if ancestor is not species_tree.root:
                            sscore = self.max_float

                    if node_link.best_score > sscore:
                        res['newick'] = MakeNewick(res,node_link,frac)
                        node_link.best_score = sscore
                        node_link.best_res = res

                else:
                    UpdateScores(lca_tables,ancestor,res)

def MiniReconcile(tree,node1,node2,node_link):
    ''' perform a speciation / duplication / loss reconciliation '''

    k = [node1,node2]
    k.sort()
    s = str(k)

    if s in tree.rec_dict:
        return tree.rec_dict[s]

    counts = {}
    counts['dup'] = 0
    counts['los'] = 0
    counts['spc'] = 0

    # call dup if children overlap:
    is_dup = False
    for ns in node1.subnodes:
        if ns in node2.subnodes:
            is_dup = True
            break

    # if there is a dup, figure out who's older:
    older_node = None
    if len(node1.subnodes) > len(node2.subnodes):
        older_node = node1
        young_node = node2
    else:
        older_node = node2
        young_node = node1

    if is_dup:
        counts['dup'] += 1

    # can we also just directly find the losses?
    llca = tree.root.FindLCA(node1,node2)

    # make this union, not amalgamation
    union = node1.leaves.copy()
    node1_leaves = node1.leaves
    node2_leaves = node2.leaves
    for k in node2_leaves:
        if not k in node1_leaves:
            union[k] = 1

    # look for losses
    if is_dup:
        loss_str, l_counts = CountLosses(llca,young_node.leaves,node_link)
    else:
        loss_str, l_counts = CountLosses(llca,union,node_link)

    counts['los'] = l_counts['los']
    counts['spc'] = l_counts['spc']

    if node_link.learn_events > 0:
        event_list = []
        if is_dup:
            event_list.extend(["[dup]: " + str(older_node.species)])
        event_list.extend(loss_str)
        event_str = event_list
        del event_list
    else:
        event_str = None

    tree.rec_dict[s] = event_str, counts
    return event_str, counts


def CountLosses(this_node,node_range,node_link):
    ''' recursively count losses '''

    species_tree = this_node.tree

    # store results we already know
    ##los_dict = species_tree.los_dict
    ##this_tuple = str((this_node,node_range))
    ##if this_tuple in los_dict:
    ##    event_str = los_dict[this_tuple][0]
    ##    counts = los_dict[this_tuple][1]
    ##    return event_str, counts

    # can i delete this node and still have the necessary leaves
    # to produce 'node_range'?
    set1 = set(node_range)
    set2 = set(this_node.leaves)
    intersect = set1.intersection(set2)
    counts = {}
    counts['los'] = 0
    counts['spc'] = 0
    len_intersect = len(intersect)

    if node_link.learn_events > 0:
        event_list = []
        if len_intersect < 1:
            # ok to lose:
            counts['los'] += 1
            event_list.extend(["[los]: " + str(this_node.species)])
        elif len_intersect != len(this_node.leaves):
            # not ok to lose:
            counts['spc'] += 1
            event_list.extend(["[spc]: " +str(this_node.species)])
            for kid in this_node.kid_nodes.keys():
                that_list, l_counts = CountLosses(kid,node_range,node_link)
                counts['los'] += l_counts['los']
                counts['spc'] += l_counts['spc']
                event_list.extend(that_list)
        event_str = event_list
    else:
        event_str = None
        if len_intersect < 1:
            # ok to lose:
            counts['los'] += 1
        elif len_intersect != len(this_node.leaves):
            # not ok to lose:
            counts['spc'] += 1
            for kid in this_node.kid_nodes.keys():
                that_list, l_counts = CountLosses(kid,node_range,node_link)
                counts['los'] += l_counts['los']
                counts['spc'] += l_counts['spc']

    ## los_dict[this_tuple] = (event_str, counts)
    return event_str, counts


def TryTransfer(parent1,parent2,node1,node2,node_link,lca_tables):
    ''' evaluate putative transfer scenarios '''

    # load in necessary variables
    model = node_link.model
    species_tree = node_link.species_tree
    ultra_tree = node_link.ultra_tree
    tree_link = node_link.tree_link
    event_guide = model.event_guide

    # grab out the lcas
    if node1.tree.unrooted:
        r1 = repr(node1.leaf_dict[parent1])
        r2 = repr(node2.leaf_dict[parent2])

        # avoid duplicating work (don't have to worry about branch
        # lengths, since we've averaged over them already.  note that
        # order matters here, since we're doing transfers ...
        ID_str = ''.join([r1,r2])
        if ID_str in node_link.tr_eval:
            return
        node_link.tr_eval[ID_str] = 1

        nl1 = node_link.global_lcas[r1]
        nl2 = node_link.global_lcas[r2]

    else:
        nl1 = node1.lca_lookups
        nl2 = node2.lca_lookups

    node1_lca_scores = nl1["scores"]
    node2_lca_scores = nl2["scores"]
    node1_lca_ranges = nl1["ranges"]
    node2_lca_ranges = nl2["ranges"]
    if node_link.learn_events > 0:
        node1_lca_events = nl1["events"]
        node2_lca_events = nl2["events"]
        node1_lca_paradox = nl1["paradox"]
        node2_lca_paradox = nl2["paradox"]

    # score their assembly
    for lca1 in node1_lca_scores:
        transfer_dict = {}
        for lca2 in node2_lca_scores:
            # if there's an ultrametric tree:
            if node_link.tree_link is not None:
                # get donor and acceptor age ranges
                donor_range = node1_lca_ranges[lca1]
                accep_range = node2_lca_ranges[lca2]
                # what's the minimum age of the donor?
                min_donor_age = None
                # case 1 of no overlap
                if donor_range[0] >= accep_range[1]:
                    continue
                # case 2 of no overlap
                if donor_range[1] <= accep_range[0]:
                    continue
                # potential overlap scenarios
                if donor_range[0] > accep_range[0]:
                    min_donor_age = donor_range[0]
                else:
                    min_donor_age = accep_range[0]
            else:
                min_donor_age = 0.0

            # for now, disallow phantom transfers within clade:
            if lca1.AreRelated(lca2):
                continue

            lca_str = str(lca1) + "-->" + str(lca2)
            species_tree.possible_hgt[lca_str] = 1
            transfer_dict[lca2] = [node2_lca_scores[lca2], min_donor_age]

            ######################################################
            # apply penalty according to how distant parents are #
            ######################################################
            if model.special == 'hr':
                host_dist = lca1.DistTo(lca2)[1]
                dist_penalty = math.pow(10.0,host_dist*node_link.hr_scaling)
                transfer_dict[lca2][0] += dist_penalty

        # find the cheapest transfer:
        min_val = node_link.max_float
        min_key = None
        min_min_donor_age = None
        for key in transfer_dict:
            if transfer_dict[key][0] < min_val:
                min_val = transfer_dict[key][0]
                min_key = key
                min_min_donor_age = transfer_dict[key][1]

        # if there were no eligible transfers, can move on:
        if min_key is None:
            continue

        # add lowest score
        sscore = node1_lca_scores[lca1]
        sscore += min_val + model.hgt_penalty

        # can back out already if the score is too high
        if lca_tables is None:
            if sscore > node_link.best_score:
                continue
        else:
            cur_lca_scores = lca_tables['scores']
            if lca1 in cur_lca_scores:
                if sscore > cur_lca_scores[lca1]:
                    continue

        res = {}
        res['node1'] = node1
        res['node2'] = node2
        res['lca1'] = lca1
        res['lca2'] = min_key
        res['parent1'] = parent1
        res['parent2'] = parent2
        res['node_link'] = node_link
        res['ancestor'] = lca1

        # pass variables to update score
        dist_penalty = 0
        frac = 0.5
        dist_ratios = None
        sscore += dist_penalty

        # can back out already if the score is too high
        if lca_tables is None:
            if sscore > node_link.best_score:
                continue
        else:
            cur_lca_scores = lca_tables['scores']
            if lca1 in cur_lca_scores:
                if sscore > cur_lca_scores[lca1]:
                    continue

        # keep track of records
        if node_link.learn_events > 0:
            all_events = []
            hgt_str = "[hgt]: " + str(lca1.species) + " --> " + str(min_key.species)
            all_events.extend([hgt_str])
            # solicit input from event guide
            if event_guide is not None:
                sscore -= EventGuideAdjust(all_events,event_guide,model)

            all_events.extend(node1_lca_events[lca1])
            all_events.extend(node2_lca_events[min_key])

            # count the number of impossible transfer scenarios
            donor_node = lca1
            accep_list = node2_lca_events[min_key]
            this_paradox = []
            paradox1 = node1_lca_paradox[lca1]
            paradox2 = node2_lca_paradox[min_key]
            this_paradox.extend(paradox1)
            this_paradox.extend(paradox2)
            for event in accep_list:
                if "hgt" in event:
                    accep_str = event.split(' --> ')[-1]
                    accep_node = species_tree.node_dict[accep_str]
                    node_relation = accep_node.IsAncestor(donor_node)
                    if node_relation == 1:
                        this_paradox.append([hgt_str,accep_list])
        else:
            this_paradox = []
            all_events = None

        # make a note on the dist ratios that these dists came
        # out of a speciation

        if dist_ratios is not None:
            dist_ratios.append('yes hgt')

        if node_link.tree_link is not None:
            final_range = [min_min_donor_age, donor_range[1]]
        else:
            final_range = None

        res['range'] = final_range
        res['paradox'] = this_paradox
        res['score'] = sscore
        res['events'] = all_events
        res['newick'] = None
        res['dist_ratios'] = dist_ratios

        # what to do when trying to root
        if lca_tables is None:
            # if you're rooting at the LUCA, make sure to give
            # infinite score to things not born at the LUCA
            if model.luca == "Root" or model.special == "hr":
                if lca1 is not species_tree.root:
                    sscore = self.max_float

            if sscore < node_link.best_score:
                res['newick'] = MakeNewick(res,node_link,frac)
                node_link.best_score = sscore
                node_link.best_res = res
        else:
            UpdateScores(lca_tables,lca1,res)


def MergeLCALookups(global_lcas,local_lcas,node_link):
    ''' merge dictionaries of lcas '''

    local_scores = local_lcas["scores"]
    local_ranges = local_lcas["ranges"]
    local_events = local_lcas["events"]
    local_newick = local_lcas["newick"]
    local_paradox = local_lcas["paradox"]

    for lca in local_scores:
        res = {}
        res['score'] = local_scores[lca]
        res['range'] = local_ranges[lca]
        res['newick'] = local_newick[lca]
        res['node_link'] = node_link
        res['paradox'] = local_paradox[lca]
        if node_link.learn_events > 0:
            res['events'] = local_events[lca]

        UpdateScores(global_lcas,lca,res)


def UpdateScores(lca_lookups,lcas_lca,res):
    ''' update scores in lca lookup tables '''

    sscore = res['score']
    ranges = res['range']
    newick = res['newick']
    node_link = res['node_link']
    if node_link.learn_events > 0:
        event0 = res['events']
    lca_scores = lca_lookups["scores"]
    lca_ranges = lca_lookups["ranges"]
    lca_newick = lca_lookups["newick"]
    lca_paradox = lca_lookups["paradox"]
    if node_link.learn_events > 0:
        lca_events = lca_lookups["events"]

    update_table = False
    if not lcas_lca in lca_scores:
        update_table = True
        update_newick = True
    elif lca_scores[lcas_lca] > sscore:
        update_table = True

    if update_table:
        lca_scores[lcas_lca] = sscore
        lca_paradox[lcas_lca] = res['paradox']
        lca_ranges[lcas_lca] = ranges
        if newick is None:
            newick = MakeNewick(res,node_link,None)
        lca_newick[lcas_lca] = newick

        # need to keep track of events when reconciling final tree or
        # modeling homologous recombination
        if node_link.learn_events > 0:
            lca_events[lcas_lca] = event0
            # can stop here if you're just tracking events for
            # homologous recombination
            if node_link.learn_events < 2:
                return

            node1 = res['node1']
            node2 = res['node2']
            parent1 = res['parent1']
            parent2 = res['parent2']
            parent1.dist_ratios[res['ancestor']] = res['dist_ratios']
            if parent1 is not parent2:
                print "parents not equivalent"
                sys.exit(1)
            # track lca assignments (useful when printing out events)
            if lcas_lca not in parent1.lca_map:
                parent1.lca_map[lcas_lca] = {}
            parent1.lca_map[lcas_lca][node1] = res['lca1']
            parent1.lca_map[lcas_lca][node2] = res['lca2']


def MakeNewick(res,node_link,frac):
    ''' generate the newick string '''

    node1 = res['node1']
    node2 = res['node2']
    lca1 = res['lca1']
    lca2 = res['lca2']
    parent1 = res['parent1']
    parent2 = res['parent2']
    ancestor = res['ancestor']

    if node1.tree.unrooted:
        leaves_1 = node1.leaf_dict[parent1]
        leaves_2 = node2.leaf_dict[parent2]
        node1_lca_newick = node_link.global_lcas[repr(leaves_1)]["newick"]
        node2_lca_newick = node_link.global_lcas[repr(leaves_2)]["newick"]
        # get branch lengths
        if parent1 not in node1.branch_dict:
            len1 = parent1.branch_dict[node1].length
        else:
            len1 = node1.branch_dict[parent1].length
        if parent2 not in node2.branch_dict:
            len2 = parent2.branch_dict[node2].length
        else:
            len2 = node2.branch_dict[parent2].length

        if frac is not None:
            len1 *= frac
            len2 *= 1-frac

    else:
        node1_lca_newick = node1.lca_lookups["newick"]
        node2_lca_newick = node2.lca_lookups["newick"]
        len1 = node1.parent_branch.length
        len2 = node2.parent_branch.length


    # write out node names
    if len(node1.branch_list) < 2:
        node1_name = node1.name + ":" + str(len1)
    else:
        node1_name = node1_lca_newick[lca1] + ":" + str(len1)

    if len(node2.branch_list) < 2:
        node2_name =  node2.name + ":" + str(len2)
    else:
        node2_name = node2_lca_newick[lca2] + ":" + str(len2)

    # handle cases where names get all mixed up:
    a_str = str(ancestor).replace('-','')
    newick = "(" + node1_name + "," + node2_name + ")" + a_str

    return newick



def EventGuideAdjust(event_list,event_guide,model):
    ''' provide score adjustment based on event guide '''
    score_adjust = 0.0
    for event in event_list:
        if event in event_guide:
            if 'spc' in event:
                score_adjust += event_guide[event]*model.spc_penalty
            elif 'los' in event:
                score_adjust += event_guide[event]*model.los_penalty
            elif 'hgt' in event:
                score_adjust += event_guide[event]*model.hgt_penalty
            elif 'dup' in event:
                score_adjust += event_guide[event]*model.dup_penalty
            else:
                print "strange event guide error."
                sys.exit(1)
    return score_adjust

