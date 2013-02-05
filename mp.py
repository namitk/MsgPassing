#!/usr/bin/python

import random
import numpy as np
from union_find import union_find
import utils
import sys

class message_passing:
    graph_file = ""
    potential_file = ""
    num_nodes = 0 # number of nodes
    num_edges = 0 # number of edges
    node_arity = [] # arity of each node
    adj_list = [] # Initially, it is the adjacency list of the given graph. Is updated as the graph is triangulated
    adj_mat = [] # Initially, it is the adjacency matrix of the given graph. Is updated as the graph is triangulated
    rev_elim_order = [] # ordering returned by max cardinality search algorithm
    max_cliques = [] # maximal cliques in the graph after triangulation
    clique_factors = [] # will contain factors corresponding to each clique
    MAP_assignments = {} # will contain the assignments of the variables which maximise the factors.
    MAP_dependencies = {} # will contain the sepset that a variables depends for MAP
    MAP_evaluated = False
    jt_adj_list = [] # adjacency list of the junction tree structure
    normalizer = 0.0
    dummy_message=False

    def __init__(self, graph_file, potential_file):
        """ Does the following: 
        1. Read graph structure. 
        2. Perform MCS to triangulate the graph
        3. Enumerate max cliques in the triangulated graph
        4. Build junction tree
        5. Read factors of the distribution from file and assign them to one of the nodes of JT
        6. Calculate normalizing factor Z of the distribution
        """
        self.graph_file = graph_file
        self.potential_file = potential_file
        self.read_graph(graph_file)        
        self.max_card_search()
        num_added, added_edges = self.triangulate()        
        print num_added,
        self.enumerate_max_cliques()
        # print self.max_cliques
        tree = self.junction_tree()
        self.read_CPT(potential_file)        
        print max([len(x) for x in self.max_cliques]), 
        if self.num_nodes==1:
            print 0
        else:
            print max([len(set(self.max_cliques[x[0]]).intersection(self.max_cliques[x[1]])) for x in tree])
        #for i in range(len(self.max_cliques)):
        #    print i, len(self.max_cliques[i]), self.max_cliques[i]

        # Calculating the normalizing factor : Do one forward pass of sum-product MP algorithm and marginalize over the root 
        sum_marginalize = True
        root = 0 # 0th clique chosen as the root by default
        dfs = utils.DFS(self.jt_adj_list, root)
        self.jt_pass(dfs, sum_marginalize)
        # Z is the message received by an empty clique from the root
        self.max_cliques.append([])
        self.normalizer = self.calc_message(root, len(self.max_cliques)-1, sum_marginalize)[1]
        self.max_cliques.pop()
        
        # reset clique_factors for the original graph
        self.clique_factors = [{} for x in self.max_cliques]
        self.MAP_assignments = {}
        self.read_CPT(potential_file)
    
    def update_graph_structures(self, v0, v1):
        """ Update the adjacency list and adjacency matrix structures with  edge v0-v1 """
        self.adj_list[v0].append(v1)
        self.adj_list[v1].append(v0)
        self.adj_mat[v0][v1]=1
        self.adj_mat[v1][v0]=1 
    
    def read_graph(self, filename):
        """ Read graph structure from file and populate the adjacency matric and list structures of graph  """
        f = open(filename,'r')
        self.num_nodes = int(f.readline().strip())
        self.num_edges = int(f.readline().strip())
        self.node_arity = [0]*self.num_nodes
        for i in range(self.num_nodes):
            tok = f.readline().strip().split(' ')
            self.node_arity[int(tok[0])] = int(tok[1])
        self.adj_list = [[] for i in range(self.num_nodes)]
        self.adj_mat = [[0]*self.num_nodes for i in range(self.num_nodes)]
        for i in range(self.num_edges):
            tok = f.readline().strip().split(' ')
            v0 = int(tok[0])
            v1 = int(tok[1])
            self.update_graph_structures(v0, v1)       

    def max_nbr(self, processed_nodes, rem_nodes):
        """ From rem_nodes, find a node that has maximum number of neighbors in processed_nodes """
        l = list(rem_nodes)
        nbrs = [processed_nodes.intersection(set(self.adj_list[l[i]])) for i in range(len(l))]
        m = -1
        index = -1
        for i in range(len(nbrs)):
            if len(nbrs[i])>m:
                m = len(nbrs[i])
                index = i
        return l[index]

    def max_card_search(self):
        """ Maximum cardinality search algorithm """
        # r = random.randint(0,self.num_nodes-1)
        r = 0
        processed_nodes = set([r])
        rem_nodes = set(range(0, r)+range(r+1,self.num_nodes))
        self.rev_elim_order = [r]
        while len(rem_nodes)>0:            
            add = self.max_nbr(processed_nodes, rem_nodes)
            self.rev_elim_order.append(add)
            processed_nodes.add(add)
            rem_nodes.remove(add)

    def triangulate(self):
        """ Traingulate graph using the elimination order found by MCS """
        added_edges = []
        for i in range(len(self.rev_elim_order)-1, 0, -1):
            processed = set(self.rev_elim_order[0:i]) 
            nbrs = set(self.adj_list[self.rev_elim_order[i]])
            intersect = processed.intersection(nbrs)
            intersect.add(self.rev_elim_order[i])
            cur_clique = list(intersect)
            for i in range(len(cur_clique)):
                for j in range(i+1, len(cur_clique)):                       
                    if self.adj_mat[cur_clique[i]][cur_clique[j]]==0:
                        added_edges.append((cur_clique[i], cur_clique[j]))                        
                        self.update_graph_structures(cur_clique[i], cur_clique[j])        
        return len(added_edges), added_edges

    def enumerate_max_cliques(self):
        """ Find all the maximal cliques in the graph after triangulation """ 
        if self.num_nodes==1:
            self.max_cliques = [[0]]
            self.jt_adj_list = [[] for x in self.max_cliques] # initialize those many lists in jt_adj_list
            self.clique_factors = [{} for x in self.max_cliques]
            return
        n = [0]
        for i in range(1, len(self.rev_elim_order)):
            processed = set(self.rev_elim_order[0:i])
            nbrs = set(self.adj_list[self.rev_elim_order[i]])
            intersect = processed.intersection(nbrs)
            n.append(len(intersect))
        n.append(0)        
        for i in range(len(n)-2, 0, -1):
            if n[i]>=n[i+1]:
                processed = set(self.rev_elim_order[0:i])
                nbrs = set(self.adj_list[self.rev_elim_order[i]])
                intersect = processed.intersection(nbrs)                
                self.max_cliques.append([self.rev_elim_order[i]]+list(intersect))
        self.jt_adj_list = [[] for x in self.max_cliques] # initialize those many lists in jt_adj_list
        self.clique_factors = [{} for x in self.max_cliques]
        #self.MAP_assignments = [{} for x in self.max_cliques]
       
    def junction_tree(self):
        """ Junction tree of the triangulated graph """
        l = len(self.max_cliques)
        graph = [[0]*l for c in self.max_cliques]
        for i in range(l):
            for j in range(i+1, l):
                weight = len(set(self.max_cliques[i]).intersection(set(self.max_cliques[j])))
                graph[i][j]=-weight
                graph[j][i]=-weight
        tree = utils.min_spanning_tree(graph)
        for edges in tree:
            self.jt_adj_list[edges[0]].append(edges[1])
            self.jt_adj_list[edges[1]].append(edges[0])
        return tree

    def read_CPT(self, filename):
        """ Read potentials corresponding to the graph """ 
        f = open(filename, 'r')
        l = [int(x) for x in f.readline().strip().split(' ')[1:]]
        while(l):
            arity = [self.node_arity[x] for x in l]
            numLines = 1
            for a in arity:
                numLines = numLines*a 
            cur = np.empty(arity)            
            for i in range(numLines):
                ip = f.readline().strip().split(' ')                
                indices = [int(x) for x in ip[0:-1]]
                cur[tuple(indices)]=float(ip[-1])
            """ Assigning potentials / factors to one of the superset cliques """
            for i in range(len(self.max_cliques)):
                if set(l) <= set(self.max_cliques[i]):
                    self.clique_factors[i][tuple(l)] = cur.copy()
            l = [int(x) for x in f.readline().strip().split(' ')[1:]]
                
    def calc_message(self, sending_clique, recv_clique, sum_marginalize): 
        """
        Calculate the message to be sent from sending clique to receiving clique
        Message is the product of all factors in the sending clique marginalized over the 
        variables not in the separating set between sending and receiving clique
        """        
        sepset = tuple(set(self.max_cliques[sending_clique]).intersection(set(self.max_cliques[recv_clique])))
        factor = self.clique_factors[sending_clique].copy()
        if sepset in factor and not self.dummy_message:
            del(factor[sepset]) # removing message that was received from the receiving clique        
        variables = ()
        for vs in factor.keys():
            variables = variables + vs        
        variables = tuple(set(variables))
        sepset=tuple(set(sepset).intersection(set(variables)))
        nonsepset = tuple(set(variables)-set(sepset))
        if not self.MAP_evaluated:
            for x in nonsepset:
                self.MAP_dependencies[x] = sepset
                if len(sepset)>0:
                    self.MAP_assignments[x]=np.zeros(tuple([self.node_arity[y] for y in sepset]))
                else:
                    self.MAP_assignments[x]=np.zeros((1,))
        potential = np.zeros(tuple([self.node_arity[x] for x in sepset]))
        assignment = [0]*len(variables)
        sepset_assignment = [0]*len(sepset)
        while True:            
            result = 1.
            for f in factor:
                cur_values = []
                for var in f:
                    cur_values.append(assignment[variables.index(var)])
                result = result * factor[f][tuple(cur_values)]

            if sum_marginalize:
                potential[tuple(sepset_assignment)] += result
            else:
                if result >= potential[tuple(sepset_assignment)]:
                    potential[tuple(sepset_assignment)] = result # PENDING: store assignment for which this maximum value was achieved  
                    if not self.MAP_evaluated:
                        for x in nonsepset:
                            self.MAP_assignments[x][tuple(sepset_assignment)] = assignment[variables.index(x)]
            """ code to iterate through all possible assignments to variables. 
            Similar approach as used in counters """
            node = 0
            while(node < len(variables)):
                assignment[node] += 1
                if assignment[node] < self.node_arity[variables[node]]:
                    break
                assignment[node]=0
                node=node+1
            # print assignment
            if node==len(variables):
                break # finished all possible assignments to sepset
            sepset_assignment = [assignment[variables.index(x)] for x in sepset]
        return sepset, potential
        
    def send_message(self, sepset, potential, sending_clique, recv_clique):
        """ Send message i.e update the clique factors """ 
        if sepset in self.clique_factors[recv_clique]:
            potential = potential*self.clique_factors[recv_clique][sepset]
        self.clique_factors[recv_clique][sepset]=potential        

    def jt_pass(self, pass_order, sum_marginalize):
        """ 
        One pass of the MP algorithm in direction as governed by pass_order
        Two passes - forwad and backward - are enough to calibrate the JT
        """
        ready_cliques = [x for x in range(len(pass_order)) if not pass_order[x]['recv']]
        finished_cliques = []
        while ready_cliques:
            sending_clique = ready_cliques[0] # found a ready clique. What if none exists? Put try catch
            for recv_clique in pass_order[sending_clique]['send']:
                sepset, potential = self.calc_message(sending_clique, recv_clique, sum_marginalize)
                self.send_message(sepset, potential, sending_clique, recv_clique)
            finished_cliques.append(ready_cliques[0])
            ready_cliques.remove(ready_cliques[0]) # else the same clique will keep getting selected. Verify this statement

            # update cliques, if any, that have now become ready due to the sent messages
            for recv_clique in pass_order[sending_clique]['send']:
                notReady = False
                for recv_from in pass_order[recv_clique]['recv']:
                    if recv_from not in finished_cliques:
                        notReady = True
                if not notReady and pass_order[recv_clique]['send']:
                    ready_cliques.append(recv_clique)
    
    def MP(self, sum_marginalize):
        """ 
        The main message passing algorithm 
        Runs sum-product algorithm if sum_marginalize is True else runs max-product algorithm
        """
        root = 0 # 0th clique chosen as the root by default
        # Forward Pass
        dfs = utils.DFS(self.jt_adj_list, root)
        self.jt_pass(dfs, sum_marginalize)
        self.MAP_evaluated=True
        # Backward Pass 
        back_dfs = [{'send':x['recv'],'recv':x['send']} for x in dfs]
        self.jt_pass(back_dfs, sum_marginalize)        
        self.MAP_evaluated=False
     
    def calc_MAP_assignments(self):
        """ Evaluates the most probable global assignment to variables """
        num_variables=len(self.MAP_assignments)
        settled=()
        while len(settled)<num_variables:
            for x in self.MAP_assignments:
                if x not in settled:
                    if len(self.MAP_dependencies[x])==0:
                        self.MAP_assignments[x]=self.MAP_assignments[x][0]
                        settled=settled+(x,)
                    else:
                        contained = True
                        for a in self.MAP_dependencies[x]:
                            if a not in settled:
                                contained = False
                        if contained:
                            index = [self.MAP_assignments[y] for y in self.MAP_dependencies[x]]
                            self.MAP_assignments[x] = self.MAP_assignments[x][tuple(index)]
                            settled=settled+(x,)
                    
    def MAP(self):
        """ Find most probable global assignment to the variables """
        self.MP(False)
        root=0
        self.max_cliques.append([])
        prob_of_asgn = self.calc_message(root, len(self.max_cliques)-1, False)[1]*1.0/self.normalizer
        self.max_cliques.pop()
        self.calc_MAP_assignments();
        asgn = self.MAP_assignments
        print asgn, prob_of_asgn

        # reset clique_factors for the original graph
        self.clique_factors = [{} for x in self.max_cliques]
        self.MAP_assignments = {}
        self.read_CPT(self.potential_file)

    def sum_marginal_query(self):
        """ Given nodes that are guaranteed to be contained in one of the cliques, 
        finds the marginal probability of this set of variables
        """
        self.MP(True)
        print "Enter nodes whose marginal probability table is required: ",
        ip = [int(x) for x in raw_input().strip().split(' ')]
        for i in range(len(self.max_cliques)):
            if set(ip) <= set(self.max_cliques[i]):            
                """ Output format required: 
                X_1		X_2		Pr(X_1,X_2)
                
                0		0		Pr(X_1=0,X_2=0)
                0		1		Pr(X_1=0,X_2=1)
                0		2		Pr(X_1=0,X_2=2)
                1		0		Pr(X_1=1,X_2=0)
                1		1		Pr(X_1=1,X_2=1)
                1		2		Pr(X_1=1,X_2=2)
                ... (for all combinations of values of nodes)
                """
                # if set(ip)==set(self.max_cliques[i]):
                #     # get potential corresponding to this ip from the factors and divide it by the normalizer. Refer else case
                #     print self.clique_factors[i]
                # else:
                self.max_cliques.append(ip)
                self.dummy_message=True
                ip,potential = self.calc_message(i, len(self.max_cliques)-1, True)
                # Table header
                s=""
                for x in ip:
                    s=s+'X_'+str(x)+'\t'
                s=s+'Pr('
                for x in ip:
                    s=s+'X_'+str(x)+','
                l = list(s)
                l[len(l)-1]=')'
                print "".join(l)
                self.dummy_message=False
                norm = np.sum(potential)
                potential = np.divide(potential, norm)
                self.max_cliques.pop()
                # print probability table by dividing potential by the sum of elements in the table
                values = [(x,) for x in range(self.node_arity[ip[0]])]
                for i in range(1,len(ip)):
                    values = [x+(y,) for x in values for y in range(self.node_arity[ip[i]])]
                for a in values:
                    s=""
                    value=potential
                    for b in a:
                        s=s+str(b)+'\t'      
                    s=s+str(potential[a])
                    print s
                break

        # reset clique_factors for the original graph
        self.clique_factors = [{} for x in self.max_cliques]
        self.MAP_assignments = {}
        self.read_CPT(self.potential_file)

if len(sys.argv) < 3:
    print "Usage: "+sys.argv[0]+" graph_file potential_file"
    sys.exit(-1)
object = message_passing(sys.argv[1], sys.argv[2])
object.MAP()
object.sum_marginal_query()
