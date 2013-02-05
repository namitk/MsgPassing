#!/usr/bin/env python

from optparse import OptionParser, Option
import os, sys, fileinput, shutil, re, string
import numpy
from numpy import *
import scipy
from scipy import stats

parser = OptionParser()
parser.add_option("-n", "--nodes", dest="nodes", type="int", default=10, help="Number of Nodes")
parser.add_option("-c", "--cardinality", dest="cardinality", type="int", default=3, help="Maximum Cardinality of Node")
parser.add_option("-p", "--potentials", dest="potentials", type="int", default=30, help="Number of Potentials")
(options, args) = parser.parse_args()

nodes = options.nodes
cardinality = options.cardinality
potentials = options.potentials
graph_log = "graph.log"
potentials_log = "potentials.log"

node_cardinality = []
node_values = []
node_list = []
edge_set = []
done = False

graph_handler = open(graph_log, 'w+')
potentials_handler = open(potentials_log, 'w+')

if cardinality < 2:
    exit()

for i in range(nodes):
    temp_cardinality = random.randint(2, cardinality+1)
    node_cardinality.append(temp_cardinality)

"""
i = 0
flag = 0
while i < edges:
    temp_node1 = random.randint(0, nodes)
    temp_node2 = random.randint(0, nodes)
    if temp_node1 == temp_node2:
        temp_node2 = (temp_node2 + random.randint(1,nodes))%nodes
    if temp_node1 > temp_node2:
        temp = temp_node1
        temp_node1 = temp_node2
        temp_node2 = temp
    for temp_edge in edge_set:
        if temp_edge[0] == temp_node1 and temp_edge[1] == temp_node2:
            flag = 1
            break
    if flag == 1:
        flag = 0
        continue
    edge_set.append((temp_node1,temp_node2))
    i=i+1
    graph_handler.write(str(temp_node1) + " " + str(temp_node2) + "\n")
"""

def pareto(mode, shape):
    return mode * pow(1.0 - random.random(),-1.0 / shape)

def increment(pos):
    global node_cardinality, node_values, node_list, done
    if pos == len(node_values):
        done = True
        return    
    node_values[len(node_values)-pos-1] = (node_values[len(node_values)-pos-1]+1)%(node_cardinality[node_list[len(node_values)-pos-1]])
    if node_values[len(node_values)-pos-1] == 0:
        increment(pos+1)

i = 0
temp_node_list = []
flag = 0
while i < potentials:
    node_list = set([])
    num_nodes = int(pareto(1.0,2)) + 1
    for n in range(num_nodes):
        node_list = node_list | set([random.randint(0,nodes)])
    if len(node_list) < num_nodes:
        continue
    else:
        i=i+1
        temp_node_list = []
        for n in node_list:
            temp_node_list.append(n)
        node_list = temp_node_list
        sort(node_list)
        for x in node_list:
            for y in node_list:
                if x == y or x > y:
                    continue
                for temp_edge in edge_set:
                    if temp_edge[0] == x and temp_edge[1] == y:
                        flag = 1
                        break
                if flag == 1:
                    flag = 0
                    continue
                edge_set.append((x,y))
        potentials_handler.write("# ")
        for n in node_list:
            potentials_handler.write(str(n) + " ")
        potentials_handler.write("\n")
        node_values = []
        for n in range(len(node_list)):
            node_values.append(0)
        done = False
        while done == False:
            for n in node_values:
                potentials_handler.write(str(n) + " ")
            potentials_handler.write(str(random.uniform(0.0,100.0)) + "\n")
            increment(0)

edges = len(edge_set)

graph_handler.write(str(nodes) + "\n" + str(edges) + "\n")

for i in range(nodes):
    graph_handler.write(str(i) + " " + str(node_cardinality[i]) + "\n")

for edge in edge_set:
    graph_handler.write(str(edge[0]) + " " + str(edge[1]) + "\n")

