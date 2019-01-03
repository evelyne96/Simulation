#!/usr/bin/env python3
import numpy as np
import random 

N =100
nr_of_edges = 2500
graph = np.zeros((N,N)).astype(int)

def init_graph_random():
    for i in range(nr_of_edges):
        i = random.randint(0, N-1)
        j = random.randint(0, N-1)
        while graph[i][j] == 1:
            i = random.randint(0, N-1)
            j = random.randint(0, N-1)
        
        graph[i][j] = 1
        graph[j][i] = 1

def calc_conn():
    graph_connections = [0 for i in range(N)]

    for i in range(N):
        for j in range(N):
            graph_connections[i] += graph[i][j]
    return graph_connections

def distribution_plot():
    graph_connections = calc_conn()
    dist_prob = np.zeros(N).astype(float)
    # print(graph_connections)

    for i in range(N):
        dist_prob[i] = graph_connections.count(i)/N

    # print(dist_prob[graph_connections[0]])

    import matplotlib.pyplot as plt
    plt.plot(dist_prob)
    plt.show()

init_graph_random()
distribution_plot()

