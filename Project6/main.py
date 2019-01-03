#!/usr/bin/env python3
import numpy as np
import random 

N_max =100
N = 10
nr_of_edges = 2500
graph = np.zeros((N_max,N_max)).astype(int)

def init_graph_random():
    for i in range(nr_of_edges):
        i = random.randint(0, N_max-1)
        j = random.randint(0, N_max-1)
        while graph[i][j] == 1:
            i = random.randint(0, N_max-1)
            j = random.randint(0, N_max-1)
        
        graph[i][j] = 1
        graph[j][i] = 1

def init_graph_scale_free(N):
    for i in range(nr_of_edges):
        if i % 40 == 0:
            print("Timestap ", i)
            N = N + 1 if N < N_max else N

        i = random.randint(0, N-1)
        j = random.randint(0, N-1)
        try_count = 1
        while graph[i][j] == 1 and try_count<10:
            i = random.randint(0, N-1)
            j = random.randint(0, N-1)
            try_count += 1
        
        graph[i][j] = 1
        graph[j][i] = 1

def calc_conn():
    graph_connections = [0 for i in range(N_max)]

    for i in range(N_max):
        for j in range(N_max):
            graph_connections[i] += graph[i][j]
    return graph_connections

def distribution_plot():
    graph_connections = calc_conn()
    dist_prob = np.zeros(N_max).astype(float)
    # print(graph_connections)

    for i in range(N_max):
        dist_prob[i] = graph_connections.count(i)/N_max

    # print(dist_prob[graph_connections[0]])

    import matplotlib.pyplot as plt
    plt.plot(dist_prob)
    plt.show()

init_graph_scale_free(N)
distribution_plot()