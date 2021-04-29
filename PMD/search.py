import heapq
import copy
from PMD import utils as utils
from PMD import analysis as ana
from PMD import nm
import numpy as np
import pickle
import os
from time import sleep

def CM_addition(CM1,CM2):
    CM = dict()
    if len(CM1)==0:
        if len(CM2)==0:
            return CM
        else:
            return CM2
    elif len(CM2) ==0:
        return CM1
    for key in CM1.keys():
        if key in CM2.keys():
            v1 = set(CM1[key])
            v2 = set(CM2[key])
            CM[key] = list(v1.union(v2))
        else:
            CM[key] = CM1[key]
    for key in CM2.keys():
        if key not in CM1.keys():
            CM[key] = CM2[key]
    return CM

class PriorityQueue:
    def __init__(self):
        self.list_elements = []
    
    def empty(self) -> bool:
        return len(self.list_elements)==0
    
    def push(self,item,priority) -> None:
        heapq.heappush(self.list_elements,(priority,item))
    
    def pop(self):
        if self.empty():
            return None
        else:
            return heapq.heappop(self.list_elements)[1]

class Node:
    def __init__(self,path):
        self.id = path[-1]
        self.path = path
        self.depth = len(path)
        self.cost = None
    
    def get_id(self):
        return self.id
    
    def get_path(self):
        return self.path

    def update_cost(self,cost):
        self.cost = cost

    def get_cost(self):
        return self.cost


class Astar:
    def __init__(self,init_id,target_id,pdb_file,path_dir):
        self.pdb = pdb_file
        self.nm = nm.nm(pdb_file,path_dir)
        self.init_id = init_id
        self.target_id = target_id
        self.explored_id = [init_id]
        self.connection_map = self.load_CM()
        self.current_path = [init_id]
        self.init_node = Node(self.current_path)
        self.current_node = Node(self.current_path)
        cost = self.evaluation_function(self.current_node)+self.heuristic_function(self.current_node)
        self.current_node.update_cost(cost)
        self.frontier = PriorityQueue()
        self.search_path = list()
        self.search_path.append(self.current_node)
        
    def get_successors(self,node,mode=0,frequency=1,force=100,time=500,velocity=False,temperature=300,cuda='0',threshold=3,window=1):
        res_id = int(node.id)
        nr = self.nm.pump(res_id,mode=mode,frequency=frequency,force=force,time=time,velocity=velocity,temperature=temperature,cuda=cuda,threshold=threshold,window=window)
        return nr

    def evaluation_function(self,node):
        '''
        Return the cost to a node
        '''
        path = node.get_path()
        id1 = path[0]
        dist = 0
        for id2 in path:
            dist = dist + 10*utils.get_dist(self.pdb,id1,id2)
            id1 = id2
        return dist

    def heuristic_function(self,node):
        '''
        Calculate the heuristic cost from node
        '''
        return 10*utils.get_dist(self.pdb,self.target_id,node.get_id())+0.1*np.random.rand()

       
    def load_CM(self):
        if os.path.exists('connection_map.pkl'):
            while True:
                try:
                    f = open('connection_map.pkl','rb')
                    connection_map = pickle.load(f)
                    f.close()
                    break
                except:
                    sleep(1)
        else:
            connection_map = dict()
        return connection_map

    def dump_CM(self,CM1):
        CM2 = self.load_CM()        
        CM = CM_addition(CM1,CM2)
        while True:
            try:
                f = open('connection_map.pkl','wb')
                pickle.dump(CM,f)
                f.close()
                break
            except:
                sleep(1)

    def show_search_path(self):
        if len(self.search_path)>0:
            for node in self.search_path:
                print(node.get_path())
    
    def execute_search(self,mode=0,frequency=1,force=100,time=500,velocity=False,temperature=300,cuda='0',threshold=3,window=1,limit=20):
        count=0
        current_id = self.current_node.get_id()
        while current_id!=self.target_id:
            self.connection_map = CM_addition(self.load_CM(),self.connection_map)
            print(self.current_node.get_cost())
            self.show_search_path()
            print('explored: ',self.explored_id)
            if current_id in self.connection_map.keys():
                successors = self.connection_map[current_id]
            else:
                if count > limit:
                    print('Search failed')
                    return None
                successors = self.get_successors(self.current_node,mode=mode,frequency=frequency,force=force,time=time,velocity=velocity,temperature=temperature,cuda=cuda,threshold=threshold,window=window)
                self.connection_map[current_id] = successors
                count = count + 1
            self.dump_CM(self.connection_map)
            if successors is not None:
                for successor in successors:
                    if successor not in self.explored_id: 
                        print('check:',successor)
                        path_temp = copy.deepcopy(self.current_path)
                        path_temp.append(successor) 
                        node_temp = Node(path_temp)
                        if node_temp.get_id() == self.target_id:
                            self.current_node = node_temp
                            self.search_path.append(node_temp)
                            return self.current_node
                        priority_temp = self.evaluation_function(node_temp) + self.heuristic_function(node_temp)
                        #print(priority_temp)
                        node_temp.update_cost(priority_temp)
                        self.frontier.push(node_temp,priority_temp)
            if self.frontier.empty():
                print('Search failed')
                return None
            else:
                self.current_node = self.frontier.pop()
                self.current_path = self.current_node.get_path()
                self.search_path.append(self.current_node)
                current_id = self.current_node.get_id()
                while(current_id in self.explored_id):
                    if self.frontier.empty():
                        print('Search failed')
                        return None
                    else:
                        self.current_node = self.frontier.pop()
                        self.current_path = self.current_node.get_path()
                        self.search_path.append(self.current_node)
                        current_id = self.current_node.get_id()
                self.explored_id.append(current_id)
        return self.current_node

                
                

