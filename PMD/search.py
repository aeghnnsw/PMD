import heapq
import copy
from PMD import utils as utils
from PMD import analysis as ana
from PMD import nm


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
    
    def get_id(self):
        return self.id
    
    def get_path(self):
        return self.path


class Astar:
    def __init__(self,init_id,target_id,pdb_file):
        self.nm = nm.nm(pdb_file)
        self.init_id = init_id
        self.target_id = target_id
        self.explored_id = [init_id]
        self.connection_map = dict()
        self.current_path = [init_id]
        self.init_node = Node(self.current_path)
        self.current_node = Node(self.current_path)
        self.frontier = PriorityQueue()
        
    def get_successors(self,node):
        res_id = node.id
        nr = self.nm.pump(res_id)
        return nr

    def evaluation_function(self,node):
        '''
        Return the cost to a node
        '''
        return len(node.get_path())

    def heuristic_function(self,node):
        '''
        Calculate the heuristic cost from node
        '''
        return 0
    
    def execute_search(self):
        current_id = self.current_node.get_id()
        while current_id!=self.target_id:
            if current_id in self.connection_map.keys():
                successors = self.connection_map[current_id]
            else:
                successors = self.get_successors(current_node)
                self.connection_map[current_id] = successors
            if successors is not None:
                for successor in successors:
                    if successor not self.explored_id:
                        path_temp = copy.deepcopy(self.current_path)
                        path_temp.append(successor) 
                        node_temp = Node(path_temp)
                        priority_temp = self.evaluation_function(node_temp) + self.heuristic_function(node_temp)
                        self.frontier.push(node_temp,priority_temp)
            if self.frontier.empty():
                print('Search failed')
                return None
            else:
                self.explored_id.append(current_id)
                self.current_node = self.frontier.pop()
                self.current_path = self.current_node.get_path()
                current_id = self.current_node.get_id()
        return current_node

                
                

