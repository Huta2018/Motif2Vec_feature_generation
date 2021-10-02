# Here I try to write a class to represent the motiff and its connections
class PrincipalComponent(object): 
    Graph = {}
    # The class "constructor"
    def __init__(self):
        self.Graph = {}
    
    # Given a number and compound symbol - create a prefix (row header)
    # Example: 3, H20  -----> H20-3-dimension
    @staticmethod
    
    def Dimension_node(self, n, prefix):
         if n == 0:
             Dimension = 'no_sharing'
         elif n ==1:
             Dimension = 'Corner'
         elif  n == 2:
             Dimension = 'Edge'
         elif n >= 3:
             Dimension = 'Face'
         return str(prefix) + '_' + str(Dimension)
     
    def ExtendNode(self, SubGraph, Vertex, Vertex_Key):
        if not self.Dimension_node(self, Vertex, Vertex_Key) in SubGraph:
            SubGraph[self.Dimension_node(self,Vertex, Vertex_Key )] =1
        else:
             SubGraph[self.Dimension_node(self, Vertex, Vertex_Key)] = SubGraph[self.Dimension_node(self, Vertex, Vertex_Key)] +1
        return SubGraph
    
    
    def AddANode(self, Vertex_Key, Connection, Connection_Key, new_list):
        Vertex_Key = new_list 
        Connection_Key = Connection_Key 
        if not Vertex_Key in self.Graph:
            self.Graph[Vertex_Key] = {}
        self.Graph[Vertex_Key] = self.ExtendNode(self.Graph[Vertex_Key], Connection, Connection_Key)
        return self.Graph
