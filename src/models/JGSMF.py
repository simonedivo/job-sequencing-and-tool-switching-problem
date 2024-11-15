import gurobipy as gp
from gurobipy import GRB
import networkx as nx

class JGSMFModel:
    def __init__(self, jobs, tools, magazine_capacity, job_tools_requirements):

        """
        jobs: list of job IDs, from 1 to N
        tools: list of tool IDs, from 1 to M
        magazine_capacity: integer, capacity of the tool magazine
        job_tools_requirements: dict, key is job ID, value is list of required tools for that job
        """

        self.jobs = jobs
        self.tools = tools
        self.magazine_capacity = magazine_capacity
        self.job_tools_requirements = job_tools_requirements
        k = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14] #placeholder for number of bins @TODO: implement dinamic bin number allocations
        self.bins = k
        self.nodes = list(range(0, len(k) + 3)) # The 3 extra nodes are H at position 0, R at position k+1, and D at position k+2
        
        self.cliques = find_cliques(jobs, job_tools_requirements, magazine_capacity)
        #print("Cliques:", self.cliques)
        self.Q = max(len(clique) for clique in self.cliques)

        # Special nodes
        self.START = 0
        self.REPO = len(k) + 1
        self.DETACHED = len(k) + 2
        self.K = len(k) # Latest bin index

        #print("START:", self.START)
        #print("REPO:", self.REPO)
        #print("DETACHED:", self.DETACHED)
        #print("K:", self.K)
        #print("K-2: ", self.K-2)
        #print("range K-1:", list(range(1, self.K)))
        #print("bins: ", self.bins)

        self.model = gp.Model("JGSMF")

        self.x = None # x[i,k]: 1 if job i is in bin k
        self.z = None # z[k]: 1 if 
        self.f = None # f[a,b,t]: 1 if tool t is transfered from a (H ∪ R ∪ B) and b (R ∪ D ∪ B)

        self.setup_variables()
        self.setup_constraints()
        self.setup_objective()


    def setup_variables(self):

        self.x = self.model.addVars(self.jobs, self.bins, vtype=GRB.BINARY, name="X")
        self.z = self.model.addVars(self.bins, vtype=GRB.BINARY, name="Z")
        self.f = self.model.addVars(self.nodes, self.nodes, self.tools, vtype=GRB.BINARY, name="F")


    def setup_constraints(self):

        for i in self.jobs:
            self.model.addConstr(gp.quicksum(self.x[i,k] for k in self.bins) == 1, name=f"Job{i}InExactlyOneBin(3b)")

        for k in self.bins:
            self.model.addConstr(gp.quicksum(self.f[k-1,k,t] for t in self.tools) == self.magazine_capacity, name=f"CapacityLimit(3c)")

        for i in self.jobs:
            required_tools = self.job_tools_requirements.get(i, [])
            for k in self.bins:
                for t in required_tools:
                    self.model.addConstr(self.x[i,k] <= self.f[k-1,k,t], name=f"AllNecessaryToolsAvailable(3d)")

        for t in self.tools:
            self.model.addConstr(self.f[self.START, 1, t] + self.f[self.START, self.REPO, t] == 1, name=f"AllToolsAvailableAtStart(3e)")

        for t in self.tools:
            self.model.addConstr(gp.quicksum(self.f[k,self.DETACHED,t] for k in self.bins) == 1, name=f"AllToolsRemovedInTheEnd(3f)")

        # Doppio f uguale alla fine? -> il doppione è un errore!
        for t in self.tools:
            for k in self.bins:
                self.model.addConstr(self.f[k-1,k,t] + self.f[self.REPO,k,t] == self.f[k,k+1,t] + self.f[k,self.DETACHED,t], name=f"FlowConservation(3g)")

        for t in self.tools:
            self.model.addConstr(self.f[self.K-2,self.K-1,t] + self.f[self.REPO,self.K-1,t] == self.f[self.K-1,self.K,t] + self.f[self.K-1,self.DETACHED,t], name=f"FlowConservation(3h)")

        for t in self.tools:
            self.model.addConstr(self.f[self.K-1,self.K,t] == self.f[self.K,self.DETACHED,t], name=f"FlowConservation(3i)")

        # in range() self.K is already K-1 -> self.K-1 within range() is K-2
        for t in self.tools:
            self.model.addConstr(self.f[self.START,self.REPO,t] + gp.quicksum(self.f[k,self.REPO,t] for k in range(1, self.K-1)) == gp.quicksum(self.f[self.REPO,k,t] for k in range(1,self.K)), name=f"FlowConservation(3j)")

        for i in self.jobs:
            required_tools = self.job_tools_requirements.get(i, [])    
            for k in self.bins:
                self.model.addConstr(self.x[i,k] <= gp.quicksum(self.f[self.REPO,k-1,t] for t in required_tools), name=f"SimmetryBreakingCut(3k)")

        #Node 0 is the start node
        for k in self.bins:
            if k > 1:
                self.model.addConstr(self.magazine_capacity * gp.quicksum(self.f[self.REPO,k-1,t] for t in self.tools) >= gp.quicksum(self.f[self.REPO,k,t] for t in self.tools), name=f"SimmetryBreakingCut(3l)")
        
        #Problem with this constraint
        for k in self.bins:
            if k > 1:
                self.model.addConstr(self.magazine_capacity * gp.quicksum(self.f[k-1,self.DETACHED,t] for t in self.tools) >= gp.quicksum(self.f[k,self.DETACHED,t] for t in self.tools), name=f"SimmetryBreakingCut(3m)")

        for t in self.tools:
            for k in self.bins:
                self.model.addConstr(gp.quicksum(self.x[i,k] for i in self.jobs if t in self.job_tools_requirements[i]) >= self.f[k,self.DETACHED,t], name=f"SimmetryBreakingCut(3n)")
        
        for k in self.bins:
            for i in self.jobs:
                required_tools = self.job_tools_requirements.get(i, [])
                Ti_cardinality = len(required_tools)
                self.model.addConstr(gp.quicksum(self.x[i,b] for b in range(1, self.K+1)) >= gp.quicksum(self.f[k-1,k,t] for t in required_tools) - Ti_cardinality + 1, name=f"SimmetryBreakingCut(3o)")

        for k in self.bins:
            for i in self.jobs:
                self.model.addConstr(self.z[k] >= self.x[i,k], name=f"TighteningConstraint(3p)")

        #Added if because it was going out of range
        for k in self.bins:
            if k < self.K:
                self.model.addConstr(self.z[k+1] <= gp.quicksum(self.f[k,self.DETACHED,t] for t in self.tools), name=f"TighteningConstraint(3q)")

        for k in self.bins:
            if k < self.K:
                self.model.addConstr(self.z[k] >= self.z[k+1], name=f"TighteningConstraint(3r)")

        for l in self.cliques:
            for k in self.bins:
                self.model.addConstr(self.z[k] >= gp.quicksum(self.x[i,k] for i in l), name=f"TighteningCliquesConstraint(3s)")

        for k in range(1, self.Q+1):
            self.model.addConstr(self.z[k] == 1, name=f"Set1BinsUntilQ(3t)")   

        #Constraints added to ensure flow correctness based on my logic of nodes (I added them)
        for t in self.tools:
            self.model.addConstr(self.f[self.START, self.DETACHED, t] == 0, name=f"NoFlowFromStartToDetached({t})")

        for t in self.tools:
            self.model.addConstr(gp.quicksum(self.f[k, self.START, t] for k in range(1, self.K+1)) == 0, name=f"NoFlowIntoStart({t})")

        for t in self.tools:
            self.model.addConstr(self.f[self.REPO, self.DETACHED, t] == 0, name=f"NoFlowFromRepoToDetached({t})")
            self.model.addConstr(self.f[self.REPO, self.START, t] == 0, name=f"NoFlowFromRepoToStart({t})")
        
        for t in self.tools:
            self.model.addConstr(gp.quicksum(self.f[self.DETACHED, k, t] for k in range(self.K+1)) == 0, name=f"NoFlowOutOfDetached({t})")


    def setup_objective(self):
        
        #with range() self.K is already K-1
        objective = gp.quicksum(self.f[k,self.DETACHED,t] for t in self.tools for k in range(1, self.K)) + gp.quicksum(self.f[k,self.REPO,t] for t in self.tools for k in range(1, self.K-1))
        self.model.setObjective(objective, GRB.MINIMIZE)

    
    def optimize(self):

        self.model.optimize()

        if self.model.status == GRB.OPTIMAL:
            print("Optimal objective:", self.model.objVal)
        elif self.model.status == GRB.INFEASIBLE:
            print("Model is infeasible")
            #self.model.computeIIS()
            #self.model.write("infeasibility_report.ilp")
        else:
            print("Optimization ended with status", self.model.status)

    
    def get_solution(self):

        if self.model.status == gp.GRB.OPTIMAL:
            job_order = sorted((k, i) for i in self.jobs for k in self.bins if self.x[i, k].x > 0.5)
            job_order = [i for k, i in job_order]
            #print("Job order:", job_order)
            return job_order
        else:
            print("No optimal solution found.")
            return None
        
def find_cliques(jobs, job_tools_requirements, magazine_capacity):

    G = nx.Graph()
    G.add_nodes_from(jobs)

    for i in jobs:
        for j in jobs:
            if i < j:
                combined_tools = set(job_tools_requirements[i]).union(set(job_tools_requirements[j]))
                if len(combined_tools) > magazine_capacity:
                    G.add_edge(i, j)

    cliques = list(nx.find_cliques(G))
    return cliques
        

# Example usage
jobs = [1, 2, 3, 4, 5, 6]
tools = [1, 2, 3, 4]
job_tools_requirements = {1: [1, 2], 2: [2], 3: [1, 2, 4], 4: [1, 2], 5: [1, 3, 4], 6: [2, 3]}
magazine_capacity = 3

model = JGSMFModel(jobs, tools, magazine_capacity, job_tools_requirements)
model.optimize()
print(model.get_solution())
