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

        self.cliques = find_cliques(self.jobs, self.job_tools_requirements, self.magazine_capacity)
        self.Q = max(len(clique) for clique in self.cliques)
        print("Q: ", self.Q)

        N = len(self.jobs)
        bins = max(5, min(self.Q, N))
        print("Bins: ", bins)

        k = list(range(1, bins + 1))
        self.bins = k
        self.nodes = list(range(0, len(k) + 3)) # The 3 extra nodes are H at position 0, R at position K+1, and D at position K+2
        self.zbins = list(range(1, len(k) + 2)) # The extra bin is at position K+1 to set at 0 in order to avoid loss of generality
        #print("Nodes: ", self.nodes)
        
        

        # Special nodes
        self.START = 0
        self.REPO = len(k) + 1
        self.DETACHED = len(k) + 2
        self.K = len(k) # Latest bin index

        #print("JOBS:", self.jobs)
        #print("TOOLS:", self.tools)
        #print("START:", self.START)
        #print("REPO:", self.REPO)
        #print("DETACHED:", self.DETACHED)
        #print("K:", self.K)
        #print("K-2: ", self.K-2)
        #print("range K-1:", list(range(1, self.K)))
        #print("bins: ", self.bins)
        
        self.model = gp.Model("JGSMF")
        self.model.setParam("OutputFlag", 0)
        self.model.setParam("TimeLimit", 3600)

        self.x = None # x[i,k]: 1 if job i is in bin k
        self.z = None # z[k]: 1 if bin k has at least one job
        self.f = None # f[a,b,t]: 1 if tool t is transfered from a (H ∪ R ∪ B) and b (R ∪ D ∪ B)

        self.setup_variables()
        self.setup_constraints()
        self.setup_objective()


    def setup_variables(self):

        self.x = self.model.addVars(self.jobs, self.bins, vtype=GRB.BINARY, name="X")
        self.z = self.model.addVars(self.zbins, vtype=GRB.BINARY, name="Z")
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
            self.model.addConstr(gp.quicksum(self.f[k,self.DETACHED,t] for k in self.bins) == 1, name=f"AllToolsWillBeRemoved(3f)")
        
        for t in self.tools:
            for k in self.bins:
                self.model.addConstr(self.f[k-1,k,t] + self.f[self.REPO,k,t] == self.f[k,k+1,t] + self.f[k,self.REPO,t] + self.f[k,self.DETACHED,t], name=f"FlowConservation(3g)")
        
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
        
        for k in self.bins:
            self.model.addConstr(self.magazine_capacity * gp.quicksum(self.f[self.REPO,k-1,t] for t in self.tools) >= gp.quicksum(self.f[self.REPO,k,t] for t in self.tools), name=f"SimmetryBreakingCut(3l)")
        """
        for k in self.bins:
            self.model.addConstr(self.magazine_capacity * gp.quicksum(self.f[k-1,self.DETACHED,t] for t in self.tools) >= gp.quicksum(self.f[k,self.DETACHED,t] for t in self.tools), name=f"SimmetryBreakingCut(3m)")
        
        for t in self.tools:
            for k in self.bins:
                self.model.addConstr(gp.quicksum(self.x[i,k] for i in self.jobs if t in self.job_tools_requirements[i]) >= self.f[k,self.DETACHED,t], name=f"SimmetryBreakingCut(3n)")
        """
        for k in self.bins:
            for i in self.jobs:
                required_tools = self.job_tools_requirements.get(i, [])
                Ti_cardinality = len(required_tools)
                self.model.addConstr(gp.quicksum(self.x[i,b] for b in range(1, k+1)) >= gp.quicksum(self.f[k-1,k,t] for t in required_tools) - Ti_cardinality + 1, name=f"SimmetryBreakingCut(3o)")
        
        self.model.addConstr(self.z[self.K+1] == 0, name=f"SetExtraBinAt0(3extra)")
        for k in self.bins:
            for i in self.jobs:
                self.model.addConstr(self.z[k] >= self.x[i,k], name=f"TighteningConstraint(3p)")
        
        for k in self.bins:
            self.model.addConstr(self.z[k+1] <= gp.quicksum(self.f[k,self.DETACHED,t] for t in self.tools), name=f"TighteningConstraint(3q)")
        
        for k in self.bins:
            self.model.addConstr(self.z[k] >= self.z[k+1], name=f"TighteningConstraint(3r)")
        
        for l in self.cliques:
            for k in self.bins:
                self.model.addConstr(self.z[k] >= gp.quicksum(self.x[i,k] for i in l), name=f"TighteningCliquesConstraint(3s)")

        for k in range(1, self.Q+1):
            self.model.addConstr(self.z[k] == 1, name=f"Set1BinsUntilQ(3t)")

            

    def setup_objective(self):
        
        #with range() self.K is already K-1
        objective = gp.quicksum(self.f[k,self.DETACHED,t] for t in self.tools for k in range(1, self.K)) + gp.quicksum(self.f[k,self.REPO,t] for t in self.tools for k in range(1, self.K-1))
        self.model.setObjective(objective, GRB.MINIMIZE)

    
    def optimize(self):

        self.model.optimize()

        if self.model.status == GRB.OPTIMAL:
            print("Optimal switches:", int(self.model.objVal))
        elif self.model.status == GRB.INFEASIBLE:
            print("Model is infeasible")
            self.model.setParam("TimeLimit", 60)
            self.model.computeIIS()
            self.model.write("infeasibility_report.ilp")
        else:
            print("Optimization ended with status", self.model.status)

        print("Detached tools: ")
        print(self.print_detached_tools())
        print("Active bins: ")
        print(self.active_bins())
        print("Jobs position: ")
        print(self.jobs_position())
        print("Jobs in bins: ")
        print(self.print_jobs_in_bins())
        print("Tools flow: ")
        print(self.tools_flow())

    
    def get_solution(self):

        if self.model.status == gp.GRB.OPTIMAL:
            job_order = sorted((k, i) for i in self.jobs for k in self.bins if self.x[i, k].x > 0.5)
            job_order = [i for k, i in job_order]
            print("Job order:", job_order)
            bins_used = set(k for i in self.jobs for k in self.bins if self.x[i, k].x > 0.5)
            num_bins_used = len(bins_used)
            print("Bins used:", bins_used)
            print("Number of bins used:", num_bins_used)
            return job_order
        else:
            print("No optimal solution found.")
            return None
        
    def print_detached_tools(self):
        if self.model.status == gp.GRB.OPTIMAL:
            for k in self.bins:
                for t in self.tools:
                    if self.f[k, self.DETACHED, t].x > 0.5:
                        print(f"Tool {t} is sent to detached from bin {k}")

    def print_jobs_in_bins(self):
        if self.model.status == gp.GRB.OPTIMAL:
            bins_jobs = {k: [] for k in self.bins}
            for k in self.bins:
                for i in self.jobs:
                    if self.x[i, k].x > 0.5:
                        bins_jobs[k].append(i)
            for k, jobs in bins_jobs.items():
                print(f"Bin {k}: Jobs {jobs}")
        else:
            print("No optimal solution found.")

    def active_bins(self):
        if self.model.status == gp.GRB.OPTIMAL:
            bins = [k for k in self.bins if self.z[k].x > 0.5]
            return bins
        else:
            print("No optimal solution found.")
            return None
        
    def jobs_position(self):
        if self.model.status == gp.GRB.OPTIMAL:
            for i in self.jobs:
                for k in self.bins:
                    if self.x[i, k].x > 0.5:
                        print(f"Job {i} is in bin {k}")
        else:
            print("No optimal solution found.")

    def tools_flow(self):
        if self.model.status == gp.GRB.OPTIMAL:
            for k in self.nodes:
                for j in self.nodes:
                    for t in self.tools:
                        if self.f[k, j, t].x > 0.5:
                            print(f"Tool {t} is transfered from bin {k} to bin {j}")
        else:
            print("No optimal solution found.")

        
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

