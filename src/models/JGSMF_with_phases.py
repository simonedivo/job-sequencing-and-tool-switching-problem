import gurobipy as gp
from gurobipy import GRB
import networkx as nx
from datetime import datetime

class JGSMFModel:
    def __init__(self, jobs, tools, magazine_capacity, job_tools_requirements, bins, time_limit, isPhase3):

        """
        jobs: list of job IDs, from 1 to N
        tools: list of tool IDs, from 1 to M
        magazine_capacity: integer, capacity of the tool magazine
        job_tools_requirements: dict, key is job ID, value is list of required tools for that job
        bins: list of bin IDs, from 1 to K
        """

        self.jobs = jobs
        self.tools = tools
        self.magazine_capacity = magazine_capacity
        self.job_tools_requirements = job_tools_requirements

        self.bins = bins
        self.nodes = list(range(0, len(bins) + 3)) # The 3 extra nodes are H at position 0, R at position K+1, and D at position K+2
        
        self.cliques = find_cliques(self.jobs, self.job_tools_requirements, self.magazine_capacity)
        self.Q = max(len(clique) for clique in self.cliques)

        # Special nodes
        self.START = 0
        self.REPO = len(bins) + 1
        self.DETACHED = len(bins) + 2
        self.K = len(bins) # Latest bin index

        #print("START:", self.START)
        #print("REPO:", self.REPO)
        #print("DETACHED:", self.DETACHED)
        #print("K:", self.K)
        #print("K-2: ", self.K-2)
        #print("range K-1:", list(range(1, self.K)))
        #print("bins: ", self.bins)
        
        self.model = gp.Model("JGSMF")
        self.model.setParam("OutputFlag", 0)
        self.model.setParam("TimeLimit", time_limit)

        self.x = None # x[i,k]: 1 if job i is in bin k
        self.z = None # z[k]: 1 if bin k has at least one job
        self.f = None # f[a,b,t]: 1 if tool t is transfered from a (H ∪ R ∪ B) and b (R ∪ D ∪ B)

        self.setup_variables()
        self.setup_constraints()
        if isPhase3:
            self.setup_phase3_constraints()
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
            self.model.addConstr(gp.quicksum(self.f[k,self.DETACHED,t] for k in self.bins) == 1, name=f"AllToolsWillBeRemoved(3f)")
        
        # Doppio f uguale alla fine? -> il doppione è un errore! Uno dei due deve essere f[k, self.REPO, t]!
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
        
        #Node 0 is the start node
        for k in self.bins:
            if k > 1:
                self.model.addConstr(self.magazine_capacity * gp.quicksum(self.f[self.REPO,k-1,t] for t in self.tools) >= gp.quicksum(self.f[self.REPO,k,t] for t in self.tools), name=f"SimmetryBreakingCut(3l)")
        
        #Problem with this constraint, added k < self.K because, otherwise, it would conflict with constraint 3i
        for k in self.bins:
            if k > 1 and k < self.K:
                self.model.addConstr(self.magazine_capacity * gp.quicksum(self.f[k-1,self.DETACHED,t] for t in self.tools) >= gp.quicksum(self.f[k,self.DETACHED,t] for t in self.tools), name=f"SimmetryBreakingCut(3m)")
        
        #for t in self.tools:
        #    for k in self.bins:
        #        self.model.addConstr(gp.quicksum(self.x[i,k] for i in self.jobs if t in self.job_tools_requirements[i]) >= self.f[k,self.DETACHED,t], name=f"SimmetryBreakingCut(3n)")
        
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


    def setup_phase3_constraints(self):

        for k in range(1, self.K-1):
            self.model.addConstr(gp.quicksum(self.f[k,self.REPO,t] + self.f[k,self.DETACHED,t] for t in self.tools) == 1, name=f"AdditionalConstraint(4a)")

        self.model.addConstr(gp.quicksum(self.f[self.K-1,self.DETACHED,t] for t in self.tools) == 1, name=f"AdditionalConstraint(4b)")

        for k in self.bins:
            self.model.addConstr(self.z[k] == 1, name=f"AdditionalConstraint(4c)")
            
    
    def setup_objective(self):
        
        #with range() self.K is already K-1
        objective = gp.quicksum(self.f[k,self.DETACHED,t] for t in self.tools for k in range(1, self.K)) + gp.quicksum(self.f[k,self.REPO,t] for t in self.tools for k in range(1, self.K-1))
        self.model.setObjective(objective, GRB.MINIMIZE)

    def optimize(self):

        self.model.optimize()
        return self.model.status

    
    def get_solution(self):

        if self.model.status == gp.GRB.OPTIMAL:
            job_order = sorted((k, i) for i in self.jobs for k in self.bins if self.x[i, k].x > 0.5)
            job_order = [i for k, i in job_order]
            #print("Job order:", job_order)
            #bins_used = set(k for i in self.jobs for k in self.bins if self.x[i, k].x > 0.5)
            #num_bins_used = len(bins_used)
            #print("Bins used:", bins_used)
            #print("Number of bins used:", num_bins_used)
            return job_order
        else:
            #print("No optimal solution found.")
            return None
        
    def count_switches(self):
        return int(self.model.ObjVal)
    
        
        
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


def phase_1(jobs, tools, magazine_capacity, job_tools_requirements, time_limit):
    cliques = find_cliques(jobs, job_tools_requirements, magazine_capacity)
    Q = max(len(clique) for clique in cliques)
    N = len(jobs)
    K1 = max(5, min(N, Q))
    T1 = time_limit
    start_time = datetime.now()

    while T1 > 0:
        print(f"Trying Phase 1 with {K1} bins...")
        elapsed_time = (datetime.now() - start_time).total_seconds()
        T1 = round(T1 - elapsed_time, 3)

        bins = list(range(1, K1+1))
        model = JGSMFModel(jobs, tools, magazine_capacity, job_tools_requirements, bins, T1, False)
        status = model.optimize()

        if status == GRB.SUBOPTIMAL or status == GRB.OPTIMAL:
            print(f"Phase 1 found a solution with {K1} bins and {model.count_switches()} switches.")
            return K1, T1, model.get_solution(), model.count_switches()
        elif status == GRB.INFEASIBLE:
            K1 += 2

    print("Time limit reached in phase 1.")
    return None, None, None, None


def phase_2(jobs, tools, magazine_capacity, job_tools_requirements, T1, K1, S1):
    K2 = K1 + 1
    T2 = (T1) / 2
    best_solution = None
    best_switches = S1
    start_time = datetime.now()

    while T2 > 0:
        print(f"Trying Phase 2 with {K2} bins...")
        elapsed_time = (datetime.now() - start_time).total_seconds()
        T2 = round(T2 - elapsed_time, 3)

        bins = list(range(1, K2+1))
        model = JGSMFModel(jobs, tools, magazine_capacity, job_tools_requirements, bins, T2, False)
        status = model.optimize()

        
        if status == GRB.SUBOPTIMAL or status == GRB.OPTIMAL:
            best_solution = model.get_solution()
            current_switches = model.count_switches()
            if current_switches < best_switches:
                best_switches = current_switches
                K2 += 1
            if current_switches >= best_switches:
                break
        elif status == GRB.INFEASIBLE:
            break
    
    if T2 <= 0:
        print("Time limit reached in phase 2")
        return K2, best_solution, best_switches
    else:
        print(f"Found less or equal switches in phase 2 with {K2} bins and {best_switches} switches.")
        return K2, best_solution, best_switches

    
    

def phase_3(jobs, tools, magazine_capacity, job_tools_requirements, T1, S3, job_order2, S2):
    K3 = S3
    T3 = (T1)/2
    best_solution = job_order2
    start_time = datetime.now()
    best_switches = S2

    while T3 > 0:
        print(f"Checking in phase 3...")
        elapsed_time = (datetime.now() - start_time).total_seconds()
        T3 = round(T3 - elapsed_time, 3)

        bins = list(range(1, K3+1))
        model = JGSMFModel(jobs, tools, magazine_capacity, job_tools_requirements, bins, T3, True)
        status = model.optimize()

        if status == GRB.INFEASIBLE:
            return K3, best_solution, best_switches
        if status == GRB.SUBOPTIMAL or status == GRB.OPTIMAL:
            K3 -= 1
            best_solution = model.get_solution()
            best_switches = model.count_switches()
    
    return None, None


def solve_with_phases(jobs, tools, magazine_capacity, job_tools_requirements, time_limit):

    K1, T1, job_order1, S1 = phase_1(jobs, tools, magazine_capacity, job_tools_requirements, time_limit)
    
    if job_order1 is None:
        return None

    K2, job_order2, S2 = phase_2(jobs, tools, magazine_capacity, job_tools_requirements, T1, K1, S1)

    K3 = min(S1, S2)
    K3, job_order3, S3 = phase_3(jobs, tools, magazine_capacity, job_tools_requirements, T1, K3, job_order2, S2)
    if job_order3 is None:
        print("No feasible solution found.")
        return None, None
    else:
        print("Solution found.")
        return job_order3, S3
    

    
# Example usage
#jobs = [1, 2, 3, 4, 5, 6, 7, 8, 9]
#tools = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
#magazine_capacity = 15
#job_tools_requirements = {1: [3, 4, 10, 11, 16, 18], 2: [5, 9, 10, 16, 17, 19, 20], 3: [3, 4, 5, 6, 12, 13, 15, 16, 17, 19], 4: [3, 7, 11, 12], 5: [4, 9, 11, 12, 15, 16], 6: [1, 2, 3, 4, 7, 9, 14, 15, 16, 19], 7: [4, 8, 12, 13, 14, 17, 18, 19], 8: [1, 3, 11, 12], 9: [5, 10, 11, 16, 18]}

#Tabela 4 number 84
#jobs = [1, 2, 3, 4, 5, 6, 7, 8, 9]
#tools = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
#magazine_capacity = 15
#job_tools_requirements = {1: [3, 4, 7, 9, 10, 17, 18, 19, 20, 22], 2: [1, 3, 4, 5, 7, 8, 14, 15, 19, 24, 25], 3: [5, 8, 9, 12, 14, 16, 19, 21, 23, 24, 25], 4: [1, 5, 10, 11, 14, 16, 17, 18, 20, 23, 24, 25], 5: [1, 2, 4, 5, 8, 10, 12, 14, 15, 16, 17, 21, 22, 23, 25], 6: [1, 4, 7, 11, 12, 13, 14, 15, 17, 18, 21, 22, 24, 25], 7: [5, 6, 7, 8, 9, 10, 11, 15, 17, 18, 19, 20, 23, 24, 25], 8: [1, 2, 3, 4, 5, 13, 14, 17, 19, 20, 21, 22, 24, 25], 9: [1, 5, 9, 11, 12, 13, 14, 15, 16, 18, 21, 22, 24, 25]}

#solution = solve_with_phases(jobs, tools, magazine_capacity, job_tools_requirements, 3600)
#print(solution)