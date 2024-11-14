import gurobipy as gp
from gurobipy import GRB

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
        k = [1, 2, 3, 4, 5, 6] #placeholder for number of bins
        self.bins = k
        self.nodes = list(range(0, len(k) + 3)) # The 3 extra nodes are H at position 0, R at position k+1, and D at position k+2
        print("nodes:", self.nodes)

        # Special nodes
        self.START = 0
        self.REPO = len(k) + 1
        self.DETACHED = len(k) + 2

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
        #print("f:", self.f)


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


    def setup_objective(self):
        
        #with range() len(self.bins) is already K-1
        objective = gp.quicksum(self.f[k,self.DETACHED,t] for t in self.tools for k in range(1, len(self.bins))) + gp.quicksum(self.f[k,self.REPO,t] for t in self.tools for k in range(1, len(self.bins)-1))
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
        

# Example usage
jobs = [1, 2, 3, 4]
tools = [1, 2]
job_tools_requirements = {1: [1, 2], 2: [2], 3: [1], 4: [1, 2]}
magazine_capacity = 2

model = JGSMFModel(jobs, tools, magazine_capacity, job_tools_requirements)
model.optimize()
print(model.get_solution())
