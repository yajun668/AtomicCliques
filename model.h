#ifndef model_h
#define model_h
#endif /* model_h */
#include "common.h"
#include "graph.h"
#include "gurobi_c++.h"
#include <map>

//store solution results
struct ACPResult
{
    string g_name="";
    string opt_status;
    double ub = -1;
    double obj = -1;
    double wallTime = 0.0;
    double readTime = 0.0; //reading graph time
    double peelingTime = 0.0; //edge peeling time
    double heurTime = 0.0; //heuristic time
    double corePeelTime = 0.0; //core peeling time
    long heurSize = 0;
    double gurobiSolveTime = 0.0;
    long num_nodes = 0; //number of universal vertics
    long num_graphs = 0; //number of graph collections
    long num_edges =0; //sum of edges across graph collections
    long new_nodes =0; //after edge peeling, number of vertices in auxiliary graph
    long new_num_edges =0;//after edge peeling, number of edges in auxiliary graph
    long num_components =0; //after edge peeling, number of components in auxiliary graph
    vector<long> bestSol;
};
typedef vector<ACPResult*> ACPRs;

void SolveBatch(string graph_list, long solver_type, ACPRs &results); //solve instances in batch
ACPResult* SolveACP_Extend(string file_name); // Atomic Clique Extended formulation with z
ACPResult* SolveACP_Monolithic(string file_name); //Monolithic Atomic Clique Formulation
ACPResult* SolveACP_EdgePeeling(string file_name); //Solve MACP using edge peeling
void ClearResults(ACPRs &results); //clear results
GRBLinExpr GRBSum(int n, GRBVar* vars, double* coeff = NULL);
GRBVar* GRBVarArray(int n, GRBModel& m_val, double lb = 0, double ub = GRB_INFINITY, char type = GRB_CONTINUOUS, double obj_coeff = 1);
vector<long> ACPHeuristic(const GSeq & gs_val);//Get Heuristic
vector<long> ACPHeuristic_edgeOverlap(const GSeq & gs_val, const VSet &v0);//Get Heuristic based on edge overlap