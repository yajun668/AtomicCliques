#include "model.h"
void SolveBatch(string graph_list, long solver_type, ACPRs &results)
{
    //this function is to solve instances in batch
    results.clear();
    string str;
    ifstream fs(graph_list);
    if (!fs)
    {
        cout << "Cannot open input file!" << endl;
        return;
    }

    while (!fs.eof())
    {
        getline(fs, str);
        if (str.size() < 1)
            continue;

        ACPResult* r = NULL;
        switch (solver_type) {
            case 1:
                r = SolveACP_Extend(str); //enhanced formulation
                break;
            case 2:
                r = SolveACP_EdgePeeling(str); //Edge Peeling + Gurobi Solver
                break;
            default:
                r = SolveACP_Monolithic(str); //monolithic formulation
                break;
        }
        results.push_back(r);
    }
    fs.close();
}

//****************************Monolithic Atomic Clique Formulation*************************
ACPResult* SolveACP_Monolithic(string file_name)
{
    // Load graph sequence:
    cout << "=====================================================================================" << endl;
    cout << "Loading "<< file_name << " ..." << endl;
    GSeq gs;
    VSet v0;
    ACPResult* r = new ACPResult(); // get results
    chrono_time_point read_time = chrono_clock::now();
    LoadGSeq(gs, v0, file_name);
    chrono::duration<double>readTime_span = chrono::duration_cast<std::chrono::duration<double> >(chrono_clock::now() - read_time);
    r->readTime = readTime_span.count();
    cout << "Loaded!" << endl << endl;
    // Get the wallclock start time
    chrono_time_point wall_time = chrono_clock::now();

    // Solve APC
    cout << "Solving the MAX AC problem using Monolithic Formulation ..." << endl;
    long n_nodes = v0.V.size();
    long p = gs.size();
    long u, v, k, t;
    Graph* gk = NULL;

    // max atomic clique
    //Gurobi env setup
    GRBEnv* env = 0;
    try {
        env = new GRBEnv();
        GRBModel m = GRBModel(*env);
        GRBVar *x = GRBVarArray(n_nodes, m, 0.0, 1.0, GRB_BINARY); // Create variables
        m.setObjective(GRBSum(n_nodes, x), GRB_MAXIMIZE); //set up objective function

        //add constraints
        vector<bool> isInGk;
        for (k = 0; k < p; k++) {
            gk = gs[k]; // for each graph
            long n1 = gk->adj.size(); // number of vertices in graph gk
            isInGk = vector<bool>(n_nodes, false); //used to check if a vertex  exist in graph Gk
            Adj::iterator i, j;// i,j->first are used to retrieve vertex set of graph gk

            //add constraint x_u + x_v <=1 for all u,v are NOT adjacent in graph G_k, k∈[p]
            for (i = gk->adj.begin(); i != gk->adj.end(); i++) {
                u = v0.v2idx[i->first]; //convert original vertex set to 0, 1,.., n-1 for Gurobi constraint index
                isInGk[u] = true;
                for (j = next(i); j != gk->adj.end(); j++) {
                    if (!gk->IsAdj(i->first, j->first)) {
                        v = v0.v2idx[j->first];
                        m.addConstr(x[u] + x[v] <= 1, "Clique constraints_" + itos(i->first) + "_" + itos(j->first));
                    }
                }
            }

            // add constraint: x_u + x_v <=1, for all u∈V(Gk),v∈V0 \V(Gk),k∈[p]
            for (i = gk->adj.begin(); i != gk->adj.end(); i++) {
                u = v0.v2idx[i->first]; //convert original vertex set to 0,1,.., n-1 for Gurobi constraint index
                for (v = 0; v < n_nodes; v++) {
                    if (!isInGk[v]) {
                        m.addConstr(x[u] + x[v] <= 1, "Atomicity constraints_" + itos(i->first) + "_" + itos(v0.V[v]));
                    }
                }
            }
        }

        //Set maximum time limit
        m.getEnv().set(GRB_DoubleParam_TimeLimit,3600);
        //Set Gurobi screen display flag: 0=switch off; 1=default
        m.getEnv().set(GRB_IntParam_OutputFlag,1);
        m.update();

        //m.write("formulation.lp");// for verification by writing constraints
        m.optimize(); // solve
        r->num_nodes = n_nodes;
        r->num_graphs = p;
        r->gurobiSolveTime = m.get(GRB_DoubleAttr_Runtime);
        r->g_name = file_name;
        r->opt_status = "Model Problematic";
        if (m.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            r->opt_status = "Model Solved";
        } else if (m.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
            r->opt_status = "Hit Time Limit";
        }

        r->ub = m.get(GRB_DoubleAttr_ObjBound);
        r->obj = m.get(GRB_DoubleAttr_ObjVal);
        for (t = 0; t < n_nodes; t++) {
            if (x[t].get(GRB_DoubleAttr_X) >= 0.6) //>=0.6 -> 1 since it's binary
                r->bestSol.push_back(v0.V[t]);
        }

        sort(r->bestSol.begin(), r->bestSol.end()); //sort the best solution
        //count the wall clock time of solving the model
        chrono::duration<double> WallTime_span = chrono::duration_cast<std::chrono::duration<double> >(
                chrono_clock::now() - wall_time);
        r->wallTime = WallTime_span.count();
        // Clean up program resources
        ClearGSeq(gs);
        delete[] x;
    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }

    delete env;

    cout << r->opt_status << endl;
    //write result into an Excel file for each graph
    ofstream fout;
    string outputFile = "max_atomic_clique_monolithic.csv";
    fout.open(outputFile.c_str(), ios::app);
    fout<<r->g_name<<","<<r->num_nodes<<","<<r->num_graphs<<","<<r->readTime<<","<<r->gurobiSolveTime<<","<<r->wallTime<<","<<r->obj<<","<<r->ub<<","<<r->opt_status<<","<<(r->ub - r->obj)/r->obj*100;
    fout<<"\n";
    fout.close();

    // Output Solution
    cout << "----------Output in Graph Collection: " <<r->g_name<<"-------------------"<< endl;
    cout << "Objective:\t" << r->obj << endl;
    cout << "Bound:\t" << r->ub << endl;
    cout << "The best solution:" << endl;
    for (t = 0; t < r->obj; t++) {
        cout << r->bestSol[t] << "-->";
    }
    cout << endl;
    cout << "Total Time : " << r->wallTime << "s " << endl;
    cout << "=====================================================================================" << endl << endl;
    return r;
}


//******************Atomic Clique Enhanced formulation***********************************************
ACPResult* SolveACP_Extend(string file_name)
{
    // Load graph sequence:
    cout << "=====================================================================================" << endl;
    cout << "Loading "<< file_name << " ..." << endl;
    GSeq gs;
    VSet v0;
    ACPResult* r = new ACPResult(); // get results
    chrono_time_point read_time = chrono_clock::now();
    LoadGSeq(gs, v0, file_name);
    chrono::duration<double>readTime_span = chrono::duration_cast<std::chrono::duration<double> >(chrono_clock::now() - read_time);
    r->readTime = readTime_span.count();
    cout << "Loaded!" << endl << endl;

    long n_nodes = v0.V.size();
    long p = gs.size();
    long u, v, k, t = -1;
    Graph* gk = NULL;
    for (long i1 =0; i1<p; i1++) {
        r->num_edges += gs[i1]->num_edges;
    }

    // Solve APC
    cout << "Solving the MAX AC problem using Extended formulation ..." << endl;
    // Get the start time
    chrono_time_point wall_time = chrono_clock::now();

    // Create variables
    //Gurobi env setup
    GRBEnv *env;
    try {
        env = new GRBEnv();
        GRBModel m = GRBModel(*env);
        GRBVar *x = GRBVarArray(n_nodes, m, 0.0, 1.0, GRB_BINARY);
        GRBVar *z = GRBVarArray(p, m, 0.0, 1.0, GRB_BINARY);
        GRBVar **y = new GRBVar *[p];

        vector<bool> isInGk(n_nodes, false); //used to check a vertex if exist in graph Gk
        Adj::iterator i1, j1;// i,j->first are used to retrieve vertex set of graph gk

        //add constraints
        vector<bool> is_y(p,false);//if the component size of a graph >=2, we create a y variable; no y variable otherwise
        for (k = 0; k < p; k++) {
            gk = gs[k];
            //add constraints: sum_y_i^k <= 1, x_u <= y_i^k, and x_u + x_v <= 1
            vector<vector<long>> componentList = gk->getComponent();
            long c_k = componentList.size(); // number of components in gk
            if (c_k >=2) //only if number of component >= 2, we need to use extend formulation; i.e., define and use y variables
            {
                GRBLinExpr sum_y = 0;
                y[k] = new GRBVar[c_k];
                is_y[k] = true;
                for (long i = 0; i < c_k; i++) {
                    y[k][i] = m.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                    sum_y += y[k][i];
                    long sizeC = componentList[i].size();
                    for (long j = 0; j < sizeC; j++) {
                        t = componentList[i][j];
                        u = v0.v2idx[t]; //convert original vertex set to 0,.., n-1 for Gurobi constraint index
                        m.addConstr(x[u] <= y[k][i]); //add constraint: x_u <= y_i^k
                        for (long l = j + 1; l < sizeC; l++) {
                            long t1 = componentList[i][l];
                            v = v0.v2idx[t1]; //convert original vertex set to 1,.., n for Gurobi constraint index
                            if (!gk->IsAdj(t, t1))
                                m.addConstr(x[u] + x[v] <= 1); //add constraint: x_u + x_v <= 1
                        }
                    }
                }
                m.addConstr(sum_y <= 1); //add constraint: sum_y_i^k <= 1
            } else {
                //if the graph itself is connected, we do NOT need y variables; i.e., we use the base version's clique constraints
                //add constraint x_u + x_v <=1 for all u,v are NOT adjacent in graph G_k, k∈[p]
                Adj::iterator i, j;// i,j->first are used to retrieve vertex set of graph gk
                for (i = gk->adj.begin(); i != gk->adj.end(); i++) {
                    u = v0.v2idx[i->first]; //convert original vertex set to 1,.., n for Gurobi constraint index
                    for (j = next(i); j != gk->adj.end(); j++) {
                        if (!gk->IsAdj(i->first, j->first)) {
                            v = v0.v2idx[j->first];
                            m.addConstr(x[u] + x[v] <= 1,
                                        "Clique constraints_" + itos(i->first) + "_" + itos(j->first));
                        }
                    }
                }
            }

            // add constraints related to z variables
            isInGk = vector<bool>(n_nodes, false); //used to check a vertex if exist in graph Gk
            for (i1 = gk->adj.begin(); i1 != gk->adj.end(); i1++) {
                u = v0.v2idx[i1->first]; //convert original vertex set to 1,.., n for Gurobi constraint index
                isInGk[u] = true;
            }

            for (u = 0; u < n_nodes; u++) {
                if (isInGk[u])
                    m.addConstr(x[u] <= z[k]); //add constraint: x_u <= z_k if u is in Gk
                else
                    m.addConstr(x[u] <= 1 - z[k]); //add constraint: x_u <= 1 - z_k if u is NOT in Gk
            }
        }


        //set up objective function
        m.setObjective(GRBSum(n_nodes, x), GRB_MAXIMIZE);

        //Set maximum time limit
        m.getEnv().set(GRB_DoubleParam_TimeLimit,3600);
        //Set Gurobi screen display flag: 0=switch off; 1=default
        m.getEnv().set(GRB_IntParam_OutputFlag,1);
        m.update();

        //m.write("Extend_formulation.lp");// for verification by writing constraints
        m.optimize(); // solve
        r->num_graphs = p;
        r->gurobiSolveTime = m.get(GRB_DoubleAttr_Runtime);
        r->g_name = file_name;
        r->opt_status = "Model Problematic";
        if (m.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            r->opt_status = "Model Solved";
        } else if (m.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
            r->opt_status = "Hit Time Limit";
        }

        r->ub = m.get(GRB_DoubleAttr_ObjBound);
        r->obj = m.get(GRB_DoubleAttr_ObjVal);
        for (t = 0; t < n_nodes; t++) {
            if (x[t].get(GRB_DoubleAttr_X) >= 0.6) //>=0.6 -> 1 since it's binary
                r->bestSol.push_back(v0.V[t]);
        }

        sort(r->bestSol.begin(), r->bestSol.end()); //sort the best solution
        //count the wall clock time of solving the model
        chrono::duration<double> WallTime_span = chrono::duration_cast<std::chrono::duration<double> >(
                chrono_clock::now() - wall_time);
        r->wallTime = WallTime_span.count();
        // Clean up program resources
        ClearGSeq(gs);
        delete[] x;
        delete[] z;
        for (k = 0; k < p; k++) {
            //if is_y[k] is true, it means that we created a y variable using new GRBVar[c_k];
            if (is_y[k]){
                delete[] y[k];
            }
        }
        delete[] y;

    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }

    delete env;

    cout << r->opt_status << endl;

    //write result into an Excel file for each graph
    ofstream fout;
    string outputFile = "max_atomic_clique_extended.csv";
    fout.open(outputFile.c_str(), ios::app);
    fout<<r->g_name<<","<<r->num_nodes<<","<<r->num_edges<<","<<r->num_graphs<<","<<r->readTime<<","<<r->gurobiSolveTime<<","<<r->wallTime<<","<<r->obj<<","<<r->ub<<","<<r->opt_status<<","<<(r->ub - r->obj)/r->obj*100;
    fout<<"\n";
    fout.close();

    // Output Solution
    cout << "----------Output in Graph Collection: " <<r->g_name<<"-------------------"<< endl;
    cout << "Objective:\t" << r->obj << endl;
    cout << "Bound:\t" << r->ub << endl;
    cout << "The best solution:" << endl;
    for (t = 0; t < r->obj; t++) {
        cout << r->bestSol[t] << "-->";
    }
    cout << endl;
    cout << "Total Time : " << r->wallTime << "s " << endl;
    cout << "=====================================================================================" << endl << endl;
    return r;
}

//Edge Peeling + Gurobi (IP based approach): After edge peeling, MACP is reduced to the max clique problem
ACPResult* SolveACP_EdgePeeling(string file_name)
{   // Load graph sequence:
    cout << "=====================================================================================" << endl;
    cout << "Loading "<< file_name << " ..." << endl;
    GSeq gs;
    VSet v0;
    ACPResult* r = new ACPResult(); // get results
    chrono_time_point read_time = chrono_clock::now();
    LoadGSeq(gs, v0, file_name);
    chrono::duration<double>readTime_span = chrono::duration_cast<std::chrono::duration<double> >(chrono_clock::now() - read_time);
    r->readTime = readTime_span.count();
    cout << "Loaded!" << endl << endl;
    // Solve APC
    cout << "Solving the MAX AC problem using Edge Peeling ..." << endl;
    // Get the start time
    chrono_time_point wall_time = chrono_clock::now();
    long p = gs.size();
    long i,j,k,u,v,t = -1;
    //count total edges before edge peeling
    for (i =0; i<p; i++) {
       r->num_edges += gs[i]->num_edges;
    }

    r->g_name = file_name;
    r->num_nodes = v0.V.size();;
    r->num_graphs = p;

    cout<<"Before Edge peeling, number of vertices and edges are: "<<r->num_nodes<<"->"<<r->num_edges<<endl;

    cout<<"Edge peeling......"<<endl;
    chrono_time_point peeling_time = chrono_clock::now();
    Graph g = EdgePeeling(v0,p); //get auxiliary graph after edge peeling
    chrono::duration<double>peelTime_span = chrono::duration_cast<std::chrono::duration<double> >(chrono_clock::now() - peeling_time);
    r->peelingTime = peelTime_span.count();

    cout<<"Edge peeling is complete, which takes "<<r->peelingTime<<" seconds!"<<endl;

    long num_verts = g.num_verts;
    //Convert original vertex index to 0, 1, 2, .., n-1 where n is the size of the graph
    vector<long> C; //store all vertex set of the graph
    Adj::iterator i1;// i1->first is used to retrieve vertex set of graph
    map<long, long> toIdx; //Convert original vertex index to 0, 1, 2, .., n-1 where n is the size of the graph; eg. toIdx[3] = 1 means that vertex 3's index is 1

    //Used for vertex index conversion
    long count =0;
    for (i1 = g.adj.begin(); i1 != g.adj.end(); i1++) {
        C.push_back(i1->first);
        toIdx[i1->first] = count;
        count++;
    }
    //find component list of the graph
    vector<vector<long>> componentList = g.getComponent();
    long num_Comp = componentList.size();
    r->num_components = num_Comp;
    r->new_nodes = g.num_verts;
    r->new_num_edges = g.num_edges;

    cout<<"After edge peeling, the auxiliary graph' #vertex, #edge, and #component are: "<<r->new_nodes<<"--"<<r->new_num_edges<<"--"<<r->num_components<<endl;

    // Create variables
    //Gurobi env setup
    GRBEnv* env = 0;
    try {
        env = new GRBEnv();
        GRBModel m = GRBModel(*env);
        GRBVar *x = GRBVarArray(num_verts, m, 0.0, 1.0, GRB_BINARY);
        GRBVar *y = GRBVarArray(num_Comp, m, 0.0, 1.0, GRB_BINARY);

        //add constraints of max clique problem with extended variables y
        GRBLinExpr sum_y = 0;
        for (i = 0; i < num_Comp; i++) {
            sum_y += y[i];
            for (j = 0; j < componentList[i].size(); j++) {
                u = componentList[i][j];
                //remember to convert original vertex index to 0, 1, 2, .., n-1
                m.addConstr(x[toIdx[u]] <= y[i]); //add constraint: x_v <= y_i for each vertex v in component i
                for (k = j + 1; k < componentList[i].size(); k++) {
                    v = componentList[i][k];
                    if (!g.IsAdj(u, v))
                        m.addConstr(x[toIdx[u]] + x[toIdx[v]] <=
                                    1); //add constraint: x_i + x_j <= 1 if i and j are NOT adjacent
                }
            }
        }
        m.addConstr(sum_y <= 1); //add constraint: sum_y_i^k <= 1

        //set up objective function
        m.setObjective(GRBSum(num_verts, x), GRB_MAXIMIZE);

        //Set maximum time limit
        m.getEnv().set(GRB_DoubleParam_TimeLimit,3600);
        //Set Gurobi screen display flag: 0=switch off; 1=default
        m.getEnv().set(GRB_IntParam_OutputFlag,1);
        m.update();

        //m.write("EdgePeeling_formulation.lp");// for verification by writing constraints
        m.optimize(); // solve
        r->gurobiSolveTime = m.get(GRB_DoubleAttr_Runtime);
        r->opt_status = "Model Problematic";
        if (m.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            r->opt_status = "Model Solved";
        } else if (m.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
            r->opt_status = "Hit Time Limit";
        }
        r->ub = m.get(GRB_DoubleAttr_ObjBound);
        r->obj = m.get(GRB_DoubleAttr_ObjVal);
        for (t = 0; t < num_verts; t++) {
            if (x[t].get(GRB_DoubleAttr_X) >= 0.6) //>=0.6 -> 1 since it's binary
                r->bestSol.push_back(C[t]);
        }
        sort(r->bestSol.begin(), r->bestSol.end()); //sort the best solution
        //count the wall clock time of solving the model
        chrono::duration<double> WallTime_span = chrono::duration_cast<std::chrono::duration<double> >(
                chrono_clock::now() - wall_time);
        r->wallTime = WallTime_span.count();
        // Clean up program resources
        ClearGSeq(gs);
        delete[] x;
        delete[] y;
        cout << r->opt_status << endl;

    }catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }
    delete env;


    //write result into an Excel file for each graph
    ofstream fout;
    string outputFile = "max_atomic_clique_edgePeeling.csv";
    fout.open(outputFile.c_str(), ios::app);
    fout<<r->g_name<<","<<r->num_nodes<<","<<r->num_edges<<","<<r->num_graphs<<","<<r->new_nodes<<","<<r->new_num_edges<<","<<r->num_components<<","<<r->readTime<<","<<r->peelingTime<<","<<r->gurobiSolveTime<<","<<r->wallTime<<","<<r->obj<<","<<r->ub<<","<<r->opt_status<<","<<(r->ub - r->obj)/r->obj*100;
    fout<<"\n";
    fout.close();

    // Output Solution
    cout << "----------Output in Graph Collection: " <<r->g_name<<"-------------------"<< endl;
    cout << "Objective:\t" << r->obj << endl;
    cout << "Bound:\t" << r->ub << endl;
    cout << "The best solution:" << endl;
    for (t = 0; t < r->obj; t++) {
        cout << r->bestSol[t] << "-->";
    }
    cout << endl;
    cout << "Total Time : " << r->wallTime << "s " << endl;
    cout << "=====================================================================================" << endl << endl;
    return r;
}


// Guribi auxiliary functions
GRBVar* GRBVarArray(int n, GRBModel& m_val, double lb, double ub,
                    char type, double obj_coeff)
{
    // delete after use
    GRBVar* vars = new GRBVar[n];

    for (int i = 0; i < n; i++)
        vars[i] = m_val.addVar(lb, ub, obj_coeff, type);

    return vars;
}
GRBLinExpr GRBSum(int n, GRBVar* vars, double* coeff)
{
    int i;
    GRBLinExpr expr = 0;

    if (coeff != NULL)
    {
        for (i = 0; i < n; i++)
            expr += (vars[i] * coeff[i]);
    }
    else
    {
        for (i = 0; i < n; i++)
            expr += vars[i];
    }
    return expr;
}


// clear results
void ClearResults(ACPRs &results)
{
    ACPRs::iterator i;
    for (i = results.begin(); i != results.end(); i++) {
        delete* i;
    }
}

//Get heuristic solution of ACP using the Algo presented in Lu et al. (2021) "Clustering
// temporal disease networks to assist clinical decision support systems in visual analytics of comorbidity progression

vector<long> ACPHeuristic(const GSeq & gs_val){
    vector<long> K, D;
    long p = gs_val.size();
    vector<bool> M(p, false); //maintain graph index containing D
    Graph* gk = NULL;
    Adj::iterator i;// i->first is used to retrieve vertex index of graph gk
    //Step 1: the for-loop is to find a common set D
    for (long k = 0; k < p; k++) {
        gk = gs_val[k];
        //find vertex set of gk and then it to tmp
        vector<long> tmp;
        for (i = gk->adj.begin(); i != gk->adj.end(); i++) {
            tmp.push_back(i->first);
        }
        sort(tmp.begin(),tmp.end());
        if (k == 0){
           D = tmp;//initialize D by assigning first graph's vertex set to D
           M[k] = true;
           continue;
        }
        vector<long> commonV = findCommonV(tmp,D);
        if (!commonV.empty()){
            D = commonV;
            M[k] = true;
        }
    }

    if (D.empty())
        return K;

    /* Step 2: find a heuristic solution k using Gurobi Solver
     Find a subset K is subset D such that K is a clique in G_i[D] for all i in M
     * */

    //Gurobi env setup
    long n_nodes = D.size();
    long u,v, j1, j2 = -1;
    GRBEnv* env = 0;
    try {
        env = new GRBEnv();
        GRBModel m = GRBModel(*env);
        GRBVar *x = GRBVarArray(n_nodes, m, 0.0, 1.0, GRB_BINARY); // Create variables
        m.setObjective(GRBSum(n_nodes, x), GRB_MAXIMIZE); //set up objective function
        //add constraint
        for (j1 = 0; j1 < n_nodes; j1++) {
            u = D[j1];
            for (j2 = j1 + 1;  j2< n_nodes; j2++) {
                v = D[j2];
                for (long k = 0; k < p; k++) {
                    if (!M[k])
                        continue;
                    gk = gs_val[k];
                    if (!gk->IsAdj(u,v)){
                        m.addConstr(x[j1] + x[j2] <= 1);
                        break;
                    }
                }
            }
        }

        //Set maximum time limit
        m.getEnv().set(GRB_DoubleParam_TimeLimit,1800); //solve the model within timelimit 30 minutes
        //Set Gurobi screen display flag: 0=switch off; 1=default
        m.getEnv().set(GRB_IntParam_OutputFlag,1);
        m.update();

        //m.write("formulation.lp");// for verification by writing constraints
        m.optimize(); // solve
        if (m.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            cout<<"Heuristic Model Solved"<<endl;
        } else if (m.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
            cout<<"Heuristic Model Hit Time Limit 30 minutes!"<<endl;
        }

        for (j1 = 0; j1 < n_nodes; j1++) {
            if (x[j1].get(GRB_DoubleAttr_X) >= 0.6) //>=0.6 -> 1 since it's binary
                K.push_back(D[j1]);
        }
        // Clean up program resources
        delete[] x;
    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }

    delete env;
    return K;
}

//heuristic based on edge overlap
vector<long> ACPHeuristic_edgeOverlap(const GSeq & gs_val, const VSet & v0_val){
    //Step 1. intersect all edges in the collection
    Graph g; //create a graph that is generated by intersecting all edges in the collection
    long num_graphs = gs_val.size();
    long num_nodes = v0_val.V.size(); //number of universal vertex
    //edge overlap checking
    long u, v =-1;
    LongSet edgeGraphSet;
    for (auto i1 = v0_val.edgeGraphSet.begin(); i1 != v0_val.edgeGraphSet.end() ; i1++) {
        u = i1->first.first;
        v = i1->first.second;
        edgeGraphSet = i1->second;
        long size_g = edgeGraphSet.size();
        if (size_g == num_graphs){
            g.AddE(u,v); //add edge uv into g
        }
    }
    g.CountVE(); //calculate number of vertices and edges
    vector<long> AC; //store atomic clique
    if (g.num_verts==0)
        return AC;

    //Step 2. find a greedy clique in g.
    vector<long> degree, V; //store each vertex's degree; vertex set of g
    for (auto it = g.adj.begin(); it != g.adj.end(); it++) {
        V.push_back(it->first);
        degree.push_back(g.adj[it->first].size());
    }

    long maxElementIndex = max_element(degree.begin(),degree.end()) - degree.begin();
    AC.push_back(V[maxElementIndex]); //initialize AC
    degree[maxElementIndex] = -1;// this vertex will NOT be considered in the future
    //find new max value and index
    long maxElement = *max_element(degree.begin(),degree.end());
    maxElementIndex = max_element(degree.begin(),degree.end()) - degree.begin();
    while(maxElement != -1){
        bool flag = true;
        for (long i = 0; i<AC.size(); i++) {
            if (!g.IsAdj(V[maxElementIndex], AC[i])){
                flag= false;
                break;
            }
        }
        if (flag)
            AC.push_back(V[maxElementIndex]);

        degree[maxElementIndex] = -1;// this vertex will NOT be considered in the future
        //find new max item and index

        maxElement = *max_element(degree.begin(),degree.end());
        maxElementIndex = max_element(degree.begin(),degree.end()) - degree.begin();
    }
    return AC;
}