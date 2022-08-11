#include <iostream>
#include "model.h"

int main()
{
    ifstream parameter_file("parameter.txt"); // set up Parameters: data file location and solver type are set up in this txt file
    string line;
    long solver_type = 1;//default
    string fileLocation;
    long count =0;
    while (getline(parameter_file, line)){
        count++;
        if (line[0]=='%')
            continue;
        if (count ==2)
            fileLocation = line; //second line is the data file location
        else if (count ==3)
            solver_type = atol(line.c_str());// 3rd line is the solver type
    }

    chdir(fileLocation.c_str()); //change working directory
    //chdir("./data/DIMACS/"); //change working directory

    ACPRs results;
    ACPRs::iterator i;

    //solve the problem using different solvers and write results into an Excel file
    // 1: enhanced formulation; 2: edge peeling + WB solver; default: monolithic formulation
    if (solver_type == 1){
        //write results into an Excel file
        string outputFile = "max_atomic_clique_extended.csv"; // Result file
        ofstream fout;
        fout.open(outputFile.c_str(), ios::out);
        //FILE OUTPUT: columns headers
        fout<<"Name,#vertex, #edges, #graphs, ReadTime, grbSolveTime, WallTime, Best Obj, Best UB, Status, MIP (%) \n";
        fout.close();
        SolveBatch("input.txt", solver_type, results); //solve model using extended formulation

    }else if (solver_type ==2){
        //First, we use edge peeling to create auxiliary graphs
        //Then, we use WB solver (Walteros and Buchanan 2020 OR paper) for the MCP on these auxiliary graphs
        WriteGraphBatch("input.txt");//generating auxiliary graphs

        //The following code provides another way to solve MCP problem directly using Gurobi solver: Edge peeling + Gurobi
        /*
        //write results into an Excel file
        string outputFile = "max_atomic_clique_edgePeeling.csv"; // Result file
        ofstream fout;
        fout.open(outputFile.c_str(), ios::out);
        //FILE OUTPUT: columns headers
        fout<<"Name,#vertex,#edges, #graphs, #newVertics, #newEdges, #components, ReadTime,PeelTime, grbSolveTime, WallTime,Best Obj, Best UB, Status, MIP (%) \n";
        fout.close();
        SolveBatch("input.txt", solver_type, results); //solve model using edge peeling
        */
    }else{
        //write results into an Excel file
        string outputFile = "max_atomic_clique_monolithic.csv"; // Result file
        ofstream fout;
        fout.open(outputFile.c_str(), ios::out);
        //FILE OUTPUT: columns headers
        fout<<"Name,#vertex,#graphs,ReadTime, grbSolveTime, WallTime,Best Obj, Best UB, Status, MIP (%) \n";
        fout.close();
        SolveBatch("input.txt", solver_type, results); //solve model using monolithic formulation
    }

    //Print results
    cout<<"Graph\tObj\tUB\tTime\tStatus"<<endl;
    for(i=results.begin();i!=results.end();i++){
        ACPResult* r = *i;
        cout<<r->g_name<<"\t"<<r->obj<<"\t"<<r->ub<<"\t"<<r->wallTime<<"\t"<<r->opt_status<<endl;
    }
    ClearResults(results);
}
