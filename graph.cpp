#include "graph.h"
#include "common.h"

//Constructor
Graph::Graph()
{
    num_verts=0;
    num_edges=0;
}


//Destructor
Graph::~Graph()
{
    adj.clear();
}

//add edge a-b
void Graph::AddE(const long & a, const long & b, bool recount_v_e)
{
    if (a != b) {
        adj[a].insert(b);
        adj[b].insert(a);
    }/* //if a =b,it is a loop; we exclude this case
    else if (adj.find(a) == adj.end()) {
        //adj[a] = LongSet();
    }*/
    //else do nothing

    if (recount_v_e) {
        CountVE();
    }
}
//count number of vertex and edges
void Graph::CountVE()
{
    num_verts = 0;
    num_edges = 0;
    map<long, LongSet>::iterator i;
    for (i = adj.begin(); i != adj.end(); i++) {
        num_edges += i->second.size();
    }
    num_verts = adj.size();
    num_edges /= 2;
}

//check if a graph contains vertex a
bool Graph::IsContainV(long a)
{
    return adj.find(a) != adj.end();
}

//check if a and b are adjacent for undirected graph
bool Graph::IsAdj(const long & a, const long & b)
{
    if (adj.find(a) == adj.end()) {
        return false;
    }

    if (adj.find(b) == adj.end()) {
        return false;
    }

    return (adj[a].find(b) != adj[a].end());
}


//Identify connected components
vector<vector<long>> Graph::getComponent(){
    vector<bool> reached(num_verts,false);
    vector<long> tempVector;
    vector<vector<long>> componentList;
    long s, t = -1;

    vector<long> C; //store all vertex set of the graph in S
    Adj::iterator i1;// i1->first is used to retrieve vertex set of graph
    map<long, long> toIdx; //Convert original vertex index to 0, 1, 2, .., n-1 where n is the size of the graph; eg. toIdx[3] = 1 means that vertex 3's index is 1

    //Used for vertex index conversion
    long count =0;
    for (i1 = adj.begin(); i1 != adj.end(); i1++) {
        C.push_back(i1->first);
        toIdx[i1->first] = count;
        count++;
    }

    //find a root node that is not reached; do BFS
    for(long r=0; r<num_verts; r++)
    {
        if(!reached[r])
        {
            reached[r]=true;
            s = C[r]; //find the vertex when index = r
            tempVector = kBFS(s,num_verts+1); //tempVector does not include root s
            tempVector.push_back(s);
            sort(tempVector.begin(),tempVector.end());
            for(long j=0;j<tempVector.size();j++)
            {
                t = toIdx[tempVector[j]]; //convert original vertex set to 0， 1,.., n-1
                reached[t]=true;
            }
            componentList.push_back(tempVector);
            tempVector.clear();
        }// if reached
    }//for r
    return componentList;
}//end getComponent


//BFS up to level-k;
vector<long> Graph::kBFS(long s, long k){
    long bigM = num_verts*10;
    vector<long> dist(num_verts,bigM);
    queue<long> Q;
    vector<long> ReachedVertices; //store component vertex starting from s

    //Convert original vertex index to 0, 1, 2, .., n-1 where n is the size of the graph
    map<long, long> toIdx; // eg. toIdx[3] = 1 means that vertex 3's index is 1
    Adj::iterator i1;// i1->first is used to retrieve vertex set of graph

    long count =0;
    for (i1 = adj.begin(); i1 != adj.end(); i1++) {
        toIdx[i1->first] = count;
        count++;
    }

    dist[toIdx[s]] = 0; //initialize root
    Q.push(s);

    long u,v,t;
    bool flag = true;
    while ( (!Q.empty()) && flag ) {
        t = Q.front();
        u = toIdx[t];//convert original vertex set to 0, 1,.., n-1
        Q.pop();
        LongSet::iterator i;// i is used to retrieve adjacent list of vertex u
        for (i = adj[t].begin(); i != adj[t].end(); i++) {
            v = toIdx[*i]; //convert original vertex set to 0, 1,.., n-1
            if(dist[v] > dist[u] + 1){
                dist[v] = dist[u] + 1;
                // BFS only up to level k
                if (dist[v] >= k+1)
                {
                    flag = false;
                    break;
                }
                else{
                    Q.push(*i);
                    ReachedVertices.push_back(*i);
                }
            }
        }
    }// while
    return ReachedVertices; ////does not include root s
}

// Function to load graph sequence
void LoadGSeq(vector<Graph*> & gs_val, VSet& v0_val, const string & file_val)
{
    //vertex index can be non-consecutive, but graph index starts from 1 to n
    int gindx = 0;
    Graph* g = NULL;

    string str;
    ifstream fs(file_val.c_str());
    if (!fs)
    {
        cout << "Cannot open instance file" <<file_val.c_str()<< "!!!"<<endl;
        return;
    }

    while (!fs.eof())
    {
        getline(fs, str);
        if (str.size() < 1 || str[0] != 'e')
            continue;
        StrV sv;
        SplitStr(sv, str); //split string vector of each line; eg: Line e 3 5 3 means a = 3, b =5,c = 3
        long a = stoi(sv[1]);
        long b = stoi(sv[2]);
        long c = stoi(sv[3]);
        if (a == b)
            continue;//exclude self loop in a graph

        //store set of graphs containing vertex v in vertexGraphSet
        //note that v0_val.vertexGraphSet[a] is a Set (which allows only unique data), so no need to check if c exists explicitly.
        v0_val.vertexGraphSet[a].insert(c);
        v0_val.vertexGraphSet[b].insert(c);
        ////store set of graphs containing edge ab; maintain the order a<b for pair {a,b}
        if (a < b)
            v0_val.edgeGraphSet[{a,b}].insert(c);
        else if (a > b)
            v0_val.edgeGraphSet[{b,a}].insert(c);

        //store universal vertex set in V, and set Universal vertex index from 0, 1, .., n-1 in v2idx
        if (v0_val.v2idx.find(a) == v0_val.v2idx.end()) {
            v0_val.V.push_back(a);
            v0_val.v2idx[a] = (long)v0_val.V.size() - 1;
        }

        if (v0_val.v2idx.find(b) == v0_val.v2idx.end()) {
            v0_val.V.push_back(b);
            v0_val.v2idx[b] = (long)v0_val.V.size() - 1;
        }


        if (c > gindx)
        {
            if (c > 1) {
                g->CountVE();
                gs_val.push_back(g);
            }//when c==1: g is still NULL
            g = new Graph();
            gindx = c;
        }
        g->AddE(a, b);

    }
    g->CountVE();
    gs_val.push_back(g); // push back the last graph
    fs.close();
}


//clear graph collections
void ClearGSeq(vector<Graph*> & gs_val)
{
    GSeq::iterator i;
    for (i = gs_val.begin(); i != gs_val.end(); i++) {
        delete* i;
    }
}

//function to do edge peeling
Graph EdgePeeling(VSet & v0_val, const long num_graphs){
    Graph g; //create an auxiliary graph
    map<pair<long ,long>, LongSet>::iterator i1;
    LongSet::iterator i2;
    map<long, LongSet>::iterator i3;
    LongSet::iterator i4;

    long num_nodes = v0_val.V.size(); //number of universal vertex

    //check if a graph contains v using boolean variables: eg, if isContainV[i][j] = true, it means that graph j+1 contains vertex i
    map<long,vector<bool>> isContainV;
    for (long i = 0; i < num_nodes; i++) {
        isContainV[v0_val.V[i]] = vector<bool>(num_graphs, false);//initialization
    }

    long k1, k2 = -1;
    for (i3 = v0_val.vertexGraphSet.begin();i3 != v0_val.vertexGraphSet.end(); i3++) {
        k1 = i3->first; //vertex index
        for (i4 = i3->second.begin(); i4 != i3->second.end(); i4++) {
            k2 = *i4; //graph index containing vertex k1
            isContainV[k1][k2-1] = true; //note that graph index starting 1 instead of 0;
        }
    }

    //edge peeling algorithm implementation
    long u, v =-1;
    LongSet edgeGraphSet;
    bool flag = false;
    for (i1 = v0_val.edgeGraphSet.begin(); i1 != v0_val.edgeGraphSet.end() ; i1++) {
        u = i1->first.first;
        v = i1->first.second;
        edgeGraphSet = i1->second;
        long size_g = edgeGraphSet.size();
        if (size_g == num_graphs){
            g.AddE(u,v); //add edge uv into g
        } else{
            flag = true;
            vector<bool> isInGSet(num_graphs, false);
            for (i2 = edgeGraphSet.begin();  i2 != edgeGraphSet.end() ; i2++) {
                isInGSet[*i2-1] = true; //note that graph index starting 1 instead of 0;
            }
            //iterate each graph that not exists in edgeGraphSet
            for (long t = 1; t < num_graphs+1; t++) {
                //if graph t does not contain edge uv
                if (!isInGSet[t-1]){
                    //if graph t does not contain vertex u or v
                    if ((!isContainV[u][t-1]) && (!isContainV[v][t-1])){
                        flag = true;
                    } else{
                        flag = false;
                        break;
                    }
                }
            }
            if (flag) {
                g.AddE(u, v); //add edge uv into g
            }
        }
    }
    g.CountVE(); //calculate number of vertices and edges
    return g;
}


 //write the auxiliary graph generated by edge peeling into a txt file with adjacent list format: graph index is 1, 2, .., n
 /* Adjacency lists of a 4-vertex, 5-edge graph
  * 4 5
  * 2 3 4
  * 1 3
  * 1 2 4
  * 1 3
 */
 void WriteGraph2AdjList(string file_name, vector<long> & auxVset){
     cout << "Writing the auxiliary graph generated by edge peeling graph "<< file_name<<" into a txt file" << " ..." << endl;
     GSeq gs;
     VSet v0;
     LoadGSeq(gs, v0, file_name); //load graph collection
     long p = gs.size();
     Graph g = EdgePeeling(v0,p); //get auxiliary graph after edge peeling
     long num_verts = g.num_verts;
     long num_edges = g.num_edges;

     //Convert original vertex index to 1, 2, .., n  where n is the size of the graph
     Adj::iterator i1;// i1->first is used to retrieve vertex set of graph
     map<long, long> toIdx; //Convert original vertex index to 1, 2, .., n where n is the size of the graph; eg. toIdx[3] = 1 means that vertex 3's index is 1
     //auxVset: store all original vertex set of the graph
     long count =1;
     for (i1 = g.adj.begin(); i1 != g.adj.end(); i1++) {
         auxVset.push_back(i1->first);
         toIdx[i1->first] = count;
         count++;
     }

     //write auxiliary graph into a txt file with adjacent list format
     long filename_size = file_name.size();
     string output_filename = file_name.substr(0,filename_size-7) + ".txt"; //Output file name is obtained by removing ".DIMACS"

     ofstream output;
     output.open(output_filename.c_str(), ios::out);
     if (!output.is_open()){
         cout << "File " << output_filename << " could not be opened!!!\n";
         return;
     }

     output << num_verts << " " << num_edges << "\n"; //write number of vertices and edges in the first line
     cout<<"print each vertex and its associated adjacency list................"<<endl;
     LongSet::iterator j;
     for (i1 = g.adj.begin(); i1 != g.adj.end(); i1++) {
         cout<<"vertex: "<<i1->first<<"----"<<endl;
         for (j = g.adj[i1->first].begin(); j != g.adj[i1->first].end(); j++) {
             if (next(j) == g.adj[i1->first].end()){
                 output << toIdx[*j]; //for last vertex of i1 neighbors, no need space
             }else{
                 output << toIdx[*j] << " ";
             }
         }
         output << endl;
     }
     output.close();
     // Clean up program resources
     ClearGSeq(gs);
}


/* write the auxiliary graph generated by edge peeling into a txt file with edge list format
 * The first line of the file must include the number vertices and edges. The vertices must be labeled from 0 to n-1.
 * %% Edge list of a 4-vertex, 4-edge graph
 *4 4
 * 0 1
 * 0 3
 * 1 2
 * 2 3
*/
void WriteGraph2EdgeList(string file_name, vector<long> & auxVset){
    cout << "Writing the auxiliary graph generated by edge peeling graph "<< file_name<<" into a txt file" << " ..." << endl;
    GSeq gs;
    VSet v0;
    LoadGSeq(gs, v0, file_name); //load graph collection
    long p = gs.size();

    long total_num_edges = 0;
    //count total edges before edge peeling
    for (long i =0; i<p; i++) {
        total_num_edges += gs[i]->num_edges;
    }
    chrono_time_point peeling_time = chrono_clock::now();
    Graph g = EdgePeeling(v0,p); //get auxiliary graph after edge peeling
    chrono::duration<double>peelTime_span = chrono::duration_cast<std::chrono::duration<double> >(chrono_clock::now() - peeling_time);

    long num_verts = g.num_verts;
    long num_edges = g.num_edges;

    long num_Comp = 0;
    //find component list of the graph
    //vector<vector<long>> componentList = g.getComponent();
    //long num_Comp = componentList.size();

    //write result into an Excel file for each graph
    ofstream fout;
    string outputFile = "graph_stat.csv";
    fout.open(outputFile.c_str(), ios::app);
    fout<<file_name<<","<<v0.V.size()<<","<<total_num_edges<<","<<p<<","<<num_verts<<","<<num_edges<<","<<num_Comp<<","<<peelTime_span.count()<<",";
    fout<<"\n";
    fout.close();

    //Convert original vertex index to 0, 1, 2, .., n-1  where n is the size of the graph
    Adj::iterator i1;// i1->first is used to retrieve vertex set of graph
    map<long, long> toIdx; //Convert original vertex index to 0, 1, 2, .., n-1 where n is the size of the graph; eg. toIdx[3] = 1 means that vertex 3's index is 1
    //auxVset: store all original vertex set of the graph
    long count = 0;
    for (i1 = g.adj.begin(); i1 != g.adj.end(); i1++) {
        auxVset.push_back(i1->first);
        toIdx[i1->first] = count;
        count++;
    }


    //write auxiliary graph into a txt file with edge list format
    long filename_size = file_name.size();
    string output_filename = "./EdgeList/"+ file_name.substr(0,filename_size-7) + ".txt"; //Output file name is obtained by removing ".DIMACS"

    ofstream output;
    output.open(output_filename.c_str(), ios::out);
    if (!output.is_open()){
        cout << "File " << output_filename << " could not be opened!!!\n";
        return;
    }
    cout<<"%% Edge list of a "<< num_verts <<"-vertex, "<< num_edges <<"-edge graph"<<endl;
    output << num_verts << " " << num_edges << "\n"; //write number of vertices and edges in the first line
    LongSet::iterator j;
    for (i1 = g.adj.begin(); i1 != g.adj.end(); i1++) {
        for (j = g.adj[i1->first].begin(); j != g.adj[i1->first].end(); j++) {
            if (toIdx[i1->first] < toIdx[*j]){
                output << toIdx[i1->first] << " "<< toIdx[*j]<< endl;
            }
        }
    }
    output.close();
    // Clean up program resources
    ClearGSeq(gs);
}

//Generate auxiliary graphs with edge list format in batch
void WriteGraphBatch(string graph_list)
{
    string str;
    ifstream fs(graph_list);

    //write graph statistics into an Excel file
    string outputFile = "graph_stat.csv"; // Result file
    ofstream fout;
    fout.open(outputFile.c_str(), ios::out);
    //FILE OUTPUT: columns headers
    fout<<"Name,#vertex,#edges, #graphs, #newVertics, #newEdges, #components,PeelTime \n";
    fout.close();

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

        vector<long> S;
        WriteGraph2EdgeList(str,S);
    }
    fs.close();

    cout<<"Finish the writing of graphs!"<<endl;
}

//preprocess a collection of graphs using k-core
void KCore(GSeq & gs,VSet & v0, long k)
{
    if (k<=1)
        return;
    long p = gs.size();
    long n_nodes = v0.V.size();
    vector<bool> toBeDel(n_nodes, false); //if there exists a vertex v in one graph such that deg_v <=k-1, then marked as true and it will be deleted in all other graphs.
    Graph* gk = NULL;
    Adj::iterator i;// i->first is used to retrieve vertex index of graph gk
    long j,u,v = -1;
    bool flag = true;
    while (flag){
        flag = false;
        for (j = 0; j < p; j++) {
            gk = gs[j];
            for (i = gk->adj.begin(); i!= gk->adj.end(); i++) {
                v = i->first;
                if (gk->adj[v].empty())
                    continue;
                else if (toBeDel[v0.v2idx[v]]){
                    flag = true;
                    LongSet tmpAdj = gk->adj[v];
                    //  remove v from v's neighbors
                    LongSet::iterator i1;
                    for (i1 = tmpAdj.begin(); i1 != tmpAdj.end(); i1++) {
                        u = *i1;//neighbor of v
                        gk->adj[u].erase(v);
                    }
                    // set v degree =0
                    gk->adj[v].clear();
                }else if (gk->adj[v].size() <= k-1){
                    //degree of v in gk is at most k-1, then 1): set toBeDel[v0.V.idx[v]] = true
                    toBeDel[v0.v2idx[v]] = true;
                    flag = true;
                    LongSet tmpAdj = gk->adj[v];
                    // 2): remove v from v's neighbors
                    LongSet::iterator i1;
                    for (i1 = tmpAdj.begin(); i1 != tmpAdj.end(); i1++) {
                        u = *i1;//neighbor of v
                        gk->adj[u].erase(v);
                    }
                    // 3):set degree v =0
                    gk->adj[v].clear();
                }
            }
        }
    }

    //remove v of degree =0 from adj in each graph
    for (j = 0; j < p; j++) {
        gk = gs[j];
        for (auto it= gk->adj.cbegin(); it!= gk->adj.cend();) {
            if (gk->adj[it->first].empty()){
                gk->adj.erase(it++);
            }else
            {
                ++it;
            }
        }
    }
    // store new universal vertex set
    LongSet V;
    for (j = 0; j < p; j++) {
        gk = gs[j];
        for (auto it= gk->adj.cbegin(); it!= gk->adj.cend(); it++) {
            V.insert(it->first); //Set does not allow duplicate values; so no need to check existence
        }
        gk->CountVE();
    }
    //clear V and v2idx
    v0.V.clear();
    v0.v2idx.clear();

    //update v0: only store vertex v of degree >=1
    for (auto it = V.cbegin();it!=V.cend();it++) {
        v0.V.push_back(*it);
        v0.v2idx[*it] = v0.V.size()-1;
    }
}
