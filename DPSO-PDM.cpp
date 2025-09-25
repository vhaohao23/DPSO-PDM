#include<bits/stdc++.h>
using namespace std;

#define rep(i,a,b,x)  for(int i=a;i<=b;i+=x)

random_device rd;   
mt19937 gen(rd());

const int N=20;
const int T=20;
const double mutationProb=0.1;

int n,m;
vector<vector<int>> E;
vector<vector<int>> P(N+1);
vector<vector<int>> Pb(N+1);
vector<int> Pg;
vector<int> Q(N+1,0);
vector<int> Qb(N+1,0);
int Qg=0;
vector<int> k;
vector<vector<int>> A;

double modularity(vector<int> l){
    double Q=0;

    for (int i=1;i<=n;i++)
        for (int j=1;j<=n;j++)
                Q+=double(A[i][j]-k[i]*k[j]/double(2*m))*double(l[i]==l[j]);
    return Q/double(2*m);
}

void LAR_rand(vector<vector<int>> &a){
    rep(u,1,n,1){
        if (!E[u].size()) continue;
        uniform_int_distribution<int> disv(0,E[u].size()-1);
        int v=E[u][disv(gen)];
        a[u].push_back(v);
        a[v].push_back(u);
    }
}

vector<int> decoding(vector<vector<int>> a){
    bool dd[n+1]={};
    vector<int> l(n+1);
    int cnt=0;

    rep(i,1,n,1)
        if (!dd[i]){
            ++cnt;
            queue<int> q;
            q.push(i);
            while (!q.empty()){
                int u=q.front();
                q.pop();
                l[u]=cnt;
                for (int v:a[u])
                    if (!dd[v]){
                        dd[v]=true;
                        q.push(v);
                    }
            }
        }
    return l;
}

void initialization(){
    rep(i,1,N,1){
        vector<vector<int>> a(n+1);
        LAR_rand(a);
        P[i]=decoding(a);
        Q[i]=modularity(P[i]);
        Pb[i]=P[i];
        Qb[i]=Q[i];
        if (Qb[i]>Qg){
            Qg=Qb[i];
            Pg=Pb[i];
        }
    } 
}

void mutation(vector<int> &Pbi){
    uniform_real_distribution<double> rand_dist(0.0, 1.0);
    bool dd[n+1]={};

    for (int k = 1; k <= n; k++) {
        double random_value = rand_dist(gen);
        
        if (random_value <= mutationProb) {
            // Find the community with largest modularity increment among neighbors
            int best_community = Pbi[k];
            double best_increment = 0;
            
            // Calculate original modularity
            double orig_mod = modularity(Pbi);
            dd[Pbi[k]] = 1;
            for (int neighbor : E[k]) {
                int l_neigh = Pbi[neighbor];
                if (dd[l_neigh]) continue;
                dd[l_neigh] = 1;
                
                // Try moving node k to community l_neigh
                int original_community = Pbi[k];
                Pbi[k] = l_neigh;
                
                // Calculate new modularity
                double new_mod = modularity(Pbi);
                
                // Calculate increment
                double increment = new_mod - orig_mod;
                
                if (increment > best_increment) {
                    best_increment = increment;
                    best_community = l_neigh;
                }
                
                // Restore original community
                Pbi[k] = original_community;
            }
            
            // Apply mutation to the best community found
            if (best_community != Pbi[k]) {
                Pbi[k] = best_community;
            }
        }
    }
}

void DPSO_PDM(){

}
int main(){
    cin>>n>>m;

    E.resize(n+1);
    k.resize(n+1);
    A.resize(n+1,vector<int>(n+1,0));
    
    int u,v;
    rep(i,1,m,1){
        cin>>u>>v;
        E[u].push_back(v);
        E[v].push_back(u);
        A[u][v]=1,A[v][u]=1;
        k[u]++,k[v]++;
    }
}