#include<bits/stdc++.h>
using namespace std;

#define rep(i,a,b,x)  for(int i=a;i<=b;i+=x)

random_device rd;   
mt19937 gen(rd());

const int N=20;
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

void DPSO_PDM(){

}
int main(){
    cin>>n>>m;
    int u,v;
    E.resize(n+1);
    k.resize(n+1);
    A.resize(n+1,vector<int>(n+1,0));
    rep(i,1,m,1){
        cin>>u>>v;
        E[u].push_back(v);
        E[v].push_back(u);
        A[u][v]=1,A[v][u]=1;
        k[u]++,k[v]++;
    }
}