#include<bits/stdc++.h>
#include<unistd.h>
using namespace std;
typedef unsigned long long ll;
ll qow(ll a,ll p) {ll ans=1;for(;p;a=a*a,p>>=1) if(p&1) ans=ans*a;return ans;}
int main(int argc,char** argv)
{
    vector<int> tmp;
    tmp.assign(100000000,1);
}