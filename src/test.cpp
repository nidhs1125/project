#include<bits/stdc++.h>
#include<unistd.h>
using namespace std;
typedef long long ll;
typedef unsigned long long ull;
#define mod 998244353ll
ll qow(ll a,ll p) {ll ans=1;for(;p;a=a*a,p>>=1) if(p&1) ans=ans*a;return ans;}
ll get_hash_value(string s)
{
    ll val=0;
    ll base=1;
    for(int i=0;i<s.length();i++){
        val=(val+base*s[i])%mod;
        base*=131;
    }
    return val;
}
int main(int argc,char** argv)
{
    // string str="AAAAAAAAAAAAAAAAAAAA";
    // string str2="AAAAAAAAAAAAAAAAAAA";
    // cout<<get_hash_value(str)<<' '<<get_hash_value(str2);
    cout<<(ull)-1<<'\n';
    for(int i=0;i<=10;i++) cout<<~i<<'\n';
}