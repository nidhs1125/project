#include<bits/stdc++.h>
#include<unistd.h>
#include <sys/time.h>
using namespace std;
typedef long long ll;
#define pb push_back
vector<string> vec1,vec2;
char symm(char c,bool is)
{
    char ret;
    if(is==0){
        ret=c;
    }
    else{
        switch(c){
            case 'A':
                ret='T';
                break;
            case 'T':
                ret='A';
                break;
            case 'C':
                ret='G';
                break;
            case 'G':
                ret='C';
                break;
            default:
                ret='N';
        }
    }
    
    return ret;
}
string cal_symm(string& s)
{
    string tmp="";
    for(int p=s.length()-1;p>=0;p--){
        tmp+=symm(s[p],1);
    }
    return tmp;
}
ll qow(ll a,ll p,ll mod) {ll ans=1;for(;p;a=a*a%mod,p>>=1) if(p&1) ans=ans*a%mod;return ans;}
signed main()
{
    string out_path="test";
    ofstream fout;
    ifstream fin;
    char ch;

    fout.open(out_path+".mine");
    fin.open(out_path+".bsc1");
    while(fin>>ch) fout<<ch;
    fin.close();
    fin.open(out_path+".bsc2");
    while(fin>>ch) fout<<ch;
    fin.close();
    fin.open(out_path+".bsc3");
    while(fin>>ch) fout<<ch;
    fout.close();
    fin.close();
    return 0;
}