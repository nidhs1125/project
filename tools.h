#ifndef tools_h
#define tools_h
#include<bits/stdc++.h>
#include"setting.h"
using namespace std;

ll qow(ll a,ll p) {ll ans=1;for(;p;a=a*a%mod,p>>=1) if(p&1) ans=ans*a%mod;return ans;}
ll inv(ll a) {return qow(a,mod-2);}


void init()
{
    hash_base[0]=invhash_base[0]=1;
    for(int i=1;i<max_str_length;i++) hash_base[i]=hash_base[i-1]*step%mod,invhash_base[i]=inv(hash_base[i]);
}

/*
calculate the k-mer of a single string s,return the pos,indexed by 0
the algorithm is baseed on lyndon decomposition
*/
int cal_k_mer(int k,string& s)
{
    if(s.length()<k) return 0;
    int n=s.size(),i=0;
    int pos=0;
    while(i<=n-k){
        pos=i;
        int j=i+1,k=i;
        while (j<n&&s[k]<=s[j]) {
            if(s[k]<s[j]) k=i;
            else k++;
            j++;
        }
        while(i<=k) i+=j-k;
    }
    return pos;

}

ll get_hash_value(string s)
{
    ll val=0;
    for(int i=0;i<s.length();i++){
        val=(val+hash_base[i]*s[i])%mod;
    }
    return val;
}

int trans(char c)
{
    int ret;
    switch(c){
        case 'A':
            ret=0;
            break;
        case 'C':
            ret=1;
            break;
        case 'G':
            ret=2;
            break;
        case 'T':
            ret=3;
            break;
        default:
            ret=4;
    }
    return ret;
}

char rtrans(int c)
{
    char ret;
    switch(c){
        case 0:
            ret='A';
            break;
        case 1:
            ret='C';
            break;
        case 2:
            ret='G';
            break;
        case 3:
            ret='T';
            break;
        default:
            ret='N';
    }
    return ret;
}

bool cmp1(Read read1,Read read2)
{
    if(!!(~read1.cid)!=!!(~read2.cid)) return read1.cid<read2.cid;//-1 go ahead
    else if(read1.val==read2.val) return read1.k_mer_pos>read2.k_mer_pos;
    else return read1.val<read2.val;
}

bool cmp2(Contig con1,Contig con2)
{
    return con1.str<con2.str;
}

bool cmp3(Read read1,Read read2)
{
    return read1.rid<read2.rid;
}

bool cmp4(Contig con1,Contig con2)
{
    return con1.cid<con2.cid;
}

#endif