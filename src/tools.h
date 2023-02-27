#ifndef tools_h
#define tools_h
#include<bits/stdc++.h>
#include"setting.h"
using namespace std;

ll qow(ll a,ll p,ll mod) {ll ans=1;for(;p;a=a*a%mod,p>>=1) if(p&1) ans=ans*a%mod;return ans;}
ll inv(ll a,ll mod) {return qow(a,mod-2,mod);}

void init()
{
    hash_base1.assign(max_str_length,1);
    hash_base2.assign(max_str_length,1);
    invhash_base1.assign(max_str_length,1);
    invhash_base2.assign(max_str_length,1);
    for(int i=1;i<max_str_length;i++){
        hash_base1[i]=hash_base1[i-1]*step1%mod1;
        invhash_base1[i]=inv(hash_base1[i],mod1);
        hash_base2[i]=hash_base2[i-1]*step2%mod2;
        invhash_base2[i]=inv(hash_base2[i],mod2);
    }
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

string cal_symm(string& s)
{
    string tmp="";
    for(int p=s.length()-1;p>=0;p--){
        tmp+=symm(s[p]);
    }
    return tmp;
}
//find the symmetric char 
char symm(char c)
{
    char ret;
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
            ret='T';//'N' is considered as 'A',so its symmetry is 'T'
    }
    return ret;
}



ull invhash(ull x,int k) {
    // http://xorshift.di.unimi.it/splitmix64.c
    x += 0x9e3779b97f4a7c15;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
    return x ^ (x >> 31);
}


// ull invhash(ull x,int k) {
//     return x;
// }

/*
ull invhash(ull key, int k)
{
    assert(k<32);
    ull mask=(1ull<<2*k)-1;
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}
*/

/*
calculate the k-mer of a single read, identified by its id

*/
//int trans(char c);
int cal_k_mer(int k,int id,ull& val)
{
    int n=vecr[id].str.size(),i=0;
    int pos=0;
    ull cur=0;
    ull rcur=0;
    int shift1=(2*k-2);
    ull mask=(1ull<<2*k)-1;
    if(vecr[id].str.find("N")!=-1){//has N
        val=-1;//max,
        return pos;
    }

    for(int i=0;i<min(n,k);i++){
        int tmp=trans(vecr[id].str[i]);
        assert(tmp>=0&&tmp<=3);
        cur=((cur<<2)|tmp)&mask;
        rcur=(rcur>>2)|(tmp^3ull)<<shift1;
    }
    //val=min(invhash(cur,k),invhash(rcur,k));
    val=invhash(cur,k);
    for(int i=1;i+k<=n;i++){
        int tmp=trans(vecr[id].str[i+k-1]);
        cur=((cur<<2)|tmp)&mask;
        rcur=(rcur>>2)|(tmp^3ull)<<shift1;
        //if(cur==rcur) continue; // skip "symmetric k-mers" as we don't know it strand
        if(invhash(cur,k)<val) val=invhash(cur,k),pos=i;
        //if(invhash(rcur,k)<val) val=invhash(rcur,k),pos=i;
    }
    return pos;
}

puu get_hash_val(string& s)
{
    ull ret1=0,ret2=0;
    for(int i=0;i<s.length();i++){
        char ch=s[i];
        if(ch=='N') ch='A';
        ret1=(ret1+hash_base1[i]*ch)%mod1;
        ret2=(ret2+hash_base2[i]*ch)%mod2;
    }
    return puu(ret1,ret2);
}

bool check(int id1,int id2,int bias)//id1->id2
{
    int tot=0;
    for(int j=0;j<min(vecr[id2].str.length(),vecr[id1].str.length()-bias);j++){
        if(vecr[id1].str[j+bias]!=vecr[id2].str[j]) tot++;
    }
    return tot<=threshold;
}

bool cmp1(int i,int j)
{
    //isrepeat = 0 or -1 is the same
    //if(vecr[i].isrepeat!=vecr[j].isrepeat) return vecr[i].isrepeat<vecr[j].isrepeat;
    if(vecr[i].val!=vecr[j].val) return vecr[i].val<vecr[j].val;
    else return vecr[i].k_mer_pos>vecr[j].k_mer_pos;
}

bool cmp2(spre& s1,spre& s2)
{
    if(s1.val1!=s2.val1) return s1.val1<s2.val1;
    else return s1.val2<s2.val2;
}

bool cmp3(Read& r1,Read& r2)
{
    return r1.isrepeat<r2.isrepeat;
}

#endif