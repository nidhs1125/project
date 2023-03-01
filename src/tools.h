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
            ret=0;//'A'
    }
    return ret;
}

int transr(char c)//for reverse str
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
            ret=3;//'T'
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
//find the symmetric char 
char symm(char c,bool is)
{
    char ret;
    if(is==0){
        switch(c){
            case 'N':
                ret='A';
                break;
            default:
                ret=c;//'N' is considered as 'A',so its symmetry is 'T'
        }
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
                ret='A';//'N' is considered as 'A',so its symmetry is 'T'
        }
    }
    
    return ret;
}



ull invhash(ull x,int k) {
    #if testflag
        return x;
    #else
        // http://xorshift.di.unimi.it/splitmix64.c
        x += 0x9e3779b97f4a7c15;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
        x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
        return x ^ (x >> 31);
    #endif
}



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
int cal_k_mer(int k,int id,ull& val,bool& iskmersymm)
{
    int n=vecr[id].str.size(),i=0;
    int pos=0;
    ull cur=0;
    ull rcur=0;
    int shift1=(2*k-2);
    ull mask=(1ull<<2*k)-1;
    iskmersymm=0;
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
    #if testflag
        val=invhash(cur,k);
        for(int i=1;i+k<=n;i++){
            int tmp=trans(vecr[id].str[i+k-1]);
            cur=((cur<<2)|tmp)&mask;
            if(invhash(cur,k)<val) val=invhash(cur,k),pos=i;
        }
        return pos;
    #else
        val=invhash(cur,k);
        if(invhash(rcur,k)<val) val=invhash(rcur,k),iskmersymm=1,pos=0+k-1;
        for(int i=1;i+k<=n;i++){
            int tmp=trans(vecr[id].str[i+k-1]);
            cur=((cur<<2)|tmp)&mask;
            rcur=(rcur>>2)|(tmp^3ull)<<shift1;
            if(cur==rcur) continue; // skip "symmetric k-mers" as we don't know it strand
            if(invhash(cur,k)<val) val=invhash(cur,k),pos=i,iskmersymm=0;
            if(invhash(rcur,k)<val) val=invhash(rcur,k),pos=i+k-1,iskmersymm=1;
        }
        return pos;
    #endif
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


bool check(int id1,int id2,int k)//id1->id2,in the k-th round
{
    int tot=0;
    int k_mer_pos1=vec_pos[k][id1],issymm1=vec_is[k][id1];
    int k_mer_pos2=vec_pos[k][id2],issymm2=vec_is[k][id2];
    int start_pos1,step1,len1,tag1;
    int start_pos2,step2,len2,tag2;

    tag1=vecr[id1].str.length();
    if(issymm1==0){
        step1=1;
        start_pos1=0;
        len1=vec_pos[k][id1];
    }
    else{
        step1=-1;
        start_pos1=vecr[id1].str.length()-1;
        len1=vecr[id1].str.length()-vec_pos[k][id1];
    }

    tag2=vecr[id2].str.length();
    if(issymm2==0){
        step2=1;
        start_pos2=0;
        len2=vec_pos[k][id2];
    }
    else{
        step2=-1;
        start_pos2=vecr[id2].str.length()-1;
        len2=vecr[id2].str.length()-vec_pos[k][id2];
    }
    
    assert(len1>=len2);
    start_pos2+=step2*(len1-len2);
    int p1=start_pos1,p2=start_pos2;
    while(p1>=0&&p1<tag1&&p2>=0&&p2<tag2){
        if(symm(vecr[id1].str[p1],issymm1)!=symm(vecr[id2].str[p2],issymm2)) tot++;
        p1+=step1;
        p2+=step2;
    }
    return tot<=threshold;
}

bool cmp1(int i,int j)//before sorting, use vecr directly
{
    //if(vecr[i].isrepeat!=vecr[j].isrepeat) return vecr[i].isrepeat<vecr[j].isrepeat;
    if(vecr[i].val!=vecr[j].val) return vecr[i].val<vecr[j].val;
    else{
        //if issymm has diff direction
        //so we need to compare with the distance between start(end) pos of k_mer and the head(tail).
        int len1,len2;
        if(vecr[i].iskmersymm==0) len1=vecr[i].k_mer_pos;
        else len1=vecr[i].str.length()-vecr[i].k_mer_pos;
        return len1>len2;
    }
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