#ifndef tools_h
#define tools_h
#include<bits/stdc++.h>
#include<unistd.h>
#include <sys/time.h>
#include"setting.h"

int trans(char c);

ll qow(ll a,ll p,ll mod) {ll ans=1;for(;p;a=a*a%mod,p>>=1) if(p&1) ans=ans*a%mod;return ans;}
ll inv(ll a,ll mod) {return qow(a,mod-2,mod);}

puu get_hash_val(string& s)
{
    ull ret1=0,ret2=0;
    for(int i=0;i<s.length();i++){
        char ch=s[i];
        ret1=(ret1+hash_base1[i]*ch)%mod1;
        ret2=(ret2+hash_base2[i]*ch)%mod2;
    }
    return puu(ret1,ret2);
}

puu get_seg_value(string& str,int l,int r)
{
    assert(r-l+1<32);
    assert(l<=r&&l>=0&&r<str.length());
    int k=r-l+1;
    ull cur=0;
    ull rcur=0;
    int shift1=2*k-2;
    ull mask=(1ull<<2*k)-1;

    for(int i=l;i<=r;i++){
        ull tmp=trans(str[i]);
        if(tmp==4) tmp=0;//regard as 'A'
        cur=(cur<<2|tmp)&mask;
        rcur=rcur>>2|(tmp^3ull)<<shift1;
    }
    return puu(cur,rcur);
}

inline void print_time()
{
    gettimeofday(&end1, NULL);
    cout<<"time consumption:"<<1000*(end1.tv_sec - start1.tv_sec) + (end1.tv_usec - start1.tv_usec)/1000<<"ms\n";
}

inline int trans(char c)
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
            ret=4;//'N'
    }
    return ret;
}

inline char trans(int c)//when compressing, use only 0-3
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

//find the symmetric char 
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

inline int cal_len(int id)
{
    return iskmersymm[id]?read_len-1-k_mer_pos[id]:k_mer_pos[id];
}

inline void mutex_init()
{
    my_mutex=PTHREAD_MUTEX_INITIALIZER;
    array_mutex=PTHREAD_MUTEX_INITIALIZER;
    my_id=0,id_upbound=vecr.size();
    
    thread_id.assign(thread_num,-1);
    thread_tmp.resize(thread_num+1);
    for(int i=0;i<=thread_num;i++) thread_tmp[i]=i;
}

inline ull invhash(ull x,int k) {
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
calculate the k-mer of a single read, identified by its id
*/
void* k_minimizer(void* arg)
{
    int k=*(int*)arg;
    while(1){
        int curid;
        pthread_mutex_lock(&my_mutex);
        curid=my_id++;
        pthread_mutex_unlock(&my_mutex);
        if(curid>=id_upbound) break;
        int id=vecr[curid];
        int i=0;
        k_mer_pos[id]=0;
        val[id]=0;
        ull cur=0;
        ull rcur=0;
        int shift1=(2*k-2);
        ull mask=(1ull<<2*k)-1;
        iskmersymm[id]=0;
        string &str=vread[id];
        assert(str.find("N")==string::npos);
        assert(read_len>=k);

        for(int i=0;i<k;i++){
            ull tmp=trans(str[i]);
            cur=(cur<<2|tmp)&mask;
            rcur=rcur>>2|(tmp^3ull)<<shift1;
        }
        val[id]=invhash(cur,k);
        if(invhash(rcur,k)<val[id]) val[id]=invhash(rcur,k),iskmersymm[id]=1,k_mer_pos[id]=0+k-1;
        for(int i=1;i+k<=read_len;i++){
            ull tmp=trans(str[i+k-1]);
            cur=(cur<<2|tmp)&mask;
            rcur=rcur>>2|(tmp^3ull)<<shift1;
            //if(cur==rcur) continue; // skip "symmetric k-mers" as we don't know it strand
            if(invhash(cur,k)<val[id]) val[id]=invhash(cur,k),k_mer_pos[id]=i,iskmersymm[id]=0;
            if(invhash(rcur,k)<val[id]) val[id]=invhash(rcur,k),k_mer_pos[id]=i+k-1,iskmersymm[id]=1;
        }
    }
    return (void*)0;
}
    

void cal_k_minimizer(int k)//cal all read in [0,rcnt],with arg=k
{
    int tmp;
    
    mutex_init();
    thread_tmp.assign(thread_num+1,k);
    assert(k<=read_len);
    k_mer_pos.assign(rcnt0,-1);
    val.assign(rcnt0,-1);
    iskmersymm.assign(rcnt0,-1);
    for(int i=0;i<thread_num;i++){
        if((tmp=pthread_create(&thread_id[i],NULL,k_minimizer,&thread_tmp[i]))!=0){
            cout<<"pthread_create ERROR\n";
            exit(0);
        }
    }
    k_minimizer(&thread_tmp[thread_num]);
    for(int i=0;i<thread_num;i++){
        //cout<<i<<'\n';
        void* tmp;
        pthread_join(thread_id[i],&tmp);
    }

}

inline int find(int x)
{
    return pre[x]==x?x:pre[x]=find(pre[x]);
}

bool check(int id1,int id2,int threshold)//id1->id2,in the k-th round
{
    int tot=0;
    int k_mer_pos1=k_mer_pos[id1],issymm1=iskmersymm[id1];
    int k_mer_pos2=k_mer_pos[id2],issymm2=iskmersymm[id2];
    int start_pos1,step1,len1;
    int start_pos2,step2,len2;

    len1=cal_len(id1);
    if(issymm1==0){
        step1=1;
        start_pos1=0;
    }
    else{
        step1=-1;
        start_pos1=read_len-1;
    }

    len2=cal_len(id2);
    if(issymm2==0){
        step2=1;
        start_pos2=0;
    }
    else{
        step2=-1;
        start_pos2=read_len-1;
    }
    
    assert(len1>=len2);
    start_pos1+=step1*(len1-len2);
    assert(len1>=len2);
    int p1=start_pos1,p2=start_pos2;
    //cout<<p1<<' '<<p2<<' '<<step1<<' '<<step2<<' '<<'\n';
    while(p1>=0&&p1<read_len&&p2>=0&&p2<read_len){
        if(symm(vread[id1][p1],issymm1)!=symm(vread[id2][p2],issymm2)) tot++;
        p1+=step1;
        p2+=step2;
    }
    //cout<<"check: "<<id1<<' '<<id2<<' '<<tot<<' '<<threshold<<' '<<vecr[id1].str<<' '<<vecr[id2].str<<'\n';
    return tot<=threshold;
}


#endif