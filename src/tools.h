#ifndef tools_h
#define tools_h
#include<bits/stdc++.h>
#include"setting.h"
//using namespace std;

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

char trans(int c)//when compressing, use only 0-3
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
        if(s[p]=='N') tmp+=s[p];//if N,remain the same to guarantee the same 
        else tmp+=symm(s[p],1);
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
                ret='A';//'N' is considered as 'T',so its symmetry is 'A'
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
void cal_k_mer(int k,int id,int& k_mer_pos,ull& val,int& iskmersymm)
{
    int i=0;
    k_mer_pos=0;
    val=0;
    ull cur=0;
    ull rcur=0;
    int shift1=(2*k-2);
    ull mask=(1ull<<2*k)-1;
    iskmersymm=0;
    if(vecr[id].str.find("N")!=-1){//has N
        val=-1;//max,
        k_mer_pos=0;
    }

    for(int i=0;i<min(read_len,k);i++){
        int tmp=trans(symm(vecr[id].str[i],0));
        int tmpr=trans(symm(vecr[id].str[i],1));
        cur=((cur<<2)|tmp)&mask;
        rcur=(rcur>>2)|(tmp^3ull)<<shift1;
    }
    val=invhash(cur,k);
    if(invhash(rcur,k)<val) val=invhash(rcur,k),iskmersymm=1,k_mer_pos=0+k-1;
    for(int i=1;i+k<=read_len;i++){
        int tmp=trans(symm(vecr[id].str[i+k-1],0));
        int tmpr=trans(symm(vecr[id].str[i+k-1],1));
        cur=((cur<<2)|tmp)&mask;
        rcur=(rcur>>2)|(tmpr)<<shift1;
        if(cur==rcur) continue; // skip "symmetric k-mers" as we don't know it strand
        if(invhash(cur,k)<val) val=invhash(cur,k),k_mer_pos=i,iskmersymm=0;
        if(invhash(rcur,k)<val) val=invhash(rcur,k),k_mer_pos=i+k-1,iskmersymm=1;
    }
}

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


bool check(int id1,int id2)//id1->id2,in the k-th round
{
    int tot=0;
    int k_mer_pos1=vecr[id1].k_mer_pos,issymm1=vecr[id1].iskmersymm;
    int k_mer_pos2=vecr[id2].k_mer_pos,issymm2=vecr[id2].iskmersymm;
    int start_pos1,step1,len1,tag1;
    int start_pos2,step2,len2,tag2;

    tag1=vecr[id1].str.length();
    len1=cal_len(id1);
    if(issymm1==0){
        step1=1;
        start_pos1=0;
    }
    else{
        step1=-1;
        start_pos1=vecr[id1].str.length()-1;
    }

    tag2=vecr[id2].str.length();
    len2=cal_len(id2);
    if(issymm2==0){
        step2=1;
        start_pos2=0;
    }
    else{
        step2=-1;
        start_pos2=vecr[id2].str.length()-1;
    }
    
    assert(len1>=len2);
    start_pos1+=step1*(len1-len2);
    assert(len1>=len2);
    int p1=start_pos1,p2=start_pos2;
    //cout<<p1<<' '<<p2<<' '<<step1<<' '<<step2<<' '<<'\n';
    while(p1>=0&&p1<tag1&&p2>=0&&p2<tag2){
        if(symm(vecr[id1].str[p1],issymm1)!=symm(vecr[id2].str[p2],issymm2)) tot++;
        p1+=step1;
        p2+=step2;
    }
    //cout<<"check: "<<id1<<' '<<id2<<' '<<tot<<' '<<threshold<<' '<<vecr[id1].str<<' '<<vecr[id2].str<<'\n';
    return tot<=threshold;
}

bool check1(int id,int pos,int is)//vecc id and pos in ans and if symm
{
    int cthre=read_len/threshold_ratio_re;
    int cnt=0;
    if(pos<0||pos>ans.length()-read_len) return 0;
    if(!is){
        for(int i=0;i<read_len;i++){
            if(vecc[id].str[i]!=ans[i+pos]) cnt++;
        }
    }
    else{
        for(int i=0;i<read_len;i++){
            if(vecc[id].str[i]!=symm(ans[pos+read_len-i-1],1)) cnt++;
        }
    }
    if(cnt<=cthre) return 1;
    else return 0;
}


int cal_len(int i)
{
    int ret=0;
    if(vecr[i].iskmersymm==0) ret=vecr[i].k_mer_pos;
    else ret=vecr[i].str.length()-1-vecr[i].k_mer_pos;
    return ret;
}

bool cmp1(int id1,int id2)
{
    if(vecr[id1].val!=vecr[id2].val) return vecr[id1].val<vecr[id2].val;
    else return cal_len(id1)>cal_len(id2);
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

bool cmp4(Read& r1,Read& r2)
{
    if(r1.cpos!=r2.cpos) return r1.cpos<r2.cpos;
    else{
        int id1,id2;
        if(r1.isrepeat==1) id1=r1.repeatid;
        else id1=r1.rid;
        if(r2.isrepeat==1) id2=r2.repeatid;
        else id2=r2.rid;
        if(id1!=id2) return id1<id2;
        else return r1.isrepeat<r2.isrepeat;
    }
}


int in_int(ifstream& fin)
{
    char ch;
    int ret=0;
    for(int i=0;i<4;i++){
        ret=ret<<8;
        fin.read(&ch,1);
        ret|=(unsigned char)ch;
    }
    return ret;
}

unsigned char in_char(ifstream& fin)
{
    char ch;
    fin.read(&ch,1);
    return (unsigned char)ch;
}

void out_int(ofstream& fout,int out)
{
    fout<<(char)((out>>24)&0xff);
    fout<<(char)((out>>16)&0xff);
    fout<<(char)((out>>8)&0xff);
    fout<<(char)((out)&0xff);
}

void out_char(ofstream& fout,char out)
{
    fout<<out;
}

int find(int x)
{
    return pre[x]==x?x:pre[x]=find(pre[x]);
}

int dsu(int id1,int id2)
{
    assert(pre[id1]==id1);
    assert(pre[id2]==id2);
    if(id1==id2) return id1;
    if(basket[id1].size()<basket[id2].size()) swap(id1,id2);
    pre[id2]=id1;
    int tmprev=vecr[id1].iskmersymm^vecr[id2].iskmersymm;
    for(int i=0;i<basket[id2].size();i++){
        basket[id1].pb(basket[id2][i]);
        basket_rev[id1].pb(basket_rev[id2][i]^tmprev);
    }
    basket[id2].resize(0);
    return id1;
}

void print_time()
{
    gettimeofday(&end1, NULL);
    cout<<"time consumption:"<<1000*(end1.tv_sec - start1.tv_sec) + (end1.tv_usec - start1.tv_usec)/1000<<"ms\n";
}



#endif