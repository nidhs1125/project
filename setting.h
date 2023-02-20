#ifndef setting_h
#define setting_h
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef unsigned long long ull;
#define pb push_back
#define max_str_length 210

//announcements
class Read;
class Contig;

void contig_gen();
void compmain(ifstream& fin,ofstream& fout);
void comp();
void decompmain(ifstream& fin,ofstream& fout);
void encoding(int l,int r);
void contig_make(int l,int r);
int cal_k_mer(int k,string& s);
ll get_hash_value(string s);
void init();
void SCS_gen();
void encode(ofstream& fout);

//consts
#define mod 998244353ll
const ll step=131;
int threshold=5;//threshold
#define pir pair<int,char>

//variables
int ccnt=0;//contig cnt
int rcnt=0;//read cnt
vector<Read> vecr;//store all reads
vector<Contig> vecc;
int hash_base[max_str_length];
int invhash_base[max_str_length];
string ans="";


//class
class Read//save a single read and k-mer
{
    public:
    Read(){
        str="",k_mer_pos=0;
        cid=-1;
    }
    Read(int cnt,string str1){
        rid=cnt;
        str=str1;
        cid=-1;
    }
    void init(int k){
        k_mer_pos=cal_k_mer(k,str);
        val=get_hash_value(str.substr(k_mer_pos,k));
    }
    string str;//read itself
    int k_mer_pos;
    ll val;//hash_value
    int rid;//read_id
    int rpos;
    int cid;//contig_id;
    vector<pir> dismatch;
};

class Contig
{
    public:
    Contig(){
        str="";
    }
    Contig(int cnt,string str1){
        cid=cnt;
        str=str1;
    }
    string str;
    int cid;
    int k_mer_pos;
    int spos;
};







#endif