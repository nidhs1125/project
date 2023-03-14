#ifndef setting_h
#define setting_h
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef unsigned long long ull;
typedef unsigned int uint;
#define pic pair<int,char>
#define pll pair<ll,ll>
#define pii pair<int,int>
#define puu pair<ull,ull>
#define pb push_back

//announcements
class Read;
class Contig;

void read_align(int k);
void compmain(ifstream& fin,string& out_path);
void comp();
void decompmain(string& in_path,ofstream& fout);
void cal_k_mer(int k,int id,int& k_mer_pos,ull& val,int& iskmersymm);
void encode(string& out_path);
ull invhash(ull key, int k);
int trans(char c);
void prework();
puu get_hash_val(string& s);
string cal_symm(string& s);
char symm(char c,bool is);
void* cal_pre(void* arg);
void* cal_nxt(void* arg);

//consts
#define testflag 1
const int maxk=32;
int threshold=2;//threshold
int thread_num=20;//
int order_preserve=0;//if order preserved;
#define max_str_length 210
const ull step1=131;
const ull step2=1331;
const ull mod1=998244353;
const ull mod2=1000000007;
const int k_num=30;//the number of read considered
int read_len=0;
const int thread_bias=1333;//


//variables
int rndk=0;//the rnd of k
int ccnt=0;//contig cnt
int rcnt=0;//read cnt
int repeatcnt=0;//repeat number of reads
vector<Read> vecr;//store all reads
vector<Contig> vecc;
string ans="";
vector<ull> hash_base1,hash_base2;
vector<ull> invhash_base1,invhash_base2;
pthread_mutex_t pre_mutex = PTHREAD_MUTEX_INITIALIZER;
vector<int> pre;//disjoint set
vector<vector<int> > basket;
//in one basket
pthread_mutex_t pthread_id_mutex = PTHREAD_MUTEX_INITIALIZER;
int pthread_id;//cur cid
vector<int> order_to_id;//cur order mapping to its id
//between basket
int bas_num;
pthread_mutex_t pre_bas_mutex = PTHREAD_MUTEX_INITIALIZER;
vector<int> nxt_bas;
vector<int> pre_bas;
vector<int> bias_bas;

//class
class Read//save a single read and k-mer
{
    public:
    Read(){
        cid=-1;
        isrepeat=0;
        str="";
    }
    Read(int cnt,string str1){
        cid=-1;
        isrepeat=0;
        rid=cnt;
        str=str1;
    }
    string str;//read itself
    int k_mer_pos;
    ull val;
    int isrev;
    int rid;//read_id
    int cpos;
    int cid;//contig_id;
    
    int pre;
    vector<pic> dismatch;

    int isrepeat;
    bool issymmrepeat;
    int repeatid;
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
    int spos;
};

struct spre
{
    int id;
    ull val1;
    ull val2;
};

struct sconmake
{
    int isrev;
    int cnt;
    int k;//the rnd
};



#endif