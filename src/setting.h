#ifndef setting_h
#define setting_h
#include<bits/stdc++.h>
#include <pthread.h>  
#include<unistd.h>
#include <sys/time.h>
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
class Edge;

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
int cal_len(int i);
void bas_align(int k);
void contig_make();
int realignment(int id,int k,int& is);
void SCS_gen();
void* cal_block(void* arg);
void realignment_prework(int k);

//consts
#define testflag 1
const int maxk=32;
int threshold=2;//threshold in read without bias
int threshold_ratio_re=5;//threadhold in realignment=readlen/threshold_ratio_re
int thread_num=19;//total thread is thread_num+1
int order_preserve=0;//if order preserved;
#define max_str_length 210
const ull step1=131;
const ull step2=1331;
const ull mod1=998244353;
const ull mod2=1000000007;
const int k_num=30;//the number of read considered
int read_len=0;
const int thread_bias=1333;//
const int max_bias=50;//max bias between reads


//variables
timeval start1, end1; 
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
vector<vector<int> > basket_rev;//record the reverse state, it is guaranteed that the first read is positive(0)
//in one basket
pthread_mutex_t pthread_id_mutex = PTHREAD_MUTEX_INITIALIZER;
int pthread_id;//cur cid
vector<int> order_to_id;//cur order mapping to its id
//between baskets
int bas_num;
pthread_mutex_t bas_edge_mutex = PTHREAD_MUTEX_INITIALIZER;
vector<vector<Edge> > bas_edge;//the first place is out_edge(only one) and the other are in_edge
vector<int> first_rev;//if the basket is the root, is it reverse

//make contig
pthread_mutex_t vis_id_mutex = PTHREAD_MUTEX_INITIALIZER;
int vis_id;

pthread_mutex_t vis_mutex = PTHREAD_MUTEX_INITIALIZER;
vector<int> vis;
//realignment
unordered_map<ull,vector<int>> mp_pos;




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
    int iskmersymm;// is k_mer reverse

    int isrev;//is reverse in the contig
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
        num=0;
        is=0;
        dismatch.clear();
    }
    Contig(int cnt,string str1=""){
        cid=cnt;
        str=str1;
        num=0;
        is=0;
        dismatch.clear();
    }
    string str;
    int cid;
    int spos;
    int num;
    int is;//if reverse in ans
    vector<pic> dismatch;
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

class Edge
{
    public:
    Edge(int a=-1,int c=read_len+max_bias,int d=0){
        to=a;
        bias=c;
        rev=d;
    }
    int to;
    int bias;
    int rev;
    bool update(Edge& e){
        if(bias>e.bias){
            *this=e;
            return 1;
        }
        return 0;
    }
};



#endif