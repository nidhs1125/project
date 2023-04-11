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
struct srealign;

//consts
#define testflag 1
const int maxk=31;
int threshold1=2;//threshold in read without bias
int threshold2=2;//threshold in read with bias
int threshold_ratio_re=5;//threadhold in realignment=readlen/threshold_ratio_re
int maxdiscnt;
int thread_num=4;//total thread is thread_num+1
int order_preserve=0;//if order preserved;
#define max_str_length 210
const ull step1=131;
const ull step2=1331;
const ull mod1=998244353;
const ull mod2=1000000007;
const int k_num=30;//the number of read considered
int read_len=0;
const int max_bias=50;//max bias between reads


//variables
timeval start1, end1; 
vector<string> vread;//to reserve all read
int ccnt=0;//contig cnt
int rcnt=0,rcnt0=0;//read cnt
int repeatcnt=0;//repeat number of reads
int ncnt;//number of read with 'N'

//read related
vector<int> vecr;//store all reads' cur id
vector<int> k_mer_pos;
vector<ull> val;
vector<int> iskmersymm;
vector<int> cpos;
vector<vector<pic>> dismatch;
vector<int> isrev;
vector<int> cid;

//mutex related
pthread_mutex_t my_mutex;
pthread_mutex_t array_mutex;
int my_id;
int id_upbound;
vector<pthread_t> thread_id;
vector<int> thread_tmp;

//cal repeat read:
vector<ull> hash_base1,hash_base2;
vector<ull> invhash_base1,invhash_base2;
vector<puu> hash_val; 
vector<int> isrepeat;
vector<int> issymmrepeat;
vector<int> repeatid;//repeatid:the read id of the same read with read i
//cal read with n
vector<int> hasn;

//read_align
vector<int> pre;//disjoint set
vector<vector<int> > basket;
vector<vector<int> > basket_rev;//record the reverse state, it is guaranteed that the first read is positive(0)
vector<int> first_rev;
//bas_align
vector<vector<Edge> > bas_edge;//the first place is out_edge(only one) and the other are in_edge
vector<int> nxt_bas;

//make contig
vector<int> vis;
vector<Contig> vecc;
string ans="";//the gene
//realignment
unordered_map<ull,vector<srealign>> mp_pos;//hash_value->vector of read_id 


//class

class Contig
{
    public:
    Contig(){
        str="";
    }
    string str;
    int spos;
};

struct srealign
{
    int id;
    int pos;
    int isrev;
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