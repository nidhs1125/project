#ifndef setting_h
#define setting_h
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef unsigned long long ull;
#define pic pair<int,char>
#define pll pair<ll,ll>
#define pii pair<int,int>
#define puu pair<ull,ull>
#define pb push_back

//announcements
class Read;
class Contig;

void read_align(int k);
void compmain(ifstream& fin,ofstream& fout);
void comp();
void decompmain(ifstream& fin,ofstream& fout);
void contig_make();
int cal_k_mer(int k,int id,ull& val);
void SCS_gen();
void encode(ofstream& fout);
ull invhash(ull key, int k);
char rtrans(int c);
int trans(char c);
void prework();
puu get_hash_val(string& s);
char symm(char c);
string cal_symm(string& s);

//consts
const int maxk=32;
int threshold=5;//threshold
#define max_str_length 210
const ull step1=131;
const ull step2=1331;
const ull mod1=998244353;
const ull mod2=1000000007;



//variables
int ccnt=0;//contig cnt
int rcnt=0;//read cnt
vector<Read> vecr;//store all reads
vector<Contig> vecc;
string ans="";
vector<ull> hash_base1,hash_base2;
vector<ull> invhash_base1,invhash_base2;
vector<vector<int> > vec_id;//order to id
vector<vector<int> > vec_order;//id to order
vector<vector<int> > vec_pos;//id to k_mer_pos
vector<int> vec_k;//k of every round


//class
class Read//save a single read and k-mer
{
    public:
    Read(){
        cid=-1;
        isrepeat=0;
        str="",k_mer_pos=0;
    }
    Read(int cnt,string str1){
        cid=-1;
        isrepeat=0;
        rid=cnt;
        str=str1;
    }
    void init(int k){
        k_mer_pos=cal_k_mer(k,rid,val);
    }
    string str;//read itself
    int k_mer_pos;//a k_mer_pos corresponding to a single k
    ull val;//hash_value of k_minimizer
    int rid;//read_id
    int cpos;
    int cid;//contig_id;
    vector<pic> dismatch;

    bool isrepeat;
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
    int k_mer_pos;
    int spos;
};

struct spre
{
    int id;
    ull val1;
    ull val2;
};





#endif