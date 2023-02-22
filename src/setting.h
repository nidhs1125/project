#ifndef setting_h
#define setting_h
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef unsigned long long ull;
#define pic pair<int,char>
#define pll pair<ll,ll>
#define pii pair<int,int>
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

//consts
const int maxk=32;
int threshold=5;//threshold
#define max_str_length 210


//variables
int ccnt=0;//contig cnt
int rcnt=0;//read cnt
vector<Read> vecr;//store all reads
vector<Contig> vecc;
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
        int n=str.length();
        k_mer_pos=cal_k_mer(k,rid,val);
        pos[k]=k_mer_pos;
    }
    string str;//read itself
    int k_mer_pos;//a k_mer_pos corresponding to a single k
    int pos[maxk+1];
    ull val;//hash_value
    int rid;//read_id
    int rpos;
    int cid;//contig_id;
    vector<pic> dismatch;
    vector<pii> nxtid;//the possible next read after this read and the corresponding k
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