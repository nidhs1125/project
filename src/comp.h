#include<bits/stdc++.h>
#include<unistd.h>
#include"tools.h"
using namespace std;




void compmain(ifstream& fin,ofstream& fout)
{
    string tmp;
    int cnt=0;
    vecr.clear();
    while(getline(fin,tmp)){
        cnt++;
        if(cnt%4==2){
            vecr.pb(Read(rcnt++,tmp));
        }
        //if(cnt>=1000000) break;//read 10000 lines
    }
    cout<<cnt<<"lines in total\n";
    //接下来需要参考minicom生成contig
    comp();
    cout<<cnt<<"lines,"<<" compression finish"<<'\n';
    encode(fout);

}


void comp()
{
    read_align(4);
    cout<<"read_align finish\n";
    contig_make();
    cout<<"contig_make finish\n";
    // for(auto &it :vecc){
    //     cout<<it.str<<'\n';
    // }
    SCS_gen();
    sort(vecc.begin(),vecc.end(),cmp4);
    cout<<"SCS_gen finish\n";
    for(auto &it:vecc){
        cout<<it.str<<' '<<it.spos<<'\n';
    }
}

void read_align(int k)
{
    cout<<"read_align start\n";
    for(int i=0;i<rcnt;i++){
        vecr[i].init(k);
    }
    vector<int> veccur;//cur vec
    for(int i=0;i<rcnt;i++) veccur.pb(i);
    sort(veccur.begin(),veccur.end(),cmp1);
    int n=rcnt;
    int i=0;
    while(i<n){
        int j=i+1;
        while(j<n&&vecr[veccur[i]].val==vecr[veccur[j]].val) j++;
        //reads in section [i,j] has the same k_mer.
        //consider the next 5 reads,to choose the best one or more
        //condition:hamming distance smaller than the threshold=5
        for(int p=i;p<j;p++){
            for(int q=p+1;q<j&&q<p+5;q++){
                int id1=veccur[p],id2=veccur[q];
                int cnt=0;//total diff
                for(int j=0;j<min(vecr[id2].str.length(),vecr[id1].str.length()+vecr[id2].k_mer_pos-vecr[id1].k_mer_pos);j++){
                    if(vecr[id2].str[j]!=vecr[id1].str[j+vecr[id2].k_mer_pos-vecr[id1].k_mer_pos]) cnt++;
                }
                if(cnt<threshold) vecr[id1].nxtid.pb(pii(id2,vecr[id1].k_mer_pos-vecr[id2].k_mer_pos));
            }
        }
        
        i=j;
    }
    
}

/*

according to the nxtid recored in read,choose the best one to continue.
choose the most frequently occurent one
if many, choose the bias min one.
if many, choose any one

*/
void contig_make()
{
    cout<<"contig_make start\n";
    vector<int> pre(rcnt,-1),nxt(rcnt,-1);
    for(int i=0;i<rcnt;i++){
        map<pii,int> mp;
        for(pii it:vecr[i].nxtid){//id and bias
            //cout<<"|||"<<i<<' '<<it.first<<' '<<it.second<<'\n';
            if(pre[it.first]!=-1) continue;
            mp[it]++;
        }
        int fre=0,bias=0,nid=-1;
        for(auto& it:mp){
            if(it.second>fre){
                fre=it.second;
                nid=it.first.first;
                bias=it.first.second;
            }
            else if(it.second==fre&&it.first.second<bias){
                nid=it.first.first;
                bias=it.first.second;
            }
        }
        if(nid!=-1){
            nxt[i]=nid;
            pre[nid]=i;
        }
    }

    vector<int> seq;

    for(int t=0;t<rcnt;t++){
        if(pre[t]==-1){//start
            seq.clear();
            int cur=t;
            while(cur!=-1){
                seq.pb(cur);
                cur=nxt[cur];
            }
            int l=0,r=seq.size()-1;
            vecc.push_back(Contig());
            Contig &con=vecc.back();
            con.cid=ccnt++;
            int L=1e9,R=-1e9;
            for(int i=l;i<=r;i++){
                L=min(L,-vecr[i].k_mer_pos);
                R=max(R,(int)vecr[i].str.length()-1-vecr[i].k_mer_pos);
            }
            con.k_mer_pos=-L;
            for(int j=L;j<=R;j++){//for every pos
                int cnt[5]={0};
                for(int i=l;i<=r;i++){//for every read
                    if(-vecr[i].k_mer_pos<=j&&(int)vecr[i].str.length()-1-vecr[i].k_mer_pos>=j){
                        cnt[trans(vecr[i].str[vecr[i].k_mer_pos+j])]++;
                    }
                }
                int mx=cnt[4],mxpos=4;
                for(int i=3;i>=0;i--) if(cnt[i]>=mx) mx=cnt[i],mxpos=i;
                con.str+=rtrans(mxpos);
            }
            for(int i=l;i<=r;i++){//for every read
                //calculate the hamming distance and setting cid
                int cnt=0;
                for(int j=0;j<vecr[i].str.length();j++){
                    if(vecr[i].str[j]!=con.str[j-vecr[i].k_mer_pos+con.k_mer_pos]) cnt++;
                }
                //if(cnt<=threshold){//join in and record the dismatch
                    vecr[i].cid=con.cid;
                    vecr[i].rpos=con.k_mer_pos-vecr[i].k_mer_pos;
                    for(int j=0;j<vecr[i].str.length();j++){
                        if(vecr[i].str[j]!=con.str[j-vecr[i].k_mer_pos+con.k_mer_pos]){
                            vecr[i].dismatch.pb(pii(j,con.str[j-vecr[i].k_mer_pos+con.k_mer_pos]));
                        }
                    }
                //}
                //else hamming distance too large ,don't join in 
            }
            //cout<<ccnt<<' '<<con.str<<'\n';
        }
    }

    
}

/*
reference to PgRC
for all contig, find the SCS

*/

void SCS_gen()
{
    cout<<"SCS make start\n";
    //length:200 -> 0
    //from back to front, calculate every suffix index for every contig.
    int n=vecc.size();
    //front alignment
    vector<vector<int> > conrnk;//rnkpos->conid
    vector<int> prernk;//the last conid->rnk
    prernk.assign(n,1);
    conrnk.clear();
    conrnk.resize(max_str_length);
    conrnk[0].resize(n);
    for(int i=0;i<n;i++) conrnk[0][i]=i;
    for(int i=1;i<max_str_length;i++){//consider suffix of length i
        //cal rnk of this rnd
        //A
        for(int j=0;j<conrnk[i-1].size();j++){//con
            int p=conrnk[i-1][j];
            if(vecc[p].str.length()>=i&&vecc[p].str[vecc[p].str.length()-i]=='A') conrnk[i].pb(p);
        }
        //C
        for(int j=0;j<conrnk[i-1].size();j++){//con
            int p=conrnk[i-1][j];
            if(vecc[p].str.length()>=i&&vecc[p].str[vecc[p].str.length()-i]=='C') conrnk[i].pb(p);
        }
        //G
        for(int j=0;j<conrnk[i-1].size();j++){//con
            int p=conrnk[i-1][j];
            if(vecc[p].str.length()>=i&&vecc[p].str[vecc[p].str.length()-i]=='G') conrnk[i].pb(p);
        }
        //T
        for(int j=0;j<conrnk[i-1].size();j++){//con
            int p=conrnk[i-1][j];
            if(vecc[p].str.length()>=i&&vecc[p].str[vecc[p].str.length()-i]=='T') conrnk[i].pb(p);
        }
        //N
        for(int j=0;j<conrnk[i-1].size();j++){//con
            int p=conrnk[i-1][j];
            if(vecc[p].str.length()>=i&&vecc[p].str[vecc[p].str.length()-i]=='N') conrnk[i].pb(p);
        }
        if(conrnk[i].size()==0) break;

    }
    cout<<"SCS make part_1 finish\n";
    //for suffix length from 200 to 0, determine the best one
    vector<int> pre(n,-1),nxt(n,-1),len(n,0);
    for(int i=max_str_length-1;i>=0;i--){
        if(conrnk[i].size()==0) continue;
        int pl=0,pr=0;
        while(pl<n&&pr<conrnk[i].size()){
            if(vecc[pl].str.length()<i){
                pl++;
                continue;
            }
            if(pre[pl]!=-1){
                pl++;
                continue;
            }
            if(nxt[conrnk[i][pr]]!=-1){
                pr++;
                continue;
            }
            if(conrnk[i][pr]==pl){//the same string
                //compare the next string.
                //to determine which increase first,pl or pr.
                if(pl+1>=n) pr++;
                else if(pr+1>=conrnk[i].size()) pl++;
                else{
                    int tmp=vecc[pl+1].str.substr(0,i).compare(vecc[conrnk[i][pr+1]].str.substr(vecc[conrnk[i][pr+1]].str.length()-i,i));
                    if(tmp<=0) pl++;//if next str still equal, then pl++ or pr++ don't matter.
                    else pr++;
                }
                continue;
            }
            int tmp=vecc[pl].str.substr(0,i).compare(vecc[conrnk[i][pr]].str.substr(vecc[conrnk[i][pr]].str.length()-i,i));
            if(tmp==0){
                pre[pl]=conrnk[i][pr];
                nxt[conrnk[i][pr]]=pl;
                len[pl]=i;
                pl++,pr++;
            }
            else if(tmp<0) pl++;
            else pr++;
        }
    }
    cout<<"SCS make part_2 finish\n";
    //make SCS
    vector<int> vis(n,0);
    int bias=0;//cur length of ans;
    for(int i=0;i<n;i++){
        if(vis[i]==1) continue;
        int pre_pos=i;
        int cur_pos=pre[pre_pos];
        int min=len[pre_pos],min_pos=pre_pos;
        vis[pre_pos]=1;
        while(cur_pos!=-1&&cur_pos!=i){
            vis[cur_pos]=1;
            if(min>len[cur_pos]) min=len[cur_pos],min_pos=cur_pos;
            pre_pos=cur_pos;
            cur_pos=pre[cur_pos];
        }
        //if form a ring start from the min_pos of overlap length
        //start from pre_pos,end at i or -1;
        int start_pos=min_pos;
        cur_pos=start_pos;
        string tmp=vecc[cur_pos].str;
        vecc[cur_pos].spos=bias;
        cur_pos=nxt[cur_pos];
        while(cur_pos!=-1&&cur_pos!=start_pos){
            vecc[cur_pos].spos=bias+tmp.length()-len[cur_pos];
            tmp+=vecc[cur_pos].str.substr(len[cur_pos]);//exclude the overlap segment
            cur_pos=nxt[cur_pos];
        }
        ans+=tmp;
    }
    //cout<<ans<<'\n';
}

void encode(ofstream& fout)
{
    cout<<"encoding...\n";
    fout<<ans.length();
    for(int i=0;i<ans.length();i+=26){
        ull tmp=0;
        ull base=1;
        for(int j=0;j+i<ans.length()&&j<26;j++){
            tmp=tmp+trans(ans[i+j])*base;
            base*=5;
        }
        fout<<tmp;
    }
    fout<<vecr.size();
    //for every read,show the pos
    sort(vecr.begin(),vecr.end(),cmp3);
    
    for(Read& it:vecr){
        fout<<vecc[it.cid].spos+it.rpos;//pos;
        //cout<<it.str<<' '<<vecc[it.cid].spos+it.rpos<<'\n';
        fout<<it.dismatch.size();
        if(it.dismatch.size()!=0){
            ull tmp=0;
            ull base=1;
            for(auto it1:it.dismatch){
                fout<<it1.first;
                tmp+=base*trans(it1.second);
                base*=5;
            }
            fout<<tmp;
        }
    }
}