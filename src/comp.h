#include<bits/stdc++.h>
#include<unistd.h>
#include"tools.h"
using namespace std;




void compmain(ifstream& fin,ofstream& fout)
{
    init();

    string tmp;
    int cnt=0;
    vecr.clear();
    while(getline(fin,tmp)){
        if(tmp[0]=='$') break;
        cnt++;
        if(cnt%4==2){
            vecr.pb(Read(rcnt++,tmp));
        }
        //if(cnt>=1000000) break;//read 10000 lines
    }
    cout<<cnt<<"lines in total\n";
    //exclude the same reads
    prework();
    //接下来需要参考minicom生成contig
    comp();
    cout<<cnt<<"lines,"<<" compression finish"<<'\n';
    encode(fout);

}

void prework()
{
    cout<<"start prework\n";
    vector<spre> vec;
    for(int i=0;i<rcnt;i++){
        puu tmp=get_hash_val(vecr[i].str);
        vec.pb(spre{i,tmp.first,tmp.second});
    }
    sort(vec.begin(),vec.end(),cmp2);
    for(int i=0;i<rcnt;i++){
        if(vecr[vec[i].id].isrepeat!=0) continue;
        vecr[vec[i].id].isrepeat=-1;
        //try to find repeat read
        int j=i+1;
        while(j<rcnt&&vec[j].val1==vec[i].val1&&vec[j].val2==vec[i].val2){
            if(vecr[vec[j].id].isrepeat!=0){
                j++;
                continue;
            }
            repeatcnt++;
            vecr[vec[j].id].isrepeat=1;
            vecr[vec[j].id].issymmrepeat=0;
            vecr[vec[j].id].repeatid=vec[i].id;
            j++;
        }
        //try to find symmetric repeat read
        string tmps=cal_symm(vecr[vec[j].id].str);
        puu tmpp=get_hash_val(tmps);
        int l=0,r=rcnt,d=0;
        while(l<r){
            d=l+r>>1;
            if(tmpp.first==vec[d].val1){
                if(tmpp.second<=vec[d].val2) r=d;
                else l=d+1;
            }
            else if(tmpp.first<vec[d].val1) r=d;
            else if(tmpp.first>vec[d].val1) l=d+1;
        }
        j=l;
        while(j<rcnt&&tmpp.first==vec[j].val1&&tmpp.second==vec[j].val2){
            if(vecr[vec[j].id].isrepeat!=0){
                j++;
                continue;
            }
            repeatcnt++;
            vecr[vec[j].id].isrepeat=1;
            vecr[vec[j].id].issymmrepeat=1;
            vecr[vec[j].id].repeatid=vec[i].id;
            j++;
        }
    }
    sort(vecr.begin(),vecr.end(),cmp3);
    rcnt-=repeatcnt;
}

void comp()
{
    //read_align(4);
    read_align(20);
    //read_align(28);
    //read_align(27);
    // for(auto &it :vecr){
    //     cout<<it.str<< ' '<<it.k_mer_pos<<'\n';
    // }
    cout<<"read_align finish\n";
    contig_make();
    cout<<"contig_make finish\n";
    // for(auto &it :vecc){
    //     cout<<it.cid<<' '<<it.str<<'\n';
    // }
    //SCS is concatenation of all contigs
    SCS_gen();
    cout<<"scs_gen finish\n";
    //cout<<ans<<'\n';
}

void read_align(int k)
{
    cout<<"read_align start\n";
    for(int i=0;i<rcnt;i++){
        vecr[i].init(k);
    }
    vector<int> order_to_id;//cur vec
    for(int i=0;i<rcnt;i++) order_to_id.pb(i);
    sort(order_to_id.begin(),order_to_id.end(),cmp1);
    // cout<<"k="<<k<<' '<<"order:\n";
    // for(auto it:order_to_id){
    //     cout<<it<<' ';
    // }
    //cout<<'\n';
    vector<int> id_to_order(rcnt,0);
    for(int i=0;i<rcnt;i++) id_to_order[order_to_id[i]]=i;
    vector<int> k_pos(rcnt,0);
    for(int i=0;i<rcnt;i++) k_pos[i]=vecr[i].k_mer_pos;
    //write down the result of this round
    vec_id.pb(order_to_id);
    vec_order.pb(id_to_order);
    vec_pos.pb(k_pos);
    vec_k.pb(k);

    /*
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
    */
    
}

/*

according to the nxtid recored in read,choose the best one to continue.


*/
void contig_make()
{
    cout<<"contig_make start\n";
    //vecb:bias between cur read and the next read
    //bias need to be small
    vector<int> pre(rcnt,-1),nxt(rcnt,-1),vecb(rcnt,max_str_length);
    for(int i=0;i<rcnt;i++){//every read
        //for every read_align, check the next 10 reads to form nxtid(id,bias)
        vector<pii> nxtid;
        for(int j=0;j<vec_k.size();j++){//for every k
            int mypos=vec_order[j][i];
            int pbias=0;
            int bias;
            for(int p=mypos+1;p<rcnt&&p<mypos+10;p++){//check the read in pos=p
                int id2=vec_id[j][p];
                int id1=i;
                bias=vec_pos[j][id1]-vec_pos[j][id2];
                if(bias<pbias) break;//impossible
                pbias=bias;
                ///cout<<"+++"<<p<<' '<<mypos<<' '<<id1<<' '<<id2<<' '<<bias<<' '<<check(id1,id2,bias)<<'\n';
                if(check(id1,id2,bias)) nxtid.pb(pii(id2,bias));
            }
        }
        // cout<<"nxtid:|||"<<i<<'\n';
        // for(auto it:nxtid) cout<<it.first<<' '<<it.second<<'\n';
        //according to nxtid, find the best one to concat
        map<pii,int> mp;
        for(pii it:nxtid){//id and bias
            //cout<<"|||"<<i<<' '<<it.first<<' '<<it.second<<'\n';
            if(pre[it.first]!=-1) continue;
            mp[it]++;
        }
        int fre=0,bias=0,nid=-1;
        //first, check the diff is less than threshold
        //choose the most frequently occurent one
        //if many, choose the bias min one.
        //if many, choose any one
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
            vecb[i]=bias;
        }
    }
    // cout<<"====bias\n";
    // for(int i=0;i<rcnt;i++){
    //     cout<<i<<' '<<pre[i]<<' '<<nxt[i]<<' '<<vecb[i]<<'\n';
    // }

    vector<int> seq;
    vector<int> vis(rcnt,0);
    for(int t=0;t<rcnt;t++){
        //the reads form a chain or a ring
        //if form a ring, find the max bias pos to start
        if(vis[t]==0){
            seq.clear();
            int cur=t;
            int mx=vecb[t],mxpos=t;
            while(pre[cur]!=-1&&pre[cur]!=t){
                cur=pre[cur];
                if(vecb[cur]>mx) mx=vecb[cur],mxpos=cur;
            }
            if(pre[cur]==-1){//chain
                while(cur!=-1){
                    seq.pb(cur);
                    vis[cur]=1;
                    cur=nxt[cur];
                }
            }
            else{//ring
                cur=nxt[mxpos];
                while(1){
                    seq.pb(cur);
                    vis[cur]=1;
                    cur=nxt[cur];
                    if(cur==mxpos) break;
                }
            }

            
            int l=0,r=seq.size()-1;
            // cout<<"seq:  ";
            // for(int i=l;i<=r;i++){
            //     cout<<seq[i]<<' ';
            // }
            // cout<<'\n';
            vecc.push_back(Contig());
            Contig &con=vecc.back();
            con.cid=ccnt++;
            int LEN=0,curb=0;
            for(int i=l;i<=r;i++){
                LEN=max(LEN,curb+(int)vecr[seq[i]].str.length());
                curb=curb+vecb[seq[i]];
            }
            for(int j=0;j<LEN;j++){//for every pos
                int cnt[5]={0};
                curb=0;
                for(int i=l;i<=r;i++){//for every read
                    if(curb<=j&&(int)vecr[seq[i]].str.length()-1+curb>=j){
                        cnt[trans(vecr[seq[i]].str[j-curb])]++;
                    }
                    curb+=vecb[seq[i]];
                }
                int mx=cnt[4],mxpos=4;
                for(int i=3;i>=0;i--) if(cnt[i]>=mx) mx=cnt[i],mxpos=i;
                con.str+=rtrans(mxpos);
            }
            curb=0;
            for(int i=l;i<=r;i++){//for every read
                //calculate the hamming distance and setting cid
                //int cnt=0;
                
                // for(int j=0;j<vecr[i].str.length();j++){
                //     if(vecr[i].str[j]!=con.str[j+curb]) cnt++;
                // }
                //if(cnt<=threshold){//join in and record the dismatch
                    vecr[seq[i]].cid=con.cid;
                    vecr[seq[i]].cpos=curb;
                    for(int j=0;j<vecr[seq[i]].str.length();j++){
                        if(vecr[seq[i]].str[j]!=con.str[j+curb]){
                            vecr[seq[i]].dismatch.pb(pii(j,vecr[seq[i]].str[j]));
                        }
                    }
                //}
                //else hamming distance too large ,don't join in 
                curb+=vecb[seq[i]];
            }
            //cout<<ccnt-1<<' '<<con.str<<'\n';
            if(ccnt<=100) cout<<ccnt<<' '<<seq.size()<<'\n';
        }
    }

    
}

void SCS_gen()
{
    cout<<"SCS_gen start\n";
    int cbias=0;
    for(int i=0;i<ccnt;i++){
        ans+=vecc[i].str;
        vecc[i].spos=cbias;
        cbias+=vecc[i].str.length();
    }
}
/*
//old version 
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


*/


void encode(ofstream& fout)
{
    cout<<"encoding...\n";
    // for(auto it:vecr) cout<<it.rid<<' ';
    // cout<<'\n';
    fout<<ans.length();
    cout<<ans.length()<<' '<<ans.length()*4/1024/26/1024<<'\n';
    cout<<rcnt<<' '<<ccnt<<'\n';
    for(int i=0;i<ans.length();i+=26){
        ull tmp=0;
        ull base=1;
        for(int j=0;j+i<ans.length()&&j<26;j++){
            tmp=tmp+trans(ans[i+j])*base;
            base*=5;
        }
        fout<<tmp;
    }
    fout<<rcnt;
    vector<int> curpos(rcnt+repeatcnt,-1);
    //for every read,show the pos
    for(int i=0;i<rcnt;i++){
        Read& it=vecr[i];
        curpos[vecr[i].rid]=i;
        fout<<(int)(vecc[it.cid].spos+it.cpos);//pos;
        //cout<<it.str<<' '<<vecc[it.cid].spos+it.cpos<<'\n';
        fout<<(char)it.dismatch.size();//at most 5
        if(it.dismatch.size()!=0){
            ull tmp=0;//5 dismatch at most when k=5,so char is enough 
            ull base=1;
            for(auto it1:it.dismatch){
                fout<<(char)it1.first;
                tmp+=base*trans(it1.second);
                base*=5;
            }
            fout<<(char)tmp;
        }
    }
    fout<<repeatcnt;
    for(int i=rcnt;i<rcnt+repeatcnt;i++){
        if(curpos[vecr[i].repeatid]==-1){
            cout<<i<<' '<<rcnt<<' '<<repeatcnt<<' '<<vecr.size()<<'\n';
            cout<<vecr[i].repeatid<<' '<<vecr[i].isrepeat<<'\n';
        }
        assert(curpos[vecr[i].repeatid]!=-1);
        //cout<<vecr[i].repeatid<<'\n';
        fout<<curpos[vecr[i].repeatid];
        fout<<(char)vecr[i].issymmrepeat;
        int cnt=0;
        for(int j=0;j<vecr[i].str.length();j++){
            if(vecr[i].str[j]=='N') cnt++;
        }
        fout<<(char)cnt;
        for(int j=0;j<vecr[i].str.length();j++){
            if(vecr[i].str[j]=='N') fout<<(char)j;
        }

    }
}