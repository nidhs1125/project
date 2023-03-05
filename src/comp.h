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
    cout<<"prework finish\n";
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
        string tmps=cal_symm(vecr[vec[i].id].str);
        puu tmpp=get_hash_val(tmps);
        //binary search
        int l=0,r=rcnt-1,d=0;
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
    #if testflag
    read_align(4);
    #else
    read_align(25);
    #endif
    //read_align(28);
    //read_align(27);
    // for(int i=0;i<rcnt;i++){
    //     auto& it=vecr[i];
    //     cout<<it.str<< ' '<<vec_pos[rndk-1][i]<<'\n';
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
    vec_k.pb(k);
    vec_val.pb(vector<ull>(rcnt,0));
    vec_is.pb(vector<int>(rcnt,0));
    vec_pos.pb(vector<int>(rcnt,0));
    for(int i=0;i<rcnt;i++){
        cal_k_mer(k,i,vec_pos[rndk][i],vec_val[rndk][i],vec_is[rndk][i]);
    }
    vector<int> order_to_id;//cur vec
    for(int i=0;i<rcnt;i++) order_to_id.pb(i);
    sort(order_to_id.begin(),order_to_id.end(),cmp1);
    #if testflag
    cout<<"k="<<rndk<<'\n';
    cout<<"order:\n";
    for(auto it:order_to_id){
        cout<<it<<' ';
    }
    cout<<'\n';
    #endif
    // cout<<"detail:\n";
    // for(int i=0;i<rcnt;i++){
    //     cout<<i<<' '<<vec_val[rndk][i]<<' '<<vec_is[rndk][i]<<' '<<vec_pos[rndk][i]<<'\n';
    // }
    // cout<<'\n';
    vector<int> id_to_order(rcnt,0);
    for(int i=0;i<rcnt;i++) id_to_order[order_to_id[i]]=i;
    //write down the result of this round
    vec_id.pb(order_to_id);
    vec_order.pb(id_to_order);
    rndk++;
}

/*

according to the nxtid recored in read,choose the best one to continue.


*/
void contig_make()
{
    cout<<"contig_make start\n";
    //vecb:bias between cur read and the next read
    //bias need to be small
    //isrev:record is reverse state between the adjacent reads 
    //veck:record the k  
    vector<int> pre(rcnt,-1),nxt(rcnt,-1),vecb(rcnt,max_str_length),isrev(rcnt,0),veck(rcnt,0);
    for(int i=0;i<rcnt;i++){//every read
        //for every read_align, check the next 10 reads to form nxtid(id,bias)
        vector<pii> nxtid;
        for(int j=0;j<rndk;j++){//for every k
            int mypos=vec_order[j][i];
            for(int p=mypos+1;p<rcnt&&p<mypos+k_num;p++){//check the read in pos=p
                int id2=vec_id[j][p];
                int id1=i;
                if(vec_val[j][id1]!=vec_val[j][id2]) break;//the k_minimizer is not same
                if(check(id1,id2,j)) nxtid.pb(pii(id2,j));
                //cout<<"+++"<<' '<<id1<<' '<<id2<<' '<<bias<<' '<<check(id1,id2,bias)<<'\n';
            }
        }
        // cout<<"nxtid:|||"<<i<<'\n';
        // for(auto it:nxtid) cout<<it.first<<' '<<it.second<<'\n';
        //according to nxtid, find the best one to concat
        map<pii,sconmake> mp;
        for(pii it:nxtid){//id and bias
            //cout<<"|||"<<i<<' '<<it.first<<' '<<it.second<<'\n';
            if(pre[it.first]!=-1) continue;
            int bias=cal_len(i,it.second)-cal_len(it.first,it.second);
            assert(bias>=0);
            if(mp.count(pii(it.first,bias))){
                mp[pii(it.first,bias)].cnt++;
            }
            else{
                mp[pii(it.first,bias)]=sconmake{vec_is[it.second][i]^vec_is[it.second][it.first],1,it.second};
            }   
        }
        int fre=0,bias=max_str_length;
        sconmake tmp=sconmake{0,-1,0};//regard tmp.cnt as nid
        //???????consider more bias!instead of just frequent only, or the read with one diff will be discard
        //first, check the diff is less than threshold
        //
        //choose the bias min one.
        //if many, choose the most frequently occurent one
        //if many, choose any one
        for(auto& it:mp){
            if(it.first.second<bias){
                fre=it.second.cnt;
                bias=it.first.second;
                tmp=sconmake{it.second.isrev,it.first.first,it.second.k};
            }
            else if(it.first.second==bias&&it.second.cnt>fre){
                fre=it.second.cnt;
                tmp=sconmake{it.second.isrev,it.first.first,it.second.k};
            }
        }
        if(tmp.cnt!=-1){
            nxt[i]=tmp.cnt;
            pre[tmp.cnt]=i;
            vecb[i]=bias;
            isrev[i]=tmp.isrev;
            veck[i]=tmp.k;
        }
    }
    // cout<<"====bias\n";
    // for(int i=0;i<rcnt;i++){
    //     if(pre[i]!=-1||nxt[i]!=-1){
    //         cout<<i<<' '<<pre[i]<<' '<<nxt[i]<<' '<<vecb[i]<<' '<<isrev[i]<<'\n';
    //     }
        
    // }

    vector<int> seq;
    vector<int> vis(rcnt,0);
    for(int t=0;t<rcnt;t++){
        //the reads form a chain
        //impossible to form a ring (only in small cases and small k)
        if(vis[t]==1) continue;
        seq.clear();
        int cur=t;
        while(pre[cur]!=-1&&pre[cur]!=t){
            cur=pre[cur];
        }
        int start=cur;
        while(cur!=-1){
            seq.pb(cur);
            vis[cur]=1;
            cur=nxt[cur];
            if(cur==start) break;
        }

        //actually make contig according to the seq
        // cout<<"seq:  ";
        // for(int i=0;i<seq.size();i++){
        //     cout<<seq[i]<<' ';
        // }
        // cout<<'\n';
        vecc.push_back(Contig());
        Contig &con=vecc.back();
        con.cid=ccnt++;
        //initialize currev, use veck to determine the state
        int LEN=0,curb=0,currev=vec_is[veck[seq[0]]][seq[0]];
        for(int i=0;i<seq.size();i++){
            LEN=max(LEN,curb+(int)vecr[seq[i]].str.length());
            curb=curb+vecb[seq[i]];
        }
        for(int j=0;j<LEN;j++){//for every pos
            int cnt[4]={0};
            curb=0;
            currev=vec_is[veck[seq[0]]][seq[0]];
            for(int i=0;i<seq.size();i++){//for every read
                if(curb<=j&&(int)vecr[seq[i]].str.length()-1+curb>=j){
                    int curc;
                    if(currev==0) curc=trans(symm(vecr[seq[i]].str[j-curb],0));
                    else curc=trans(symm(vecr[seq[i]].str[vecr[seq[i]].str.length()-1-(j-curb)],1));
                    //cout<<"debug: "<<j<<' '<<i<<' '<<curb<<' '<<currev<<' '<<curc<<'\n';
                    cnt[curc]++;
                }
                curb+=vecb[seq[i]];
                currev=currev^isrev[seq[i]];
            }
            int mx=cnt[3],mxpos=3;
            for(int i=2;i>=0;i--) if(cnt[i]>=mx) mx=cnt[i],mxpos=i;
            con.str+=trans(mxpos);
        }
        curb=0;
        currev=vec_is[veck[seq[0]]][seq[0]];
        for(int i=0;i<seq.size();i++){//for every read
            //calculate the hamming distance and setting cid
            //int cnt=0;
            
            // for(int j=0;j<vecr[i].str.length();j++){
            //     if(vecr[i].str[j]!=con.str[j+curb]) cnt++;
            // }
            //if(cnt<=threshold){//join in and record the dismatch
                vecr[seq[i]].cid=con.cid;
                vecr[seq[i]].cpos=curb;
                vecr[seq[i]].isrev=currev;
                for(int j=0;j<vecr[seq[i]].str.length();j++){
                    if(currev==0){
                        if(vecr[seq[i]].str[j]!=con.str[j+curb]){
                            vecr[seq[i]].dismatch.pb(pii(j,vecr[seq[i]].str[j]));
                        }
                    }
                    else{//reverse
                        if(vecr[seq[i]].str[j]!=symm(con.str[vecr[seq[i]].str.length()-1-j+curb],1)){//contig has no 'N'
                            vecr[seq[i]].dismatch.pb(pii(j,vecr[seq[i]].str[j]));
                        }
                    }
                    
                }
            //}
            //else hamming distance too large ,don't join in 
            curb+=vecb[seq[i]];
            currev^=isrev[i];
        }
        #if testflag
        cout<<ccnt-1<<' '<<con.str<<' '<<LEN<<'\n';
        #endif
        //if(ccnt<=100) cout<<ccnt<<' '<<seq.size()<<'\n';
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
    out_int(fout,ans.length());
    cout<<rcnt<<' '<<ccnt<<'\n';
    //cout<<ans<<'\n';
    for(int i=0;i<ans.length();i+=16){
        uint tmp=0;
        uint base=0x40000000;
        for(int j=0;j<16;j++){//don't has 'N'
            if(j+i<ans.length()) tmp=tmp+trans(ans[i+j])*base;
            base>>=2;
        }
        out_int(fout,tmp);
    }
    out_int(fout,rcnt);
    vector<int> curpos(rcnt+repeatcnt,-1);
    sort(vecr.begin(),vecr.begin()+rcnt,cmp4);
    sort(vecr.begin()+rcnt,vecr.end(),cmp4);//sort separately
    //for every read,show the pos
    for(int i=0;i<rcnt;i++){
        Read& it=vecr[i];
        curpos[vecr[i].rid]=i;
        out_int(fout,(int)(vecc[it.cid].spos+it.cpos));//pos;
        out_char(fout,(char)vecr[i].str.length());//len
        //here can use bias...
        out_char(fout,(char)(vecr[i].isrev));
        //cout<<it.str<<' '<<vecc[it.cid].spos+it.cpos<<'\n';
        out_char(fout,(char)it.dismatch.size());//at most 5
        if(it.dismatch.size()!=0){
            for(auto it1:it.dismatch){
                out_char(fout,(char)it1.first);
                out_char(fout,(char)it1.second);
            }
        }
        //cout<<"+++"<<i<<' '<<(int)(vecc[it.cid].spos+it.cpos)<<' '<<(int)vecr[i].str.length()<<' '<<(int)(vecr[i].isrev)<<' '<<(int)it.dismatch.size()<<'\n';
    }
    out_int(fout,repeatcnt);
    for(int i=rcnt;i<rcnt+repeatcnt;i++){
        if(curpos[vecr[i].repeatid]==-1){
            cout<<i<<' '<<rcnt<<' '<<repeatcnt<<' '<<vecr.size()<<'\n';
            cout<<vecr[i].repeatid<<' '<<vecr[i].isrepeat<<'\n';
        }
        assert(curpos[vecr[i].repeatid]!=-1);
        //cout<<vecr[i].repeatid<<'\n';
        out_int(fout,curpos[vecr[i].repeatid]);
        out_char(fout,(char)vecr[i].issymmrepeat);
        if(order_preserve==1){
            out_int(fout,vecr[i].rid);//the rid of read,if need order preserved
        }
    }
}