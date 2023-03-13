#include<bits/stdc++.h>
#include<unistd.h>
#include"tools.h"
using namespace std;




void compmain(ifstream& fin,string& out_path)
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
    read_len=vecr[0].str.length();
    //exclude the same reads
    prework();
    cout<<"prework finish\n";
    //接下来需要参考minicom生成contig
    comp();
    cout<<cnt<<"lines,"<<" compression finish"<<'\n';
    encode(out_path);

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
    #if testflag
    // int k=4;
    // int tmp;
    // pthread_t pid_i=thread_bias-1;
    // if((tmp=pthread_create(&pid_i,NULL,read_align,(void*)k))!=0){
    //     cout<<"pthread_create ERROR\n";
    //     exit(0);
    // }
    read_align(4);
    #else
    read_align(25);
    read_align(27);
    read_align(29);
    #endif
    // for(int i=0;i<rcnt;i++){
    //     auto& it=vecr[i];
    //     cout<<it.str<< ' '<<vec_pos[rndk-1][i]<<'\n';
    // }

    cal_nxt();
    cout<<"read_align finish\n";
    contig_make();


}


void* read_align(void* arg)
{
    
    cout<<"read_align start\n";
    int k=(int)arg;
    for(int i=0;i<rcnt;i++){
        cal_k_mer(k,i,vecr[i].k_mer_pos,vecr[i].val,vecr[i].isrev);
    }
    pre.resize(rcnt);
    for(int i=0;i<rcnt;i++) pre[i]=i;
    int tmp;
    for(int i=0;i<rcnt;i++){
        pthread_t pid_i=i+thread_bias;
        if((tmp=pthread_create(&pid_i,NULL,cal_pre,(void*)i))!=0){
            cout<<"pthread_create ERROR\n";
            exit(0);
        }
    }
    for(int i=0;i<rcnt;i++){
        pthread_t pid_i=i+thread_bias;
        pthread_join(pid_i,NULL);
    }
    for(int i=0;i<rcnt;i++) find(i);//Path compression
    basket.resize(rcnt);
    for(int i=0;i<rcnt;i++){
        basket[pre[i]].pb(i);
    }
}

void* cal_pre(void* arg)
{
    int pos=(int)arg;
    for(int i=pos-1;i>=0&&i>=pos-300&&vecr[i].val==vecr[pos].val&&cal_len(i)==cal_len(pos);i--){
        if(check(i,pos)){
            pre[pos]=i;
            break;
        }
    }
}

//according to the basket,calculate the relation between root of basket
//or union basket
void cal_nxt()
{
    nxt_bas.assign(rcnt,-1);
    //...
}

void contig_make()
{
    cout<<"contig_make start\n";
    //vecb:bias between cur read and the next read
    //bias need to be small
    //isrev:record is reverse state between the adjacent reads 
    //veck:record the k  
    vector<int> pre(rcnt,-1),nxt(rcnt,-1),vecb(rcnt,max_str_length),isrev(rcnt,0),veck(rcnt,0);
    vector<int> seq;
    vector<int> vis(rcnt,0);
    for(int t=0;t<rcnt;t++){
        //the reads form a chain
        //impossible to form a ring (only in small cases and small k)
        if(vis[t]==1) continue;
        seq.clear();
        int cur=t;
        while(cur!=-1){
            seq.pb(cur);
            vis[cur]=1;
            cur=nxt_bas[cur];
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
        int LEN=0,curb=0;
        for(int i=0;i<seq.size();i++){
            LEN=max(LEN,curb+(int)vecr[seq[i]].str.length());
            curb=curb+bias_bas[seq[i]];
        }
        for(int j=0;j<LEN;j++){//for every pos
            int cnt[4]={0};
            curb=0;
            for(int i=0;i<seq.size();i++){//for every basket
                for(int m=0;m<basket[i].size();m++){
                    if(curb<=j&&(int)vecr[seq[i]].str.length()-1+curb>=j){
                        int id=basket[seq[i]][m];
                        int curc;
                        if(vecr[id].isrev==0) curc=trans(symm(vecr[id].str[j-curb],0));
                        else curc=trans(symm(vecr[id].str[vecr[id].str.length()-1-(j-curb)],1));
                        //cout<<"debug: "<<j<<' '<<i<<' '<<curb<<' '<<currev<<' '<<curc<<'\n';
                        cnt[curc]++;
                    }
                }
                curb+=vecb[seq[i]];
            }
            int mx=cnt[3],mxpos=3;
            for(int i=2;i>=0;i--) if(cnt[i]>=mx) mx=cnt[i],mxpos=i;
            con.str+=trans(mxpos);
        }
        curb=0;
        for(int i=0;i<seq.size();i++){//for every basket
            for(int m=0;m<basket[i].size();m++){
                int id=basket[seq[i]][m];
                for(int j=0;j<vecr[id].str.length();j++){
                    if(vecr[id].isrev==0){
                        if(vecr[id].str[j]!=con.str[j+curb]){
                            vecr[id].dismatch.pb(pii(j,vecr[id].str[j]));
                        }
                    }
                    else{//reverse
                        if(vecr[id].str[j]!=symm(con.str[vecr[id].str.length()-1-j+curb],1)){//contig has no 'N'
                            vecr[id].dismatch.pb(pii(j,vecr[id].str[j]));
                            // if(threshold==0){//other character must be same when thereshold=0
                            //     if(vecr[id].str[j]!='N'){
                            //         cout<<"ERROR\n";
                            //         cout<<"seq:  ";
                            //         for(int i=0;i<seq.size();i++){
                            //             cout<<id<<' ';
                            //         }
                            //         for(int i=0;i<seq.size();i++){
                            //             cout<<vecr[id].str<<'\n';
                            //         }
                            //         cout<<'\n';
                            //         cout<<i<<' '<<j<<' '<<id<<' '<<'\n';
                            //         assert(vecr[id].str[j]=='N');
                            //     }
                            // }
                        }
                    }
                }
                
            }
            curb+=vecb[seq[i]];
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

void encode(string& out_path)
{
    cout<<"encoding...\n";
    // for(auto it:vecr) cout<<it.rid<<' ';
    // cout<<'\n';
    //according to out_path, determine the out_path
    //has 3 paths in total: contig, readpos(bias) and the mismatch
    ofstream fout1,fout2,fout3;
    fout1.open(out_path+".mine1");
    fout2.open(out_path+".mine2");
    fout3.open(out_path+".mine3");

    //part 1
    out_int(fout1,ans.length());
    out_int(fout1,rcnt+repeatcnt);
    cout<<rcnt<<' '<<ccnt<<'\n';
    //cout<<ans<<'\n';
    for(int i=0;i<ans.length();i+=16){
        uint tmp=0;
        uint base=0x40000000;
        for(int j=0;j<16;j++){//don't has 'N'
            if(j+i<ans.length()) tmp=tmp+trans(ans[i+j])*base;
            base>>=2;
        }
        out_int(fout1,tmp);
    }

    //part 2
    
    vector<int> curpos(rcnt+repeatcnt,-1);
    for(int i=0;i<rcnt;i++){
        curpos[vecr[i].rid]=i;
        vecr[i].cpos+=vecc[vecr[i].cid].spos;//pos in contig +bias of contig -> pos in ans
    }
    for(int i=rcnt;i<rcnt+repeatcnt;i++){
        vecr[i].cpos=vecr[curpos[vecr[i].repeatid]].cpos;
    }
    sort(vecr.begin(),vecr.end(),cmp4);
    int ppos=0,cpos;

    int dismatchcnt=0;
    for(int i=0;i<rcnt+repeatcnt;i++){
        Read& it=vecr[i];
        cpos=vecr[i].cpos;
        out_char(fout2,cpos-ppos);//bias
        out_char(fout2,it.isrepeat==1);//if repeat
        if(it.isrepeat==1){
            out_char(fout2,it.issymmrepeat);//if symm
        }
        else{
            out_char(fout2,(char)(it.isrev));//if rev
            out_char(fout3,(char)it.dismatch.size());//at most 5
            //assert(it.dismatch.size()<=threshold);
            int nnum=0;
            for(int j=0;j<it.str.length();j++) if(it.str[j]=='N') nnum++;
            if(threshold==0) assert(nnum==it.dismatch.size());
            dismatchcnt+=it.dismatch.size()*2+1;
            if(it.dismatch.size()!=0){
                for(auto it1:it.dismatch){
                    out_char(fout3,(char)it1.first);
                    out_char(fout3,(char)it1.second);
                }
            }
        }
        ppos=cpos;
    }
    cout<<dismatchcnt<<' '<<repeatcnt<<'\n';
}

/*





*/