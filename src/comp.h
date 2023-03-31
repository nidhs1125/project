#include<bits/stdc++.h>
#include<unistd.h>
#include"tools.h"

//using namespace std;


void compmain(ifstream& fin,string& out_path)
{
    init();
    print_time();
    string tmp;
    int cnt=0;
    vecr.clear();
    while(getline(fin,tmp)){
        if(tmp[0]=='$') break;
        cnt++;
        if(cnt%4==2){
            vecr.pb(Read(rcnt++,tmp));
        }
        // if(cnt%1000000==0){
        //     cout<<"cnt= "<<cnt<<'\n';
        // }
        //if(cnt>=30000000) break;
    }
    cout<<cnt<<"lines in total\n";
    read_len=vecr[0].str.length();
    print_time();
    //exclude the same reads
    prework();
    cout<<"prework finish\n";
    print_time();
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
    pre.resize(rcnt);
    basket.resize(rcnt);
    basket_rev.resize(rcnt);
    for(int i=0;i<rcnt;i++) pre[i]=i,basket[i].pb(i),basket_rev[i].pb(0);


    //part 1, merging the reads without bias
    first_rev.assign(rcnt,0);
    #if testflag
    read_align(4);
    #else
    read_align(25);
    read_align(27);
    read_align(29);
    // read_align(31);
    // read_align(23);
    
    #endif

    cout<<"read_align finish"<<'\n';
    print_time();

    #if testflag
    for(int i=0;i<rcnt;i++){
        if(find(i)==i){
            cout<<"basket: "<<i<<'\n';
            for(int j=0;j<basket[i].size();j++){
                cout<<vecr[basket[i][j]].str<<' '<<vecr[basket[i][j]].k_mer_pos<<' '<<vecr[basket[i][j]].iskmersymm<<'\n';
            }
        }
    }
    #endif
    //test basket:
    {
        
        // for(int i=0;i<rcnt;i++){
        //     if(find(i)==i){
        //         if(basket[i].size()<=10&&basket[i].size()>=5){
        //             cout<<"small:"<<basket[i].size()<<'\n';
        //             for(int j=0;j<basket[i].size();j++){
        //                 int id=basket[i][j];
        //                 if(vecr[id].iskmersymm){
        //                     cout<<cal_symm(vecr[id].str)<<'\n';
        //                 }
        //                 else cout<<vecr[id].str<<'\n';
        //             }
        //             break;
        //         }
        //     }
        // }

        // for(int i=0;i<rcnt;i++){
        //     if(find(i)==i){
        //         if(basket[i].size()>=30&&basket[i].size()<=50){
        //             cout<<"mid:"<<basket[i].size()<<'\n';
        //             for(int j=0;j<basket[i].size();j++){
        //                 int id=basket[i][j];
        //                 if(vecr[id].iskmersymm){
        //                     cout<<cal_symm(vecr[id].str)<<'\n';
        //                 }
        //                 else cout<<vecr[id].str<<'\n';
        //             }
        //             break;
        //         }
        //     }
        // }

        // for(int i=0;i<rcnt;i++){
        //     if(find(i)==i){
        //         if(basket[i].size()>=100){
        //             cout<<"large:"<<basket[i].size()<<'\n';
        //             for(int j=0;j<basket[i].size();j++){
        //                 int id=basket[i][j];
        //                 if(vecr[id].iskmersymm){
        //                     cout<<cal_symm(vecr[id].str)<<'\n';
        //                 }
        //                 else cout<<vecr[id].str<<'\n';
        //             }
        //             break;
        //         }
        //     }
        // }
        // int mx=0;
        // map<int,int> bas_size;
        // for(int i=0;i<rcnt;i++){
        //     if(find(i)==i){
        //         mx=max(mx,(int)basket[i].size());
        //         bas_size[basket[i].size()]++;
        //     }
        // }
        // cout<<"max basket size= "<<mx<<'\n';
        // cout<<"basket_size and number of baskets:\n";
        // for(auto it:bas_size){
        //     cout<<it.first<<' '<<it.second<<'\n';
        // }


    }

    //part 2,merging baskets with bias
    bas_num=0;
    order_to_id.resize(rcnt);
    for(int i=0;i<rcnt;i++){
        if(find(i)==i){
            order_to_id[bas_num++]=i;
        }
    }
    bas_edge.resize(rcnt);
    first_rev.assign(rcnt,0);
    for(int i=0;i<rcnt;i++){
        bas_edge[i].resize(1);
        bas_edge[i][0]=Edge(i);
    }
    //firtly,find the best nxt read, for all k
    //constraint: distance<=threshold and has same k-minimizer
    #if testflag
    bas_align(4);
    #else
    bas_align(25);
    bas_align(27);
    bas_align(29);
    #endif
    //set edge and first_rev
    for(int i=0;i<rcnt;i++){
        if(bas_edge[i][0].bias<=max_bias){
            Edge& tmp=bas_edge[i][0];
            bas_edge[tmp.to].pb(Edge(i,-tmp.bias,tmp.rev));
            first_rev[tmp.to]=first_rev[i]^tmp.rev;
        }
    }
    cout<<"bas_align finish\n";
    print_time();
    //then concatenate the read and make contig
    vis.assign(rcnt,0);
    contig_make();
    cout<<"contig_make finish\n";
    print_time();
    SCS_gen();
}


void read_align(int k)
{
    cout<<"read_align start, k="<<k<<'\n';
    order_to_id.resize(rcnt);
    //this part can add multi-thread
    for(int i=0;i<rcnt;i++){
        order_to_id[i]=i;
        if(find(i)==i){
            cal_k_mer(k,i,vecr[i].k_mer_pos,vecr[i].val,vecr[i].iskmersymm);
            #if testflag
            cout<<i<<' '<<vecr[i].str<<' '<<vecr[i].k_mer_pos<<' '<<vecr[i].iskmersymm<<' '<<vecr[i].val<<'\n';
            #endif
        }
        else vecr[i].val=UINT64_MAX,vecr[i].k_mer_pos=0;//not root,discard
    }
    sort(order_to_id.begin(),order_to_id.end(),cmp1);
    pthread_id=0;
    int tmp;
    vector<pthread_t> thread_id(thread_num);
    vector<int> vectmp(thread_num+1);
    for(int i=0;i<thread_num;i++){
        //cout<<i<<'\n';
        vectmp[i]=i;
        if((tmp=pthread_create(&thread_id[i],NULL,cal_pre,(void*)&vectmp[i]))!=0){
            cout<<"pthread_create ERROR\n";
            exit(0);
        }
    }
    vectmp[thread_num]=thread_num;
    cal_pre((void*)&vectmp[thread_num]);
    //cout<<"read_align :thread_create finish\n";
    for(int i=0;i<thread_num;i++){
        //cout<<i<<'\n';
        void* tmp;
        pthread_join(thread_id[i],&tmp);
    }
    
}

void* cal_pre(void* arg)
{
    int cnt=0;
    int thread_id=*(int*)arg;
    while(1){
        int my_id;
        pthread_mutex_lock(&pthread_id_mutex);
        my_id=pthread_id++;
        pthread_mutex_unlock(&pthread_id_mutex);
        cnt++;
        int pos=order_to_id[my_id];
        //cout<<"cal_pre: "<<my_id<<' '<<pos<<'\n';
        if(my_id>=rcnt||vecr[pos].val==UINT64_MAX) break;
        for(int i=my_id-1;i>=0&&i>=my_id-300;i--){
            int ppos=order_to_id[i];
            if(vecr[pos].val!=vecr[ppos].val) break;
            if(cal_len(pos)!=cal_len(ppos)) break;
            //printf("+++ %d %d %d\n",pos,ppos,check(pos,ppos));
            if(check(pos,ppos)){
                //merge two baskets,need lock, the total complexity is O(nlog^2n)
                pthread_mutex_lock(&pre_mutex);
                int id1=find(pos),id2=find(ppos);
                pre[id1]=pre[id2]=dsu(id1,id2);
                first_rev[pos]=vecr[pos].iskmersymm;//if multi,then the first_rev is cur rev state
                pthread_mutex_unlock(&pre_mutex);
                
                break;
            }
        }
    }
    //printf("thread_id=%d,cnt=%d\n",thread_id,cnt);
    return (void*)0;
}

//according to the basket,calculate the relation between root of basket
//or union basket
void bas_align(int k)
{
    cout<<"bas_align start,k="<<k<<'\n';
    order_to_id.resize(bas_num);
    //this part can add multi-thread
    for(int i=0;i<bas_num;i++){
        int id=order_to_id[i];
        assert(find(id)==id);
        cal_k_mer(k,id,vecr[id].k_mer_pos,vecr[id].val,vecr[id].iskmersymm);
    }
    sort(order_to_id.begin(),order_to_id.begin()+bas_num,cmp1);
    #if testflag
    cout<<"order_to_id:\n";
    for(int i=0;i<bas_num;i++) cout<<order_to_id[i]<<' ';
    cout<<'\n';
    #endif
    pthread_id=0;
    int tmp;
    vector<pthread_t> thread_id(thread_num);
    vector<int> vectmp(thread_num+1);
    for(int i=0;i<thread_num;i++){
        vectmp[i]=i;
        if((tmp=pthread_create(&thread_id[i],NULL,cal_nxt,&vectmp[i]))!=0){
            cout<<"pthread_create ERROR\n";
            exit(0);
        }
    }
    vectmp[thread_num]=thread_num;
    cal_nxt(&vectmp[thread_num]);
    for(int i=0;i<thread_num;i++){
        void* tmp;
        pthread_join(thread_id[i],&tmp);
    }
    
}

void* cal_nxt(void* arg)
{
    int id=*(int*) arg;
    while(1){
        int my_id;
        pthread_mutex_lock(&pthread_id_mutex);
        my_id=pthread_id++;
        pthread_mutex_unlock(&pthread_id_mutex);
        if(my_id>=bas_num) break;
        int pos=order_to_id[my_id];
        //printf("cal_nxt: %d %d\n",my_id,pos);
        assert(vecr[pos].val!=UINT64_MAX);
        for(int i=my_id+1;i<bas_num&&i<=my_id+300;i++){
            int npos=order_to_id[i];
            if(vecr[pos].val!=vecr[npos].val) break;
            if(cal_len(pos)==cal_len(npos)) continue;//don't consider the baskets without bias
            if(cal_len(pos)-cal_len(npos)>max_bias) break;
            //if(pre_bas[npos]!=-1) continue;
            if(check(pos,npos)){
                //merge two baskets,need lock, the total complexity is O(nlog^2n)
                pthread_mutex_lock(&bas_edge_mutex);
                Edge tmpe=Edge(npos,cal_len(pos)-cal_len(npos),vecr[pos].iskmersymm^vecr[npos].iskmersymm);
                int tmp=bas_edge[pos][0].update(tmpe);
                if(tmp) first_rev[pos]=vecr[pos].iskmersymm;
                //else ,bias too large ,discard
                pthread_mutex_unlock(&bas_edge_mutex);//can't merge,but this can occur only few cases
                break;
                
            }
        }
    }
    return (void*)0;
}

void contig_make()
{
    cout<<"contig_make start\n";
    vis_id=0;
    vis.assign(rcnt,0);
    int tmp;
    vector<pthread_t> thread_id(thread_num);
    vector<int> vectmp(thread_num+1);
    for(int i=0;i<thread_num;i++){
        vectmp[i]=i;
        if((tmp=pthread_create(&thread_id[i],NULL,cal_block,&vectmp[i]))!=0){
            cout<<"pthread_create ERROR\n";
            exit(0);
        }
    }
    vectmp[thread_num]=thread_num;
    cal_block(&vectmp[thread_num]);
    for(int i=0;i<thread_num;i++){
        void* tmp;
        pthread_join(thread_id[i],&tmp);
    }

    
}

void dfs(vector<Edge>& vec,int rt,int cbias,int crev)
{
    if(vis[rt]==1) return;
    vec.pb(Edge(rt,cbias,crev));
    vis[rt]=1;
    for(Edge& it:bas_edge[rt]){
        if(it.bias<=max_bias) dfs(vec,it.to,cbias+it.bias,crev^it.rev);
    }
}

void* cal_block(void* arg)
{
    int id=*(int*)arg;
    while(1){
        int my_id;
        int contig_id;
        vector<Edge> seq;
        pthread_mutex_lock(&vis_id_mutex);
        while(vis_id<rcnt){
            if(find(vis_id)==vis_id&&vis[vis_id]==0) break;
            vis_id++;
        }
        if(vis_id>=rcnt){
            pthread_mutex_unlock(&vis_id_mutex);
            break;
        }
        my_id=vis_id++;
        dfs(seq,my_id,0,first_rev[my_id]);//the whole block should be taged
        contig_id=ccnt;
        vecc.pb(Contig(ccnt));
        ccnt++;
        pthread_mutex_unlock(&vis_id_mutex);


        //for this block, make a contig
        int L=1e9,R=-1e9;
        Contig& con=vecc[contig_id];
        for(Edge it:seq){
            L=min(L,it.bias);
            R=max(R,it.bias+read_len);
        }
        printf("cal_block: %d %d %d\n",id,L,R);
        for(int j=L;j<R;j++){//for every pos
            int cnt[4]={0,0,0,0};
            for(int i=0;i<seq.size();i++){//for every basket
                for(int m=0;m<basket[seq[i].to].size();m++){//for every read 
                    int id=basket[seq[i].to][m];
                    int tmprev=basket_rev[seq[i].to][m];
                    if(seq[i].bias<=j&&(int)vecr[id].str.length()-1+seq[i].bias>=j){
                        int curc;
                        if(seq[i].rev^tmprev==0) curc=trans(symm(vecr[id].str[j-seq[i].bias],0));
                        else curc=trans(symm(vecr[id].str[vecr[id].str.length()-1-(j-seq[i].bias)],1));
                        //cout<<"debug: "<<j<<' '<<i<<' '<<curb<<' '<<currev<<' '<<curc<<'\n';
                        cnt[curc]++;
                    }
                }
            }
            //cout<<"++++"<<j<<' '<<cnt[0]<<' '<<cnt[1]<<' '<<cnt[2]<<' '<<cnt[3]<<'\n';
            int mx=cnt[3],mxpos=3;
            for(int i=2;i>=0;i--) if(cnt[i]>=mx) mx=cnt[i],mxpos=i;
            con.str+=trans(mxpos);
        }

        for(int i=0;i<seq.size();i++){//for every basket
            for(int m=0;m<basket[seq[i].to].size();m++){//for every read
                int id=basket[seq[i].to][m];
                int tmprev=basket_rev[seq[i].to][m];
                vecr[id].cpos=seq[i].bias-L;
                vecr[id].isrev=tmprev^seq[i].rev;
                vecr[id].cid=con.cid;
                for(int j=0;j<vecr[id].str.length();j++){
                    if(seq[i].rev^tmprev==0){
                        if(vecr[id].str[j]!=con.str[j+seq[i].bias-L]){
                            vecr[id].dismatch.pb(pii(j,vecr[id].str[j]));
                        }
                    }
                    else{//reverse
                        if(vecr[id].str[j]!=symm(con.str[vecr[id].str.length()-1-j+seq[i].bias-L],1)){//contig has no 'N'
                            vecr[id].dismatch.pb(pii(j,vecr[id].str[j]));
                        }
                    }
                }
                con.num++;
            }
        }
    }
    return (void*)0;
}


void SCS_gen()
{
    cout<<"SCS_gen start\n";
    sort(vecc.begin(),vecc.end(),cmp5);
    int cbias=0;
    for(int i=0;i<ccnt;i++){
        if(vecc[i].num!=1){//if more than one, add 
            ans+=vecc[i].str;
            vecc[i].spos=cbias;
            cbias+=vecc[i].str.length();
        }
        else{//try to realign
            realignment(i);
            if(vecc[i].is==0){//add to the end of ans(SCS)
                ans+=vecc[i].str;
                vecc[i].spos=cbias;
                cbias+=vecc[i].str.length();
            }
        }
        
    }
    sort(vecc.begin(),vecc.end(),cmp6);
}


void realignment(int id)//cur vecc id
{

}


//old version
void encode(string& out_path)
{
    cout<<"encoding...\n";

    #if testflag
    cout<<ans<<'\n';
    for(int i=0;i<rcnt;i++){
        cout<<i<<' '<<vecr[i].str<<' '<<vecc[vecr[i].cid].spos<<' '<<vecr[i].cpos<<' '<<vecr[i].isrev<<' '<<vecr[i].dismatch.size()<<'\n';
        for(auto it:vecr[i].dismatch){
            cout<<it.first<<' '<<it.second<<'\n';
        }
    }
    #endif
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
    out_int(fout1,read_len);
    cout<<rcnt<<' '<<bas_num<<' '<<ccnt<<'\n';
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
            if(threshold==0&&nnum!=it.dismatch.size()){
                cout<<"ERROR!\n";
                cout<<i<<' '<<it.str<<' '<<ans.substr(cpos,read_len)<<' '<<it.isrev<<'\n';
                cout<<it.dismatch.size()<<'\n';
                for(auto it1:it.dismatch){
                    cout<<it1.first<<' '<<it1.second<<'\n';
                }
                assert(nnum==it.dismatch.size());
            }
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

    fout1.close();
    fout2.close();
    fout3.close();
    system("./bsc e test.mine1 test.bsc1");
    system("./bsc e test.mine2 test.bsc2");
    system("./bsc e test.mine3 test.bsc3");
}


/*
//new version
void encode(string& out_path)
{
    cout<<"encoding...\n";

    #if testflag
    cout<<ans<<'\n';
    for(int i=0;i<rcnt;i++){
        cout<<i<<' '<<vecr[i].str<<' '<<vecc[vecr[i].cid].spos<<' '<<vecr[i].cpos<<' '<<vecr[i].isrev<<' '<<vecr[i].dismatch.size()<<'\n';
        for(auto it:vecr[i].dismatch){
            cout<<it.first<<' '<<it.second<<'\n';
        }
    }
    #endif
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
    out_int(fout1,rcnt);
    out_int(fout1,repeatcnt);
    out_int(fout1,read_len);
    cout<<rcnt<<' '<<bas_num<<' '<<ccnt<<'\n';
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
    //for(int i=0;i<ans.length();i++) out_char(fout1,ans[i]);
    //part 2
    
    vector<int> curpos(rcnt+repeatcnt,-1);
    for(int i=0;i<rcnt;i++){
        curpos[vecr[i].rid]=i;
        vecr[i].cpos+=vecc[vecr[i].cid].spos;//pos in contig +bias of contig -> pos in ans
    }
    for(int i=rcnt;i<rcnt+repeatcnt;i++){
        vecr[i].repeatid=curpos[vecr[i].repeatid];//change repeadid to record the id
    }
    int ppos=0,cpos;

    int dismatchcnt=0;
    for(int i=0;i<rcnt;i++){
        Read& it=vecr[i];
        cpos=vecr[i].cpos;
        out_char(fout2,cpos-ppos);//bias
        out_char(fout2,(char)(it.isrev));//if rev
        out_char(fout3,(char)it.dismatch.size());//at most 5
        //assert(it.dismatch.size()<=threshold);
        dismatchcnt+=it.dismatch.size()*2+1;
        if(it.dismatch.size()!=0){
            for(auto it1:it.dismatch){
                out_char(fout3,(char)it1.first);
                out_char(fout3,(char)it1.second);
            }
        }
        ppos=cpos;
    }
    for(int i=rcnt;i<rcnt+repeatcnt;i++){
        out_int(fout2,vecr[i].repeatid);
    }
    cout<<dismatchcnt<<' '<<repeatcnt<<'\n';


    fout1.close();
    fout2.close();
    fout3.close();
    system("./bsc e test.mine1 test.bsc1");
    system("./bsc e test.mine2 test.bsc2");
    system("./bsc e test.mine3 test.bsc3");
}
*/

/*





*/