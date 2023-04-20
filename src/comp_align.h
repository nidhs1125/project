#include<bits/stdc++.h>
#include<unistd.h>
#include"tools.h"

void read_align(int k);
void* cal_pre(void* arg);
int dsu(int id1,int id2);
void bas_align(int k);
void* cal_nxt(void* arg);
int bas_dsu(int id1,int id2);
//read_align and absket align
void align()
{
    //read_align
    pre.resize(rcnt0);
    basket.resize(rcnt0);
    basket_rev.resize(rcnt0);
    //first_rev.resize(rcnt0);
    vecr.clear();
    for(int i=0;i<rcnt0;i++){
        if(isrepeat[i]==0&&hasn[i]==0){
            pre[i]=i,basket[i].pb(i),basket_rev[i].pb(0);
            vecr.pb(i);
        }
        else pre[i]=-1e9;
    }
    
    assert(vecr.size()==rcnt&&rcnt==rcnt0-repeatcnt-ncnt);
    #if testflag
    read_align(4);
    //read_align(3);
    #else
    read_align(31);
    read_align(29);
    read_align(27);
    read_align(25);
    read_align(23);
    #endif
    print_time();

    //basket align
    vecr.clear();
    block.resize(rcnt0);
    block_rev.resize(rcnt0);
    block_bias.resize(rcnt0);
    bas_pre.resize(rcnt0,-1);
    block_pos.resize(rcnt0,-1);
    for(int i=0;i<rcnt0;i++){//the first edge of every read is outedge(bas_edge[i][0]),others are inedge
        if(isrepeat[i]==0&&hasn[i]==0&&find(pre,i)==i){
            vecr.pb(i);
            block[i].pb(i);
            block_rev[i].pb(0);
            block_bias[i].pb(0);
            block_pos[i]=0;
            bas_pre[i]=i;
        }
    }
    
    #if testflag
    bas_align(4);
    // cout<<"block info:------\n";
    // for(int i=0;i<vecr.size();i++){
    //     int t=vecr[i];
    //     if(block[t].size()>1){
    //         cout<<"block_root="<<t<<'\n';
    //         for(int j=0;j<block[t].size();j++){
    //             cout<<block[t][j]<<' '<<block_rev[t][j]<<' '<<block_bias[t][j]<<' '<<vread[block[t][j]]<<'\n';

    //         }
    //     }
    // }
    // cout<<"end----------\n";
    bas_align(3);
    #else
    bas_align(31);
    bas_align(29);
    bas_align(27);
    bas_align(25);
    bas_align(23);
    #endif


    
    cout<<"bas_align finish\n";
    print_time();
}


void read_align(int k)
{
    cout<<"read_align start,k="<<k<<'\n';
    cal_k_minimizer(k);

    sort(vecr.begin(),vecr.end(),[](int id1,int id2){
        if(val[id1]!=val[id2]) return val[id1]<val[id2];
        else return cal_len(id1)>cal_len(id2);
    });
    int tmp;
    mutex_init();
    for(int i=0;i<thread_num;i++){
        if((tmp=pthread_create(&thread_id[i],NULL,cal_pre,&thread_tmp[i]))!=0){
            cout<<"pthread_create ERROR\n";
            exit(0);
        }
    }
    cal_pre(&thread_tmp[thread_num]);
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
    //cout<<"cal_pre,thread_id="<<thread_id<<'\n';
    while(1){
        int curid;
        while(my_id<vecr.size()){
            if(find(pre,vecr[my_id])!=vecr[my_id]){
                my_id++;
                continue;
            }
            else if(find(pre,vecr[my_id])==vecr[my_id]) break;
            my_id++;
        }
        if(my_id>=vecr.size()){
            break;
        }
        curid=my_id++;
        int pos=vecr[curid];
        if(curid>=vecr.size()) break;
        for(int i=curid-1;i>=0&&i>=curid-300;i--){
            int ppos=vecr[i];
            if(val[pos]!=val[ppos]) break;
            if(cal_len(pos)!=cal_len(ppos)) break;
            if(find(pre,ppos)!=ppos) continue;
            if(check(pos,ppos,threshold1)){
                //merge two baskets,need lock, the total complexity is O(nlog^2n)
                //merge root
                pre[pos]=pre[ppos]=dsu(pos,ppos);
                //first_rev[pos]=iskmersymm[pos];//if multi,then the first_rev is cur rev state
                
                break;
            }
        }
    }
    return (void*)0;
}

int dsu(int id1,int id2)
{
    assert(pre[id1]==id1);
    assert(pre[id2]==id2);
    if(id1==id2) return id1;
    if(basket[id1].size()<basket[id2].size()) swap(id1,id2);
    int tmprev=iskmersymm[id1]^iskmersymm[id2];
    for(int i=0;i<basket[id2].size();i++){
        basket[id1].pb(basket[id2][i]);
        basket_rev[id1].pb(basket_rev[id2][i]^tmprev);
    }
    basket[id2].clear();
    basket_rev[id2].clear();
    return id1;
}


void bas_align(int k)
{
    cout<<"bas_align start,k="<<k<<'\n';
    curk=k;
    cal_k_minimizer(k);


    sort(vecr.begin(),vecr.end(),[](int id1,int id2){
        if(val[id1]!=val[id2]) return val[id1]<val[id2];
        else return cal_len(id1)>cal_len(id2);
    });
    
    nxt_bas.assign(vecr.size(),1e9);
    int tpre=1e9;
    for(int i=vecr.size()-2;i>=0;i--){
        if(val[vecr[i]]!=val[vecr[i+1]]) tpre=i+1;
        else if(cal_len(vecr[i])!=cal_len(vecr[i+1])) tpre=i+1;
        nxt_bas[i]=tpre;
    }

    int tmp;
    mutex_init();

    for(int i=0;i<thread_num;i++){
        if((tmp=pthread_create(&thread_id[i],NULL,cal_nxt,&thread_tmp[i]))!=0){
            cout<<"pthread_create ERROR\n";
            exit(0);
        }
    }
    cal_nxt(&thread_tmp[thread_num]);
    for(int i=0;i<thread_num;i++){
        void* tmp;
        pthread_join(thread_id[i],&tmp);
    }
    
}

void* cal_nxt(void* arg)
{
    int id=*(int*) arg;
    while(1){
        int curid;
        // while(my_id<vecr.size()){
        //     if(find(bas_pre,vecr[my_id])!=vecr[my_id]){
        //         my_id++;
        //         continue;
        //     }
        //     else if(find(bas_pre,vecr[my_id])==vecr[my_id]) break;
        //     my_id++;
        // }
        if(my_id>=vecr.size()){
            break;
        }
        curid=my_id++;
        //printf("pthread_id=%d,curid=%d\n",id,curid);
        int pos=vecr[curid];
        for(int i=nxt_bas[curid];i<vecr.size()&&i<=nxt_bas[curid]+300;i++){
            int npos=vecr[i];
            if(val[pos]!=val[npos]) break;
            if(cal_len(pos)-cal_len(npos)>max_bias) break;
            assert(cal_len(pos)>cal_len(npos));
            int tmp=find(bas_pre,npos);
            //if(tmp!=npos) continue;//not the root
            if(check(pos,npos,threshold2)){
                //merge two baskets,need lock, the total complexity is O(nlog^2n)
                bas_pre[find(bas_pre,pos)]=bas_pre[find(bas_pre,npos)]=bas_dsu(pos,npos);
                break;
                
            }
        }
    }
    return (void*)0;
}

//merge two blocks
// basic on two reads(maybe not root)
//the two reads id1 and id2 can match and has their bias and 
int bas_dsu(int id1,int id2)
{
    
    int k=curk;
    int r1=find(bas_pre,id1),r2=find(bas_pre,id2);
    if(r1==r2) return r1;
    //cout<<"bas_dsu:"<<id1<<' '<<id2<<' '<<r1<<' '<<r2<<'\n';
    //cout<<"block_pos:"<<block_pos[id1]<<' '<<block_pos[id2]<<'\n';
    if(block[r1].size()<block[r2].size()) swap(r1,r2),swap(id1,id2);
    int rev1=iskmersymm[id1]^block_rev[r1][block_pos[id1]],rev2=iskmersymm[id2]^block_rev[r2][block_pos[id2]];
    int len1=cal_len(id1)+qow(-1,rev1,2)*block_bias[r1][block_pos[id1]],len2=cal_len(id2)+qow(-1,rev2,2)*block_bias[r2][block_pos[id2]];
    //cout<<rev1<<' '<<rev2<<' '<<len1<<' '<<len2<<' '<<block_rev[r1][block_pos[id1]]<<' '<<block_rev[r2][block_pos[id2]]<<'\n';
    int tmprev=rev1^rev2;
    int tmpdir=qow(-1,tmprev,2);
    int basic_bias=(len1-len2)*qow(-1,rev1,2);
    //cout<<tmprev<<' '<<tmpdir<<' '<<basic_bias<<'\n';
    for(int i=0;i<block[r2].size();i++){
        block[r1].pb(block[r2][i]);
        block_rev[r1].pb(block_rev[r2][i]^tmprev);
        block_bias[r1].pb(block_bias[r2][i]*tmpdir+basic_bias);
        block_pos[block[r2][i]]=block[r1].size()-1;
    }
    block[r2].clear();
    block_rev[r2].clear();
    block_bias[r2].clear();
    //cout<<"bas_dsu finish,"<<r1<<'\n';
    return r1;

}